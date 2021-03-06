/*  Lattice Boltzmann sample, written in C
 *
 *  Copyright (C) 2018 Yuankun Fu
 *  E-mail: fu121@purdue.edu
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
*/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "lb.h"
#include "assert.h"

/*------------------------Two-steps alg modifcation------------------------*/
#include "boundaries.h"
extern PressureBCData* pressureBoundary;
extern int blk_size;
extern double *myrho1, *myrho2;
extern int iT, count;
extern int thread_block;
extern int NUM_THREADS;


extern long long global_CM[NUM_EVENTS];
/********************* ******************** ***********************/


// #define DEBUG_PRINT

/* D2Q9 lattice constants                                        */
/*****************************************************************/

  // lattice weights
static const double t[9] = { 4./9., 1./9., 1./9., 1./9., 1./9.,
                             1./36., 1./36., 1./36., 1./36. };
  // lattice velocities
static const int c[9][2] = {
    {0,0},
    {1,0}, {0,1}, {-1,0}, {0,-1},
    {1,1}, {-1,1}, {-1,-1}, {1,-1}
};


/* struct Node, methods                                          */
/*****************************************************************/

  // initialize a node to default value
void constructNode(Node* node) {
    int iPop;
    for (iPop=0; iPop<9; ++iPop) {
        node->fPop[iPop] = 0.;
    }
    node->dynamics = 0;
}

  // initialize a node to its local equilibrium term
void iniEquilibrium(Node* node, double rho, double ux, double uy) {
    int iPop;
    double uSqr = ux*ux + uy*uy;
    for (iPop=0; iPop<9; ++iPop) {
        node->fPop[iPop] =
            computeEquilibrium(iPop, rho, ux, uy, uSqr);
    }
}


/* struct Simulation, methods                                    */
/*****************************************************************/

  // allocate memory for a full simulation ("constructor")
void constructSim(Simulation* sim, int lx, int ly) {
    sim->lx = lx;
    sim->ly = ly;

    sim->memoryChunk    =
        (Node*) calloc((lx+2)*(ly+2), sizeof(Node));
    sim->tmpMemoryChunk =
        (Node*) calloc((lx+2)*(ly+2), sizeof(Node));

    sim->lattice        = (Node**) calloc(lx+2, sizeof(Node*));
    sim->tmpLattice     = (Node**) calloc(lx+2, sizeof(Node*));

    int iX, iY;
    for (iX=0; iX<lx+2; ++iX) {
        sim->lattice[iX] = sim->memoryChunk + iX*(ly+2);
        sim->tmpLattice[iX] = sim->tmpMemoryChunk + iX*(ly+2);
        for (iY=0; iY<ly+2; ++iY) {
            constructNode(&(sim->lattice[iX][iY]));
            constructNode(&(sim->tmpLattice[iX][iY]));
        }
    }
}

  // free the memory for the simulation ("destructor")
void destructSim(Simulation* sim) {
    free(sim->lattice);
    free(sim->memoryChunk);
    free(sim->tmpLattice);
    free(sim->tmpMemoryChunk);
    // free(sim->tmpMemoryChunk2);
}

  // specify the dynamics for a given lattice site
void setDynamics(Simulation* sim, int iX, int iY, Dynamics *dyn) {
    sim->lattice[iX][iY].dynamics = dyn;
    sim->tmpLattice[iX][iY].dynamics = dyn;
}

  // apply collision step to a lattice node (and simulate
  //   virtual dispatch)
inline static void collideNode(Node* node) {
    node->dynamics->dynamicsFun(node->fPop,
                                node->dynamics->selfData);
}

/* some free helper functions                                    */
/*****************************************************************/

  // compute density and velocity from the f's
void computeMacros(double* f, double* rho, double* ux, double* uy) {
    double upperLine  = f[2] + f[5] + f[6];
    double mediumLine = f[0] + f[1] + f[3];
    double lowerLine  = f[4] + f[7] + f[8];
    *rho = upperLine + mediumLine + lowerLine;
    *ux  = (f[1] + f[5] + f[8] - (f[3] + f[6] + f[7]))/(*rho);
    *uy  = (upperLine - lowerLine)/(*rho);
}

  // compute local equilibrium from rho and u
double computeEquilibrium(int iPop, double rho,
                          double ux, double uy, double uSqr)
{
    double c_u = c[iPop][0]*ux + c[iPop][1]*uy;
    return rho * t[iPop] * (
               1. + 3.*c_u + 4.5*c_u*c_u - 1.5*uSqr
           );
}

  // bgk collision term
void bgk(double* fPop, void* selfData) {
    double omega = *((double*)selfData);
    double rho, ux, uy;
    computeMacros(fPop, &rho, &ux, &uy);
    double uSqr = ux*ux+uy*uy;
    int iPop;

    for(iPop=0; iPop<9; ++iPop) {
        fPop[iPop] *= (1-omega);
        fPop[iPop] += omega * computeEquilibrium (
                                  iPop, rho, ux, uy, uSqr );
    }
}


/*-------------- Original LBM -----------------------------------*/
/*****************************************************************/
  // apply collision step to all lattice nodes
void collide(Simulation* sim) {
    int iX, iY, iPop;

#ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
#endif

    for (iX=1; iX<=sim->lx; ++iX) {
        for (iY=1; iY<=sim->ly; ++iY) {
            collideNode(&(sim->lattice[iX][iY]));
        }
    }

#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        int my_rank = 0;
        int thread_count = 1;

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
#endif

}

void collide_openmp(Simulation* sim) {
    int iX, iY, iPop;

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: global_CM)
{
    #ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
    #endif

    #pragma omp for private(iX, iY) schedule(static, thread_block)
    for (iX=1; iX<=sim->lx; ++iX)
        for (iY=1; iY<=sim->ly; ++iY)
            collideNode(&(sim->lattice[iX][iY]));

    #ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          int thread_count = omp_get_num_threads();
        #else
          int my_rank = 0;
          int thread_count = 1;
        #endif

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
    #endif
}
#else
    printf("No OPENMP used");
#endif

}

  // apply propagation step with help of temporary memory
void propagate(Simulation* sim) {
    int iX, iY, iPop;
    int lx = sim->lx;
    int ly = sim->ly;
    int nextX, nextY;

#ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
#endif

    for (iX=1; iX<=lx; ++iX) {
        for (iY=1; iY<=ly; ++iY) {
            for (iPop=0; iPop<9; ++iPop) {
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }
        }
    }
    // exchange lattice and tmplattice
    Node** swapLattice = sim->lattice;
    sim->lattice = sim->tmpLattice;
    sim->tmpLattice = swapLattice;

#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        int my_rank = 0;
        int thread_count = 1;

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
#endif

}

void propagate_openmp(Simulation* sim) {
    int iX, iY, iPop;
    int lx = sim->lx;
    int ly = sim->ly;
    int nextX, nextY;

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: global_CM)
{
    #ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
    #endif

    #pragma omp for private(iX, iY, nextX, nextY) schedule(static, thread_block)
    for (iX=1; iX<=lx; ++iX)
        for (iY=1; iY<=ly; ++iY)
            for (iPop=0; iPop<9; ++iPop) {
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }


    #ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          int thread_count = omp_get_num_threads();
        #else
          int my_rank = 0;
          int thread_count = 1;
        #endif

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
    #endif
}
#else
    printf("No OPENMP used");
#endif
    // exchange lattice and tmplattice
    Node** swapLattice = sim->lattice;
    sim->lattice = sim->tmpLattice;
    sim->tmpLattice = swapLattice;
}


void finalize_stream(Simulation* sim){
    Node** swapLattice = sim->lattice;
    sim->lattice = sim->tmpLattice;
    sim->tmpLattice = swapLattice;
}


/*----------------- Fused LBM -----------------------------------*/
/*****************************************************************/
void collide_with_stream(Simulation* sim) {
    int iX, iY, iPop;
    int nextX, nextY;

#ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
#endif

    for (iX=1; iX<=sim->lx; ++iX) {
        for (iY=1; iY<=sim->ly; ++iY) {

            collideNode(&(sim->lattice[iX][iY]));
            // streamming imediately once we got the updated f
            for (iPop=0; iPop<9; ++iPop) {
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }
        }
    }

#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        int my_rank = 0;
        int thread_count = 1;

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
#endif
}

void collide_with_stream_openmp(Simulation* sim) {
    int iX, iY, iPop;
    int nextX, nextY;

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: global_CM)
{
    #ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
    #endif

    #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, thread_block)
    for (iX=1; iX<=sim->lx; ++iX) {
        for (iY=1; iY<=sim->ly; ++iY) {
            collideNode(&(sim->lattice[iX][iY]));
            // streamming imediately once we got the updated f
            for (iPop=0; iPop<9; ++iPop) {
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }
        }
    }

    #ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        int my_rank = omp_get_thread_num();
        int thread_count = omp_get_num_threads();

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
    #endif
}
#else
    printf("No OPENMP used");
#endif

}

/*-------- Two Steps Line LBM -----------------------------------*/
/*****************************************************************/
// immediately compute iX-1, iY-1
void collide_tight(Simulation* sim) {
    int iX, iY, index, iPop;
    int lx=sim->lx, ly=sim->ly;
    double ux1, uy1, ux2, uy2;
    int nextX, nextY;

#ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
#endif

    for (iX=1; iX<=sim->lx; ++iX) {
        for (iY=1; iY<=sim->ly; ++iY) {

            // step1: collision on this line y
            collideNode(&(sim->lattice[iX][iY]));

            // step 2: stream from line x-1 to x
            for (iPop=0; iPop<9; ++iPop) {
                //iPop=below[index];
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }

            if(iX>1 && iY>1){

#ifdef ZGB
                //save rho
                if( iX==(lx-1) ){
                    //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                    computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho2[iY-1], &ux2, &uy2);
                }
                if( iX==lx ){
                    computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho1[iY-1], &ux1, &uy1);
                }
#endif
                // step 3: second collision on line x-1, y-1
                // should be based on the result of first stream
                // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
                collideNode(&(sim->tmpLattice[iX-1][iY-1]));

                // another branch for iX=sim->lx-1 and iY=sim-lx-2

                // step 4: second stream from  line y-1
                for (iPop=0; iPop<9; ++iPop) {
                    nextX = iX-1 + c[iPop][0];
                    nextY = iY-1 + c[iPop][1];
                    sim->lattice[nextX][nextY].fPop[iPop] =
                        sim->tmpLattice[iX-1][iY-1].fPop[iPop];
                }
            }
        }// end of iY loop
    }// end of iX loop

    //Line iX=1~lx-1, y=ly need to compute one more time
    // iY=sim->ly;

    for(iX=1; iX<sim->lx; ++iX){
        collideNode(&(sim->tmpLattice[iX][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX][ly].fPop[iPop];
        }
    }

    //Line iY=1~ly, iX=lx need to compute one more time
    // iX=sim->lx;
    //simple optimize
    iY=1;
    collideNode(&(sim->tmpLattice[lx][iY]));
    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][iY].fPop[iPop];
    }

    for (iY=2; iY<sim->ly; ++iY){

#ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
#endif
        collideNode(&(sim->tmpLattice[lx][iY]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    }


    //compute lx, ly point
    collideNode(&(sim->tmpLattice[lx][ly]));

    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = ly + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][ly].fPop[iPop];
    }

#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          int thread_count = omp_get_num_threads();
        #else
          int my_rank = 0;
          int thread_count = 1;
        #endif

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
#endif

}// end of func

/*----------Two Steps Line LBM (worse performance) --------------*/
/*****************************************************************/
//go through whole line, then compute lower line
void collide_tight2(Simulation* sim) {
    int iX, iY, index, iPop;
    int lx=sim->lx, ly=sim->ly;
    double ux1, uy1, ux2, uy2;
    int nextX, nextY;

    for (iX=1; iX<=sim->lx; ++iX) {

#ifdef _OPENMP
#pragma omp parallel for default(shared) \
    private(iY, iPop, nextX, nextY) \
    schedule(static, thread_block)
#endif
        for (iY=1; iY<=sim->ly; ++iY) {
            // step1: collision on this line y
            collideNode(&(sim->lattice[iX][iY]));

            // step 2: stream from line x-1 to x
            for (iPop=0; iPop<9; ++iPop) {
                //iPop=below[index];
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }
        }// end of iY loop

#ifdef _OPENMP
#pragma omp parallel for default(shared) \
    private(iY, iPop, nextX, nextY) \
    schedule(static, thread_block)
#endif
        for(iY=1; iY<=sim->ly; ++iY){
            if(iX>1 && iY>1){

#ifdef ZGB
                //save rho
                if( iX==(lx-1) ){
                    //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                    computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho2[iY-1], &ux2, &uy2);
                }
                if( iX==lx ){
                    computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho1[iY-1], &ux1, &uy1);
                }
#endif
                // step 3: second collision on line x-1, y-1
                // should be based on the result of first stream
                // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
                collideNode(&(sim->tmpLattice[iX-1][iY-1]));

                // another branch for iX=sim->lx-1 and iY=sim-lx-2

                // step 4: second stream from  line y-1
                for (iPop=0; iPop<9; ++iPop) {
                    nextX = iX-1 + c[iPop][0];
                    nextY = iY-1 + c[iPop][1];
                    sim->lattice[nextX][nextY].fPop[iPop] =
                        sim->tmpLattice[iX-1][iY-1].fPop[iPop];
                }
            }
        }

    }// end of iX loop

    //Line iX=1~lx-1, y=ly need to compute one more time
    // iY=sim->ly;
#ifdef _OPENMP
#pragma omp parallel for default(shared) \
    private(iX, iPop, nextX, nextY) \
    schedule(static, thread_block)
#endif
    for(iX=1; iX<sim->lx; ++iX){
        collideNode(&(sim->tmpLattice[iX][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX][ly].fPop[iPop];
        }
    }

    //Line iY=1~ly, iX=lx need to compute one more time
    // iX=sim->lx;
    //simple optimize
    iY=1;
    collideNode(&(sim->tmpLattice[lx][iY]));
    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][iY].fPop[iPop];
    }

#ifdef _OPENMP
#pragma omp parallel for default(shared) \
    private(iY, iPop, nextX, nextY) \
    schedule(static, thread_block)
#endif
    for (iY=2; iY<sim->ly; ++iY){

#ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
#endif
        collideNode(&(sim->tmpLattice[lx][iY]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    }


    //compute lx, ly point
    collideNode(&(sim->tmpLattice[lx][ly]));

    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = ly + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][ly].fPop[iPop];
    }

}// end of func

void collide_tight_openmp(Simulation* sim) {
    int iX, iY, index, iPop;
    int lx=sim->lx, ly=sim->ly;
    double ux1, uy1, ux2, uy2;
    int nextX, nextY;

    //compute each thread upper boundary line at iX=thread_block 1st c+s

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: global_CM)
{
    #ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
    #endif


    #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
    for(iX=thread_block; iX<=sim->lx; iX+=thread_block){
        for (iY=1; iY<=sim->ly; ++iY) {

            #ifdef DEBUG_PRINT
            #ifdef _OPENMP
            int my_rank = omp_get_thread_num();
            printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
            fflush(stdout);
            #endif
            #endif

            // step1: collision on this line y
            collideNode(&(sim->lattice[iX][iY]));

            // step 2: stream from line x-1 to x
            for (iPop=0; iPop<9; ++iPop) {
                //iPop=below[index];
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }
        }
    }

    #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static, thread_block)
    for (iX=1; iX<=sim->lx; ++iX) {
        for (iY=1; iY<=sim->ly; ++iY) {

            if(iX%thread_block != 0){

                #ifdef DEBUG_PRINT
                #ifdef _OPENMP
                int my_rank = omp_get_thread_num();
                printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
                fflush(stdout);
                #endif
                #endif

                // step1: collision on this line y
                collideNode(&(sim->lattice[iX][iY]));

                // step 2: stream from line x-1 to x
                for (iPop=0; iPop<9; ++iPop) {
                    //iPop=below[index];
                    nextX = iX + c[iPop][0];
                    nextY = iY + c[iPop][1];
                    sim->tmpLattice[nextX][nextY].fPop[iPop] =
                        sim->lattice[iX][iY].fPop[iPop];
                }
            }

            if(iX>1 && iY>1){

#ifdef ZGB
                //save rho
                if( iX==(lx-1) ){
                    //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                    computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho2[iY-1], &ux2, &uy2);
                }
                if( iX==lx ){
                    computeMacros(sim->tmpLattice[iX-1][iY-1].fPop, &myrho1[iY-1], &ux1, &uy1);
                }
#endif
                if( (iX-1)%thread_block != 0){

                    #ifdef DEBUG_PRINT
                    #ifdef _OPENMP
                        int my_rank = omp_get_thread_num();
                        printf("T%d: 2nd c+s on iX=%d, iY=%d will compute iX-1=%d, iY-1=%d\n", my_rank, iX, iY, iX-1, iY-1);
                        fflush(stdout);
                    #endif
                    #endif

                    // step 3: second collision on line x-1, y-1
                    // should be based on the result of first stream
                    // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
                    collideNode(&(sim->tmpLattice[iX-1][iY-1]));

                    // another branch for iX=sim->lx-1 and iY=sim-lx-2

                    // step 4: second stream from  line y-1
                    for (iPop=0; iPop<9; ++iPop) {
                        nextX = iX-1 + c[iPop][0];
                        nextY = iY-1 + c[iPop][1];
                        sim->lattice[nextX][nextY].fPop[iPop] =
                            sim->tmpLattice[iX-1][iY-1].fPop[iPop];
                    }
                }
            }
        }// end of iY loop
    }// end of iX loop

    //compute thread boundary line at iX=thread_block 2nd c+s
    //NOTICE: 1~ly-1 !!! use tmpLattice !!!
    #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
    for(iX=thread_block; iX<sim->lx; iX+=thread_block){
        for (iY=1; iY<=(sim->ly-1); ++iY) {

            #ifdef DEBUG_PRINT
            #ifdef _OPENMP
                int my_rank = omp_get_thread_num();
                printf("T%d: 2nd c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
                fflush(stdout);
            #endif
            #endif

            // step1: collision on this line y
            collideNode(&(sim->tmpLattice[iX][iY]));

            // step 2: stream from line x-1 to x
            for (iPop=0; iPop<9; ++iPop) {
                //iPop=below[index];
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->lattice[nextX][nextY].fPop[iPop] =
                    sim->tmpLattice[iX][iY].fPop[iPop];
            }
        }
    }


    //Line iX=1~lx-1, y=ly need to compute one more time
    // iY=sim->ly;
    #pragma omp for private(iX, iPop, nextX, nextY) schedule(static, thread_block)
    for(iX=1; iX<sim->lx; ++iX){
        collideNode(&(sim->tmpLattice[iX][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX][ly].fPop[iPop];
        }
    }

#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          int thread_count = omp_get_num_threads();
        #else
          int my_rank = 0;
          int thread_count = 1;
        #endif

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
#endif

}
#else
    printf("No OPENMP used");
#endif
 
    // compute
    // Line iY=1~ly, iX=lx need to compute one more time
    // iX=sim->lx;
    // #pragma omp critical
    // {
        iY=1;
        collideNode(&(sim->tmpLattice[lx][iY]));
        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    // }

    // #pragma omp for private(iY, iPop, nextX, nextY) schedule(static, thread_block)
    for (iY=2; iY<sim->ly; ++iY){

#ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
#endif
        collideNode(&(sim->tmpLattice[lx][iY]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    }

    // #pragma omp critical
    // {
        //compute lx, ly point
        collideNode(&(sim->tmpLattice[lx][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][ly].fPop[iPop];
        }
    // }

}// end of func


/*----------- Two Steps Blocking LBM ----------------------------*/
/*****************************************************************/
void collide_tight_block(Simulation* sim) {
    unsigned int iX, iY, iPop;
    // unsigned int blk_size = 32;
    int iix, iiy;
    int lx=sim->lx, ly=sim->ly;
    double ux1, uy1, ux2, uy2;
    int nextX, nextY;

#ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
#endif

    for (iX=1; iX<=sim->lx; iX+=blk_size) {
        for (iY=1; iY<=sim->ly; iY+=blk_size) {
            for(iix = 0; iix < blk_size; iix++){
                for(iiy = 0; iiy < blk_size; iiy++){

                    // step1: collision on this line y
                    collideNode(&(sim->lattice[iX+iix][iY+iiy]));

                    // step 2: stream from line x-1 to x
                    for (iPop=0; iPop<9; ++iPop) {
                        //iPop=below[index];
                        nextX = iX+iix + c[iPop][0];
                        nextY = iY+iiy + c[iPop][1];
                        sim->tmpLattice[nextX][nextY].fPop[iPop] =
                            sim->lattice[iX+iix][iY+iiy].fPop[iPop];
                    }

                    if(iX+iix>1 && iY+iiy>1){
#ifdef ZGB
                        //save rho
                        if( (iX+iix)==(lx-1) ){
                            //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                            computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho2[iY+iiy-1], &ux2, &uy2);
                        }
                        if( (iX+iix)==lx ){
                            computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho1[iY+iiy-1], &ux1, &uy1);
                        }
#endif
                        // step 3: second collision on line x-1, y-1
                        // should be based on the result of first stream
                        // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
                        collideNode(&(sim->tmpLattice[iX+iix-1][iY+iiy-1]));

                        // another branch for iX=sim->lx-1 and iY=sim-lx-2

                        // step 4: second stream from  line y-1
                        for (iPop=0; iPop<9; ++iPop) {
                            nextX = iX + iix -1 + c[iPop][0];
                            nextY = iY + iiy -1 + c[iPop][1];
                            sim->lattice[nextX][nextY].fPop[iPop] =
                                sim->tmpLattice[iX+iix-1][iY+ iiy -1].fPop[iPop];
                        }
                    }
                }
            }
        }// end of iY loop
    }// end of iX loop

    //Line iX=1~lx-1, y=ly need to compute one more time
    // iY=sim->ly;
    for(iX=1; iX<sim->lx; ++iX){
        collideNode(&(sim->tmpLattice[iX][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX][ly].fPop[iPop];
        }
    }

    //Line iY=1~ly, iX=lx need to compute one more time
    // iX=sim->lx;
    //simple optimize
    iY=1;
    collideNode(&(sim->tmpLattice[lx][iY]));
    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][iY].fPop[iPop];
    }

    for (iY=2; iY<sim->ly; ++iY){

#ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
#endif
        collideNode(&(sim->tmpLattice[lx][iY]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    }

    //compute lx, ly point
    collideNode(&(sim->tmpLattice[lx][ly]));

    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = ly + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][ly].fPop[iPop];
    }

#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          int thread_count = omp_get_num_threads();
        #else
          int my_rank = 0;
          int thread_count = 1;
        #endif

        int i;
        for(i=0; i<NUM_EVENTS; i++){
/*             printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);*/
             /*fflush(stdout);*/
            global_CM[i] += value_CM[i];
        }
#endif

}

void collide_tight_block_openmp(Simulation* sim) {
    unsigned int iX, iY, iPop;
    // unsigned int blk_size = 32;
    int iix, iiy;
    int lx=sim->lx, ly=sim->ly;
    double ux1, uy1, ux2, uy2;
    int nextX, nextY;

    //compute each thread upper boundary line at iX=thread_block 1st c+s
#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: global_CM)
{
    #ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
    #endif

    #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
    for(iX=thread_block; iX<=sim->lx; iX+=thread_block){
        for (iY=1; iY<=sim->ly; ++iY) {

            #ifdef DEBUG_PRINT
            int my_rank = omp_get_thread_num();
            printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
            fflush(stdout);
            #endif

            // step1: collision on this line y
            collideNode(&(sim->lattice[iX][iY]));

            // step 2: stream from line x-1 to x
            for (iPop=0; iPop<9; ++iPop) {
                //iPop=below[index];
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }
        }
    }


    // int schedule_thread_chunk = sim->lx/blk_size/NUM_THREADS;
    // printf("sim->lx=%d, blk_size=%d, NUM_THREADS=%d, \n", sim->lx, blk_size, NUM_THREADS);
    // fflush(stdout);
    #pragma omp for private(iix, iiy, iX, iY, iPop, nextX, nextY) schedule(static, thread_block/blk_size)
    for (iX=1; iX<=sim->lx; iX+=blk_size) {
        for (iY=1; iY<=sim->ly; iY+=blk_size) {
            for(iix = 0; iix < blk_size; iix++){
                for(iiy = 0; iiy < blk_size; iiy++){

                    if((iX+iix)%thread_block != 0){

                        #ifdef DEBUG_PRINT
                        #ifdef _OPENMP
                        int my_rank = omp_get_thread_num();
                        printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX+iix, iY+iiy);
                        fflush(stdout);
                        #endif
                        #endif
                        // step1: collision on this line y
                        collideNode(&(sim->lattice[iX+iix][iY+iiy]));

                        // step 2: stream from line x-1 to x
                        for (iPop=0; iPop<9; ++iPop) {
                            //iPop=below[index];
                            nextX = iX+iix + c[iPop][0];
                            nextY = iY+iiy + c[iPop][1];
                            sim->tmpLattice[nextX][nextY].fPop[iPop] =
                                sim->lattice[iX+iix][iY+iiy].fPop[iPop];
                        }
                    }

                    if(iX+iix>1 && iY+iiy>1){
#ifdef ZGB
                        //save rho
                        if( (iX+iix)==(lx-1) ){
                            //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                            computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho2[iY+iiy-1], &ux2, &uy2);
                        }
                        if( (iX+iix)==lx ){
                            computeMacros(sim->tmpLattice[iX+iix-1][iY+iiy-1].fPop, &myrho1[iY+iiy-1], &ux1, &uy1);
                        }
#endif
                        if( (iX+iix-1)%thread_block != 0){

                            #ifdef DEBUG_PRINT
                                int my_rank = omp_get_thread_num();
                                printf("T%d: 2nd c+s on iX=%d, iY=%d will compute iX-1=%d, iY-1=%d\n",
                                    my_rank, iX+iix, iY+iiy, iX+iix-1, iY+iiy-1);
                                fflush(stdout);
                            #endif
                            // step 3: second collision on line x-1, y-1
                            // should be based on the result of first stream
                            // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
                            collideNode(&(sim->tmpLattice[iX+iix-1][iY+iiy-1]));

                            // another branch for iX=sim->lx-1 and iY=sim-lx-2

                            // step 4: second stream from  line y-1
                            for (iPop=0; iPop<9; ++iPop) {
                                nextX = iX + iix -1 + c[iPop][0];
                                nextY = iY + iiy -1 + c[iPop][1];
                                sim->lattice[nextX][nextY].fPop[iPop] =
                                    sim->tmpLattice[iX+iix-1][iY+ iiy -1].fPop[iPop];
                            }
                        }
                    }
                }
            }
        }// end of iY loop
    }// end of iX loop

    //compute thread boundary line at iX=thread_block 2nd c+s
    //NOTICE: 1~ly-1 !!! use tmpLattice !!!
    #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
    for(iX=thread_block; iX<sim->lx; iX+=thread_block){
        for (iY=1; iY<=(sim->ly-1); ++iY) {

            #ifdef DEBUG_PRINT
                int my_rank = omp_get_thread_num();
                printf("T%d: 2nd c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
                fflush(stdout);
            #endif

            // step1: collision on this line y
            collideNode(&(sim->tmpLattice[iX][iY]));

            // step 2: stream from line x-1 to x
            for (iPop=0; iPop<9; ++iPop) {
                //iPop=below[index];
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->lattice[nextX][nextY].fPop[iPop] =
                    sim->tmpLattice[iX][iY].fPop[iPop];
            }
        }
    }

    //Line iX=1~lx-1, y=ly need to compute one more time
    // iY=sim->ly;
    #pragma omp for private(iX, iPop, nextX, nextY) schedule(static, thread_block)
    for(iX=1; iX<sim->lx; ++iX){
        collideNode(&(sim->tmpLattice[iX][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX][ly].fPop[iPop];
        }
    }
#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

          int my_rank = omp_get_thread_num();
          int thread_count = omp_get_num_threads();

        int i;
        for(i=0; i<NUM_EVENTS; i++){
             /*printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);*/
             /*fflush(stdout);*/
            global_CM[i] += value_CM[i];
        }
#endif

}
#else
    printf("No OPENMP used");
#endif
 
    //Line iY=1~ly, iX=lx need to compute one more time
    // iX=sim->lx;
    //simple optimize
    // #pragma omp critical
    // {
        iY=1;
        collideNode(&(sim->tmpLattice[lx][iY]));
        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    // }

    // #pragma omp for private(iY, iPop, nextX, nextY) schedule(static, thread_block)
    for (iY=2; iY<sim->ly; ++iY){

#ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
#endif
        collideNode(&(sim->tmpLattice[lx][iY]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    }

    //compute lx, ly point
    // #pragma omp critical
    // {
        collideNode(&(sim->tmpLattice[lx][ly]));
        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][ly].fPop[iPop];
        }
    // }


}

/*----------- Two steps Panel LBM -------------------------------*/
/*****************************************************************/
// Traversing on x axis
void collide_tight_panel_ix(Simulation* sim) {
    unsigned int iX, iY, iPop;
    // unsigned int blk_size = 32;
    int iix;
    int lx=sim->lx, ly=sim->ly;
    double ux1, uy1, ux2, uy2;
    int nextX, nextY;

    for (iX=1; iX<=sim->lx; iX+=blk_size) {
        for(iix = 0; iix< blk_size; iix++){
            for (iY=1; iY<=sim->ly; ++iY) {

                // step1: collision on this line y
                collideNode(&(sim->lattice[iX+iix][iY]));

                // step 2: stream from line x-1 to x
                for (iPop=0; iPop<9; ++iPop) {
                    //iPop=below[index];
                    nextX = iX+iix + c[iPop][0];
                    nextY = iY + c[iPop][1];
                    sim->tmpLattice[nextX][nextY].fPop[iPop] =
                        sim->lattice[iX+iix][iY].fPop[iPop];
                }

                if(iX+iix>1 && iY>1){
#ifdef ZGB
                    //save rho
                    if( (iX+iix)==(lx-1) ){
                        //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                        computeMacros(sim->tmpLattice[iX+iix-1][iY-1].fPop, &myrho2[iY-1], &ux2, &uy2);
                    }
                    if( (iX+iix)==(lx) ){
                        computeMacros(sim->tmpLattice[iX+iix-1][iY-1].fPop, &myrho1[iY-1], &ux1, &uy1);
                    }
#endif
                    // step 3: second collision on line x-1, y-1
                    // should be based on the result of first stream
                    // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
                    collideNode(&(sim->tmpLattice[iX+iix-1][iY-1]));

                    // another branch for iX=sim->lx-1 and iY=sim-lx-2

                    // step 4: second stream from  line y-1
                    for (iPop=0; iPop<9; ++iPop) {
                        nextX = iX +iix -1 + c[iPop][0];
                        nextY = iY      -1 + c[iPop][1];
                        sim->lattice[nextX][nextY].fPop[iPop] =
                            sim->tmpLattice[iX+iix-1][iY-1].fPop[iPop];
                    }
                }
            }// end of iY loop
        }// end of iix loop
    }// end of iX loop

    //Line iX=1~lx-1, y=ly need to compute one more time
    // iY=sim->ly;
    for(iX=1; iX<sim->lx; ++iX){
        collideNode(&(sim->tmpLattice[iX][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX][ly].fPop[iPop];
        }
    }

    // Line iY=1~ly, iX=lx need to compute one more time
    // iX=sim->lx;
    //simple optimize
    iY=1;
    collideNode(&(sim->tmpLattice[lx][iY]));
    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][iY].fPop[iPop];
    }

    for (iY=2; iY<sim->ly; ++iY){

#ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
#endif
        collideNode(&(sim->tmpLattice[lx][iY]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    }

    //compute lx, ly point
    collideNode(&(sim->tmpLattice[lx][ly]));

    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = ly + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][ly].fPop[iPop];
    }
}

// Traversing on y axis
void collide_tight_panel_iy(Simulation* sim) {
    unsigned int iX, iY, iPop;
    // unsigned int blk_size = 32;
    int iiy;
    int lx=sim->lx, ly=sim->ly;
    double ux1, uy1, ux2, uy2;
    int nextX, nextY;

#ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
#endif

    for (iY=1; iY<=sim->ly; iY+=blk_size) {
        for (iX=1; iX<=sim->lx; iX++) {
            for(iiy = 0; iiy < blk_size; iiy ++){

                // printf("1st c+s: iX=%d, iY=%d\n", iX, iY+iiy);
                // fflush(stdout);
                // step1: collision on this line y
                collideNode(&(sim->lattice[iX][iY+iiy]));

                // step 2: stream from line x-1 to x
                for (iPop=0; iPop<9; ++iPop) {
                    //iPop=below[index];
                    nextX = iX     + c[iPop][0];
                    nextY = iY+iiy + c[iPop][1];
                    sim->tmpLattice[nextX][nextY].fPop[iPop] =
                        sim->lattice[iX][iY+iiy].fPop[iPop];
                }


                if(iX>1 && iY+iiy>1){
#ifdef ZGB
                    //save rho
                    if( iX==(lx-1) ){
                        //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                        computeMacros(sim->tmpLattice[iX-1][iY+iiy-1].fPop, &myrho2[iY+iiy-1], &ux2, &uy2);
                    }
                    if( iX==(lx) ){
                        computeMacros(sim->tmpLattice[iX-1][iY+iiy-1].fPop, &myrho1[iY+iiy-1], &ux1, &uy1);
                    }
#endif
                    // printf("2nd c+s: iX=%d, iY=%d\n", iX-1, iY+iiy-1);
                    // fflush(stdout);
                    // step 3: second collision on line x-1, y-1
                    // should be based on the result of first stream
                    // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
                    collideNode(&(sim->tmpLattice[iX-1][iY+iiy-1]));

                    // another branch for iX=sim->lx-1 and iY=sim-lx-2

                    // step 4: second stream from  line y-1
                    for (iPop=0; iPop<9; ++iPop) {
                        nextX = iX       -1 + c[iPop][0];
                        nextY = iY + iiy -1 + c[iPop][1];
                        sim->lattice[nextX][nextY].fPop[iPop] =
                            sim->tmpLattice[iX-1][iY+ iiy -1].fPop[iPop];
                    }
                }
            }
        }// end of iX loop
    }// end of iY loop

    //Line iX=1~lx-1, y=ly need to compute one more time
    // iY=sim->ly;
    for(iX=1; iX<sim->lx; ++iX){
        collideNode(&(sim->tmpLattice[iX][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX][ly].fPop[iPop];
        }
    }

    // Line iY=1~ly, iX=lx need to compute one more time
    // iX=sim->lx;
    //simple optimize
    iY=1;
    collideNode(&(sim->tmpLattice[lx][iY]));
    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = iY + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][iY].fPop[iPop];
    }

    for (iY=2; iY<sim->ly; ++iY){

#ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
#endif
        collideNode(&(sim->tmpLattice[lx][iY]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    }

    //compute lx, ly point
    collideNode(&(sim->tmpLattice[lx][ly]));

    for (iPop=0; iPop<9; ++iPop) {
        nextX = lx + c[iPop][0];
        nextY = ly + c[iPop][1];
        sim->lattice[nextX][nextY].fPop[iPop] =
            sim->tmpLattice[lx][ly].fPop[iPop];
    }

#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          int thread_count = omp_get_num_threads();
        #else
          int my_rank = 0;
          int thread_count = 1;
        #endif

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
#endif

}

void collide_tight_panel_iy_openmp(Simulation* sim) {
    unsigned int iX, iY, iPop;
    // unsigned int blk_size = 32;
    int iiy;
    int lx=sim->lx, ly=sim->ly;
    double ux1, uy1, ux2, uy2;
    int nextX, nextY;

#ifdef _OPENMP
#pragma omp parallel default(shared) reduction(+: global_CM)
{
    #ifdef ADDPAPI
        long long value_CM[NUM_EVENTS];
        int retval;
        retval = PAPI_start_counters(EventSet, NUM_EVENTS);
    #endif

    //compute each thread right boundary line at iY=thread_block 1st c+s
    #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
    for(iY=thread_block; iY<=sim->ly; iY+=thread_block){
        for (iX=1; iX<=sim->lx; ++iX) {

            #ifdef DEBUG_PRINT
            #ifdef _OPENMP
                int my_rank = omp_get_thread_num();
                printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
                fflush(stdout);
            #endif
            #endif

            // step1: collision on this line y
            collideNode(&(sim->lattice[iX][iY]));

            // step 2: stream from line x-1 to x
            for (iPop=0; iPop<9; ++iPop) {
                //iPop=below[index];
                nextX = iX + c[iPop][0];
                nextY = iY + c[iPop][1];
                sim->tmpLattice[nextX][nextY].fPop[iPop] =
                    sim->lattice[iX][iY].fPop[iPop];
            }
        }
    }

    // int schedule_thread_chunk = sim->lx/blk_size/NUM_THREADS;
    // printf("sim->lx=%d, blk_size=%d, NUM_THREADS=%d, \n", sim->lx, blk_size, NUM_THREADS);
    // fflush(stdout);
    #pragma omp for private(iiy, iX, iY, iPop, nextX, nextY) schedule(static, thread_block/blk_size)
    for (iY=1; iY<=sim->ly; iY+=blk_size) {
        for (iX=1; iX<=sim->lx; iX++) {
            for(iiy = 0; iiy < blk_size; iiy ++){

                if((iY+iiy)%thread_block != 0){

                    #ifdef DEBUG_PRINT
                    #ifdef _OPENMP
                        int my_rank = omp_get_thread_num();
                        printf("T%d: 1st c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
                        fflush(stdout);
                    #endif
                    #endif

                    // step1: collision on this line y
                    collideNode(&(sim->lattice[iX][iY+iiy]));

                    // step 2: stream from line x-1 to x
                    for (iPop=0; iPop<9; ++iPop) {
                        //iPop=below[index];
                        nextX = iX     + c[iPop][0];
                        nextY = iY+iiy + c[iPop][1];
                        sim->tmpLattice[nextX][nextY].fPop[iPop] =
                            sim->lattice[iX][iY+iiy].fPop[iPop];
                    }
                }

                if(iX>1 && iY+iiy>1){
#ifdef ZGB
                    //save rho
                    if( iX==(lx-1) ){
                        //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                        computeMacros(sim->tmpLattice[iX-1][iY+iiy-1].fPop, &myrho2[iY+iiy-1], &ux2, &uy2);
                    }
                    if( iX==(lx) ){
                        computeMacros(sim->tmpLattice[iX-1][iY+iiy-1].fPop, &myrho1[iY+iiy-1], &ux1, &uy1);
                    }
#endif
                    if((iY+iiy-1)%thread_block != 0){

                        #ifdef DEBUG_PRINT
                        #ifdef _OPENMP
                            int my_rank = omp_get_thread_num();
                            printf("T%d: 2nd c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
                            fflush(stdout);
                        #endif
                        #endif

                        // step 3: second collision on line x-1, y-1
                        // should be based on the result of first stream
                        // how to get velocity from direction 6 and 5(need 1 offset in x direction too)?
                        collideNode(&(sim->tmpLattice[iX-1][iY+iiy-1]));

                        // another branch for iX=sim->lx-1 and iY=sim-lx-2

                        // step 4: second stream from  line y-1
                        for (iPop=0; iPop<9; ++iPop) {
                            nextX = iX       -1 + c[iPop][0];
                            nextY = iY + iiy -1 + c[iPop][1];
                            sim->lattice[nextX][nextY].fPop[iPop] =
                                sim->tmpLattice[iX-1][iY+ iiy -1].fPop[iPop];
                        }
                    }
                }
            } // end of iiy
        }// end of iX loop
    }// end of iY loop

    //compute thread boundary line at iY=thread_block~ly-thread_block 2nd c+s
    //NOTICE: 1~lx !!! use tmpLattice !!!
    #pragma omp for private(iX, iY, iPop, nextX, nextY) schedule(static)
    for(iY=thread_block; iY<sim->ly; iY+=thread_block){
        for (iX=1; iX<=sim->lx; ++iX) {
#ifdef ZGB
            //save rho
            if( iX==(lx-1) ){
                //store rho from column iX=lx-2, iY=2~ly-1 need to be computed; iY=1, ly also computed but not used
                computeMacros(sim->tmpLattice[iX-1][iY].fPop, &myrho2[iY], &ux2, &uy2);
            }
            if( iX==(lx) ){
                computeMacros(sim->tmpLattice[iX-1][iY].fPop, &myrho1[iY], &ux1, &uy1);
            }
#endif
            #ifdef DEBUG_PRINT
            #ifdef _OPENMP
                int my_rank = omp_get_thread_num();
                printf("T%d: 2nd c+s on iX=%d, iY=%d\n", my_rank, iX, iY);
                fflush(stdout);
            #endif
            #endif

            if(iX!=sim->lx){
                // step1: collision on this line y
                collideNode(&(sim->tmpLattice[iX][iY]));

                // step 2: stream from line x-1 to x
                for (iPop=0; iPop<9; ++iPop) {
                    //iPop=below[index];
                    nextX = iX + c[iPop][0];
                    nextY = iY + c[iPop][1];
                    sim->lattice[nextX][nextY].fPop[iPop] =
                        sim->tmpLattice[iX][iY].fPop[iPop];
                }
            }
        }
    }

    //Line iX=1~lx-1, y=ly need to compute one more time
    // iY=sim->ly;
    #pragma omp for private(iX, iPop, nextX, nextY) schedule(static, thread_block)
    for(iX=1; iX<sim->lx; ++iX){
        collideNode(&(sim->tmpLattice[iX][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = iX + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[iX][ly].fPop[iPop];
        }
    }

#ifdef ADDPAPI
        retval=PAPI_stop_counters(value_CM, NUM_EVENTS);

        #ifdef _OPENMP
          int my_rank = omp_get_thread_num();
          int thread_count = omp_get_num_threads();
        #else
          int my_rank = 0;
          int thread_count = 1;
        #endif

        int i;
        for(i=0; i<NUM_EVENTS; i++){
            // printf("T%d: event[%d]=%lld\n", my_rank, i, value_CM[i]);
            // fflush(stdout);
            global_CM[i] += value_CM[i];
        }
#endif

}//end pragma parallel
#else
    printf("No OPENMP used");
#endif
 
    // Line iY=1~ly, iX=lx need to compute one more time
    // iX=sim->lx;
    //simple optimize
    // #pragma omp critical
    // {
        iY=1;
        collideNode(&(sim->tmpLattice[lx][iY]));
        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    // }

    // #pragma omp for private(iY, iPop, nextX, nextY) schedule(static, thread_block)
    for (iY=2; iY<sim->ly; ++iY){

#ifdef ZGB
        //Compute a second order extrapolation on the right boundary
        pressureBoundary[iY].rho = 4./3.* myrho1[iY] - 1./3.* myrho2[iY];
        pressureBoundary[iY].uPar = 0.;
#endif
        collideNode(&(sim->tmpLattice[lx][iY]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = iY + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][iY].fPop[iPop];
        }
    }

    //compute lx, ly point
    // #pragma omp critical
    // {
        collideNode(&(sim->tmpLattice[lx][ly]));

        for (iPop=0; iPop<9; ++iPop) {
            nextX = lx + c[iPop][0];
            nextY = ly + c[iPop][1];
            sim->lattice[nextX][nextY].fPop[iPop] =
                sim->tmpLattice[lx][ly].fPop[iPop];
        }
    // }

}

  // implement periodic boundary conditions (to be called after
  //   the propagation step)
void makePeriodic(Simulation* sim) {
    int lx = sim->lx;
    int ly = sim->ly;
    Node** lat = sim->lattice;

    int iX, iY;
    for (iX=1; iX<=lx; ++iX) {
        //1
        lat[iX][ly].fPop[4] = lat[iX][0].fPop[4];
        lat[iX][ly].fPop[7] = lat[iX][0].fPop[7];
        lat[iX][ly].fPop[8] = lat[iX][0].fPop[8];

        //2
        lat[iX][1].fPop[2] = lat[iX][ly+1].fPop[2];
        lat[iX][1].fPop[5] = lat[iX][ly+1].fPop[5];
        lat[iX][1].fPop[6] = lat[iX][ly+1].fPop[6];
    }

    for (iY=1; iY<=ly; ++iY) {
        //3
        lat[1][iY].fPop[1] = lat[lx+1][iY].fPop[1];
        lat[1][iY].fPop[5] = lat[lx+1][iY].fPop[5];
        lat[1][iY].fPop[8] = lat[lx+1][iY].fPop[8];

        //4
        lat[lx][iY].fPop[3] = lat[0][iY].fPop[3];
        lat[lx][iY].fPop[6] = lat[0][iY].fPop[6];
        lat[lx][iY].fPop[7] = lat[0][iY].fPop[7];
    }

    lat[1][1].fPop[5]   = lat[lx+1][ly+1].fPop[5];
    lat[lx][1].fPop[6]  = lat[0][ly+1].fPop[6];
    lat[lx][ly].fPop[7] = lat[0][0].fPop[7];
    lat[1][ly].fPop[8]  = lat[lx+1][0].fPop[8];
}

/*---- Output for generating Paraview video ---------------------*/
/*****************************************************************/
  // save the velocity field (norm) to disk
void saveVel(Simulation* sim, char fName[]) {
    FILE* oFile = fopen(fName, "w");
    int iX, iY;
    double ux, uy, uNorm, rho;
    double tmp=1.0;

    // compute norm to verify results & generate figure for MATLAB
#ifdef OUTPUT_MATLAB
     for (iY=1; iY<=sim->ly; ++iY) {
        for (iX=1; iX<=sim->lx; ++iX) {
            computeMacros(sim->lattice[iX][iY].fPop, &rho, &ux, &uy);
            uNorm = sqrt(ux*ux+uy*uy);
            fprintf(oFile, "%f ", uNorm);
        }
        fprintf(oFile, "\n");
     }
#endif

    // save vx, vy for paraview+catalyst
    for (iY=1; iY<=sim->ly; ++iY) {
       for (iX=1; iX<=sim->lx; ++iX) {
           computeMacros(sim->lattice[iX][iY].fPop, &rho, &ux, &uy);
           fprintf(oFile, "%f,%f,%f\n", ux, tmp, uy);
       }
    }

    fclose(oFile);
}

  // save one lattice population to disk
void saveF(Simulation* sim, int iPop, char fName[]) {
    FILE* oFile = fopen(fName, "w");
    int iX, iY;
    double ux, uy, uNorm, rho;

    for (iY=1; iY<=sim->ly; ++iY) {
       for (iX=1; iX<=sim->lx; ++iX) {
           double f = sim->lattice[iX][iY].fPop[iPop];
           fprintf(oFile, "%f ", f);
       }
       fprintf(oFile, "\n");
    }

    fclose(oFile);
}
