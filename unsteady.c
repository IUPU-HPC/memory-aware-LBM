/*  Lattice Boltzmann sample, written in C
 *
 *  Copyright (C) 2006 Jonas Latt
 *  Address: Rue General Dufour 24,  1211 Geneva 4, Switzerland
 *  E-mail: Jonas.Latt@cui.unige.ch
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

/* unsteady.c:
 * This example examines an unsteady flow past a cylinder placed in a channel.
 * The cylinder is offset somewhat from the center of the flow to make the
 * steady-state symmetrical flow unstable. At the inlet, a Poiseuille profile is
 * imposed on the velocity, where as the outlet implements an outflow condition:
 * grad_x u = 0. At Reynolds numbers around 100, an unstable periodic pattern is
 * created, the Karman vortex street.
 */

#include "lb.h"
#include "boundaries.h"
#include <stdio.h>
#include <stdlib.h>
#include "eval_tools.h"

#define NO_SAVE

  // These constants define the flow geometry and are commented in
  //   the function setConstants()
int lx, ly;
int obst_x, obst_y, obst_r;
double uMax, Re, nu, omega;
int maxT, tSave;
int NUM_THREADS;

/********************* added by Yuankun Fu ***********************/
int blk_size;
double *myrho1, *myrho2;
int iT=0, count;
int thread_block;
// #define TIGHT2

/*#ifdef ADDPAPI*/
  // int EventSet[] = {PAPI_TLB_DM, PAPI_L1_DCM, PAPI_L2_DCM, PAPI_L3_DCM};
  long long global_CM[NUM_EVENTS];
  // int EventSet[] = {PAPI_L1_DCM, PAPI_L2_DCM, PAPI_TLB_DM};
/*#endif*/
/********************* ******************** ***********************/

  // The dynamics that are to be specified on different regions of
  //   the flow: LBGK in the bulk, bounce-back on the cylinder,
  //   regularized boundary condition on the four domain boundaries.
  //   Alternatively, the Zou/He boundary condition can be used:
  //   replace upperRegularized by upperZouHe etc.
Dynamics       bulkDynamics         = { &bgk, (void*)&omega };
Dynamics       bounceBackDynamics   = { &bounceBack, 0 };

  // The velocity to be imposed on the upper/lower boundary: u=0
VelocityBCData zeroVelocityBoundary = { &bulkDynamics, 0., 0. };

Dynamics upperBoundary = { &upperRegularized,
                           (void*)&zeroVelocityBoundary };
Dynamics lowerBoundary = { &lowerRegularized,
                           (void*)&zeroVelocityBoundary };
Dynamics* leftBoundary;   // Those two objects are initialized in the
Dynamics* rightBoundary;  //   function iniData()

  // These arrays contain the velocities that are to be imposed on
  //   the inlet (poiseuilleBoundary) and the outlet
  //   (pressureBoundary) of the channel
VelocityBCData* poiseuilleBoundary;
PressureBCData* pressureBoundary;

  // The main object containing the simulation
Simulation sim;

void setConstants(int argc, char *argv[]) {
    /*lx       = 250;     // channel length
    ly       = 50;      // channel height
    obst_r = lx/10+1;   // radius of the cylinder
    obst_x = lx/5;
    obst_y = ly/2;*/

    lx = atoi(argv[1]);     //channel lenghth
    ly = atoi(argv[2]);     //channel width
    blk_size = atoi(argv[3]);

    /*obst_r = 8;   // radius of the cylinder*/
    obst_r = ly/10+1;   // radius of the cylinder
    obst_x = lx/4;      // position of the cylinder; the cylinder is
    obst_y = ly/2;      // offset from the center to break symmetry

    uMax  = 0.02;       // maximum velocity of the Poiseuille inflow
    Re    = 100;        // Reynolds number
    nu    = uMax * 2.*obst_r / Re;  // kinematic fluid viscosity
    omega = 1. / (3*nu+1./2.);      // relaxation parameter

    maxT   = atoi(getenv("STOP"));       // total number of iterations
    tSave  = 100;          // frequency of periodic saves to disk

    printf("\nlx=%d, ly=%d, omega=%f, blk_size=%d\n\n", lx, ly, omega, blk_size);
}

  // Memory allocation and default initialisation of the simulation
  //   and the left/right boundaries
void iniData() {
    int iX, iY;
    poiseuilleBoundary =
        (VelocityBCData*) calloc(ly+2, sizeof(VelocityBCData));
    pressureBoundary =
        (PressureBCData*) calloc(ly+2, sizeof(PressureBCData));

    leftBoundary  = (Dynamics*) calloc(ly+2, sizeof(Dynamics));
    rightBoundary = (Dynamics*) calloc(ly+2, sizeof(Dynamics));

#ifdef ZGB
    //add by Yuankun
    myrho1 = (double *) calloc(ly, sizeof(double));
    myrho2 = (double *) calloc(ly, sizeof(double));
#endif

    for (iY=2; iY<=ly-1; ++iY) {
        poiseuilleBoundary[iY].bulkDynamics   = &bulkDynamics;
        poiseuilleBoundary[iY].uy             = 0.;
        pressureBoundary[iY].bulkDynamics = &bulkDynamics;
        pressureBoundary[iY].rho          = 1.;
        pressureBoundary[iY].uPar         = 0.;

        leftBoundary[iY].dynamicsFun = &leftRegularized;
        leftBoundary[iY].selfData
            = (void*)&poiseuilleBoundary[iY];
        rightBoundary[iY].dynamicsFun = &rightPressureRegularized;
        rightBoundary[iY].selfData
            = (void*)&pressureBoundary[iY];
    }
}

  // De-allocation of the memory
void freeData() {
    free(rightBoundary);
    free(leftBoundary);
    free(poiseuilleBoundary);
    free(pressureBoundary);

    //add by Yuankun
    free(myrho1);
    free(myrho2);
}

  // compute parabolic Poiseuille profile
double computePoiseuille(int iY) {
    double y = (double)(iY-1);
    double L = (double)(ly-1);
    return 4.*uMax / (L*L) * (L*y-y*y);
}

  // Specify the geometry of the simulation
void iniGeometry() {
    int iX, iY;
    for(iX=1; iX<=lx; ++iX) {
        for(iY=1; iY<=ly; ++iY) {
              // All bulk nodes are initialized at equilibrium with constant
              // density and a velocity determined by a y-dependend Poiseuille
              // profile.
            double uPoiseuille = computePoiseuille(iY);
            iniEquilibrium(&sim.lattice[iX][iY], 1., uPoiseuille, 0.);
              // on the obstacle, set bounce back dynamics
            if ( (iX-obst_x)*(iX-obst_x) +
                 (iY-obst_y)*(iY-obst_y) <= obst_r*obst_r )
            {
                setDynamics(&sim, iX, iY, &bounceBackDynamics);
            }
            //   //add by Yuankun
            // else if ( (iY==1) || (iY==ly) ){
            //     setDynamics(&sim, iX, iY, &bounceBackDynamics);
            // }
              // elsewhere, use lbgk dynamics
            else {
                setDynamics(&sim, iX, iY, &bulkDynamics);
            }
        }
    }

      // upper and lower boundary: u=0
    for (iX=1; iX<=lx; ++iX) {
        setDynamics(&sim, iX, 1, &lowerBoundary);
        setDynamics(&sim, iX, ly, &upperBoundary);
    }

      // left boundary: Poiseuille profile, constant through time
      // right boundary: initially Poiseuille profile, then outlet
      //   condition grad_x u = 0
    for (iY=2; iY<=ly-1; ++iY) {
        double uPoiseuille = computePoiseuille(iY);
        poiseuilleBoundary[iY].ux   = uPoiseuille;
        setDynamics(&sim, 1, iY, &leftBoundary[iY]);
        setDynamics(&sim, lx, iY, &rightBoundary[iY]);
    }
}

  // Compute a second order extrapolation on the right boundary to
  // ensure a zero-gradient boundary condition on the pressure.
  // This must be recomputed at every time step. The velocity is
  // constrained to be perpendicular to the outflow surface.
void updateZeroGradientBoundary() {
    int iY;
    double rho1, ux1, uy1, rho2, ux2, uy2;
    for (iY=2; iY<=ly-1; ++iY) {
        computeMacros(sim.lattice[lx-1][iY].fPop, &rho1, &ux1, &uy1);
        computeMacros(sim.lattice[lx-2][iY].fPop, &rho2, &ux2, &uy2);
        pressureBoundary[iY].rho = 4./3.*rho1 - 1./3.*rho2;
        pressureBoundary[iY].uPar = 0.; //uy=0

        // if(iT==1 || iT==2){
        //     printf("in ZGB, iT=%d, rho2[%d]=%f, rho1[%d]=%f\n", iT, iY, rho2, iY, rho1);
        //     fflush(stdout);
        // }
    }
}

static const double t[9] = { 4./9., 1./9., 1./9., 1./9., 1./9.,
                             1./36., 1./36., 1./36., 1./36. };
  // lattice velocities
static const int c[9][2] = {
    {0,0},
    {1,0}, {0,1}, {-1,0}, {0,-1},
    {1,1}, {-1,1}, {-1,-1}, {1,-1}
};

//OpenMP test hello
void Hello(void){

#ifdef _OPENMP
  int my_rank = omp_get_thread_num();
  int thread_count = omp_get_num_threads();
#else
  int my_rank = 0;
  int thread_count = 1;
#endif

  printf("Hello from thread %d of %d\n", my_rank, thread_count);
}



int main(int argc, char *argv[]) {
    int i;
    char case_name[80], filename[80];

    void (*constructSim_func)(Simulation*, int, int) = NULL;
    void (*collision_func)(Simulation *) = NULL;
    void (*stream_func)(Simulation *) = NULL;

      // initialisation of a lx * ly simulation
    setConstants(argc, argv);
    int interval = maxT/10;
    // int interval = 1;


#ifdef _OPENMP
    NUM_THREADS = atoi(getenv("OMP_NUM_THREADS"));
    // chunk_size=(lx<ly?lx:ly)/NUM_THREADS;
    // chunk_size = atoi(argv[4]);
    // thread_block = lx/NUM_THREADS;
    thread_block = atoi(argv[4]);
    printf("thread_block=%d\n", thread_block);
    omp_set_num_threads(NUM_THREADS);
#endif

#ifdef _OPENMP
#pragma omp parallel
#endif
    Hello();

#if defined(COMBINE)
        sprintf(case_name, "combine");
        #ifdef ADDPAPI
        sprintf(case_name, "combine_papi");
        #endif
        constructSim_func = &constructSim;
        collision_func=&collide_with_stream;

  #ifdef _OPENMP
        sprintf(case_name, "combine_openmp");
        #ifdef ADDPAPI
        sprintf(case_name, "combine_openmp_papi");
        #endif
        collision_func=&collide_with_stream_openmp;
  #endif
        stream_func=&finalize_stream;

#elif defined(QUICKTEST)
    #define TWOSTEPS
        sprintf(case_name, "quicktest");
        constructSim_func = &constructSim;
        collision_func=&collide_with_stream_twice;
        stream_func=&finalize_stream;

#elif defined(TIGHT)
    #define TWOSTEPS
        sprintf(case_name, "tight");
        #ifdef ADDPAPI
        sprintf(case_name, "tight_papi");
        #endif
        constructSim_func = &constructSim;
        collision_func=&collide_tight;

    #ifdef _OPENMP
        sprintf(case_name, "tight_openmp");
        #ifdef ADDPAPI
        sprintf(case_name, "tight_openmp_papi");
        #endif
        collision_func=&collide_tight_openmp;
    #endif

    #ifdef TIGHT2
        sprintf(case_name, "tight2");
        collision_func=&collide_tight2;
    #endif
        // no need to swap!
        stream_func=&finalize_stream;

#elif defined(TIGHT_BLOCK)
    #define TWOSTEPS
        sprintf(case_name, "tight_block");
        #ifdef ADDPAPI
        sprintf(case_name, "tight_block_papi");
        #endif
        constructSim_func = &constructSim_blk;
        collision_func=&collide_tight_block;
    #ifdef _OPENMP
        sprintf(case_name, "tight_block_openmp");
        #ifdef ADDPAPI
        sprintf(case_name, "tight_block_openmp_papi");
        #endif
        collision_func=&collide_tight_block_openmp;
    #endif
        // no need to swap!
        stream_func=&finalize_stream;

#elif defined(PANEL)
    #define TWOSTEPS
        sprintf(case_name, "panel");
        #ifdef ADDPAPI
        sprintf(case_name, "panel_papi");
        #endif
        constructSim_func = &constructSim;
        // collision_func=&collide_tight_panel_ix;
        collision_func=&collide_tight_panel_iy;
    #ifdef _OPENMP
        sprintf(case_name, "panel_openmp");
        #ifdef ADDPAPI
        sprintf(case_name, "panel_openmp_papi");
        #endif
        collision_func=&collide_tight_panel_iy_openmp;
    #endif
        // no need to swap!
        stream_func=&finalize_stream;

#else
        sprintf(case_name, "origin");
        #ifdef ADDPAPI
        sprintf(case_name, "origin_papi");
        #endif
        constructSim_func = &constructSim;
        collision_func=&collide;
        stream_func=&propagate;
    #ifdef _OPENMP
        sprintf(case_name, "origin_openmp");
        #ifdef ADDPAPI
        sprintf(case_name, "origin_openmp_papi");
        #endif
        collision_func=&collide_openmp;
        stream_func=&propagate_openmp;
    #endif

#endif

    printf("use case %s\n", case_name);
    fflush(stdout);

    // allocate spaces for boundries
    iniData();
     printf("Pass iniData\n");
    // fflush(stdout);

    // allocae space for nodes and lattice
    // set configurations
    constructSim_func(&sim, lx, ly);
    printf("Pass constructSim\n");
    // fflush(stdout);

    // bounce back dynamics in obstacles, else use lbgk(bulk dynamics)
    iniGeometry();
    printf("Pass iniGeometry, obst_x=%d, obst_y=%d, obst_r=%d\n",
             obst_x, obst_y, obst_r);
    fflush(stdout);

    double t_start, t_end, t0, t1, t_collid=0.0, t_zgb=0.0, t_stream=0.0;
    t_start = get_cur_time();

#ifdef ADDPAPI
  int retval = PAPI_library_init (PAPI_VER_CURRENT);
  if (retval != PAPI_VER_CURRENT) {
    printf("PAPI_library_init error!\n");
    exit(1);
  }

  int thread_count;
  #ifdef _OPENMP
    if (PAPI_thread_init((unsigned long (*)(void))(omp_get_thread_num)) != PAPI_OK){
      printf("PAPI_thread_init error!\n");
      exit(1);
    }

    thread_count = NUM_THREADS;
  #else
    thread_count = 1;
  #endif
  // long long realtime[2];

  // realtime[0] = 0;
  // realtime[1] = 0;

  for(i=0; i<NUM_EVENTS; i++)
      global_CM[i] = 0;
  // retval = PAPI_create_eventset(EventSet);
  // if (PAPI_create_eventset(&EventSet) != PAPI_OK){
  //   printf("PAPI library create error!\n");
  //   exit(1);
  // }

#endif

      // the main loop over time steps
    // int iT, count;
    count =0;
#ifndef TWOSTEPS
    for (iT=0; iT<maxT; iT+=1) {
#else
    for (iT=0; iT<maxT; iT+=2) {
#endif
        if (iT%interval==0) {
            printf("t=%d\n", iT);
            fflush(stdout);
        }

        //save before computing
/*#ifndef NO_SAVE*/
        /*if (iT%tSave==0) {*/
            /*printf("iT=%d, save, count=%d\n", iT, count);*/
            /*fflush(stdout);*/
            /*sprintf(filename, "vel_%s_%d.dat", case_name, count);*/
            /*saveVel(&sim, filename);*/
            /*count++;*/
        /*}*/
/*#endif*/

#ifdef ZGB
        t0 = get_cur_time();
        // on the right boundary, outlet condition grad_x u = 0
        updateZeroGradientBoundary();
        t1 = get_cur_time();
        t_zgb += t1-t0;
#endif

        //save after ZGB
/*#ifndef NO_SAVE*/
        /*if (iT%tSave==0) {*/
            /*printf("iT=%d, save, count=%d\n", iT, count);*/
            /*fflush(stdout);*/
            /*sprintf(filename, "vel_%s_zgb_%d.dat", case_name, iT);*/
            /*saveVel(&sim, filename);*/
        /*}*/
/*#endif*/

        // step 3,4,5:compute rho, u, get f^eq, then update f.
        t0 = get_cur_time();
        collision_func(&sim);
        t1 = get_cur_time();
        t_collid += t1-t0;


//           // the data are written to disk after collision, to be
//           //   that the macroscopic variables are computed
//           //   correctly on the boundaries

/* #ifndef NO_SAVE*/
         /*if (iT%tSave==0) {*/
             /*sprintf(filename, "vel_%s_collide_%d.dat", case_name, iT);*/
             /*saveVel(&sim, filename);*/
         /*}*/

 /*#endif*/

        // step 2: streamming step
#if defined(TIGHT)
        //do nothing
#elif defined(TIGHT_BLOCK)
        //do nothing
#elif defined(PANEL)
        //do nothing
#else
        t0 = get_cur_time();
        stream_func(&sim);
        t1 = get_cur_time();
        t_stream += t1-t0;
#endif

        //save after stream
#ifndef NO_SAVE
        if (iT%tSave==0) {
            printf("iT=%d, save, count=%d\n", iT, count);
            fflush(stdout);
            sprintf(filename, "vel_%s_stream_%d.dat", case_name, count);
            saveVel(&sim, filename);
            count++;
        }
#endif

          // By default: periodic boundary conditions. In this case,
          //   this is important, because the upper and lower
          //   boundaries are horizontally periodic, so that no
          //   special corner nodes are needed.
// #if defined(TIGHT)
//         //do nothing
// #else
        // makePeriodic(&sim);
// #endif

          // the data are written to disk after streaming and Periodic, to be
          //   that the macroscopic variables are computed
          //   correctly on the boundaries

/*
        updateZeroGradientBoundary();
        collision_func(&sim);
        stream_func(&sim);
        makePeriodic(&sim);
        */
    }

#ifdef ADDPAPI
  /* Collect the data into the variables passed in */
    for(i=0; i<NUM_EVENTS; i++){
      //Event[%d]
      printf("Global_E[%02d]: %lld\n", i, global_CM[i]);
    }
    fflush(stdout);

    printf("AE: ");
    for(i=0; i<NUM_EVENTS; i++){
    //Avg_Event[%d]
    printf("%lld ", global_CM[i]/thread_count);
    }
    printf("\n");
    fflush(stdout);

#endif

    //save final data
#ifndef NO_SAVE
    sprintf(filename, "vel_%s_stream_%d.dat", case_name, count);
    saveVel(&sim, filename);
#endif

    t_end = get_cur_time();
    printf("time eclapsed in %d steps, t_total=%.3f, t_collid=%.3f, t_stream=%.3f, t_zgb=%.3f\n",
      maxT, t_end-t_start, t_collid, t_stream, t_zgb);



    destructSim(&sim);
    freeData();
}
