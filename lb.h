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

/* lb.h: Assembly of various functions for simulating a lattice
 *       Boltzmann dynamics with BGK collision term on a two-
 *       dimenional D2Q9 lattice. The code is generic: every
 *       lattice site can have a different dynamics.
 *       Other collision terms than BGK can easily be added.
 */

#ifndef LB_H
#define LB_H

#include "stddef.h"

#define ZGB

#ifdef _OPENMP
#include <omp.h>
#endif

#define NUM_EVENTS 4 
#ifdef ADDPAPI
#include "papi.h"
static int EventSet[] = {PAPI_L1_TCM, PAPI_L2_TCM, PAPI_L3_TCM,
                          PAPI_TLB_DM, //PAPI_TLB_IM,
                          //PAPI_L2_ICA, PAPI_L3_ICA,
                          // PAPI_L2_DCA, PAPI_L3_DCA,
                          //PAPI_L1_DCM, PAPI_L1_ICM,
                          //PAPI_L2_DCM, PAPI_L2_ICM
                          };
#endif

/* struct Dynamics                                               */
/*****************************************************************/
/*   emulation of a class that defines the dynamics of a lattice
 *   site. The function dynamicsFun contains the algorithm of the
 *   collision, and selfData points to the function arguments
 *   (for example the local viscosity)
 */
typedef struct {
    void   (*dynamicsFun) (double* fPop, void* selfData);
    void* selfData;
} Dynamics;


/* struct Node                                                   */
/*****************************************************************/
/*  a lattice node, containing the data (distribution functions)
 *  for a D2Q9 lattice, and a pointer to the local dynamics, i.e.
 *  the collision term. Two "methods" are added to construct and
 *  initialize a node.
 */
typedef struct {
    double    fPop[9];
    Dynamics* dynamics;
} Node;

void constructNode(Node* node);
void iniEquilibrium(Node* node, double rho, double ux, double uy);


/* struct Simulation                                             */
/*****************************************************************/
/*  a full D2Q9 LB simulation, containing the allocated memory
 *  for the lattice. Some "methods" are added to initiate the
 *  dynamics.
 */
typedef struct {
    int lx, ly;               // lx*ly lattice
    Node*  memoryChunk;       // contiguous raw memory for buf1
    Node*  tmpMemoryChunk;    // contiguous raw tmp memory for buf2
    Node** lattice;           // lattice, points to raw memory
    Node** tmpLattice;        // tmp lattice, points to raw memory
} Simulation;

/* different algorithm for collide&streaming                     */
/*****************************************************************/
void updateZeroGradientBoundary();

/* Original LBM*/
void collide(Simulation* sim);
void collide_openmp(Simulation* sim);
void propagate(Simulation* sim);
void propagate_openmp(Simulation* sim);
void finalize_stream(Simulation* sim);

/* Fused LBM*/
void collide_with_stream(Simulation* sim);
void collide_with_stream_openmp(Simulation* sim);

/* Two-steps Line LBM*/
void collide_tight(Simulation* sim);
void collide_tight2(Simulation* sim);
void collide_tight_openmp(Simulation* sim);

/* Two-steps Blocking LBM*/
void collide_tight_block(Simulation* sim);
void collide_tight_block_openmp(Simulation* sim);

/* Two-steps Panel LBM*/
void collide_tight_panel_ix(Simulation* sim);
void collide_tight_panel_iy(Simulation* sim);
void collide_tight_panel_iy_openmp(Simulation* sim);

/* Initialize sturcture sim*/
void constructSim(Simulation* sim, int lx, int ly);
void destructSim(Simulation* sim);
void setDynamics(Simulation* sim, int iX, int iY, Dynamics* dyn);

void makePeriodic(Simulation* sim);
void saveVel(Simulation* sim, char fName[]);
void saveF(Simulation* sim, int iPop, char fName[]);

/* some free helper functions                                    */
/*****************************************************************/

  // compute density and velocity from the f's
void computeMacros(double* f, double* rho, double* ux, double* uy);

  // compute local equilibrium from rho and u
double computeEquilibrium(int iPop, double rho,
                          double ux, double uy, double uSqr);
  // bgk collision term
void bgk(double* fPop, void* selfData);

#endif
