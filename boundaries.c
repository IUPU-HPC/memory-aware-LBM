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

#include "boundaries.h"
#include "lb.h"

/********************* added by Yuankun Fu ***********************/


/* D2Q9 lattice constants                                        */
/*****************************************************************/

  // opposite directions, for bounce back implementation
static const int oppositeOf[9] = { 0, 3, 4, 1, 2, 7, 8, 5, 6 };

  // lattice weights
static const double t[9] = { 4./9., 1./9., 1./9., 1./9., 1./9.,
                             1./36., 1./36., 1./36., 1./36. };
  // lattice velocities
static const int c[9][2] = {
    {0,0},
    {1,0}, {0,1}, {-1,0}, {0,-1},
    {1,1}, {-1,1}, {-1,-1}, {1,-1}
};


/* Bounce back                                                   */
/*****************************************************************/

void bounceBack(double* fPop, void* selfData) {
    static double fTmp[9];
    int iPop;

    for (iPop=0; iPop<9; ++iPop) {
        fTmp[iPop] = fPop[oppositeOf[iPop]];
    }

    for (iPop=0; iPop<9; ++iPop) {
        fPop[iPop] = fTmp[iPop];
    }
}

/* Helper functions: compute rho from u and vice versa           */
/*****************************************************************/

  /* Compute density on wall from bulk information on
     upper boundary. */
inline static double upperRho(double* fPop, double uy) {
    return 1./(1.+uy) * (
        fPop[0] + fPop[3] + fPop[1] + 2*(fPop[6]+fPop[2]+fPop[5])
    );
}

  /* Compute uy on wall from bulk information on
     upper boundary. */
inline static double upperU(double* fPop, double rho) {
    return -1. + 1./rho * (
        fPop[0] + fPop[3] + fPop[1] + 2*(fPop[6]+fPop[2]+fPop[5])
    );
}

  /* Compute density on wall from bulk information on
     lower boundary. */
inline static double lowerRho(double* fPop, double uy) {
    return 1./(1.-uy) * (
        fPop[0] + fPop[3] + fPop[1] + 2*(fPop[8]+fPop[4]+fPop[7])
    );
}

  /* Compute uy on wall from bulk information on
     lower boundary. */
inline static double lowerU(double* fPop, double rho) {
    return 1. - 1./rho * (
        fPop[0] + fPop[3] + fPop[1] + 2*(fPop[8]+fPop[4]+fPop[7])
    );
}

  /* Compute density on wall from bulk information on
     right boundary. */
inline static double rightRho(double* fPop, double ux) {
    return 1./(1.+ux) * (
        fPop[0] + fPop[2] + fPop[4] + 2*(fPop[1]+fPop[5]+fPop[8])
    );
}

  /* Compute ux on wall from bulk information on
     right boundary. */
inline static double rightU(double* fPop, double rho) {
    return -1. + 1./rho * (
        fPop[0] + fPop[2] + fPop[4] + 2*(fPop[1]+fPop[5]+fPop[8])
    );
}

  /* Compute density on wall from bulk information on
     left boundary. */
inline static double leftRho(double* fPop, double ux) {
    return 1./(1.-ux) * (
        fPop[0] + fPop[2] + fPop[4] + 2*(fPop[3]+fPop[6]+fPop[7])
    );
}

  /* Compute ux on wall from bulk information on
     left boundary. */
inline static double leftU(double* fPop, double rho) {
    return 1. - 1./rho * (
        fPop[0] + fPop[2] + fPop[4] + 2*(fPop[3]+fPop[6]+fPop[7])
    );
}

/* Zou/He helper functions and boundary implmenetations          */
/*****************************************************************/

inline static void completeUpper(double* fPop,
                                 double ux, double uy, double rho)
{
    fPop[7] = fPop[5] + 0.5 * (fPop[1]-fPop[3])
                      - rho*uy/6 - rho*ux/2.;
    fPop[8] = fPop[6] + 0.5 *(fPop[3]-fPop[1])
                      - rho*uy/6 + rho*ux/2.;
    fPop[4] = fPop[2] - 2./3.*rho*uy;
}

inline static void completeLower(double* fPop,
                                 double ux, double uy, double rho)
{
    fPop[6] = fPop[8] + 0.5 *(fPop[1]-fPop[3])
                      + rho*uy/6 - rho*ux/2.;
    fPop[5] = fPop[7] + 0.5 *(fPop[3]-fPop[1])
                      + rho*uy/6 + rho*ux/2.;
    fPop[2] = fPop[4] + 2./3.*rho*uy;
}

inline static void completeRight(double* fPop,
                                 double ux, double uy, double rho)
{
    fPop[6] = fPop[8] + 0.5 *(fPop[4]-fPop[2])
                      - rho*ux/6 + rho*uy/2.;
    fPop[7] = fPop[5] + 0.5 *(fPop[2]-fPop[4])
                      - rho*ux/6 - rho*uy/2.;
    fPop[3] = fPop[1] - 2./3.*rho*ux;
}

inline static void completeLeft(double* fPop,
                                double ux, double uy, double rho)
{
    fPop[5] = fPop[7] + 0.5 *(fPop[4]-fPop[2])
                      + rho*ux/6. + rho*uy/2.;
    fPop[8] = fPop[6] + 0.5 *(fPop[2]-fPop[4])
                      + rho*ux/6. - rho*uy/2.;
    fPop[1] = fPop[3] + 2./3.*rho*ux;
}

void upperZouHe(double* fPop, void* selfData) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    double rho = upperRho(fPop, data->uy);
    completeUpper(fPop, data->ux, data->uy, rho);
    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void lowerZouHe(double* fPop, void* selfData) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    double rho = lowerRho(fPop, data->uy);
    completeLower(fPop, data->ux, data->uy, rho);
    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void leftZouHe(double* fPop, void* selfData) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    double rho = leftRho(fPop, data->ux);
    completeLeft(fPop, data->ux, data->uy, rho);
    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void rightZouHe(double* fPop, void* selfData) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    double rho = rightRho(fPop, data->ux);
    completeRight(fPop, data->ux, data->uy, rho);
    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void upperPressureZouHe(double* fPop, void* selfData) {
    PressureBCData* data = (PressureBCData*) selfData;
    double uy = upperU(fPop, data->rho);
    completeUpper(fPop, data->uPar, uy, data->rho);
    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void lowerPressureZouHe(double* fPop, void* selfData) {
    PressureBCData* data = (PressureBCData*) selfData;
    double uy = lowerU(fPop, data->rho);
    completeLower(fPop, data->uPar, uy, data->rho);
    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void leftPressureZouHe(double* fPop, void* selfData) {
    PressureBCData* data = (PressureBCData*) selfData;
    double ux = leftU(fPop, data->rho);
    completeLeft(fPop, ux, data->uPar, data->rho);
    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void rightPressureZouHe(double* fPop, void* selfData) {
    PressureBCData* data = (PressureBCData*) selfData;
    double ux = rightU(fPop, data->rho);
    completeRight(fPop, ux, data->uPar, data->rho);
    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}


/* Regularized helper functions and boundary implmenetations     */
/*****************************************************************/

inline static void splitEqNeq (
        double* f, double* fEq, double* fNeq,
        double rho, double ux, double uy)
{
    int iPop;
    for (iPop=0; iPop<9; ++iPop) {
        fEq[iPop]
            = computeEquilibrium(iPop, rho, ux, uy, ux*ux+uy*uy);
        fNeq[iPop] = f[iPop]-fEq[iPop];
    }
}

inline static void regularizedF (
        double* f, double* fEq,
        double neqPixx, double neqPiyy, double neqPixy )
{
    int iPop;
    for (iPop=0; iPop<9; ++iPop) {
        f[iPop] = fEq[iPop] + 9./2. * t[iPop] *
            ( (c[iPop][0]*c[iPop][0]-1./3.)*neqPixx +
              (c[iPop][1]*c[iPop][1]-1./3.)*neqPiyy +
              2.*c[iPop][0]*c[iPop][1]*neqPixy );
    }
}

  // compute non-equilibrium stress tensor on
  // upper boundary
inline static void upperNeqPi(double* fNeq, double* neqPixx,
                              double* neqPiyy, double* neqPixy)
{
    *neqPixx = fNeq[1] + fNeq[3] + 2.*(fNeq[5]+fNeq[6]);
    *neqPiyy = 2. * (fNeq[6]+fNeq[2]+fNeq[5]);
    *neqPixy = 2. * (fNeq[5] - fNeq[6]);
}

  // compute non-equilibrium stress tensor on
  // lower boundary
inline static void lowerNeqPi(double* fNeq, double* neqPixx,
                              double* neqPiyy, double* neqPixy)
{
    *neqPixx = fNeq[1] + fNeq[3] + 2.*(fNeq[7]+fNeq[8]);
    *neqPiyy = 2. * (fNeq[4]+fNeq[7]+fNeq[8]);
    *neqPixy = 2. * (fNeq[7] - fNeq[8]);
}

  // compute non-equilibrium stress tensor on
  // right boundary
inline static void rightNeqPi(double* fNeq, double* neqPixx,
                              double* neqPiyy, double* neqPixy)
{
    *neqPixx = 2. * (fNeq[1]+fNeq[5]+fNeq[8]);
    *neqPiyy = fNeq[2] + fNeq[4] + 2.*(fNeq[5]+fNeq[8]);
    *neqPixy = 2. * (fNeq[5] - fNeq[8]);
}

  // compute non-equilibrium stress tensors on
  // left boundary
inline static void leftNeqPi(double* fNeq, double* neqPixx,
                             double* neqPiyy, double* neqPixy)
{
    *neqPixx = 2. * (fNeq[3]+fNeq[6]+fNeq[7]);
    *neqPiyy = fNeq[2] + fNeq[4] + 2.*(fNeq[6]+fNeq[7]);
    *neqPixy = 2. * (fNeq[7] - fNeq[6]);
}

void upperRegularized(double* fPop, void* selfData) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    double rho = upperRho(fPop, data->uy);
    double neqPixx, neqPiyy, neqPixy;
    double fEq[9], fNeq[9];

    splitEqNeq(fPop, fEq, fNeq, rho, data->ux, data->uy);
    upperNeqPi(fNeq, &neqPixx, &neqPiyy, &neqPixy);
    regularizedF(fPop, fEq, neqPixx, neqPiyy, neqPixy);

    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void lowerRegularized(double* fPop, void* selfData) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    double rho = lowerRho(fPop, data->uy);
    double neqPixx, neqPiyy, neqPixy;
    double fEq[9], fNeq[9];

    splitEqNeq(fPop, fEq, fNeq, rho, data->ux, data->uy);
    lowerNeqPi(fNeq, &neqPixx, &neqPiyy, &neqPixy);
    regularizedF(fPop, fEq, neqPixx, neqPiyy, neqPixy);

    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void leftRegularized(double* fPop, void* selfData) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    double rho = leftRho(fPop, data->ux);
    double neqPixx, neqPiyy, neqPixy;
    double fEq[9], fNeq[9];

    splitEqNeq(fPop, fEq, fNeq, rho, data->ux, data->uy);
    leftNeqPi(fNeq, &neqPixx, &neqPiyy, &neqPixy);
    regularizedF(fPop, fEq, neqPixx, neqPiyy, neqPixy);

    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void rightRegularized(double* fPop, void* selfData) {
    VelocityBCData* data = (VelocityBCData*) selfData;
    double rho = rightRho(fPop, data->ux);
    double neqPixx, neqPiyy, neqPixy;
    double fEq[9], fNeq[9];

    splitEqNeq(fPop, fEq, fNeq, rho, data->ux, data->uy);
    rightNeqPi(fNeq, &neqPixx, &neqPiyy, &neqPixy);
    regularizedF(fPop, fEq, neqPixx, neqPiyy, neqPixy);

    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void upperPressureRegularized(double* fPop, void* selfData) {
    PressureBCData* data = (PressureBCData*) selfData;
    double uy  = upperU(fPop, data->rho);
    double neqPixx, neqPiyy, neqPixy;
    double fEq[9], fNeq[9];

    splitEqNeq(fPop, fEq, fNeq, data->rho, data->uPar, uy);
    upperNeqPi(fNeq, &neqPixx, &neqPiyy, &neqPixy);
    regularizedF(fPop, fEq, neqPixx, neqPiyy, neqPixy);

    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void lowerPressureRegularized(double* fPop, void* selfData) {
    PressureBCData* data = (PressureBCData*) selfData;
    double uy  = lowerU(fPop, data->rho);
    double neqPixx, neqPiyy, neqPixy;
    double fEq[9], fNeq[9];

    splitEqNeq(fPop, fEq, fNeq, data->rho, data->uPar, uy);
    lowerNeqPi(fNeq, &neqPixx, &neqPiyy, &neqPixy);
    regularizedF(fPop, fEq, neqPixx, neqPiyy, neqPixy);

    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void leftPressureRegularized(double* fPop, void* selfData) {
    PressureBCData* data = (PressureBCData*) selfData;
    double ux = leftU(fPop, data->rho);
    double neqPixx, neqPiyy, neqPixy;
    double fEq[9], fNeq[9];

    splitEqNeq(fPop, fEq, fNeq, data->rho, ux, data->uPar);
    leftNeqPi(fNeq, &neqPixx, &neqPiyy, &neqPixy);
    regularizedF(fPop, fEq, neqPixx, neqPiyy, neqPixy);

    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

void rightPressureRegularized(double* fPop, void* selfData) {
    PressureBCData* data = (PressureBCData*) selfData;
    double ux = rightU(fPop, data->rho);
    double neqPixx, neqPiyy, neqPixy;
    double fEq[9], fNeq[9];

    splitEqNeq(fPop, fEq, fNeq, data->rho, ux, data->uPar);
    rightNeqPi(fNeq, &neqPixx, &neqPiyy, &neqPixy);
    regularizedF(fPop, fEq, neqPixx, neqPiyy, neqPixy);

    data->bulkDynamics->dynamicsFun (
            fPop, data->bulkDynamics->selfData );
}

