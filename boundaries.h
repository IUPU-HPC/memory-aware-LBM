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

/* boundaries.h: Assembly of various collision terms for the
 *               lattice boundaries. They can be used as
 *               members of the struct "Dynamics" and replace
 *               the bgk collision on boundaries. Only straight
 *               boundaries, and no corners are implemented
 *               currently. The implementation includes Zou/He
 *               boundaries, and regularized boundaries computed
 *               from the macroscopic quantities rho, u and pi
 *               (=stress tensor). Regularized boundaries tend
 *               to be more stable than Zou/He. Both Dirichlet
 *               velocity boundaries and pressure boundaries are
 *               implemented.
 *               Feel free to add your own boundary conditions.
 */

#ifndef BOUNDARIES_H
#define BOUNDARIES_H

#include "lb.h"

/* struct VelocityBCData                                         */
/*****************************************************************/
/*   This struct is used as a parameter to the velocity boundary
 *   conditions. It takes place of the variable "selfData" in the
 *   struct "Dynamics". It contains the value of the Dirichlet
 *   boundary, and a pointer to the underlying bulk dynamics, e.g.
 *   bgk.
 */
typedef struct {
    Dynamics* bulkDynamics;
    double    ux, uy;
} VelocityBCData;

/* struct PressureBCData                                         */
/*****************************************************************/
/*   This struct is used as a parameter to the pressure boundary
 *   conditions. It takes place of the variable "selfData" in the
 *   struct "Dynamics". It contains the value of the velocity
 *   tangential to the boundary, the value of the prescribed
 *   pressure, and a pointer to the underlying bulk dynamics, e.g.
 *   bgk.
 */
typedef struct{
    Dynamics* bulkDynamics;
    double    rho, uPar;
} PressureBCData;

/* All the implemented boundaries...                             */
/*****************************************************************/

  // bounce back no-slip condition. "selfData" is a null pointer,
  //   given that the function need no parameter.
void bounceBack(double* fPop, void* selfData);

  // ZouHe velocity boundaries on upper, lower, left and right
  //   boundaries
void upperZouHe(double* fPop, void* selfData);
void lowerZouHe(double* fPop, void* selfData);
void leftZouHe (double* fPop, void* selfData);
void rightZouHe(double* fPop, void* selfData);

  // ZouHe pressure boundaries on upper, lower, left and right
  //   boundaries
void upperPressureZouHe(double* fPop, void* selfData);
void lowerPressureZouHe(double* fPop, void* selfData);
void leftPressureZouHe (double* fPop, void* selfData);
void rightPressureZouHe(double* fPop, void* selfData);

  // Regularized velocity boundaries on upper, lower, left and
  // right boundaries
void upperRegularized(double* fPop, void* selfData);
void lowerRegularized(double* fPop, void* selfData);
void leftRegularized (double* fPop, void* selfData);
void rightRegularized(double* fPop, void* selfData);

  // Regularized pressure boundaries on upper, lower, left and
  // right boundaries
void upperPressureRegularized(double* fPop, void* selfData);
void lowerPressureRegularized(double* fPop, void* selfData);
void leftPressureRegularized (double* fPop, void* selfData);
void rightPressureRegularized(double* fPop, void* selfData);

#endif
