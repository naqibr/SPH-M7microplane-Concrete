/* 
 * This program is part of the following research paper:
 * 
 * [1] M. N. Rahimi and G. Moutsanidis, "SPH modeling of concrete failure using the M7 microplane model" International Journal of Mechanical Sciences
 * 
 * It also implements the Total Lagrangian Smoothed Particle Hydrodynamics formulations presented in:
 * 
 * [2] M. N. Rahimi and G. Moutsanidis, “A smoothed particle hydrodynamics approach for phase field modeling of brittle fracture,” 
 *     Computer Methods in Applied Mechanics and Engineering, Aug. 2022. doi: https://doi.org/10.1016/j.cma.2022.115191
 * 
 * [3] M. N. Rahimi and G. Moutsanidis, “Modeling dynamic brittle fracture in functionally graded materials using hyperbolic phase field and smoothed particle hydrodynamics,” 
 *     Computer Methods in Applied Mechanics and Engineering, Nov. 2022. doi: https://doi.org/10.1016/j.cma.2022.115642
 * 
 * [4] M. N. Rahimi and G. Moutsanidis, “An SPH-based FSI framework for phase-field modeling of brittle fracture under extreme hydrodynamic events,” 
 *     Engineering with Computers, Aug. 2023. doi: https://doi.org/10.1007/s00366-023-01857-0
 * 
 * [5] M. N. Rahimi and G. Moutsanidis, “IGA-SPH: Coupling isogeometric analysis with smoothed particle hydrodynamics for air-blast–structure interaction,” 
 *     Engineering with Computers, May 2024. doi: https://doi.org/10.1007/s00366-024-01978-0
 * 
 * If this code is helpfull, whether partially or entirely, in your research, please cite the relevant paper(s).
 * 
 * ============================================================================
 */
#pragma once
#ifndef __Particle_H__
#define __Particle_H__

#include "Vector3D.h"
#include "Matrix3D.h"
#define MAXNEIGHBOUR2D 80
#define MAXNEIGHBOUR3D 160

class Particle3D
{
public:
	int type;						//Type of particle = 0: normal particle to be solved, =!0 particles with special conditions
	double vol;						//Particle volume
	double mass;					//Particle mass

	Vec3d disp;						//Particle displacement u
	Vec3d pos;						//Particle current position x
	Vec3d posi;						//Particle initial position X
	Vec3d vel;						//Particle current Velocity
	Vec3d accl;						//Particle current Acceleration = dv/dt
	Vec3d artviscforce;				//Artificial viscousity
	Vec3d iforce;					//Internal force
	Vec3d contforce;				//Contact force

	Mat3d FF;						//Deformation Gradient FF=dx/dX
	double jacob;					//Determinant of defromation gradient
	Mat3d FFdot;					//Rate of defromation gradient
	Mat3d FFinv;					//Inverse of Deformation Gradient F^-1
	Mat3d LL;						//Displacement Gradient LL=du/dX
	Mat3d EE;						//Strain Tensor
	Mat3d EEdot;					//Strain Rate Tensor
	double maxPstrain;				//Maximum principle strain
	double minPstrain;				//Minimum principle strain
	double wSum;					//Sum of kernel
	int nNum;						//Number of neighbor particles for particle i
	int id;							//ID of particle i
	int pid[MAXNEIGHBOUR3D];		//ID of the neighbour particle j
	double pQ[MAXNEIGHBOUR3D];		//kernel for neighbour particle (xi-xj)
	Vec3d pG[MAXNEIGHBOUR3D];		//kernel gradient for neighbour particle (xi-xj)

	Vec3d distVec[MAXNEIGHBOUR3D];	//Current distance vector xi-xj
	double dist[MAXNEIGHBOUR3D];	//Lenght of current distance |xi-xj|
	Vec3d idistVec[MAXNEIGHBOUR3D];	//Initital distance vector Xi-Xj
	double idist[MAXNEIGHBOUR3D];	//Lenght of initial distance |Xi-Xj|
	Vec3d velDif[MAXNEIGHBOUR3D];	//Velocity difference vi-vj
	Vec3d DispDif[MAXNEIGHBOUR3D];	//Displacement difference ui-uj

	//Variables used for MicroPlane model
	Mat3d epsdot;
	Mat3d strainInc;
	Mat3d eps;
	Mat3d epsOld;
	Mat3d sigma;
	Mat3d sigmaold;
	double deps[6];
	double epsN[6];
  	double sigmaN[6];
  	double sigmaNP1[6];
  	double stateN[189];
  	double stateNP1[189];
  	double microDamageN[37];
  	double microDamageNP1[37];
	double internEnergyN;
    double inelasEnergyN;
    double internEnergyNP1;
    double inelasEnergyNP1;
    double weightedDamageNP1;
	double MaxPstrain;
	double MinPstrain;
	
	Particle3D* next;
};


class Cell3D
{
public:
	Particle3D* head;
};

#endif
