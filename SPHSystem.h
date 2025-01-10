/* 
 * This program is part of the following research paper:
 * 
 * [1] M. N. Rahimi and G. Moutsanidis, "Modeling concrete failure with smoothed particle hydrodynamics using the Microplane (M7) constitutive model" (Journal)
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
#ifndef __SPHSYSTEM_H__
#define __SPHSYSTEM_H__

#include <sys/stat.h>  // For mkdir on Unix-like systems
#include <sys/types.h>
#include <string.h>
#include <fstream>
#include <omp.h>
#include <math.h>
#include <ctime>
#include <iostream>
#include "Vector3D.h"
#include "particle.h"

#ifndef PI
#define PI 3.14159265358979323846
#endif
#ifdef _WIN32
	#include <direct.h>    // For _mkdir on Windows
    #define mkdir(x) _mkdir(x)
#endif

using namespace std;

class SPHSystem3D
{
public:
	SPHSystem3D(int argc, char *argv[]);
	~SPHSystem3D();

	// FUNCTIONS
	//************************************************************************
	// Generate initial particle distrubution
	void initParticle();
	// Add one particle at position pos with velocity vel and type tpe
	void addSingleParticle(Vec3d pos, Vec3d vel, int tpe);
	// Calculate cell position of the particle located at pos
	Vec3i calcCellPos(Vec3d pos);
	// Calculate hash of the cell located with position pos
	uint calcCellHash(Vec3i pos);
	// Build search Grids
	void buildGrid();
	// Calculate kernels with grid search
	void calcKernel();
	// Correct kernel and its gradient
	void gradientCorrection();

	// One step of time integration
	void OneStepCalculations();

	// Write results in paraview format.
	void SaveResults();

	// getter
	uint getNumParticle() { return numParticle; }
	Particle3D *getParticles() { return particles; }
	Cell3D *getCells() { return cells; }

	// kernel function and its derivative
	double CubicSplineKernel(double qq)
	{
		if (qq >= 2.0e0) return 0.0e0;
		else if (qq >= 1.0e0) return (0.25e0 / (PI * kernel * kernel * kernel)) * pow(2.0e0 - qq, 3);
		else if (qq >= 0.0e0 - INF) return (1.0e0 / (PI * kernel * kernel * kernel)) * (1.0e0 - 1.5 * qq * qq + 0.75 * qq * qq * qq);
		else
		{
			printf("Error in calculating the Cubic Spline kernel, q is less than 0.0");
			exit(0);
		}
	}
	Vec3d CubicSplineGradient(double qq, double dist, Vec3d distVec)
	{
		if (qq >= 2.0e0) return Vec3d(0.0, 0.0, 0.0);
		else if (qq >= 1.0e0) return (0.25e0 / (PI * kernel * kernel * kernel)) * (-3.0e0 * pow(2.0e0 - qq, 2)) / (dist * kernel) * distVec;
		else if (qq >= 0.0e0 - INF) return (1.0e0 / (PI * kernel * kernel * kernel)) * (-3.0e0 * qq + 9.0e0 / 4.0e0 * qq * qq) / (dist * kernel) * distVec;
		else
		{
			printf("Error in calculating the the gradient of Cubic kernel, q is less than 0.0:");
			exit(0);
		}
	}
	
	double maxTime;		  	// Maximum time of simulation [s]
	double timeStep;		// Time step used for simulation

private:
	int thrd_num; 			// Number of parallel threads
	double k_kernel; 		// Smoothing lenght coefficient for kernel
	int numParticle; 		// Total Number of particles
	

	int iTimestep;			// Current analysis step
	int maxTimeStep;		// Maximum Timestep for simulation
	double currenttime; 	// Current time

	Mat3d Ident;	  	 	// Identity matrix

	uint Type;				// Size effect type -> 1: Type 1, 2: Type 2
	double D;
	// Geometry Properties
	double Lx;				// Lenght of plate in X direction
	double Ly;				// Lenght of plate in Y direction
	double Lz;				// Lenght of plate in Z direction

	int print_step; 		// Step interval for printing informations
	int save_step;			// Step interval for saving data

	double CompressionRate; 		// Rate of compression
	double maxCompression; 	// Total percent Compression
	int RampTimeStep;				// Total time step to increase Compresssion rate from 0 to CompressionRate
	double compression;				//Amount of compression

	double Wstrain;	  				// Total strain energy
	double Wfracture; 				// Total fracture energy
	double twself;					// Kernel value at particles own position twself = CubicSplineKernel(0)
	// Variables for domain limits
	double maxX, minX, maxY, minY, maxZ, minZ;
	// Material Properties
	double young;					//Elastic Modulus
	double dens;					//Density, for M7 model and static analyses we take dens=1 to minimize dynamic effects.
	double poisson;					//Poisson's ratio

	double ArtViscCoeff;			//Constant for Artificial Viscousity
	double kernel;	 				//Smoothing Lenght, h. (kernel=1.33 dx)
	double partDist; 				//Initial Particle distance
	int numCrossecParticles;		//Number of particles at the cross section where stress is calculated
	int *ListCrossecParticles;		//List of particles at the cross section where stress is calculated

	Particle3D *particles;			//List of particles

	//Variables for faster neighbour search
	Vec3i gridSize;					//Total number of cells in x, y, and z directions.
	double cellSize;				//Size of one cell
	uint totCell;					//Total number of celss
	Cell3D *cells;					//List of cells

	// MicroPlane paramters
	double k_1;
	double k_2;
	double k_3;
	double k_4;
	double k_5;
	double MPMqn[6][37];
	double MPMql[6][37];
	double MPMqm[6][37];
	double MPMw[37];

	//Particle-to-Particle contact
	void InitContParticles();							//Generates contact particles
	Vec3d *ContPos;										//Current Position list of particles of the contacting object(S)
	Vec3d *ContPosi;									//Initial Position list of particles of the contacting object(S)
	Vec3d *ContForce;									//List of Contact force of particles of the contacting object(S)
	int *BodyPartList;									//List of all the SPH particles that potentially will be in contact with contacting object(s)
	
	/*Here I assume 3 contacting objects and gnerate lists for 3+1*/
	int BodyPartStart[4]={0,0,0,0};		//Adresses of the first SPH particle contacting with nth=1,2,3 contacting object(s)
	int ContPartStart[4]={0,0,0,0};		//Adresses of the first contact particle of the nth=1,2,3 contacting object(s)
	int numContParticles=0;				//Total number of contact particles in all the contacting objects
	int numListParticles=0;				//Total number of SPH particles that potentially will be in contact with contacting object(s)
	double width;						//Width of contacting objects
	double height;						//height of contacting objects
	double dr0=0.0;						//Extra contact range, total contact range is partDist+dr0
	double Kp=0.0;						//Contact potential
};

#endif