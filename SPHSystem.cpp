/* 
 * This program is part of the following research paper:
 * 
 * [1] M. N. Rahimi and G. Moutsanidis, "SPH modeling of concrete failure using the M7 microplane model" International Journal of Mechanical Sciences, doi: https://doi.org/10.1016/j.ijmecsci.2025.110378
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
#include "SPHSystem.h"

char DataFileName[100];
std::ofstream DataFile;

// FORTRAN Subroutine to set the micro plane M7 parameters
extern "C"
{
	void setsystem_(double MPMqn[6][37], double MPMql[6][37], double MPMqm[6][37], double MPMw[37]);
}

// Fortran Subroutine for micro plane M7 constitutive model to calculate new stress values
extern "C"
{
	void m7fmaterial_(const double *dt, double *strainInc, double *epsOld,
					  double *stressOld, double *stateOld, const double *enerInternOld,
					  const double *enerInelasOld, double *damageOld, double *stressNew,
					  double *stateNew, double *enerInternNew, double *enerInelasNew,
					  double *damageNew, double *weightedDamage, double *MaxPstrain, double *MinPstrain,
					  const double *young, const double *poisson, const double *k_1,
					  const double *k_2, const double *k_3, const double *k_4, const double *k_5,
					  double MPMqn[6][37], double MPMql[6][37], double MPMqm[6][37], double MPMw[37]);
}

bool createdir(const std::string& path) {
    #ifdef _WIN32
        return _mkdir(path.c_str()) == 0;
    #else
        return mkdir(path.c_str(), 0777) == 0;
    #endif
}

// Creates the SPHSystem3D object, initializes the variables, and runs the pre-time-integration routines.
SPHSystem3D::SPHSystem3D(int argc, char *argv[])
{
	thrd_num = 1;
	omp_set_dynamic(0);
	if (argc >= 2)
	{
		for (int i = 1; i < argc; i = i + 2)
		{
			if (strcmp(argv[i], "-np") == 0)
			{
				thrd_num = atoi(argv[i + 1]);
			}
		}
	}
	if (thrd_num <= 0) thrd_num = 1;
	omp_set_num_threads(thrd_num);
	std::cout << "\n\nNumber of parallel threads used: " << thrd_num << "\n";

	young = 30.1730e9 / 1.0e6; // m7fmaterial_ subroutine requires Elastic modulus in units of MPa
	dens = 1.0e0;			   // Set to 1 to reduce dynamic effects
	poisson = 0.18e0;
	k_1 = 105.0e-6;
	k_2 = 110.0e0;
	k_3 = 20.0e0;
	k_4 = 30.0e0;
	k_5 = 8.0e-4;

	ArtViscCoeff = 5.0;
	k_kernel = 2.0001;
	timeStep = 1.0e-10;

	// Geometry Properties
	Lx = 100.0e-3;
	Ly = 100.0e-3;
	Lz = 50.0e-3;
	partDist = 2.5e-3;
	CompressionRate = (1.0 / 100.0 * Lz) / (1000000.0 * timeStep); // Adjusted so that the prism is compressed 1\% in 1000000 steps
	
	// time optionse
	print_step = 10000;
	save_step = print_step;
	RampTimeStep = save_step;
	maxTime = (1.0 / 100.0 * Lz) / (CompressionRate) + 2.0 * save_step * timeStep;
	kernel = partDist * 1.333333e0;
	cellSize = kernel * k_kernel;

	twself = CubicSplineKernel(0.0e0);

	initParticle();

	currenttime = 0.0e0;
	compression = 0.0;
	maxTimeStep = ceilf(1.0 * maxTime / timeStep);
	iTimestep = 0;

	Ident.setZero();
	Ident.e[0][0] = 1.0e0;
	Ident.e[1][1] = 1.0e0;
	Ident.e[2][2] = 1.0e0;

	std::cout << "\n\n================================";
	std::printf("\n=> SPHSystem variables:");
	std::printf("\n dx\t= %.3g", partDist);
	std::printf("\n Total SPH Particles\t= %d", numParticle);
	std::printf("\n Simulation Timestep:\t= %.3g", timeStep);
	std::printf("\n E\t= %.3g", young);
	std::printf("\n Pr\t= %.2g", poisson);
	std::printf("\n dens\t= %.2g", dens);
	std::printf("\n k_1\t= %.3g", k_1);
	std::printf("\n k_2\t= %.3g", k_2);
	std::printf("\n k_3\t= %.3g", k_3);
	std::printf("\n k_4\t= %.3g", k_4);
	std::printf("\n k_5\t= %.3g", k_5);
	std::printf("\n Art. Visc. Coef.\t= %.3g", ArtViscCoeff);
	std::printf("\n Comp.\t= %.3g [m/s]", CompressionRate);
	std::cout << "\n================================\n";
	std::cout.flush();

    createdir("output");
    createdir("output/ParaviewFiles");
	char dir0[100];
	std::sprintf(dir0, "output/ParaviewFiles");
	createdir(dir0);
	// Create file for stress strain curve
	std::sprintf(DataFileName, "output/_Stress_Strain.txt");
	DataFile.open(DataFileName, std::ios::out | std::ios::trunc);
	DataFile.precision(11);
	DataFile.setf(std::ios::fixed);
	DataFile.setf(std::ios::showpoint);
	DataFile.close();
	buildGrid();
	calcKernel();
	gradientCorrection();

	free(cells);
	setsystem_(MPMqn, MPMql, MPMqm, MPMw);
	printf(" >> System initialized.\n");
}

void SPHSystem3D::initParticle()
{
	numParticle = 0;
	numCrossecParticles = 0;

	Vec3d pos, vel(0.0e0, 0.0e0, 0.0e0);
	double i, j, k;

	maxX = -1.0e20;
	minX = 1.0e20;
	maxY = -1.0e20;
	minY = 1.0e20;
	maxZ = -1.0e20;
	minZ = 1.0e20;
	int tpe;

	//Dummy count
	for (k = partDist / 2.0e0; k <= Lz; k += partDist)
	{
		if (k > maxZ)
			maxZ = k;
		if (k < minZ)
			minZ = k;

		for (j = partDist / 2.0e0; j <= Ly; j += partDist)
		{
			if (j > maxY)
				maxY = j;
			if (j < minY)
				minY = j;

			for (i = partDist / 2.0e0; i <= Lx; i += partDist)
			{
				if (i > maxX)
					maxX = i;
				if (i < minX)
					minX = i;
				numParticle++;
			}
		}
	}
	particles = (Particle3D *)malloc(sizeof(Particle3D) * numParticle);
	numParticle = 0;

	// Generate particles for the right half of the domain
	for (k = partDist / 2.0e0; k <= Lz; k += partDist)
	{
		for (j = partDist / 2.0e0; j <= Ly; j += partDist)
		{
			for (i = partDist / 2.0e0; i <= Lx; i += partDist)
			{
				tpe = 0; // no constrain
				if (k >= Lz - partDist)
								tpe = 3; // Apply Compression in Z direction
				if (k < 0.99 * partDist)
								tpe = 1; // Constrain displacement in Z direction
				pos.e[0] = i;
				pos.e[1] = j;
				pos.e[2] = k;
				
				addSingleParticle(pos, vel, tpe);
			}
		}
	}
	gridSize.e[0] = ceil((maxX - minX) / cellSize) + 1;
	gridSize.e[1] = ceil((maxY - minY) / cellSize) + 1;
	gridSize.e[2] = ceil((maxZ - minZ) / cellSize) + 1;

	totCell = gridSize.e[0] * gridSize.e[1] * gridSize.e[2];

	cells = (Cell3D *)malloc(sizeof(Cell3D) * totCell);

	// List of particles on the cross section with normal in Z direction, used for Force calculations
	numCrossecParticles=0;
	Particle3D *p;
	//Dummy count
	for (int k = 0; k < numParticle; k++)
	{
		p = &(particles[k]);
		if ((p->pos.e[2] < Lz -3.0 * partDist) && (p->pos.e[2] >= Lz - 4.0* partDist)) numCrossecParticles++;
	}
	ListCrossecParticles = new int[numCrossecParticles];
	numCrossecParticles = 0;
	for (int k = 0; k < numParticle; k++)
	{
		p = &(particles[k]);
		if ((p->pos.e[2] < Lz -3.0 * partDist) && (p->pos.e[2] >= Lz - 4.0* partDist)) {
			ListCrossecParticles[numCrossecParticles] = k;
			numCrossecParticles++;
		}
	}
	std::cout << "\n >> Domain discretized.";
	std::cout.flush();
}

void SPHSystem3D::buildGrid()
{
	for (uint i = 0; i < totCell; i++)
		cells[i].head = NULL;
	Particle3D *p;
	Vec3i cellPos;
	int hash;

	for (int i = 0; i < numParticle; i++)
	{
		p = &(particles[i]);

		cellPos = calcCellPos(p->pos);
		hash = calcCellHash(cellPos);
		p->id = i;
		if (cells[hash].head == NULL)
		{
			p->next = NULL;
			cells[hash].head = p;
		}
		else
		{
			p->next = cells[hash].head;
			cells[hash].head = p;
		}
	}
	std::cout << "\n >> Grids built.";
	std::cout.flush();
}

void SPHSystem3D::calcKernel()
{
	Particle3D *p;
	Particle3D *np;
	Vec3i cellPos;
	Vec3i nearPos;
	uint hash;

	for (int k = 0; k < numParticle; k++)
	{
		p = &(particles[k]);

		cellPos = calcCellPos(p->pos);

		int nid = 0;
		p->nNum = 0;
		p->wSum = twself;
		for (int s = -1; s <= 1; s++)
		{
			for (int i = -1; i <= 1; i++)
			{
				for (int j = -1; j <= 1; j++)
				{
					nearPos.e[0] = cellPos.e[0] + i;
					nearPos.e[1] = cellPos.e[1] + j;
					nearPos.e[2] = cellPos.e[2] + s;

					hash = calcCellHash(nearPos);

					if (hash == 0xffffffff)
						continue;

					np = cells[hash].head;

					while (np != NULL)
					{
						Vec3d distVec = p->posi - np->posi;
						double disti = distVec.Length();

						if (disti < INF || disti >= kernel * k_kernel)
						{
							np = np->next;
							continue;
						}

						double qq = disti / kernel;

						p->pid[nid] = np->id;

						p->pQ[nid] = CubicSplineKernel(qq);
						p->pG[nid] = CubicSplineGradient(qq, disti, distVec);

						p->nNum = p->nNum + 1;
						p->wSum = p->wSum + p->pQ[nid];

						p->idistVec[nid] = 0.0 - distVec;
						p->idist[nid] = disti;

						p->distVec[nid] = 0.0 - distVec;
						p->dist[nid] = disti;

						nid++;
						np = np->next;
					}
				}
			}
		}
		p->vol = 1.0e0 / p->wSum;
		p->mass = dens * p->vol;
	}

	std::cout << "\n >> Kernels calculated.";
	std::cout.flush();
}

void SPHSystem3D::gradientCorrection()
{
	Particle3D *p;
	Particle3D *np;

	Mat3d Ccor;
	Mat3d TempMat;

	for (int k = 0; k < numParticle; k++)
	{
		p = &(particles[k]);

		Ccor.setZero();
		double wsum0 = 0.0;
		for (int i = 0; i < p->nNum; i++)
		{

			np = &(particles[p->pid[i]]);

			TempMat = (p->idistVec[i].Dyadic(p->pG[i])).transpose() * np->vol;
			Ccor = Ccor + TempMat;
			wsum0 = wsum0 + np->vol * p->pQ[i];
		}
		wsum0 = wsum0 + p->vol * twself;

		if (fabs(Ccor.det()) < INF)
		{
			printf("ERROR! in correcting the kernel derivative, determinant is zero. determinant: %g", Ccor.det());
			exit(0);
		}

		Mat3d Cinv = Ccor.inverse();

		double pwself = twself / wsum0;
		p->wSum = pwself;
		for (int i = 0; i < p->nNum; i++)
		{
			p->pQ[i] = p->pQ[i] / wsum0;
			p->wSum = p->wSum + p->pQ[i];
			p->pG[i] = Cinv.DotVec(p->pG[i]);
		}
	}

	std::cout << "\n >> Gradient corrected.\n";
	std::cout.flush();
}

void SPHSystem3D::OneStepCalculations()
{
	Particle3D *p;
	Particle3D *np;
	Mat3d dFt;

	double dCompRate;
	dCompRate = CompressionRate;
	// Ramp compression from 0 to CompressionRate in RampTimeStep steps
	if (iTimestep < RampTimeStep)
		dCompRate = CompressionRate * iTimestep / RampTimeStep;

	compression = compression + dCompRate * timeStep;
#pragma omp parallel
	{
		//Stress calculations for each particle
#pragma omp for private(p, np) schedule(static)
		for (int k = 0; k < numParticle; k++)
		{
			p = &(particles[k]);
			p->LL.setZero();
			p->FFdot.setZero();
			p->epsOld = p->eps;

			for (int i = 0; i < p->nNum; i++)
			{

				np = &(particles[p->pid[i]]);

				p->distVec[i] = np->pos - p->pos;
				p->dist[i] = p->distVec[i].Length();

				p->velDif[i] = np->vel - p->vel;
				p->DispDif[i] = np->disp - p->disp;

				p->LL = p->LL  + p->DispDif[i].Dyadic(p->pG[i]) * np->vol;
				p->FFdot = p->FFdot + p->velDif[i].Dyadic(p->pG[i]) * np->vol;
			}
			p->FF = p->LL + Ident;
			p->jacob = p->FF.det();
			p->FFinv = p->FF.inverse();

			Mat3d l = p->FFdot.DotMat(p->FFinv);
			p->EEdot = 0.5 * (l + l.transpose());
			p->strainInc = p->EEdot * timeStep;
			// voigt notation
			p->epsN[0] = p->epsOld.e[0][0];
			p->epsN[1] = p->epsOld.e[1][1];
			p->epsN[2] = p->epsOld.e[2][2];
			p->epsN[3] = p->epsOld.e[1][2];
			p->epsN[4] = p->epsOld.e[0][2];
			p->epsN[5] = p->epsOld.e[0][1];

			p->deps[0] = p->strainInc.e[0][0];
			p->deps[1] = p->strainInc.e[1][1];
			p->deps[2] = p->strainInc.e[2][2];
			p->deps[3] = p->strainInc.e[1][2];
			p->deps[4] = p->strainInc.e[0][2];
			p->deps[5] = p->strainInc.e[0][1];
			//Change stress from Pa to MPa for m7fmaterial_ subroutine input
			p->sigmaN[0] = p->sigmaold.e[0][0] / 1.0e6;
			p->sigmaN[1] = p->sigmaold.e[1][1] / 1.0e6;
			p->sigmaN[2] = p->sigmaold.e[2][2] / 1.0e6;
			p->sigmaN[3] = p->sigmaold.e[1][2] / 1.0e6;
			p->sigmaN[4] = p->sigmaold.e[0][2] / 1.0e6;
			p->sigmaN[5] = p->sigmaold.e[0][1] / 1.0e6;
			// call the Fortran subroutine for m7 constitutive model 
			m7fmaterial_(&timeStep, p->deps, p->epsN, p->sigmaN, p->stateN, &p->internEnergyN,
						 &p->inelasEnergyN, p->microDamageN, p->sigmaNP1, p->stateNP1, &p->internEnergyNP1,
						 &p->inelasEnergyNP1, p->microDamageNP1, &p->weightedDamageNP1, &p->MaxPstrain, &p->MinPstrain,
						 &young, &poisson, &k_1, &k_2, &k_3, &k_4, &k_5, MPMqn, MPMql, MPMqm, MPMw);

			for (int i1 = 0; i1 < 189; i1++)
			{
				p->stateN[i1] = p->stateNP1[i1];
			}
			for (int i1 = 0; i1 < 37; i1++)
			{
				p->microDamageN[i1] = p->microDamageNP1[i1];
			}
			p->internEnergyN = p->internEnergyNP1;
			p->inelasEnergyN = p->inelasEnergyNP1;
			//Change stress output of m7fmaterial_ subroutine from MPa to Pa
			p->sigma.e[0][0] = p->sigmaNP1[0] * 1.0e6;
			p->sigma.e[0][1] = p->sigmaNP1[5] * 1.0e6;
			p->sigma.e[0][2] = p->sigmaNP1[4] * 1.0e6;
			p->sigma.e[1][0] = p->sigmaNP1[5] * 1.0e6;
			p->sigma.e[1][1] = p->sigmaNP1[1] * 1.0e6;
			p->sigma.e[1][2] = p->sigmaNP1[3] * 1.0e6;
			p->sigma.e[2][0] = p->sigmaNP1[4] * 1.0e6;
			p->sigma.e[2][1] = p->sigmaNP1[3] * 1.0e6;
			p->sigma.e[2][2] = p->sigmaNP1[2] * 1.0e6;

			p->eps = p->epsOld + p->strainInc;
			p->epsOld = p->eps;
			p->sigmaold = p->sigma;
			p->PK = p->jacob * p->FFinv.DotMat(p->sigmaold);
		}

#pragma omp for private(p, np, dFt) schedule(static)
		for (int k = 0; k < numParticle; k++)
		{
			p = &(particles[k]);
			p->iforce.setZero();
			p->artviscforce.setZero();
			for (int i = 0; i < p->nNum; i++)
			{
				np = &(particles[p->pid[i]]);
				dFt = p->PK + np->PK;
				p->iforce = p->iforce + dFt.DotVec(p->pG[i]) * np->mass;
				
				//Artificial Viscousity
				double tt = p->velDif[i].DotVec(p->distVec[i]);
				if (tt < 0.0){
					p->artviscforce = p->artviscforce + (tt / (0.0001*kernel*kernel+p->dist[i] * p->dist[i])) * p->pG[i] * np->mass;
				}
			}
			p->iforce = p->iforce / (dens * dens);
			p->artviscforce = ArtViscCoeff * kernel * p->jacob * p->FFinv.DotVec(p->artviscforce) * 3200; //3200 is the approximate speed of sound in concretes

			p->accl = p->iforce + p->artviscforce;
			p->vel = p->vel + p->accl * timeStep;
			if (p->type == 1)
				p->vel.e[2] = 0.0;
			else if (p->type == 3)
				p->vel.e[2] = 0.0 - dCompRate;
			p->disp = p->disp + p->vel * timeStep;
			p->pos = p->posi + p->disp;
		}
	}

	//Export results when needed
	if (iTimestep % save_step < 1)
	{
		SaveResults();
	}

	if (iTimestep % print_step < 1)
	{

		std::printf("\nTimestep: %i, time: %.3e = %.2f %% of %.3e seconds", iTimestep, currenttime, currenttime / maxTime * 100.0e0, maxTime);
	}
	iTimestep++;
	currenttime = currenttime + timeStep;
	std::cout.flush();
}

Vec3i SPHSystem3D::calcCellPos(Vec3d pos)
{
	Vec3i res;

	res.e[0] = floor((pos.e[0] - minX) / cellSize);
	res.e[1] = floor((pos.e[1] - minY) / cellSize);
	res.e[2] = floor((pos.e[2] - minZ) / cellSize);
	return res;
}

uint SPHSystem3D::calcCellHash(Vec3i pos)
{
	uint hash = 0;
	if (pos.e[0] < 0 || pos.e[0] > gridSize.e[0] || pos.e[1] < 0 || pos.e[1] > gridSize.e[1] || pos.e[2] < 0 || pos.e[2] > gridSize.e[2])
	{
		return 0xffffffff;
	}

	hash = pos.e[2] * gridSize.e[0] * gridSize.e[1] + pos.e[1] * gridSize.e[0] + pos.e[0];
	if (hash > totCell)
	{
		std::printf("ERROR! in calcCellHash\n");
		exit(0);
	}
	return hash;
}

void SPHSystem3D::SaveResults()
{
	Particle3D *p;
	double CompStress = 0.0;		//Global Compressive stress
	double Fz=0.0;					//Total force in Z direction
	double Area=0.0;				//Cross-section area
	double fzp;
	for (int j = 0; j < numCrossecParticles; j++)
	{
		p = &(particles[ListCrossecParticles[j]]);
		fzp = (p->sigma.e[2][2] + p->sigma.e[2][1] + p->sigma.e[2][0]) * partDist * partDist;
		Fz -= fzp;
		Area += partDist * partDist;
	}
	CompStress = Fz / Area;
	DataFile.open(DataFileName, std::ios::app);
	DataFile << compression / Lz << "\t" << CompStress << endl;
	DataFile.close();

	// MESH:
	char filename2[200];
	std::sprintf(filename2, "output/ParaviewFiles/SPHM7_%i.vtk", iTimestep);

	std::ofstream paramesh(filename2);
	paramesh << "# vtk DataFile Version 5.1";
	paramesh << "\nvtk output";
	paramesh << "\nASCII";
	paramesh << "\nDATASET POLYDATA";

	paramesh << "\nPOINTS " << numParticle << " float\n";
	for (int i = 0; i < numParticle; i++)
	{
		p = &(particles[i]);
		paramesh << p->posi.e[0] << " " << p->posi.e[1] << " " << p->posi.e[2] << " ";
	}
	paramesh << "\n\nVERTICES 2 " << numParticle;
	paramesh << "\nOFFSETS vtktypeint64\n0 " << numParticle;
	paramesh << "\nCONNECTIVITY vtktypeint64\n";
	for (int i = 0; i < numParticle; i++)
	{
		paramesh << i << " ";
	}
	paramesh << "\nPOINT_DATA " << numParticle << "\n";

	paramesh << "\n\nFIELD Transport_fields 1\n";
	paramesh << "Displacement 3 " << numParticle << " double \n";
	for (int i = 0; i < numParticle; i++)
	{
		p = &(particles[i]);
		paramesh << p->disp.e[0] << " " << p->disp.e[1] << " " << p->disp.e[2] << " ";
	}

	paramesh << "\n\nFIELD Transport_fields 1\n";
	paramesh << "\nDamage 1 " << numParticle << " double \n";
	for (int i = 0; i < numParticle; i++)
	{
		p = &(particles[i]);
		paramesh << p->weightedDamageNP1 << " ";
	}

	paramesh << "\n\nFIELD Transport_fields 1\n";
	paramesh << "\nPrincipleStrainMax 1 " << numParticle << " double \n";
	for (int i = 0; i < numParticle; i++)
	{
		p = &(particles[i]);
		paramesh << p->MaxPstrain << " ";
	}

	paramesh.close();
}

SPHSystem3D::~SPHSystem3D()
{
	free(particles);
	free(cells);
}

void SPHSystem3D::addSingleParticle(Vec3d pos, Vec3d vel, int type)
{
	Particle3D *p = &(particles[numParticle]);
	p->posi = pos;
	p->accl.setZero();
	p->iforce.setZero();
	p->disp.setZero();
	p->type = type;
	p->FF.setZero();

	p->pos = p->posi;
	p->vel = vel;

	for (int i = 0; i < 6; i++)
	{
		p->deps[i] = 0.0e0;
		p->epsN[i] = 0.0e0;
		p->sigmaN[i] = 0.0e0;
		p->sigmaNP1[i] = 0.0e0;
	}
	for (int i = 0; i < 189; i++)
	{
		p->stateN[i] = 0.0e0;
		p->stateNP1[i] = 0.0e0;
	}
	for (int i = 0; i < 37; i++)
	{
		p->microDamageN[i] = 0.0e0;
		p->microDamageNP1[i] = 0.0e0;
	}

	p->internEnergyN = 0.0e0;
	p->inelasEnergyN = 0.0e0;
	p->internEnergyNP1 = 0.0e0;
	p->inelasEnergyNP1 = 0.0e0;

	p->next = NULL;
	numParticle++;
}
