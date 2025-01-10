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
#include "SPHSystem.h"

SPHSystem3D *sph;

int main(int argc, char *argv[])
{
	sph = new SPHSystem3D(argc, argv);
	std::cout << "\n >> Time integration started:\n";
	time_t startTime;
    time(&startTime);
	bool estimated30Seconds = false;
    bool estimated60Seconds = false;
    for (double tt = 0.0; tt <= sph->maxTime; tt = tt + sph->timeStep) {
		sph->OneStepCalculations();
		if (!estimated30Seconds) {
            time_t endTime;
    		time(&endTime);
            double elapsedSeconds = difftime(endTime, startTime);

            if (elapsedSeconds >= 30.0) {
                double estimatedCompletionTime = elapsedSeconds / (tt+sph->timeStep) * sph->maxTime;
				int days = static_cast<int>(estimatedCompletionTime  / (24.0*3600.0));
				int hours = static_cast<int>((estimatedCompletionTime-24.0*3600.0*days) / 3600.0);
				int minutes = static_cast<int>((estimatedCompletionTime -24.0*3600.0*days - hours * 3600.0) / 60.0);
				int seconds = static_cast<int>(estimatedCompletionTime -24.0*3600.0*days - hours * 3600.0 - minutes * 60.0);
				std::cout << "\n >> Finishes in: " <<days<<" Day(s) "<< hours << " Hour(s) "<< minutes << " Minute(s) " <<seconds << " Second(s)";
				estimated30Seconds = true;
            }
        }

        if (!estimated60Seconds) {
            time_t endTime;
    		time(&endTime);
            double elapsedSeconds = difftime(endTime, startTime);

            if (elapsedSeconds >= 60.0) {
				double estimatedCompletionTime = elapsedSeconds / (tt+sph->timeStep) * sph->maxTime;
				int days = static_cast<int>(estimatedCompletionTime  / (24.0*3600.0));
				int hours = static_cast<int>((estimatedCompletionTime-24.0*3600.0*days) / 3600.0);
				int minutes = static_cast<int>((estimatedCompletionTime -24.0*3600.0*days - hours * 3600.0) / 60.0);
				int seconds = static_cast<int>(estimatedCompletionTime -24.0*3600.0*days - hours * 3600.0 - minutes * 60.0);
				std::cout << "\n >> Finishes in: " <<days<<" Day(s) "<< hours << " Hour(s) "<< minutes << " Minute(s) " <<seconds << " Second(s)";
				estimated60Seconds = true;
            }
        }
	}

	free(sph);
	return 0;
}