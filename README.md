# SPH-M7 Microplane Concrete Modeling Repository

## Overview
This is the repository for the source code accompanying the paper "**Modeling concrete failure with smoothed particle hydrodynamics using the Microplane (M7) constitutive model.** M. N. Rahimi & G. Moutsanidis". This repository provides an open-access implementation of the SPH-M7 Microplane algorithm tailored for modeling concrete materials, including all scripts and data necessary for reproducing the results for selected cases presented in the paper.

The repository is designed to promote transparency, reproducibility, and further development by the research community.

## Features
- Implementation of the SPH-M7 algorithm.
- Implementation of a particle-to-particle contact algorithm for SPH.
- Material model specifically for concrete, including failure and dynamic response modeling.
- Includes the M7 constitutive model based on the paper "Caner, F.C., and Bažant, Z.P. (2013). \`\`Microplane model M7 for plain concrete: I. formulation." *ASCE J. of Engrg. Mechanics* 139 (12), Dec., 1714--1723. The original source can be downloaded [here](http://www.civil.northwestern.edu/people/bazant/m7-coding/m7_cyc_schell_v1.f).

## Installation
1. Clone the repository:
   ```bash
   git clone https://github.com/naqibr/SPH-M7-Concrete.git
   cd SPH-M7-Concrete
   ```

2. Install dependencies using your package manager or manually as needed.

3. Build the code:
   ```bash
   make
   ```
This will generate the executable `Concrete`.

4. Run the executable from the build directory:
   ```bash
   ./Concrete -np [number of parallel threads] (for compression of short prism)
   ```
   or
   ```bash
   ./Concrete -np [number of parallel threads] -type [1 or 2] -D [D] (for Type 1 and Type 2 size effects)
   ```
   Forexample
   ```bash
   ./Concrete -np 4 
   ```
   or
   ```bash
   ./Concrete -np 4 -type 1 -D 40.0e-3
   ```
This will run the executable with 4 parallel threads. The default thread number is 1.

5. Output files will be generated in the `output/` directory.

### Visualization
Use ParaView to visualize the `.vtk` files generated during the simulation. Sometimes to visualize the generated particles, one needs to change the "Display (GeometryRepresentation)" option to "Point Gaussian".

### Examples
The repository includes two branches each tailored for one of the selected examples in the paper:
- **main**: This branch generates results for the case "4.1. Short rectangular prism under uniaxial compression" reported in the article.
- **Type12_SizeEffect**: This branch generates results for the case "4.7. Type 1 and Type 2 size effect" reported in the article. This branch also covers the application of the contact algorithm presented in Section 2.3. (Particle-to-particle contact for SPH) of the paper.

## File Descriptions
- **Makefile**: Contains commands to compile the codes into an executable.
- **Main.cpp**: The main file of the executable that creates the SPH object as `sph = new SPHSystem3D(argc, argv);`.
- **SPHSystem.h**: Contains the class `SPHSystem3D`.
- **SPHSystem.cpp**: Contains the function definitions for the class `SPHSystem3D`.
- **particle.h**: Contains the class `Particle3D`.
- **Matrix3D.h**: Contains the class `Mat3d`, a local class for 3x3 matrices with necessary matrix arithmetic operations.
- **Vector3D.h**: Contains the class `Vec3d`, a local class for 3D vectors with necessary vector arithmetic operations.
- **Vec3DMatrix3D.cpp**: Contains some arithmetic operations related to `Vec3d` and `Mat3d` classes.
- **m7fmaterial.f**: A cleaned-up version of the file provided [here](http://www.civil.northwestern.edu/people/bazant/m7-coding/m7_cyc_schell_v1.f), containing FORTRAN subroutines for initializing the M7 constitutive model parameters and evaluating the model itself.
- **params.inc**: Contains the parameters used in `m7fmaterial.f`.

## How the Code Works
1. **Initialization**:
   - The object `sph` is created as `sph = new SPHSystem3D(argc, argv);`.
   - The object is initialized in `SPHSystem3D(argc, argv)` in `SPHSystem.cpp`, where system variables are set up, particles are generated, neighbor search is performed, SPH kernels and their derivatives are corrected, and variables for the M7 constitutive model are initialized.

2. **Time Integration**:
   - Time integration begins, and `sph->OneStepCalculations();` is called for each time step.
   - Inside `OneStepCalculations`, the following steps occur:
     - Deformation gradient and strain rate are calculated.
     - `m7fmaterial_` is called to evaluate new Cauchy stress tensors for each particle.
     - Artificial viscosity forces are calculated.
     - Contact forces are calculated if necessary (Type12_SizeEffect branach only).
     - Results are exported every required time step.

## Citation
If this code is helpfull, whether partially or entirely, in your research, please cite the relevant paper(s):

### Main Paper
Paper Title: **Modeling concrete failure with smoothed particle hydrodynamics using the Microplane (M7) constitutive model**  
Authors: **M. Naqib Rahimi, Georgios Moutsanidis**  
Journal: **[Journal Name]**  
DOI: **[Insert DOI here]**  

### Related Papers
This program also implements the Total Lagrangian Smoothed Particle Hydrodynamics formulations presented in the following papers:

1. **A smoothed particle hydrodynamics approach for phase field modeling of brittle fracture,**  
   *M. N. Rahimi and G. Moutsanidis,*  
   *Computer Methods in Applied Mechanics and Engineering,* Aug. 2022.  
   DOI: [https://doi.org/10.1016/j.cma.2022.115191](https://doi.org/10.1016/j.cma.2022.115191)  

3. **Modeling dynamic brittle fracture in functionally graded materials using hyperbolic phase field and smoothed particle hydrodynamics,**  
   *M. N. Rahimi and G. Moutsanidis,*  
   *Computer Methods in Applied Mechanics and Engineering,* Nov. 2022.  
   DOI: [https://doi.org/10.1016/j.cma.2022.115642](https://doi.org/10.1016/j.cma.2022.115642)

4. **An SPH-based FSI framework for phase-field modeling of brittle fracture under extreme hydrodynamic events,**  
   *M. N. Rahimi and G. Moutsanidis,*  
   *Engineering with Computers,* Aug. 2023.  
   DOI: [https://doi.org/10.1007/s00366-023-01857-0](https://doi.org/10.1007/s00366-023-01857-0)

5. **IGA-SPH: Coupling isogeometric analysis with smoothed particle hydrodynamics for air-blast–structure interaction,**  
   *M. N. Rahimi and G. Moutsanidis,*  
   *Engineering with Computers,* May 2024.  
   DOI: [https://doi.org/10.1007/s00366-024-01978-0](https://doi.org/10.1007/s00366-024-01978-0)  

---

### BibTeX
```bibtex
@article{RahimiSPHM7,
  title={Modeling concrete failure with smoothed particle hydrodynamics using the Microplane (M7) constitutive model},
  author={Rahimi, M. Naqib and Moutsanidis, Georgios},
  journal={Journal Name},
  year={2025},
  volume={XX},
  pages={XX--XX},
  doi={DOI}
}

@article{RahimiSPHPhaseField,
  title={A smoothed particle hydrodynamics approach for phase field modeling of brittle fracture},
  author={Rahimi, M. Naqib and Moutsanidis, Georgios},
  publisher = {Elsevier},
  month = {8},
  year = {2022},
  journal = {\textbf{Computer Methods in Applied Mechanics and Engineering}},
  doi = {https://doi.org/10.1016/j.cma.2022.115191}
}

@article{RahimiDynamicBrittleFracture,
  title={Modeling dynamic brittle fracture in functionally graded materials using hyperbolic phase field and smoothed particle hydrodynamics},
  author={Rahimi, M. Naqib and Moutsanidis, Georgios},
  journal = {\textbf{Computer Methods in Applied Mechanics and Engineering}},
  month = {11},
  year = {2022},
  doi = {https://doi.org/10.1016/j.cma.2022.115642}
}

@article{RahimiSPHFsiFramework,
  title={An SPH-based FSI framework for phase-field modeling of brittle fracture under extreme hydrodynamic events},
  author={Rahimi, M. Naqib and Moutsanidis, Georgios},
  month = {8},
  year = {2023},
  journal = {\textbf{Engineering with Computers}},
  doi = {https://doi.org/10.1007/s00366-023-01857-0}
}

@article{RahimiIGASPH,
  title={IGA-SPH: Coupling isogeometric analysis with smoothed particle hydrodynamics for air-blast–structure interaction},
  author={Rahimi, M. Naqib and Moutsanidis, Georgios},
  journal={\textbf{Engineering with Computers}},
  pages={1--22},
  year={2024},
  month={05},
  publisher={Springer},
  doi = {https://doi.org/10.1007/s00366-024-01978-0}
}
```

## License
This project is licensed under the MIT License. See the `LICENSE` file for details.

## Contact
For questions or support, please contact:
- **Naqib Rahimi**: Naqib.rahimy123@gmail.com
- **GitHub Issues**: Use the Issues tab for bug reports or feature requests.
---

We hope you find this repository useful! Happy modeling!

