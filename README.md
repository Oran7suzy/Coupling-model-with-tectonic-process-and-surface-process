# Coupling-model-with-tectonic-process-and-surface-process
This is couping model with tectonic process and surface process. The model describes deep lithospheric deformation based on viscoelastic-plastic rheology. At the surface, erosion and deposition processes such as river incision and hillslope diffusion were resolved using the FastScape.

The repository includes:
coupled_model.for — main Fortran source code containing UMAT, UMESHMOTION, and TEMPARRAY_1 subroutines.
job-1.inp — example Abaqus input file.
surface_inp.py — Python preprocessing script for generating surface-process input files.

---
# 1. Requirements
<Abaqus>-tested on 2020/2021/2022
<Fortran compiler>-GNU Fortran)
<Python 3.x>-for preprocessing)
<FastScape>-Version 2.8.4 (Reference: https://fastscape.org/fastscapelib-fortran/#_fastscape_set_dt)

---
# 2. Subroutine Notes and Required Settings
**2.1 TEMPARRAY_1 Subroutine**
BASE_PATH — directory containing surface-process input files  
PATH1-9 — filenames for required datasets

These files include:
- 'mo_pre.txt'
- 'model_param.txt'
- 'ada_ele.txt'
- 'u_slip.txt'
- 'hh.txt'
- 'drainage_area.txt'
- 'total_erosion.txt'
- 'erosion_rate.txt'
- 'catchment.txt'
- 'basement.txt'

Use the preprocessing script:
This script outputs all required files into the folder specified by `BASE_PATH`.

**2.2 UMAT Subroutine**
You must define 9 parameters:
1. Young’s modulus  
2. Poisson’s ratio  
3. Density  
4. Viscosity  
5. Cohesion  
6. Cohesion drop  
7. Internal friction angle  
8. Plastic modulus  
9. (Optional) Additional softening/strength parameters

NOTE: Although parameters can be hard-coded in coupled_model.for, it is strongly recommended to define them in the Abaqus using:
This allows assigning different material properties to different regions.

**2.3 UMESHMOTION Subroutine**
You must define:
- Erosion coefficient (K)
- River incision exponents M and N (power-law)
- Hillslope diffusion coefficient
- Sedimentation coefficient
- Rainfall rate
- Boundary conditions (fixed or mobile nodes)
- Time step scaling

These parameters follow FastScape formulations.  
Detailed explanation:  
https://fastscape.org/fastscapelib-fortran/#_fastscape_set_dt

You may define spatially variable fot rainfall, erodibility, and diffusion coefficient  
These can also be generated using surface_inp.py.

---
# 3. Running the Model
Step 1 - Prepare DEM or initial inp file from abaqus  
Step 2 - Run the surface_inp.py to get file required surface-process input files (2.1)
Step 3 - Set all parameters for coupled_model.for.
Step 4 - Run Abaqus job

---
# 5. Data Availability
All model input files, source code, and preprocessing scripts are included in this repository.  
Users may modify parameters in the Fortran code to reproduce or extend the simulations.

---
# 6. Citation
If you use this model, please cite this repository (Zenodo DOI will appear here after release):

---
# 7. Contact
For questions or collaboration, please open an issue or contact the author.
