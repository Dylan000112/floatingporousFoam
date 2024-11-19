# FloatingporousFoam
![Image loading](./doc/pics.png)
**FloatingporousFoam** is a library for simulating porous structure's motion under wave and current in OpenFOAM

## Usage
Case have been tested for OpenFOAM v2206 only.

Requires `olaFlow`: https://github.com/phicau/olaFlow.

For the moored-floating porous sturctures, requires `foamMooring`: https://gitlab.com/hfchen20/foamMooring/-/tree/master/.

## Installation
1. Copy the `floatingporousFoam` floder to `$WM_PROJECT_USER_DIR`.
2. Run `./Allwmake`.
3. Include the following lines in your `contorlDict`:
   
   `libs
      (
        poroussixDoFRigidBodyMotion
      );`

4. Run your case with `porousolaDyMFlow`

## Publications
Yiyong Dong, Weikai Tan, Hao Chen, Jing Yuan; Numerical modeling of wave interaction with a porous floating structure consisting of uniform spheres. Physics of Fluids 1 August 2024; 36 (8): 087133. https://doi.org/10.1063/5.0222161

## Acknowledgements
**OpenFOAM** is free, open source software for computational fluid dynamics (CFD), developed primarily by CFD Direct, on behalf of the OpenFOAM Foundation. (https://openfoam.org/)

**olaFlow** is an open source project developed within the OpenFOAM® framework as a continuation of the work in Higuera et al. for simulating wave and porous structure interactions in the coastal and offshore fields. (https://github.com/phicau/olaFlow)

**foamMooring** is a mooring restraints library for rigid body motions in OpenFOAM®, developed by Haifei Chen. (https://gitlab.com/hfchen20/foamMooring/-/tree/master/).
