# FloatingporousFoam

A specialized OpenFOAM library for simulating the interaction between waves, currents, and moving porous structures.

<p align="center">
  <img src="./doc/DeformationMesh_Mooring.png" alt="Moored Floating Porous Structure" width="45%">
  <img src="./doc/Resolved.png" alt="Fully Resolved Porous Structure" width="45%">
</p>

## Description

**FloatingporousFoam** is an open-source solver developed to address the complex hydrodynamics of porous structures in marine environments. It extends the capabilities of standard OpenFOAM solvers to handle:

*   **Porous Structure Motion:** Simulation of porous bodies moving under wave and current loads.
*   **Fully Resolved Models:** High-fidelity simulation of flow through and around porous media.
*   **Mooring Systems:** Coupling with mooring dynamics (via foamMooring).

This library is particularly useful for coastal engineering applications, such as floating breakwaters, aquaculture cages, and other permeable marine structures.

## Compatibility

*   **OpenFOAM Version:** Tested specifically on **OpenFOAM v2206**.
*   **System:** Linux (Standard OpenFOAM environment).

## Dependencies

Before installing, ensure the following libraries are present in your environment:

1.  **olaFlow** (Required for wave generation and absorption)
    *   Source: [github.com/phicau/olaFlow](https://github.com/phicau/olaFlow)
2.  **foamMooring** (Required for moored cases only)
    *   Source: [gitlab.com/hfchen20/foamMooring](https://gitlab.com/hfchen20/foamMooring/-/tree/master/)

## Installation

1.  **Download the source code** to your user directory:
    ```bash
    cp -r floatingporousFoam $WM_PROJECT_USER_DIR/
    ```

2.  **Compile the library and solvers**:
    ```bash
    cd $WM_PROJECT_USER_DIR/floatingporousFoam
    ./Allwmake
    ```

## Usage & Configuration

### 1. ControlDict Setup
Depending on your simulation type, you must load the appropriate libraries in `system/controlDict`.

*   **For Dynamic Mesh cases:**
    ```cpp
    libs
    (
        "libporoussixDoFRigidBodyMotion.so"
    );
    ```

*   **For Overset Mesh cases:**
    ```cpp
    libs
    (
        "liboversetPadding.so"
        "libporoussixDoFRigidBodyMotion.so"
    );
    ```

### 2. Running the Solver
Run the simulation using the custom solver:
```bash
porousolaDyMFlow
```

## Publications

If you use this code in your research, please cite the following paper:

> **Yiyong Dong, Weikai Tan, Hao Chen, Jing Yuan.** (2024). "Numerical modeling of wave interaction with a porous floating structure consisting of uniform spheres." *Physics of Fluids*, 36(8), 087133. [DOI: 10.1063/5.0222161](https://doi.org/10.1063/5.0222161)

**BibTeX:**
```bibtex
@article{Dong2024,
  author = {Dong, Yiyong and Tan, Weikai and Chen, Hao and Yuan, Jing},
  title = {Numerical modeling of wave interaction with a porous floating structure consisting of uniform spheres},
  journal = {Physics of Fluids},
  volume = {36},
  number = {8},
  pages = {087133},
  year = {2024},
  doi = {10.1063/5.0222161}
}
```

## Acknowledgements

*   **[OpenFOAM](https://openfoam.org/)**: The open source CFD toolbox.
*   **[olaFlow](https://github.com/phicau/olaFlow)**: Developed by Higuera et al. for simulating wave and porous structure interactions.
*   **[foamMooring](https://gitlab.com/hfchen20/foamMooring/-/tree/master/)**: Developed by Haifei Chen for mooring restraints in OpenFOAM.