# FloatingporousFoam

![OpenFOAM Version](https://img.shields.io/badge/OpenFOAM-v2206-blueviolet.svg)
![Status](https://img.shields.io/badge/Maintenance-Active-green.svg)

<p align="center">
  <img src="./doc/DeformationMesh_Mooring.png" alt="Moored Floating Porous Structure" width="45%">
  <img src="./doc/Resolved.png" alt="Fully Resolved Porous Structure" width="45%">
</p>

<p align="center">
  <b>A specialized OpenFOAM library for simulating the interaction between <br> waves, currents, and moving porous structures.</b>
</p>

---

## 📋 Table of Contents
- [Description](#-description)
- [Features](#-features)
- [Prerequisites](#-prerequisites)
- [Installation](#-installation)
- [Usage & Configuration](#-usage--configuration)
- [Publications](#-publications)
- [Acknowledgements](#-acknowledgements)

---

## 📖 Description

**FloatingporousFoam** is an open-source solver developed to address the complex hydrodynamics of porous structures in marine environments. It extends the capabilities of standard OpenFOAM solvers to handle porous structure motion, high-fidelity flow resolution, and coupling with mooring dynamics.

This library is particularly useful for **coastal engineering applications**, such as:
*   Floating porous breakwaters
*   Aquaculture cages
*   Permeable marine barriers
*   Porous membranes

## ✨ Features
## ✨ Features
*   **Porous Structure Motion**: Simulation of porous bodies moving under wave and current loads.
*   **Versatile Mesh Support**: Seamlessly supports both **Dynamic Mesh** and **Overset Mesh** techniques for handling complex motions.
*   **Mooring Systems**: Coupling with mooring dynamics (via `foamMooring`).
*   **Fully Resolved Models**: High-fidelity simulation of flow through and around porous media.
---

## 🛠 Prerequisites

Before installing, ensure the following libraries are present in your environment.

| Dependency | Version | Link | Notes |
| :--- | :--- | :--- | :--- |
| **OpenFOAM** | `v2206` | [OpenFOAM.com](https://www.openfoam.com/) | **Strictly tested on v2206**. |
| **olaFlow** | Latest | [GitHub](https://github.com/phicau/olaFlow) | Required for wave generation/absorption. |
| **foamMooring**| Master | [GitLab](https://gitlab.com/hfchen20/foamMooring) | Required only for moored cases. |

> [!NOTE]
> Ensure your environment is running Linux with a standard OpenFOAM setup.

---

## 🚀 Installation

### 1. Download Source Code
Clone this repository or copy the source folder to your OpenFOAM user directory (`$WM_PROJECT_USER_DIR`).

```bash
# Option A: Git Clone (Recommended)
cd $WM_PROJECT_USER_DIR
git clone https://github.com/Dylan000112/floatingporousFoam.git

# Option B: Manual Copy
cp -r floatingporousFoam $WM_PROJECT_USER_DIR/
```

### 2. Compile
Navigate to the directory and execute the make script.

```bash
cd $WM_PROJECT_USER_DIR/floatingporousFoam

# Clean previous builds (optional)
wclean

# Compile library and solvers
./Allwmake
```

---

## 💻 Usage & Configuration

### 1. ControlDict Setup
You must load the appropriate libraries in `system/controlDict` depending on your simulation strategy.

#### A. For Dynamic Mesh Cases
Use this configuration for standard mesh deformation:
```cpp
libs
(
    "libporoussixDoFRigidBodyMotion.so"
);
```

#### B. For Overset Mesh Cases
Use this configuration when using overset (chimera) meshes:
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

---

## 📚 Publications

If you use this code in your research, please cite the following paper:

> **Yiyong Dong, Weikai Tan, Hao Chen, Jing Yuan.** (2024).  
> "Numerical modeling of wave interaction with a porous floating structure consisting of uniform spheres."  
> *Physics of Fluids*, 36(8), 087133.  
> [DOI: 10.1063/5.0222161](https://doi.org/10.1063/5.0222161)

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

---

## 👏 Acknowledgements

*   **[OpenFOAM](https://openfoam.org/)**: The open source CFD toolbox.
*   **[olaFlow](https://github.com/phicau/olaFlow)**: Developed by *Higuera et al.* for simulating wave and porous structure interactions.
*   **[foamMooring](https://gitlab.com/hfchen20/foamMooring/-/tree/master/)**: Developed by *Haifei Chen* for mooring restraints in OpenFOAM.

---

<p align="right">(<a href="#top">back to top</a>)</p>