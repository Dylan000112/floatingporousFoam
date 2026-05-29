# FloatingporousFoam: Tutorials and Validation Cases

This repository provides OpenFOAM case templates for modeling the hydrodynamic responses of Porous Floating Breakwater (PFB).

## 2D Macroscopic Methods
* **`PFB_2D_oversetMesh` & `PFB_2D_DeformationMesh`**: 2D simulations of floating porous breakwater under wave action using a macroscopic force model (Dong et al., 2024).
* **`PFB_2D_Mooring`**: Coupled fluid-structure interaction simulations utilizing the `foamMooring` module for restoring force calculations.
* **`Force_validation`**: Verification case ensuring numerical consistency between macroscopic resistance forces and integrated pressure gradients.

## 3D Macroscopic Methods
* **`Irregular_wave_3D_OversetMesh` & `Irregular_wave_3D_DeformationMesh`**: 3D dynamic response simulations of PFBs in multi-directional irregular waves, with oversetmesh recommended for large-amplitude motions.
* **`WC_EmptyFlume_3D`**: 3D numerical wave tank benchmark for synthesizing and validating multi-directional irregular wave fields.

## 3D Micro-scale Resolved LES
* **`PFB_3D_Resolved_LES`**: 3D high-fidelity pore-scale resolved Large Eddy Simulation for calibrating and validating macroscopic force formulations.

