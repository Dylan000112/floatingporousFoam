### Core Source Modules

* **`porousforce`**: Core library for porous force calculation. Solves the macroscopic porous resistance models and performs accurate volumetric integration of fluid forces within the porous media zone.
* **`poroussixDoFRigidMotion` & `inverDistancePadding`**: Motion coupling and mesh adaptation modules. Implements the coupling of integrated volumetric porous forces with 6-DoF rigid body dynamics, and provides underlying spatial displacement and interpolation adaptation for Dynamic Mesh and Overset Mesh methodologies.
* **`porousolaFlow`**: Hydrodynamic solver for moving porous bodies, specifically designed to solve the transient hydrodynamic responses of floating porous structures.
