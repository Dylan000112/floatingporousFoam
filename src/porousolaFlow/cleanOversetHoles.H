/*---------------------------------------------------------------------------*\
  File
    cleanOversetHoles.H

  Description
    Utility to clamp flow variables in 'HOLE' cells for Overset meshes.
    Prevents spurious velocities and pressure spikes in inactive regions.

    - Target field: volScalarField "cellTypes" (value ~ 2.0)
    - Actions: Clamps U -> 0, p_rgh -> 0, alpha.water -> 0
    - Parallel: Fully MPI compatible with global reduction

\*---------------------------------------------------------------------------*/

#ifndef cleanOversetHoles_H
#define cleanOversetHoles_H

#include "fvCFD.H"
#include "volFields.H"
#include "DynamicList.H"

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void cleanOversetHoles(const fvMesh& mesh)
{
    // Constants for hole identification
    static const word cellTypesName("cellTypes");
    static const scalar holeValue(2.0);
    static const scalar tolerance(0.1);

    // 1. Guard: Check if cellTypes exists
    if (!mesh.foundObject<volScalarField>(cellTypesName))
    {
        return;
    }

    // 2. Access raw field data
    const auto& cellTypes = mesh.lookupObject<volScalarField>(cellTypesName);
    const scalarField& cTypes = cellTypes.primitiveField();

    // 3. Acquire pointers to target fields (if they exist)
    // Using pointers allows for optional handling of fields without exceptions
    volVectorField* UPtr = 
        mesh.foundObject<volVectorField>("U") 
        ? &const_cast<volVectorField&>(mesh.lookupObject<volVectorField>("U")) 
        : nullptr;

    volScalarField* pPtr = 
        mesh.foundObject<volScalarField>("p_rgh") 
        ? &const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("p_rgh")) 
        : nullptr;

    volScalarField* alphaPtr = 
        mesh.foundObject<volScalarField>("alpha.water") 
        ? &const_cast<volScalarField&>(mesh.lookupObject<volScalarField>("alpha.water")) 
        : nullptr;

    // 4. Processing Loop
    label nHolesLocal = 0;

    // Direct primitive field access for performance in the loop
    vectorField* UData = UPtr ? &UPtr->primitiveFieldRef() : nullptr;
    scalarField* pData = pPtr ? &pPtr->primitiveFieldRef() : nullptr;
    scalarField* alphaData = alphaPtr ? &alphaPtr->primitiveFieldRef() : nullptr;

    forAll(cTypes, i)
    {
        // Check for HOLE type with tolerance for floating point safety
        if (mag(cTypes[i] - holeValue) < tolerance)
        {
            if (UData)     (*UData)[i] = vector::zero;
            if (pData)     (*pData)[i] = 0.0;
            if (alphaData) (*alphaData)[i] = 0.0;

            nHolesLocal++;
        }
    }

    // 5. Parallel Reduction & Logging
    // Calculate global sum across all processors
    label nHolesGlobal = nHolesLocal;
    reduce(nHolesGlobal, sumOp<label>());

    // Only log on master if holes were actually found and cleaned
    if (nHolesGlobal > 0 && Pstream::master())
    {
        // Construct list of cleaned fields for clean output
        DynamicList<word> cleanedFields;
        if (UData) cleanedFields.append("U");
        if (pData) cleanedFields.append("p_rgh");
        if (alphaData) cleanedFields.append("alpha.water");

        Info<< "    [HoleClean] Clamped " << nHolesGlobal 
            << " hole cells. Fields: " << cleanedFields << endl;
    }
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

#endif

// ************************************************************************* //
