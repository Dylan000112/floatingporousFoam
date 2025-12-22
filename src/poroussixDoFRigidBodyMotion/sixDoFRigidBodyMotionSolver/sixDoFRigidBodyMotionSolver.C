/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "sixDoFRigidBodyMotionSolver.H"
#include "addToRunTimeSelectionTable.H"
#include "polyMesh.H"
#include "pointPatchDist.H"
#include "pointConstraints.H"
#include "uniformDimensionedFields.H"
#include "porousforces.H"
#include "mathematicalConstants.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
    defineTypeNameAndDebug(sixDoFRigidBodyMotionSolver, 0);

    addToRunTimeSelectionTable
    (
        motionSolver,
        sixDoFRigidBodyMotionSolver,
        dictionary
    );
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionSolver::sixDoFRigidBodyMotionSolver
(
    const polyMesh& mesh,
    const IOdictionary& dict
)
:

    displacementMotionSolver(mesh, dict, typeName),
    motion_
    (
        coeffDict(),
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            mesh.time().timeName(),
            "uniform",
            mesh
        ).typeHeaderOk<IOdictionary>(true)
      ? IOdictionary
        (
            IOobject 
            (
                "sixDoFRigidBodyMotionState",
                mesh.time().timeName(),
                "uniform",
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE,
                false
            )
        )
      : coeffDict(),
        mesh.time()
    ),
    patches_(coeffDict().get<wordRes>("patches")),
    patchSet_(mesh.boundaryMesh().patchSet(patches_)),
    di_(coeffDict().get<scalar>("innerDistance")),
    do_(coeffDict().get<scalar>("outerDistance")),
    cramp_(coeffDict().get<scalar>("cramp")),
    test_(coeffDict().getOrDefault("test", false)),
    rhoInf_(1.0),
    rhoName_(coeffDict().getOrDefault<word>("rho", "rho")),
    scale_
    (
        IOobject
        (
            "motionScale",
            mesh.time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE,
            false
        ),
        pointMesh::New(mesh),
        dimensionedScalar(dimless, Zero)
    ),
    curTimeIndex_(-1),
    cOfGdisplacement_(coeffDict().getOrDefault<word>("cOfGdisplacement", "none"))
{
    if (rhoName_ == "rhoInf")
    {
        coeffDict().readEntry("rhoInf", rhoInf_);
    }

/**************************************************************************************************/
    // Calculate scaling factor everywhere

    // ... 在 sixDoFRigidBodyMotionSolver 构造函数中 ...
    
    {
        const pointMesh& pMesh = pointMesh::New(mesh);

        // 1. 初始化距离场为 "无穷大" (GREAT)
        // 这一步对 Overset 至关重要，保证了没被计算的点（背景网格）距离永远是无限大
        pointPatchDist pDist(pMesh, patchSet_, points0());        
        pDist.primitiveFieldRef() = Foam::GREAT; 

        // 2. 读取用户开关 "overset"
        // 默认为 false，即默认当作普通动网格处理
        const bool isOverset = coeffDict().getOrDefault("overset", false);

        const pointField& pts = points0();
        scalarField& pDistField = pDist.primitiveFieldRef();

        // 预计算几何参数
        const vector cofr = motion_.centreOfRotation();
        // 防止除以0，加一个极小值检查
        const scalar Px_sq = sqr(motion_.Px() > SMALL ? motion_.Px()/2.0 : 1.0);
        const scalar Py_sq = sqr(motion_.Py() > SMALL ? motion_.Py()/2.0 : 1.0);
        const scalar Pz_sq = sqr(motion_.Pz() > SMALL ? motion_.Pz()/2.0 : 1.0);

        // =========================================================
        // 分支逻辑：Overset (拓扑隔离) vs DynamicMesh (全场变形)
        // =========================================================

        if (isOverset)
        {
            // --- 模式 A: Overset 模式 ---
            
            // 必须提供 Zone 的名字
            if (!coeffDict().found("porousZone"))
            {
                FatalErrorInFunction 
                    << "You set 'overset on', so you must specify 'porousZone' name!" 
                    << exit(FatalError);
            }
            
            word zoneName = coeffDict().get<word>("porousZone");
            
            Info << "SixDoF: Overset mode ON. limiting motion to zone: " << zoneName << endl;

            const label zoneID = mesh.pointZones().findZoneID(zoneName);
            if (zoneID == -1)
            {
                FatalErrorInFunction 
                    << "PointZone '" << zoneName << "' not found in polyMesh." 
                    << " Please run topoSet."
                    << exit(FatalError);
            }

            const labelList& zonePoints = mesh.pointZones()[zoneID];

            // 只循环 Zone 里的点
            forAll(zonePoints, i)
            {
                label pointI = zonePoints[i];
                vector d = pts[pointI] - cofr;
                pDistField[pointI] = sqrt(
                    (d.x()*d.x())/Px_sq + (d.y()*d.y())/Py_sq + (d.z()*d.z())/Pz_sq
                );
            }
        }
        else
        {
            // --- 模式 B: 普通动网格模式 (变形) ---
            
            Info << "SixDoF: Overset mode OFF. Applying global deformation." << endl;

            // 循环全场所有点
            forAll(pts, pointI)
            {
                vector d = pts[pointI] - cofr;
                pDistField[pointI] = sqrt(
                    (d.x()*d.x())/Px_sq + (d.y()*d.y())/Py_sq + (d.z()*d.z())/Pz_sq
                );
            }
        }

        // 3. 计算 Scale (通用)
        // Overset模式下：背景网格 pDist=GREAT -> scale=0 -> 不动
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    (do_ - pDist.primitiveField())/(do_ - di_),
                    scalar(0)
                ),
                scalar(1)
            );



        // Convert the scale function to a cosine
        scale_.primitiveFieldRef() =
            min
            (
                max
                (
                    0.5
                  - 0.5
                   *cos(scale_.primitiveField()
                   *Foam::constant::mathematical::pi),
                    scalar(0)
                ),
                scalar(1)
            );
        pointConstraints::New(pMesh).constrain(scale_);
        scale_.write();
    }
}


/**************************************************************************************************/


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::sixDoFRigidBodyMotion&
Foam::sixDoFRigidBodyMotionSolver::motion() const
{
    return motion_;
}


Foam::tmp<Foam::pointField>
Foam::sixDoFRigidBodyMotionSolver::curPoints() const
{
    tmp<pointField> newPoints
    (
        points0() + pointDisplacement_.primitiveField()
    );

    if (!moveAllCells())
    {
        tmp<pointField> ttransformedPts(new pointField(mesh().points()));
        pointField& transformedPts = ttransformedPts.ref();

        UIndirectList<point>(transformedPts, pointIDs()) =
            pointField(newPoints.ref(), pointIDs());

        return ttransformedPts;
    }

    return newPoints;
}


void Foam::sixDoFRigidBodyMotionSolver::solve() 
{
    
    const Time& t = mesh().time();

    if (mesh().nPoints() != points0().size()) 
    {
        FatalErrorInFunction
            << "The number of points in the mesh seems to have changed." << endl
            << "In constant/polyMesh there are " << points0().size()
            << " points; in the current mesh there are " << mesh().nPoints()
            << " points." << exit(FatalError);
    }


    // Store the motion state at the beginning of the time-step
    bool firstIter = false;
    if (curTimeIndex_ != this->db().time().timeIndex())

    {
        motion_.newTime();
        curTimeIndex_ = this->db().time().timeIndex();
        firstIter = true;
    }

    dimensionedVector g("g", dimAcceleration, Zero);

    if (db().time().foundObject<uniformDimensionedVectorField>("g")) 
    {
        g = db().time().lookupObject<uniformDimensionedVectorField>("g");
    }
    else
    {
        coeffDict().readIfPresent("g", g);
    }


    const scalar ramp = 1.0;

    if (test_)   
    {
        motion_.update
        (
            firstIter,
            ramp*(motion_.mass()*g.value()),
            ramp*(motion_.mass()*(motion_.momentArm() ^ g.value())),
            t.deltaTValue(),
            t.deltaT0Value()
        );
          
    }
            
    else
    {
        dictionary forcesDict;

        forcesDict.add("type", functionObjects::porousforces::typeName);
        forcesDict.add("patches", patches_);
        forcesDict.add("rhoInf", rhoInf_);
        forcesDict.add("rho", rhoName_);
        forcesDict.add("CofR", motion_.centreOfRotation());
        forcesDict.add("myrhoporo", motion_.myrhoporo());
        forcesDict.add("Linear_vel", motion_.v());
        forcesDict.add("Angular_vel", motion_.omega());
        vector oldPos = motion_.centreOfRotation();        
        functionObjects::porousforces f("forces", db(), forcesDict);  

        f.calcForcesMoments();          
 
        Info<< "Total porous force = " << f.forceEff() << nl << endl;
        Info<< "Total porous moment = " << f.momentEff() << nl << endl; 


        Info<< "PFB velocity = " << motion_.v() << nl << endl; 
        Info<< "cramp = " << cramp_ << nl << endl; 	
        	
        	
        Info<< "cramp*motion_.v() = " << cramp_*motion_.v() << nl << endl; 
        Info<< "f.forceEff() + cramp_*motion_.v() = " << f.forceEff() + cramp_*motion_.v() << nl << endl; 
        
        	
        motion_.update
        (
            firstIter, 
            ramp*(f.forceEff() - cramp_*motion_.v() + motion_.Rigidmass()*g.value()),           
            
            ramp
           *(
               f.momentEff()
             + motion_.Rigidmass()*(motion_.momentArm() ^ g.value())   
            ),           
            
            t.deltaTValue(),  
            t.deltaT0Value()  
        );
        
        if (cOfGdisplacement_ != "none")
        {
            if
            (
                db().time().foundObject<uniformDimensionedVectorField>
                (
                    cOfGdisplacement_
                )
            )
            {
                auto& disp =
                    db().time().lookupObjectRef<uniformDimensionedVectorField>
                    (
                        cOfGdisplacement_
                    );

                disp += (motion_.centreOfRotation() - oldPos);
            }
        }
    }

    // Update the displacements
    pointDisplacement_.primitiveFieldRef() =
        motion_.transform(points0(), scale_) - points0();
    
    // Displacement has changed. Update boundary conditions
    pointConstraints::New
    (
        pointDisplacement_.mesh()
    ).constrainDisplacement(pointDisplacement_);
    
    
    if (mesh().time().writeTime())
    {
        this->writeObject(IOstreamOption(), true);
    }
}


bool Foam::sixDoFRigidBodyMotionSolver::writeObject
(
    IOstreamOption streamOpt,
    const bool valid
) const
{
    IOdictionary dict
    (
        IOobject
        (
            "sixDoFRigidBodyMotionState",
            mesh().time().timeName(),
            "uniform",
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE,
            false
        )
    );

    motion_.state().write(dict);
    
    
    return dict.regIOobject::write();
    
}


bool Foam::sixDoFRigidBodyMotionSolver::read()
{
    if (displacementMotionSolver::read())
    {
        motion_.read(coeffDict());

        return true;
    }

    return false;
}


// ************************************************************************* //
