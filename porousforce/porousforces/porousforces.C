/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2015-2022 OpenCFD Ltd.
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

#include "porousforces.H"
#include "fvcGrad.H"
#include "porosityModel.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"
#include "cartesianCS.H"
#include "addToRunTimeSelectionTable.H"
//#include "math.H"
#include "fvCFD.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(porousforces, 0);
    addToRunTimeSelectionTable(functionObject, porousforces, dictionary);
}
}


// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void Foam::functionObjects::porousforces::setCoordinateSystem 
(
    const dictionary& dict,
    const word& e3Name,
    const word& e1Name
)
{
    coordSysPtr_.clear();

    point origin(Zero);
    if (dict.readIfPresent<point>("CofR", origin))
    {
        const vector e3 = e3Name == word::null ?
            vector(0, 0, 1) : dict.get<vector>(e3Name);
        const vector e1 = e1Name == word::null ?
            vector(1, 0, 0) : dict.get<vector>(e1Name);

        coordSysPtr_.reset(new coordSystem::cartesian(origin, e3, e1));
    }
    else
    {
        // The 'coordinateSystem' sub-dictionary is optional,
        // but enforce use of a cartesian system if not found.

        if (dict.found(coordinateSystem::typeName_()))
        {
            // New() for access to indirect (global) coordinate system
            coordSysPtr_ =
                coordinateSystem::New
                (
                    obr_,
                    dict,
                    coordinateSystem::typeName_()
                );
        }
        else
        {
            coordSysPtr_.reset(new coordSystem::cartesian(dict));
        }
    }
}


Foam::volVectorField& Foam::functionObjects::porousforces::force() 
{
    auto* forcePtr = mesh_.getObjectPtr<volVectorField>(scopedName("force"));

    if (!forcePtr)
    {
        forcePtr = new volVectorField
        (
            IOobject
            (
                scopedName("force"),
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimForce, Zero)
        );

        mesh_.objectRegistry::store(forcePtr);
    }

    return *forcePtr;
}


Foam::volVectorField& Foam::functionObjects::porousforces::moment() 
{
    auto* momentPtr = mesh_.getObjectPtr<volVectorField>(scopedName("moment"));

    if (!momentPtr)
    {
        momentPtr = new volVectorField
        (
            IOobject
            (
                scopedName("moment"),
                time_.timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedVector(dimForce*dimLength, Zero)
        );

        mesh_.objectRegistry::store(momentPtr);
    }

    return *momentPtr;
}


void Foam::functionObjects::porousforces::initialise() 
{
    if (initialised_)
    {
        return;
    }

    if (directForceDensity_)
    {
        if (!foundObject<volVectorField>(fDName_))
        {
            FatalErrorInFunction
                << "Could not find " << fDName_ << " in database"
                << exit(FatalError);
        }
    }
    else
    {
        if
        (
            !foundObject<volVectorField>(UName_)
         || !foundObject<volScalarField>(pName_)
        )
        {
            FatalErrorInFunction
                << "Could not find U: " << UName_
                << " or p:" << pName_ << " in database"
                << exit(FatalError);
        }

        if (rhoName_ != "rhoInf" && !foundObject<volScalarField>(rhoName_))
        {
            FatalErrorInFunction
                << "Could not find rho:" << rhoName_ << " in database"
                << exit(FatalError);
        }
    }

    initialised_ = true;
}


void Foam::functionObjects::porousforces::reset()  
{
    sumPatchForcesP_ = Zero;
    sumPatchForcesV_ = Zero;
    sumPatchMomentsP_ = Zero;
    sumPatchMomentsV_ = Zero;


    sumInternalForces_ = Zero;
    sumInternalMoments_ = Zero;

    auto& force = this->force();
    auto& moment = this->moment();
    force == dimensionedVector(force.dimensions(), Zero);
    moment == dimensionedVector(moment.dimensions(), Zero);
}


Foam::tmp<Foam::volSymmTensorField>
Foam::functionObjects::porousforces::devRhoReff() const
{
    typedef compressible::turbulenceModel cmpTurbModel;
    typedef incompressible::turbulenceModel icoTurbModel;

    if (foundObject<cmpTurbModel>(cmpTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<cmpTurbModel>(cmpTurbModel::propertiesName);

        return turb.devRhoReff();
    }
    else if (foundObject<icoTurbModel>(icoTurbModel::propertiesName))
    {
        const auto& turb =
            lookupObject<icoTurbModel>(icoTurbModel::propertiesName);

        return rho()*turb.devReff();
    }
    else if (foundObject<fluidThermo>(fluidThermo::dictName))
    {
        const auto& thermo = lookupObject<fluidThermo>(fluidThermo::dictName);

        const auto& U = lookupObject<volVectorField>(UName_);

        return -thermo.mu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        const auto& U = lookupObject<volVectorField>(UName_);

        return -rho()*laminarT.nu()*dev(twoSymm(fvc::grad(U)));
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
            lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        const auto& U = lookupObject<volVectorField>(UName_);

        return -rho()*nu*dev(twoSymm(fvc::grad(U)));
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for viscous stress calculation"
            << exit(FatalError);

        return volSymmTensorField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::porousforces::mu() const
{
    if (foundObject<fluidThermo>(basicThermo::dictName))
    {
        const auto& thermo = lookupObject<fluidThermo>(basicThermo::dictName);

        return thermo.mu();
    }
    else if (foundObject<transportModel>("transportProperties"))
    {
        const auto& laminarT =
            lookupObject<transportModel>("transportProperties");

        return rho()*laminarT.nu();
    }
    else if (foundObject<dictionary>("transportProperties"))
    {
        const auto& transportProperties =
             lookupObject<dictionary>("transportProperties");

        const dimensionedScalar nu("nu", dimViscosity, transportProperties);

        return rho()*nu;
    }
    else
    {
        FatalErrorInFunction
            << "No valid model for dynamic viscosity calculation"
            << exit(FatalError);

        return volScalarField::null();
    }
}


Foam::tmp<Foam::volScalarField> Foam::functionObjects::porousforces::rho() const 
{
    if (rhoName_ == "rhoInf")
    {
        return tmp<volScalarField>::New
        (
            IOobject
            (
                "rho",
                mesh_.time().timeName(),
                mesh_
            ),
            mesh_,
            dimensionedScalar(dimDensity, rhoRef_)
        );
    }

    return (lookupObject<volScalarField>(rhoName_));
}


Foam::scalar Foam::functionObjects::porousforces::rho(const volScalarField& p) const
{
    if (p.dimensions() == dimPressure)
    {
        return 1;
    }

    if (rhoName_ != "rhoInf")
    {
        FatalErrorInFunction
            << "Dynamic pressure is expected but kinematic is provided."
            << exit(FatalError);
    }

    return rhoRef_;
}


void Foam::functionObjects::porousforces::addToPatchFields 
(
    const label patchi,
    const vectorField& Md,
    const vectorField& fP,
    const vectorField& fV
)
{
    sumPatchForcesP_ += sum(fP);
    sumPatchForcesV_ += sum(fV);
    force().boundaryFieldRef()[patchi] += fP + fV;

    const vectorField mP(Md^fP);
    const vectorField mV(Md^fV);

    sumPatchMomentsP_ += sum(mP);
    sumPatchMomentsV_ += sum(mV);
    moment().boundaryFieldRef()[patchi] += mP + mV;
}


void Foam::functionObjects::porousforces::addToInternalField  
(
    const labelList& cellIDs,
    const vectorField& Md,
    const vectorField& f
)
{
    auto& force = this->force();
    auto& moment = this->moment();

    forAll(cellIDs, i)
    {
        const label celli = cellIDs[i];

        sumInternalForces_ += f[i];
        force[celli] += f[i];

        const vector m(Md[i]^f[i]);
        sumInternalMoments_ += m;
        moment[celli] = m;
    }
}


void Foam::functionObjects::porousforces::createIntegratedDataFiles() 
{
    if (!forceFilePtr_.valid())
    {
        forceFilePtr_ = createFile("force");
        writeIntegratedDataFileHeader("Force", forceFilePtr_());
    }

    if (!momentFilePtr_.valid())
    {
        momentFilePtr_ = createFile("moment");
        writeIntegratedDataFileHeader("Moment", momentFilePtr_());
    }   
}


void Foam::functionObjects::porousforces::writeIntegratedDataFileHeader  
(
    const word& header,
    OFstream& os
) const
{
    const auto& coordSys = coordSysPtr_();
    const auto vecDesc = [](const word& root)->string
    {
        return root + "_x " + root + "_y " + root + "_z";
    };
    writeHeader(os, header);
    writeHeaderValue(os, "CofR", coordSys.origin());
    writeHeader(os, "");
    writeCommented(os, "Time");
    writeTabbed(os, vecDesc("Total"));
    writeTabbed(os, vecDesc("Pressureforces"));
    writeTabbed(os, vecDesc("Gravity"));
    writeTabbed(os, vecDesc("Porousforces"));

    os  << endl;
}


void Foam::functionObjects::porousforces::writeIntegratedDataFiles() 
{
    const auto& coordSys = coordSysPtr_();

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchForcesP_),
        coordSys.localVector(sumPatchForcesV_),
        coordSys.localVector(sumInternalForces_),
        forceFilePtr_()
    );

    writeIntegratedDataFile
    (
        coordSys.localVector(sumPatchMomentsP_),
        coordSys.localVector(sumPatchMomentsV_),
        coordSys.localVector(sumInternalMoments_),
        momentFilePtr_()
    );
}


void Foam::functionObjects::porousforces::writeIntegratedDataFile  
(
    const vector& pres,
    const vector& vis,
    const vector& internal,
    OFstream& os
) const
{
    writeCurrentTime(os);

    writeValue(os, pres + vis + internal);
    writeValue(os, pres);
    writeValue(os, vis);
    writeValue(os, internal);


    os  << endl;
}


void Foam::functionObjects::porousforces::logIntegratedData  
(
    const string& descriptor,
    const vector& pres,
    const vector& vis,
    const vector& internal
) const
{
    if (!log)
    {
        return;
    }

    Log << "    Sum of " << descriptor.c_str() << nl
        << "        Total    : " << (pres + vis + internal) << nl
        << "        Pressureforces : " << pres << nl
        << "        Gravity  : " << vis << nl
        << "        Porousforces   : " << internal << nl;

}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::functionObjects::porousforces::porousforces
(
    const word& name,
    const Time& runTime,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, runTime, dict),
    writeFile(mesh_, name),
    sumPatchForcesP_(Zero),
    sumPatchForcesV_(Zero),
    sumPatchMomentsP_(Zero),
    sumPatchMomentsV_(Zero),
    sumInternalForces_(Zero),
    sumInternalMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSysPtr_(nullptr),
    patchSet_(),
    rhoRef_(VGREAT),
    pRef_(0),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    
    //add  
    BuoyancyFilePtr_(),
    GravityFilePtr_(),
    DragforceFilePtr_(),
    alphaName_("alpha.water"),                             
    //add
    
    fDName_("fD"),
    directForceDensity_(false),
    //porosity(false),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


Foam::functionObjects::porousforces::porousforces
(
    const word& name,
    const objectRegistry& obr,
    const dictionary& dict,
    bool readFields
)
:
    fvMeshFunctionObject(name, obr, dict),
    writeFile(mesh_, name),
    sumPatchForcesP_(Zero),
    sumPatchForcesV_(Zero),
    sumPatchMomentsP_(Zero),
    sumPatchMomentsV_(Zero),
    sumInternalForces_(Zero),
    sumInternalMoments_(Zero),
    forceFilePtr_(),
    momentFilePtr_(),
    coordSysPtr_(nullptr),
    patchSet_(),
    rhoRef_(VGREAT),
    pRef_(0),
    pName_("p"),
    UName_("U"),
    rhoName_("rho"),
    
    //add  
    BuoyancyFilePtr_(),                                      
    GravityFilePtr_(),
    DragforceFilePtr_(),                                     
    alphaName_("alpha.water"),                          
    //add
    
    fDName_("fD"),
    directForceDensity_(false),
    //porosity_(false),
    writeFields_(false),
    initialised_(false)
{
    if (readFields)
    {
        read(dict);
        setCoordinateSystem(dict);
        Log << endl;
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::porousforces::read(const dictionary& dict)   
{
    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict))
    {
        return false;
    }

    initialised_ = false;

    Info<< type() << " " << name() << ":" << endl;

    patchSet_ =
        mesh_.boundaryMesh().patchSet
        (
            dict.get<wordRes>("patches")
        );

    dict.readIfPresent("directForceDensity", directForceDensity_);
    if (directForceDensity_)
    {
        // Optional name entry for fD
        if (dict.readIfPresent<word>("fD", fDName_))
        {
            Info<< "    fD: " << fDName_ << endl;
        }
    }
    else
    {
        // Optional field name entries
        if (dict.readIfPresent<word>("p", pName_))     
        {
            Info<< "    p: " << pName_ << endl;
        }
        if (dict.readIfPresent<word>("U", UName_))
        {
            Info<< "    U: " << UName_ << endl;
        }
        if (dict.readIfPresent<word>("rho", rhoName_))
        {
            Info<< "    rho: " << rhoName_ << endl;
        }
        
        
        //add
        if (dict.readIfPresent<word>("alpha.water", alphaName_))  
        {                                                         
            Info<< "    alpha: " << alphaName_ << endl;           
        }
        if (dict.readIfPresent<scalar>("myrhoporo", myrhoporo_))  
        {                                                         
            Info<< "    myrhoporo = " << myrhoporo_ << endl;           
        }
        if (dict.readIfPresent("Linear_vel", Linear_vel_))  
        {                                                         
            Info<< "    Linear_vel = " << Linear_vel_ << endl;           
        }
        if (dict.readIfPresent("Angular_vel", Angular_vel_))  
        {                                                         
            Info<< "    Angular_vel = " << Angular_vel_ << endl;           
        }
                                                         
        //add
        
    
        // Reference density needed for incompressible calculations
        if (rhoName_ == "rhoInf")
        {
            rhoRef_ = dict.getCheck<scalar>("rhoInf", scalarMinMax::ge(SMALL));
            Info<< "    Freestream density (rhoInf) set to " << rhoRef_ << endl;
        }

        // Reference pressure, 0 by default
        if (dict.readIfPresent<scalar>("pRef", pRef_))
        {
            Info<< "    Reference pressure (pRef) set to " << pRef_ << endl;
        }
    }
/*
    dict.readIfPresent("porosity", porosity_);
    if (porosity_)
    {
        Info<< "    Including porosity effects" << endl;
    }
    else
    {
        Info<< "    Not including porosity effects" << endl;
    }
*/
    writeFields_ = dict.getOrDefault("writeFields", false);
    if (writeFields_)
    {
        Info<< "    Fields will be written" << endl;
    }


    return true;
}


void Foam::functionObjects::porousforces::calcForcesMoments()  
{
    initialise();

    reset();

    const point& origin = coordSysPtr_->origin();

	//Rigid body
	/*
    if (directForceDensity_)
    {
        const auto& fD = lookupObject<volVectorField>(fDName_);

        const auto& Sfb = mesh_.Sf().boundaryField();

        for (const label patchi : patchSet_)
        {
            const vectorField& d = mesh_.C().boundaryField()[patchi];

            const vectorField Md(d - origin);

            const scalarField sA(mag(Sfb[patchi]));

            // Pressure force = surfaceUnitNormal*(surfaceNormal & forceDensity)
            const vectorField fP
            (
                Sfb[patchi]/sA
               *(
                    Sfb[patchi] & fD.boundaryField()[patchi]
                )
            );

            // Viscous force (total force minus pressure fP)
            const vectorField fV(sA*fD.boundaryField()[patchi] - fP);

            addToPatchFields(patchi, Md, fP, fV);
        }
    }
    else
    {

    	{
	    const auto& p = lookupObject<volScalarField>(pName_);

		const auto& Sfb = mesh_.Sf().boundaryField();

		tmp<volSymmTensorField> tdevRhoReff = devRhoReff();
		const auto& devRhoReffb = tdevRhoReff().boundaryField();

		// Scale pRef by density for incompressible simulations
		const scalar rhoRef = rho(p);
		const scalar pRef = pRef_/rhoRef;

		for (const label patchi : patchSet_)
		{
		    const vectorField& d = mesh_.C().boundaryField()[patchi];

		    const vectorField Md(d - origin);

		    //pressure force
		    const vectorField fP
		    (
		        rhoRef*Sfb[patchi]*(p.boundaryField()[patchi] - pRef)
		    );

		    //viscous force
		    const vectorField fV(Sfb[patchi] & devRhoReffb[patchi]);

		    addToPatchFields(patchi, Md, fP, fV);
		}
    	}
    		
	            
    }

    if (porosity_)
    {
        const auto& U = lookupObject<volVectorField>(UName_);
        const volScalarField rho(this->rho());
        const volScalarField mu(this->mu());

        const auto models = obr_.lookupClass<porosityModel>();

        if (models.empty())
        {
            WarningInFunction
                << "Porosity effects requested, "
                << "but no porosity models found in the database"
                << endl;
        }

        forAllConstIters(models, iter)
        {
            // Non-const access required if mesh is changing
            auto& pm = const_cast<porosityModel&>(*iter());

            const vectorField fPTot(pm.force(U, rho, mu));

            const labelList& cellZoneIDs = pm.cellZoneIDs();

            for (const label zonei : cellZoneIDs)
            {
                const cellZone& cZone = mesh_.cellZones()[zonei];

                const vectorField d(mesh_.C(), cZone);
                const vectorField fP(fPTot, cZone);
                const vectorField Md(d - origin);

                addToInternalField(cZone, Md, fP);
            }
        }
    }
 	        
	*/

	
	//Porous structure
    //read from each cells
    const auto& alpha = lookupObject<volScalarField>(alphaName_);
    const auto& U = lookupObject<volVectorField>(UName_);     
    //const volScalarField rho(this->rho());   
    volVectorField acc = fvc::ddt(U);   // a = du/dt        
    const volScalarField mu(this->mu());  
    const auto& porosityC = lookupObject<volScalarField>("porosity");  
    const auto& D50C = lookupObject<volScalarField>("D50Field");
    const auto& alphaAC = lookupObject<volScalarField>("aPorField");
    const auto& betaC = lookupObject<volScalarField>("bPorField"); 
    const auto& cC = lookupObject<volScalarField>("cPorField"); 
    const auto& rhoC = lookupObject<volScalarField>("rho"); 
    const auto& KCC = lookupObject<volScalarField>("KCPorField");
    const scalar dens = myrhoporo_;	 //Gravity input: porous structure's density
    const auto& P = lookupObject<volScalarField>(pName_); 
    const volVectorField gradP = Foam::fvc::grad(P);               
        
    //ZoneID
	const label ZoneID = mesh_.cellZones().findZoneID("porous");
	const labelList& cells = mesh_.cellZones()[ZoneID];

    //initialize
	//scalar sumBuoyancy(0);
	scalar sumGravity(0);
	vector sumPressureforce;   sumPressureforce = Zero;
	vector sumDragforce;	sumDragforce = Zero;    
	vector sumDragforceA;	sumDragforceA = Zero;	
	vector sumDragforceB;	sumDragforceB = Zero;	
	vector sumDragforceC;	sumDragforceC = Zero;	
	volVectorField UR = U*0;	//ini relative U
	volScalarField F_U = P*0;	//F*U
	volVectorField F_Drag = U*0;	//Drag force
	volVectorField F_Pressure = U*0;	//Pressure force	   	      
    //Md
	const vectorField d(mesh_.C(), cells); //every cell's center
	const vectorField Md(d - origin);      //vector point->origin
	vector m;

		
    //Cells
	forAll(cells, i)
	{   
		//initialize
		const scalar poro = porosityC[cells[i]]; 
		const scalar D50 = D50C[cells[i]]; 
		const scalar alphaA = alphaAC[cells[i]]; 
		const scalar beta = betaC[cells[i]]; 
		const scalar c = cC[cells[i]]; 
		const scalar rou = rhoC[cells[i]];
		const scalar KC = KCC[cells[i]];  	       
	       
		//UR
		UR[cells[i]][0] = U[cells[i]][0] - Linear_vel_[0] - Angular_vel_[1] * Md[i][2] + Angular_vel_[2] * Md[i][1] ;
		UR[cells[i]][1] = U[cells[i]][1] - Linear_vel_[1] - Angular_vel_[2] * Md[i][0] + Angular_vel_[0] * Md[i][2] ;
		UR[cells[i]][2] = U[cells[i]][2] - Linear_vel_[2] - Angular_vel_[0] * Md[i][1] + Angular_vel_[1] * Md[i][0] ;	 	
	
		//Drag force
			//A
			vector Dragforce_A = 
				alpha[cells[i]]* alphaA * (1-poro)*(1-poro)/ poro/poro * 
				mu[cells[i]] / D50/D50 * UR[cells[i]]  * mesh_.V()[cells[i]];			
			//B
				scalar Dragforce_B1 =
					beta * (1+0 * 7.5/KC) * (1-poro) /poro/poro *
					rou / D50 * mesh_.V()[cells[i]]; 	
				vector Dragforce_B2 =
					UR[cells[i]] *  mag(UR[cells[i]]);    //U*|U| 
			vector Dragforce_B =  alpha[cells[i]] * Dragforce_B1 * Dragforce_B2;  							
	
			//C
			vector Dragforce_C =
				alpha[cells[i]] * c * acc[cells[i]] * mesh_.V()[cells[i]]/poro* rou;	//du/dt	 		
		
			vector Dragforce = Dragforce_A + Dragforce_B + Dragforce_C;
			sumDragforce = sumDragforce + Dragforce;
			sumDragforceA = sumDragforceA + Dragforce_A;
			sumDragforceB = sumDragforceB + Dragforce_B;
			sumDragforceC = sumDragforceC + Dragforce_C;
		
		//Pressure force
		    vector Pressureforce =  - (1 - poro) * gradP[cells[i]] * mesh_.V()[cells[i]] * alpha[cells[i]];
			sumPressureforce = sumPressureforce + Pressureforce;				
	        	       
	    //Buoyancy        
			//scalar Buoyancy = alpha[cells[i]]*rou*mesh_.V()[cells[i]]*(1-poro)*9.8;
			//sumBuoyancy += Buoyancy;

		
		//Gravity
			scalar Gravity = - dens * mesh_.V()[cells[i]] * (1-poro) * 9.8;
			sumGravity += Gravity;				
		
		//sumInternalForces
			sumInternalForces_ = sumDragforce;		    		    
			sumPatchForcesV_[2] =  sumGravity;
			sumPatchForcesP_ = sumPressureforce;

		//sumInternalmoment
			vector addforce = Dragforce_A + Dragforce_B + Dragforce_C - (1-poro) * gradP[cells[i]] * mesh_.V()[cells[i]];
			addforce[2] = addforce[2]  + Gravity;
				m[1] = -Md[i][0] * addforce[2];
				m[2] = -Md[i][1] * addforce[0];
				m[0] = -Md[i][2] * addforce[1];
			sumInternalMoments_ += m;
		
		
		//FU - energy dissipation
			F_U[cells[i]] =  U[cells[i]][0] * (Dragforce[0] ) 
					+ U[cells[i]][1] * (Dragforce[1] )
					+ U[cells[i]][2] * (Dragforce[2] );					
					
		//F Drag
			F_Drag[cells[i]] = Dragforce;
			
		//F Pressure
			F_Pressure[cells[i]] = Pressureforce;
			
	}
        
	/*
	Info<< "sumDragforceA = " << sumDragforceA << nl << endl;
	Info<< "sumDragforceB = " << sumDragforceB << nl << endl;
	Info<< "sumDragforceC = " << sumDragforceC << nl << endl;
	*/
                              
    /*  
    //write
    if (!BuoyancyFilePtr_.valid())
    {
        BuoyancyFilePtr_ = createFile("sumDragforceA");
    }
    writeCurrentTime(BuoyancyFilePtr_());       
    writeValue(BuoyancyFilePtr_(),  sumDragforceA);
    BuoyancyFilePtr_()  << endl;
        
    if (!GravityFilePtr_.valid())
    {
        GravityFilePtr_ = createFile("sumDragforceB");
    }
    writeCurrentTime(GravityFilePtr_());    
    writeValue(GravityFilePtr_(),  sumDragforceB);
    GravityFilePtr_()  << endl;
        
    if (!DragforceFilePtr_.valid())
    {
        DragforceFilePtr_ = createFile("sumDragforceC");
    }
    writeCurrentTime(DragforceFilePtr_());
    writeValue(DragforceFilePtr_(),  sumDragforceC);
    DragforceFilePtr_()  << endl;
	*/        

	reduce(sumPatchForcesP_, sumOp<vector>());
	reduce(sumPatchForcesV_, sumOp<vector>());
	reduce(sumPatchMomentsP_, sumOp<vector>());
	reduce(sumPatchMomentsV_, sumOp<vector>());
	reduce(sumInternalForces_, sumOp<vector>());
	reduce(sumInternalMoments_, sumOp<vector>());
    
  	//Output
  	Info<< "Dragforce = " << sumInternalForces_ << nl << endl;
	Info<< "Pressureforce = " << sumPatchForcesP_ << nl << endl;	    
	Info<< "Gravity = " << sumPatchForcesV_ << nl << endl;

    
	//volVectorField URR = lookupObject<volVectorField>("UR");  
	//UR.write("UR");
	
	//Output —— FU  Drag  Pressure  
	if( mesh_.time().writeTime() ) 
        {
            Info<< "Output —— FU  FDrag  Fpressure  UR:" << nl << endl;
	    volScalarField FU
	    		 (IOobject
	    		 	(  
	    		 	"FU",
	    		 	mesh_.time().timeName(),
	    		 	mesh_,
	    		 	IOobject::NO_READ,
	    		 	IOobject::NO_WRITE
	    		 	),
	    		 mesh_,
	    		 dimensionedScalar("zero",dimensionSet(1,-1,-2,0,0,0,0),Foam::scalar(0)),
	    		 P.boundaryField().types()
	    		 );
	    FU = F_U;
	    FU.write();	 
	    
	    volVectorField FDrag
	    		 (IOobject
	    		 	(  
	    		 	"FDrag",
	    		 	mesh_.time().timeName(),
	    		 	mesh_,
	    		 	IOobject::NO_READ,
	    		 	IOobject::NO_WRITE
	    		 	),
	    		 mesh_,
	    		 dimensionedVector("zero",dimensionSet(0,1,-1,0,0,0,0),Foam::vector(0,0,0)),
	    		 U.boundaryField().types()
	    		 );
	    FDrag = F_Drag;
	    FDrag.write();	 
	    
	    volVectorField FPressure
	    		 (IOobject
	    		 	(  
	    		 	"FPressure",
	    		 	mesh_.time().timeName(),
	    		 	mesh_,
	    		 	IOobject::NO_READ,
	    		 	IOobject::NO_WRITE
	    		 	),
	    		 mesh_,
	    		 dimensionedVector("zero",dimensionSet(0,1,-1,0,0,0,0),Foam::vector(0,0,0)),
	    		 U.boundaryField().types()
	    		 );
	    FPressure = F_Pressure;
	    FPressure.write();	
	    
	    volVectorField U_R
	    		 (IOobject
	    		 	(  
	    		 	"U_R",
	    		 	mesh_.time().timeName(),
	    		 	mesh_,
	    		 	IOobject::NO_READ,
	    		 	IOobject::NO_WRITE
	    		 	),
	    		 mesh_,
	    		 dimensionedVector("zero",dimensionSet(0,1,-1,0,0,0,0),Foam::vector(0,0,0)),
	    		 U.boundaryField().types()
	    		 );
	    U_R = UR;
	    U_R.write();		 
	} 	   
	
	//transport UR    
    const Time& t = mesh_.time();
    if (t.foundObject<IOdictionary>("motion"))
    {
        const dictionary& motionDict =
        	     t.lookupObject<IOdictionary>("motion");
        	    
		//update for the motion solver
		dictionary updateDbDictionary = motionDict;
		//updateDbDictionary.set("UR", UR);
		updateDbDictionary.set("Md", Md);
		updateDbDictionary.set("Linear_vel_", Linear_vel_);
		updateDbDictionary.set("Angular_vel_", Angular_vel_);
		const_cast <dictionary&>(motionDict) = updateDbDictionary;
    }

} 
    
	/*
	//add - return Ur
	Foam::volVectorField  Foam::functionObjects::porousforces::UR_() const
	{
    	volVectorField U = lookupObject<volVectorField>(UName_);
    	return U;
	}*/


Foam::vector Foam::functionObjects::porousforces::forceEff() const
{
    return sumPatchForcesP_ + sumPatchForcesV_ + sumInternalForces_;
}


Foam::vector Foam::functionObjects::porousforces::momentEff() const
{
    return sumPatchMomentsP_ + sumPatchMomentsV_ + sumInternalMoments_;
}


bool Foam::functionObjects::porousforces::execute()
{
    calcForcesMoments();

    Log << type() << " " << name() << " write:" << nl;

    const auto& coordSys = coordSysPtr_();

    const auto localFp(coordSys.localVector(sumPatchForcesP_));
    const auto localFv(coordSys.localVector(sumPatchForcesV_));
    const auto localFi(coordSys.localVector(sumInternalForces_));

    logIntegratedData("forces", localFp, localFv, localFi);

    const auto localMp(coordSys.localVector(sumPatchMomentsP_));
    const auto localMv(coordSys.localVector(sumPatchMomentsV_));
    const auto localMi(coordSys.localVector(sumInternalMoments_));

    logIntegratedData("moments", localMp, localMv, localMi);

    setResult("pressureForce", localFp);
    setResult("viscousForce", localFv);
    setResult("internalForce", localFi);
    setResult("pressureMoment", localMp);
    setResult("viscousMoment", localMv);
    setResult("internalMoment", localMi);

    return true;
}


bool Foam::functionObjects::porousforces::write()
{
    if (writeToFile())
    {
        Log << "    writing force and moment files." << endl;

        createIntegratedDataFiles();
        writeIntegratedDataFiles();
    }

    if (writeFields_)
    {
        Log << "    writing force and moment fields." << endl;

        force().write();
        moment().write();
    }

    Log << endl;

    return true;
}


// ************************************************************************* //
