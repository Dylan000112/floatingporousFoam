/*---------------------------------------------------------------------------*\
  Porous Forces Function Object - Source (Final Fixed v2)
\*---------------------------------------------------------------------------*/

#include "porousforces.H"
#include "fvcGrad.H"
#include "fvcDdt.H"
#include "coordinateSystem.H"
#include "addToRunTimeSelectionTable.H"
#include "fvCFD.H"
#include "turbulentTransportModel.H"
#include "turbulentFluidThermoModel.H"

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(porousforces, 0);
    addToRunTimeSelectionTable(functionObject, porousforces, dictionary);
}
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
    sumDragForce_(Zero), sumDragForceA_(Zero), sumDragForceB_(Zero), sumDragForceC_(Zero),
    sumPressureForce_(Zero), sumGravityForce_(Zero), sumTotalMoment_(Zero),
    coordSysPtr_(nullptr),
    myrhoporo_(0), Linear_vel_(Zero), Angular_vel_(Zero),
    pName_("p"), UName_("U"), rhoName_("rho"), alphaName_("alpha.water"),
    writeFields_(false), initialised_(false)
{ if (readFields) read(dict); }

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
    sumDragForce_(Zero), sumDragForceA_(Zero), sumDragForceB_(Zero), sumDragForceC_(Zero),
    sumPressureForce_(Zero), sumGravityForce_(Zero), sumTotalMoment_(Zero),
    coordSysPtr_(nullptr),
    myrhoporo_(0), Linear_vel_(Zero), Angular_vel_(Zero),
    pName_("p"), UName_("U"), rhoName_("rho"), alphaName_("alpha.water"),
    writeFields_(false), initialised_(false)
{ if (readFields) read(dict); }

// --- 持久化场管理实现 ---
template<class GeoField>
GeoField& Foam::functionObjects::porousforces::getOrMakeField
(
    const word& name, 
    const dimensionSet& dims
)
{
    auto* ptr = mesh_.getObjectPtr<GeoField>(name);
    if (!ptr)
    {
        ptr = new GeoField
        (
            IOobject(name, mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::AUTO_WRITE),
            mesh_,
            dimensioned<typename GeoField::value_type>(dims, Zero)
        );
        mesh_.objectRegistry::store(ptr);
    }
    return *ptr;
}

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool Foam::functionObjects::porousforces::read(const dictionary& dict)
{
    if (!fvMeshFunctionObject::read(dict) || !writeFile::read(dict)) return false;

    initialised_ = false;

    dict.readIfPresent("p", pName_);
    dict.readIfPresent("U", UName_);
    dict.readIfPresent("rho", rhoName_);
    dict.readIfPresent("alpha.water", alphaName_);

    myrhoporo_   = dict.get<scalar>("myrhoporo");
    Linear_vel_  = dict.get<vector>("Linear_vel");
    Angular_vel_ = dict.get<vector>("Angular_vel");

    if (dict.found("coordinateSystem"))
        coordSysPtr_ = coordinateSystem::New(obr_, dict, "coordinateSystem");
    else
    {
        point origin = dict.getOrDefault<point>("CofR", dict.getOrDefault<point>("origin", Zero));
        coordSysPtr_.reset(new coordinateSystem("coord", origin, vector(1,0,0), vector(0,0,1)));
    }

    writeFields_ = dict.getOrDefault("writeFields", false);
    return true;
}

void Foam::functionObjects::porousforces::resetSum()
{
    sumDragForce_ = sumDragForceA_ = sumDragForceB_ = sumDragForceC_ = Zero;
    sumPressureForce_ = sumGravityForce_ = sumTotalMoment_ = Zero;
}

void Foam::functionObjects::porousforces::calcForcesMoments()
{
    resetSum();

    const auto& p        = lookupObject<volScalarField>(pName_);
    const auto& U        = lookupObject<volVectorField>(UName_);
    const auto& alpha    = lookupObject<volScalarField>(alphaName_);
    const auto& rhoCell  = lookupObject<volScalarField>(rhoName_);
    const auto& porosity = lookupObject<volScalarField>("porosity");

    tmp<volScalarField> tmu;
    if (foundObject<fluidThermo>(basicThermo::dictName))
        tmu = lookupObject<fluidThermo>(basicThermo::dictName).mu();
    else if (foundObject<transportModel>("transportProperties"))
        tmu = rhoCell * lookupObject<transportModel>("transportProperties").nu();
    else
        tmu = lookupObject<volScalarField>("mu");
    const auto& mu = tmu();

    const volVectorField gradP = fvc::grad(p);
    const volVectorField acc   = fvc::ddt(U);
    
    const auto& D50Field   = lookupObject<volScalarField>("D50Field");
    const auto& aPorField  = lookupObject<volScalarField>("aPorField");
    const auto& bPorField  = lookupObject<volScalarField>("bPorField");
    const auto& cPorField  = lookupObject<volScalarField>("cPorField");
    const auto& KCPorField = lookupObject<volScalarField>("KCPorField");

    const label zoneID = mesh_.cellZones().findZoneID("porous");
    const point& origin = coordSysPtr_->origin();

    // 持久化场初始化
    auto& FU = getOrMakeField<volScalarField>("FU", dimensionSet(1, -1, -3, 0, 0, 0, 0));
    auto& FDrag = getOrMakeField<volVectorField>("FDrag", dimensionSet(1, -2, -2, 0, 0, 0, 0));
    auto& FPressure = getOrMakeField<volVectorField>("FPressure", dimensionSet(1, -2, -2, 0, 0, 0, 0));
    auto& URField = getOrMakeField<volVectorField>("U_R", dimensionSet(0, 1, -1, 0, 0, 0, 0));

    // 使用显式赋值重置场
    FU = dimensionedScalar("0", FU.dimensions(), 0.0);
    FDrag = dimensionedVector("0", FDrag.dimensions(), Zero);
    FPressure = dimensionedVector("0", FPressure.dimensions(), Zero);
    URField = dimensionedVector("0", URField.dimensions(), Zero);

    // --- 修复分支类型冲突 Bug ---
    if (zoneID != -1)
    {
        const labelList& cells = mesh_.cellZones()[zoneID];
        vectorField MdField(cells.size());

        forAll(cells, i)
        {
            const label celli = cells[i];
            const scalar V = mesh_.V()[celli];
            const scalar poro = porosity[celli];
            const vector Md = mesh_.C()[celli] - origin;
            MdField[i] = Md;

            const vector UR = U[celli] - (Linear_vel_ + (Angular_vel_ ^ Md));
            URField[celli] = UR;

            vector fA = alpha[celli] * aPorField[celli] * sqr(1.0-poro)/pow3(poro) * mu[celli] / sqr(D50Field[celli]) * UR * V;
            vector fB = alpha[celli] * (bPorField[celli] * (1.0 + 7.5/KCPorField[celli]) * (1.0-poro)/pow3(poro) * rhoCell[celli] / D50Field[celli] * V) * UR * mag(UR);
            vector fC = alpha[celli] * cPorField[celli] * acc[celli] * V / poro * rhoCell[celli];

            const vector cellDrag = fA + fB + fC;
            const vector cellPressure = -(1.0 - poro) * gradP[celli] * V * alpha[celli];
            const vector cellGravity(0, 0, -myrhoporo_ * V * (1.0 - poro) * 9.81);

            sumDragForceA_ += fA; sumDragForceB_ += fB; sumDragForceC_ += fC;
            sumDragForce_  += cellDrag;
            sumPressureForce_ += cellPressure;
            sumGravityForce_  += cellGravity;
            sumTotalMoment_   += (Md ^ (cellDrag + cellPressure + cellGravity));

            FU[celli] = (U[celli] & cellDrag) / V;
            FDrag[celli] = cellDrag / V;
            FPressure[celli] = cellPressure / V;
        }

        // 运动字典交互
        if (mesh_.time().foundObject<IOdictionary>("motion"))
        {
            IOdictionary& mDict = const_cast<IOdictionary&>(mesh_.time().lookupObject<IOdictionary>("motion"));
            mDict.set("Md", MdField);
            mDict.set("sumDragForce", vector(sumDragForce_));
            mDict.set("sumTotalMoment", vector(sumTotalMoment_));
            mDict.set("Linear_vel_", vector(Linear_vel_));
            mDict.set("Angular_vel_", vector(Angular_vel_));
        }
    }

    // 并行规约
    reduce(sumDragForce_, sumOp<vector>());
    reduce(sumDragForceA_, sumOp<vector>());
    reduce(sumDragForceB_, sumOp<vector>());
    reduce(sumDragForceC_, sumOp<vector>());
    reduce(sumPressureForce_, sumOp<vector>());
    reduce(sumGravityForce_, sumOp<vector>());
    reduce(sumTotalMoment_, sumOp<vector>());
}

Foam::vector Foam::functionObjects::porousforces::forceEff() const
{ return sumDragForce_ + sumPressureForce_ + sumGravityForce_; }

Foam::vector Foam::functionObjects::porousforces::momentEff() const
{ return sumTotalMoment_; }

bool Foam::functionObjects::porousforces::execute()
{
    calcForcesMoments();
    const auto& coordSys = *coordSysPtr_;
    setResult("dragForce",     coordSys.localVector(sumDragForce_));
    setResult("pressureForce", coordSys.localVector(sumPressureForce_));
    setResult("gravityForce",  coordSys.localVector(sumGravityForce_));
    setResult("totalMoment",   coordSys.localVector(sumTotalMoment_));
    return true;
}

bool Foam::functionObjects::porousforces::write()
{
    if (writeToFile())
    {
        if (!forceFilePtr_.valid())
        {
            forceFilePtr_ = createFile("forces");
            writeHeader(forceFilePtr_(), "Porous Forces (Total, DragA, DragB, DragC, Pressure, Gravity)");
        }
        writeCurrentTime(forceFilePtr_());
        writeValue(forceFilePtr_(), sumDragForce_);
        writeValue(forceFilePtr_(), sumDragForceA_);
        writeValue(forceFilePtr_(), sumDragForceB_);
        writeValue(forceFilePtr_(), sumDragForceC_);
        writeValue(forceFilePtr_(), sumPressureForce_);
        writeValue(forceFilePtr_(), sumGravityForce_);
        forceFilePtr_() << endl;

        if (!momentFilePtr_.valid())
        {
            momentFilePtr_ = createFile("moments");
            writeHeader(momentFilePtr_(), "Total Moment");
        }
        writeCurrentTime(momentFilePtr_());
        writeValue(momentFilePtr_(), sumTotalMoment_);
        momentFilePtr_() << endl;
    }

    if (mesh_.time().writeTime())
    {
        mesh_.lookupObject<volScalarField>("FU").write();
        mesh_.lookupObject<volVectorField>("FDrag").write();
        mesh_.lookupObject<volVectorField>("FPressure").write();
        mesh_.lookupObject<volVectorField>("U_R").write();
    }

    return true;
}

void Foam::functionObjects::porousforces::writeHeader(OFstream& os, const word& title) const
{
    os  << "# " << title << endl;
    os  << "# Origin: " << coordSysPtr_->origin() << endl;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

