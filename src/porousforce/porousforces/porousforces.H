
/*---------------------------------------------------------------------------*\
  Porous Forces Function Object - Header (Final Fixed v2)
\*---------------------------------------------------------------------------*/

#ifndef Foam_functionObjects_porousforces_H
#define Foam_functionObjects_porousforces_H

#include "fvMeshFunctionObject.H"
#include "writeFile.H"
#include "coordinateSystem.H"
#include "volFieldsFwd.H"

namespace Foam
{
namespace functionObjects
{

class porousforces
:
    public fvMeshFunctionObject,
    public writeFile
{
protected:

    // --- 累计统计受力 (全局规约值) ---
    vector sumDragForce_;
    vector sumDragForceA_;
    vector sumDragForceB_;
    vector sumDragForceC_;
    vector sumPressureForce_;
    vector sumGravityForce_;
    vector sumTotalMoment_;

    // --- 文件输出流 ---
    autoPtr<OFstream> forceFilePtr_;
    autoPtr<OFstream> momentFilePtr_;

    // --- 配置参数 (匹配旧算例键名) ---
    autoPtr<coordinateSystem> coordSysPtr_;

    word pName_;
    word UName_;
    word rhoName_;
    word alphaName_;

    scalar myrhoporo_;
    vector Linear_vel_;
    vector Angular_vel_;

    bool writeFields_;
    bool initialised_;

    // --- 内部持久化场管理 ---
    template<class GeoField>
    GeoField& getOrMakeField(const word& name, const dimensionSet& dims);

    void resetSum();
    void writeHeader(OFstream& os, const word& title) const;

public:

    //- Runtime type information
    TypeName("porousforces");

    // --- 构造函数 ---
    porousforces
    (
        const word& name,
        const Time& runTime,
        const dictionary& dict,
        const bool readFields = true
    );

    porousforces
    (
        const word& name,
        const objectRegistry& obr,
        const dictionary& dict,
        const bool readFields = true
    );

    porousforces(const porousforces&) = delete;
    void operator=(const porousforces&) = delete;

    virtual ~porousforces() = default;

    // --- 核心成员函数 ---
    virtual bool read(const dictionary& dict);
    virtual void calcForcesMoments();
    virtual vector forceEff() const;
    virtual vector momentEff() const;
    virtual bool execute();
    virtual bool write();
};

}
}

#endif
