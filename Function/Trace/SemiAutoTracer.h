#ifndef SEMIAUTOTRACER_H
#define  SEMIAUTOTRACER_H
#include "Function\ineuronprocessobject.h"
#include "ngtypes\basetypes.h"
#include "ngtypes\volume.h"
#include "Function\volumealgo.h"
#include "ngtypes\tree.h"
#include "ngtypes\ParamPack.h"
class SemiAutoTracer;
NG_SMART_POINTER_TYPEDEF(SemiAutoTracer, NGSemiAutoTracer);
class SemiAutoTracer :
    public INeuronProcessObject
{
public:
    static NGSemiAutoTracer New(){ return NGSemiAutoTracer(new SemiAutoTracer()); }
    SemiAutoTracer();
    ~SemiAutoTracer();
    virtual ProcStatPointer Update();
    void SetParam(NGParamPack arg){ paramPack = arg; }
   void GetSemiTraceLine(VectorVec5d &arg)const{ arg = traceLine_; }
   const VectorVec5d& GetSemiTraceLine()const{ return traceLine_; }
   void SetBeginEnd(const Vec3d &bP, const Vec3d &eP){ beginPt_ = bP; endPt_ = eP; }
   void SetGlobalOffset(int xO, int yO, int zO){ xOffSet_ = xO; yOffSet_ = yO; zOffSet_ = zO; }
protected:
	Vec3d CalcNormDir(const Vec3d& p1, const Vec3d &p2);
	void CalcPrinDirAndDistList(const VectorVec3d &ptLine, Vec3d &mainDirection);
    virtual IDataPointer ReleaseData();
    virtual ConstIDataPointer GetOutput();
private:
    int xOffSet_, yOffSet_, zOffSet_;
    Vec3d beginPt_, endPt_;
    VectorVec5d traceLine_;
    std::shared_ptr<const Volume<unsigned short> > origImgPointer;
    NGParamPack paramPack;
};

#endif