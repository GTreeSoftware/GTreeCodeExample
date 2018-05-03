#ifndef TRACEUTIL_H
#define TRACEUTIL_H
#include "../../ngtypes/basetypes.h"
template<typename T> class Volume;

struct TraceUtil{
    static void GetGradientVectorFlowForTrace(const Volume<double> &sphereRayWet, Volume<double> &sphereRayGradient);
    static void SmoothGradientCurves(const std::vector<double> &ttk1, std::vector<double> &S_diffttk);
    static void ExtractSubRegionOP(const Vec3d&, Volume<NGCHAR>&, VectorVec3d&);//
    static void ExtractSubRegionAboveThrev(const VectorVec3d&, VectorVec3d&, Volume<NGCHAR>&);//extract area of which value is greater than threshold
};
template<class T>
struct VectorDecorator{
    VectorDecorator(std::vector<T>& p):m_Vec(p){}
    std::vector<T> &m_Vec;
    inline VectorDecorator& operator << (const T &t){
        m_Vec.push_back(t); return *this;
    }
    inline VectorDecorator& operator , (const T &t){
        m_Vec.push_back(t); return *this;
    }
};

#endif // TRACEUTIL_H
