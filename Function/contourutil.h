
#ifndef CONTOURUTIL_H
#define CONTOURUTIL_H
#include <vector>
#include "../ngtypes/basetypes.h"
template<typename T> class Volume;
typedef Volume<NGCHAR> CVolume;
typedef Volume<unsigned short> SVolume;
typedef Volume<double> DVolume;

struct ContourUtil
{
    /**/
    static void CalAreaVolume(const VectorVec3d &dataPP, const int Theta, const int Phi,
        VectorVec3d &Patch, double &area, double &cellVol);

    static void AddAreaVolume(const Vec3d &center, const Vec3d &A, const Vec3d &B, const Vec3d &C,
        Vec3d &AB, Vec3d &AC, Vec3d &OA, Vec3d &OB, Vec3d &OC,
        double &area, double &cellVol);

    /*߱*/
    static void GetRayLimit(const DVolume &sphereRayWet,
        const double constrictionThrev,
        std::vector<std::vector<double> > &ray_limit);
    /**/
    static void CalculateOneRayLimit(const std::vector<double> &ray,
        const double constriction_threv,
        int &one_ray_limit);


    static void CalculateSphereOneNode(const CVolume &locOrigImg,
        const double value, 
        const double x, 
        const double y, 
        const double z, 
        double &segmentDense,
        double& segmentWet);

    static void GetGradientVectorFlow(const DVolume &Sphere_XX, DVolume &Uz);

    static void SmoothGradientCurves(const std::vector<double> &ttk, std::vector<double> &S_diffttk);

    static void GetBoundaryBack(const std::vector<double> &outer_shell, 
        const double threv, 
                                std::vector<double> &boundary_back);

    //2015-6-8
    static void CalculateSphereOneNode(const SVolume &locOrigImg, const double value,
                                       const double x, const double y, const double z,
                                       double &segmentDense, double &segmentWet);
};



#endif
