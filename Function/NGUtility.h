#ifndef NGUTILITY
#define  NGUTILITY
#include "../ngtypes/basetypes.h"
#include "../ngtypes/volume.h"
#pragma warning(disable : 4996)
namespace NGUtility{
    int sign(double);
    void GetPartVectorVec3d(const VectorVec3d&, int minId, int maxId, VectorVec3d&);
    void GetPartVectorVec3d(const VectorVec4d&, int minId, int maxId, VectorVec3d&);
    void GetPartVectorVec3d(const VectorVec5d&, int minId, int maxId, VectorVec3d&);
    void GetReversePartVectorVec3d(const VectorVec3d&, int minId, int maxId, VectorVec3d&);
    void GetReversePartVectorVec3d(const VectorVec4d&, int minId, int maxId, VectorVec3d&);
    void GetReversePartVectorVec3d(const VectorVec5d&, int minId, int maxId, VectorVec3d&);
    bool IsEqualDouble(double a, double b);
    void GetPartVectorVec3dStep(const VectorVec3d& orig, int minId, int maxId, int step, VectorVec3d& dst);
    void GetPartVectorVec3dStep(const VectorVec5d& orig, int minId, int maxId, int step, VectorVec3d& dst);
    void GetPartVectorVec3i(const VectorVec5d&, int minId, int maxId, VectorVec3i&);
    void Vec3d2Vec3i(const Vec3d&, Vec3i&);
    void Vec3i2Vec3d(const Vec3i&, Vec3d&);
    Vec5d MakeVec5d(double a1, double a2, double a3, double a4 = 0, double a5 = 0);
    Vec5d MakeVec5d(Vec3d, double a4 = 0, double a5 = 0);
    void WriteVectorVec3i(const char*, const char*, const VectorVec3i&);
    void WriteVectorVec3d(const char*, const char*, const VectorVec3d&);
    void WriteVectorVec3d(const char*, const char*, const VectorVec3i&);
    void WriteVectorVec5d(const char*, const char*, const VectorVec5d&);
    void WriteVectorVectorVec5d(const char *path, const std::vector<VectorVec5d> &arg);
    void Principald(const VectorVec3d &dataL, Vec3d &x1);
    Vec6d MakeVec6d(double a1, double a2, double a3, double a4/*=0*/, double a5/*=0*/, double a6/*=0*/);
    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, std::vector<double> &rayNodeWet);
    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, VecXd &rayNodeWet);
    VecXd WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg);
    double WeighRayValue(const Vec3d &rayNode, const Volume<unsigned short> &locOrigImg);

    bool IsPointInRange(const Vec3d&, int, int, int, int, int, int);
    bool IsPointInRange(const Vec5d&, int, int, int, int, int, int);
    void CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir);

    template<class T>
    struct VectorDecorator{
        VectorDecorator(std::vector<T>& p) :m_Vec(p){}
        std::vector<T> &m_Vec;
        inline VectorDecorator& operator << (const T &t){
            m_Vec.push_back(t); return *this;
        }
        inline VectorDecorator& operator , (const T &t){
            m_Vec.push_back(t); return *this;
        }
    };

    template<typename T>
    int Round(T val)
    {
        return int(val + 0.5 * (val > 0 ? 1 : -1));
    }
   

    template<class T>
    void VectorAddSwapBack(std::vector<T>& arg, T& item)
    {
        arg.push_back(T());
        arg.back().swap(item);
    }

    template<class T>
    bool CompareVector(std::vector<T>& lhs, std::vector<T> &rhs){
        if (lhs.size() != rhs.size()) return false;
        std::sort(lhs.begin(), lhs.end()); std::sort(rhs.begin(), rhs.end());
        for (size_t i = 0; i < lhs.size(); ++i) {
            if (lhs[i] != rhs[i]) return false;
        }
        return true;
    }
    template<>
    bool CompareVector<double>(std::vector<double>& lhs, std::vector<double> &rhs);

    template<class T>
    bool IsInVector(const std::vector<T>& lhs, const T& arg){
        if (lhs.empty()) return false;
        for (size_t i = 0; i < lhs.size(); ++i) {
            if (lhs[i] == arg) return true;
        }
        return false;
    }
    template<>
    bool IsInVector<double>(const std::vector<double>& lhs, const double& arg);

}


#endif // !UTILITY
