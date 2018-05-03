#include "NGUtility.h"
#include "contourutil.h"

namespace NGUtility{
    int sign(double a)
    {
        if (a < 0) return -1;
        else if (a > 0) return 1;
        return 0;
    }

    //************************************
    // Method:    CalcPtCurveDirection
    // FullName:  NGUtility::CalcPtCurveDirection
    // Access:    public 
    // Returns:   void
    // Qualifier:
    // Parameter: const VectorVec3d & ptCurve
    // Parameter: Vec3d & fir
    //************************************
    void CalcPtCurveDirection(const VectorVec3d &ptCurve, Vec3d &fir)
    {
        int nx = int(ptCurve.size());
        Vec3d tmpFir;
        tmpFir.setZero();
        double tmp(0.0);
        Vec3d tmpVec;
        for (int i = 0; i < nx - 1; ++i){
            tmpVec = ptCurve[i + 1] - ptCurve[i];
            tmp = tmpVec.norm();
            tmpFir += tmp * tmpVec;
        }
        fir = tmpFir.normalized();
    }

    void GetPartVectorVec3d(const VectorVec3d &orig, int minId, int maxId, VectorVec3d &dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetPartVectorVec3d.\n");
            //LOG(ERROR) << "error size in GetPartVectorVec3d.";
            return;
        }
        dst.resize(maxId - minId + 1);
        for (int i = minId; i <= maxId; ++i){
            dst[i - minId] = orig[i];
        }
    }

    void GetPartVectorVec3d(const VectorVec4d &orig, int minId, int maxId, VectorVec3d &dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetPartVectorVec3d.\n");
            //LOG(ERROR) << "error size in GetPartVectorVec3d.";
            return;
        }
        dst.resize(maxId - minId + 1);
        for (int i = minId; i <= maxId; ++i){
            dst[i - minId] = orig[i].block(0, 0, 3, 1);
        }
    }

    void GetPartVectorVec3d(const VectorVec5d &orig, int minId, int maxId, VectorVec3d &dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetPartVectorVec3d.\n");
            //LOG(ERROR) << "error size in GetPartVectorVec3d.";
            return;
        }
        dst.resize(maxId - minId + 1);
        Vec3d tmp;
        for (int i = minId; i <= maxId; ++i){
            tmp << orig[i](0), orig[i](1), orig[i](2);
            dst[i - minId] = tmp;
        }
    }

    void GetReversePartVectorVec3d(const VectorVec3d &orig, int minId, int maxId, VectorVec3d &dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetReversePartVectorVec3d.\n");
            //LOG(ERROR) << "error size in GetReversePartVectorVec3d.";
            return;
        }
        int len = maxId - minId + 1;
        dst.resize(len);
        for (int i = 0; i < len; ++i){
            dst[i] = orig[maxId - i];
        }
    }

    void GetReversePartVectorVec3d(const VectorVec4d &orig, int minId, int maxId, VectorVec3d &dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetReversePartVectorVec3d.\n");
            //LOG(ERROR) << "error size in GetReversePartVectorVec3d.";
            return;
        }
        int len = maxId - minId + 1;
        dst.resize(len);
        for (int i = 0; i < len; ++i){
            dst[i] = orig[maxId - i].block(0, 0, 3, 1);
        }
    }

    void GetReversePartVectorVec3d(const VectorVec5d &orig, int minId, int maxId, VectorVec3d &dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetReversePartVectorVec3d.\n");
            //LOG(ERROR) << "error size in GetReversePartVectorVec3d.";
            return;
        }
        int len = maxId - minId + 1;
        dst.resize(len);
        for (int i = 0; i < len; ++i){
            dst[i] = orig[maxId - i].block(0, 0, 3, 1);
        }
    }

    bool IsEqualDouble(double a, double b)
    {
        return std::abs(a - b) < 0.0001;
    }

    void GetPartVectorVec3dStep(const VectorVec3d& orig, int minId, int maxId, int step, VectorVec3d& dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetPartVectorVec3dStep.\n");
            //LOG(ERROR) << "error size in GetPartVectorVec3dStep.";
            return;
        }
        //dst.resize(maxId - minId + 1);
        for (int i = minId; i <= maxId; i += step){
            dst.push_back(orig[i]);
        }
    }

    void GetPartVectorVec3dStep(const VectorVec5d& orig, int minId, int maxId, int step, VectorVec3d& dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetPartVectorVec3dStep.\n");
            //LOG(ERROR) << "error size in GetPartVectorVec3dStep.";
            return;
        }
        //dst.resize(maxId - minId + 1);
        Vec3d tmp;
        for (int i = minId; i <= maxId; i += step){
            tmp << orig[i](0), orig[i](1), orig[i](2);
            dst.push_back(tmp);
        }
    }

    void GetPartVectorVec3i(const VectorVec5d& orig, int minId, int maxId, VectorVec3i& dst)
    {
        if (minId > maxId || maxId >= (int)orig.size() || minId < 0) {
            dst.clear();
            printf("error size in GetPartVectorVec3dStep.\n");
            //LOG(ERROR) << "error size in GetPartVectorVec3dStep.";
            return;
        }
        dst.resize(maxId - minId + 1);
        Vec3i tmp;
        for (int i = minId; i <= maxId; ++i){
            tmp << Round(orig[i](0)), Round(orig[i](1)), Round(orig[i](2));
            dst[i - minId] = tmp;
        }
    }
    void Vec3d2Vec3i(const Vec3d& orig, Vec3i& dst)
    {
        dst << Round(orig(0)), Round(orig(1)), Round(orig(2));
    }

    void Vec3i2Vec3d(const Vec3i& orig, Vec3d& dst)
    {
        dst << orig(0), orig(1), orig(2);
    }

    Vec5d MakeVec5d(double a1, double a2, double a3, double a4/*=0*/, double a5/*=0*/)
    {
        Vec5d tmp; tmp << a1, a2, a3, a4, a5; return tmp;
    }

    Vec5d MakeVec5d(Vec3d arg, double a4/*=0*/, double a5/*=0*/)
    {
        Vec5d tmp; tmp << arg(0), arg(1), arg(2), a4, a5; return tmp;
    }

    Vec6d MakeVec6d(double a1, double a2, double a3, double a4/*=0*/, double a5/*=0*/, double a6/*=0*/)
    {
        Vec6d tmp; tmp << a1, a2, a3, a4, a5, a6; return tmp;
    }


    void WriteVectorVec3i(const char* path, const char* param, const VectorVec3i& arg)
    {
        FILE* fp;
        if ((fp = fopen(path, "w")) == NULL) {
            printf("WriteVectorVec3d error.\n");
            return;
        }
        int id = 0;
        if (strcmp(param, "dot") == 0) {
            for (auto &it : arg) {
                fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, double(it(0)), double(it(1)), double(it(2)) );
            }
        }
        else if (strcmp(param, "line") == 0) {
            fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, double(arg[0](0)), double(arg[0](1)), double(arg[0](2)));
            int pid = 0;
            for (auto &it : arg) {
                fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", ++id, double(it(0)), double(it(1)), double(it(2)), ++pid);
            }
        }
        fclose(fp);
    }


    void WriteVectorVec3d(const char *path, const char *param, const VectorVec3d &arg)
    {
        FILE* fp;
        if ((fp = fopen(path, "w")) == NULL) {
            printf("WriteVectorVec3d error.\n");
            return;
        }
        int id = 0;
        if (strcmp(param, "dot") == 0) {
            for (auto &it : arg) {
                fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, it(0), it(1), it(2));
            }
        }
        else if (strcmp(param, "line") == 0) {
            fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, arg[0](0), arg[0](1), arg[0](2));
            int pid = 0;
            for (auto &it : arg) {
                fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", ++id, it(0), it(1), it(2), ++pid);
            }
        }
        fclose(fp);
    }

    void WriteVectorVec3d(const char *path, const char *param, const VectorVec3i &arg)
    {
        FILE* fp;
        if ((fp = fopen(path, "w")) == NULL) {
            printf("WriteVectorVec3d error.\n");
            return;
        }
        int id = 0;
        if (strcmp(param, "dot") == 0) {
            for (auto &it : arg) {
                fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, double(it(0)), double(it(1)), double(it(2)));
            }
        }
        else if (strcmp(param, "line") == 0) {
            fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, double(arg[0](0)), double(arg[0](1)), double(arg[0](2)));
            int pid = 0;
            for (auto &it : arg) {
                fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", ++id, double(it(0)), double(it(1)), double(it(2)), ++pid);
            }
        }
        fclose(fp);
    }

    void WriteVectorVec5d(const char *path, const char *param, const VectorVec5d &arg)
    {
        FILE* fp;
        if ((fp = fopen(path, "w")) == NULL) {
            printf("WriteVectorVec3d error.\n");
            return;
        }
        int id = 0;
        if (strcmp(param, "dot") == 0) {
            for (auto &it : arg) {
                fprintf(fp, "%d 1 %lf %lf %lf %lf -1\n", ++id, it(0), it(1), it(2), it(3));
            }
        }
        else if (strcmp(param, "line") == 0) {
            fprintf(fp, "%d 1 %lf %lf %lf 1.0 -1\n", ++id, arg[0](0), arg[0](1), arg[0](2));
            int pid = 0;
            for (auto &it : arg) {
                fprintf(fp, "%d 1 %lf %lf %lf 1.0 %d\n", ++id, it(0), it(1), it(2), it(3), ++pid);
            }
        }
        fclose(fp);
    }

    void WriteVectorVectorVec5d(const char *path, const std::vector<VectorVec5d> &arg)
    {
        FILE* fp;
        if ((fp = fopen(path, "w")) == NULL) {
            printf("WriteVectorVec3d error.\n");
            return;
        }
        int globalIndex = 1;
        for (size_t i = 0; i < arg.size(); ++i){
            //One curve
            const VectorVec5d &localCurve = arg[i];
            if (localCurve.empty()) continue;
            fprintf(fp, "%d %d %lf %lf %lf %lf -1\n", globalIndex, 2, localCurve[0](0),
                localCurve[0](1), localCurve[0](2), 1.0);
            ++globalIndex;
            for (size_t k = 1; k < localCurve.size(); ++k){
                fprintf(fp, "%d %d %lf %lf %lf %lf %d\n", globalIndex, 2, localCurve[k](0),
                    localCurve[k](1), localCurve[k](2), 1.0, globalIndex - 1);
                ++globalIndex;
            }
        }
        fclose(fp);
        fclose(fp);
    }

    void Principald(const VectorVec3d &dataL, Vec3d &x1)
    {
        int nx = int(dataL.size());
        Vec3d nn;
        nn.setZero();
        double tmp(0.0);
        Vec3d tmpVec;
        for (int i = 0; i < nx - 1; ++i){
            tmpVec = dataL[i + 1] - dataL[i];
            tmp = tmpVec.norm();
            nn += tmp * tmpVec;
        }
        x1 = nn.normalized();
    }

    VecXd WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg)
    {
        typedef double spheredist;
        int nxss = int(rayNode.size());
        VecXd rayNodeWet; rayNodeWet.setZero(nxss);

        spheredist x, y, z;
        spheredist segmentDense, segmentWet;//, w1;
        for (int i = 0; i < nxss; ++i){
            x = rayNode[i](0);
            y = rayNode[i](1);
            z = rayNode[i](2);
            segmentDense = segmentWet = 0.0;
            ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
            rayNodeWet(i) = segmentDense / (segmentWet + 0.0001);
        }
        return rayNodeWet;
    }

    double WeighRayValue(const Vec3d &rayNode, const Volume<unsigned short> &locOrigImg)
    {
        typedef double spheredist;
        spheredist x, y, z;
        spheredist segmentDense, segmentWet;//, w1;
        x = rayNode(0);
        y = rayNode(1);
        z = rayNode(2);
        segmentDense = segmentWet = 0.0;
        ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
        return segmentDense / (segmentWet + 0.0001);
    }

    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg,
        std::vector<double> &rayNodeWet)
    {
        typedef double spheredist;
        int nxss = int(rayNode.size());
        rayNodeWet.clear();

        //зјБъ
        typedef double spheredist;
        spheredist x, y, z;
        spheredist segmentDense, segmentWet;//, w1;
        for (int i = 0; i < nxss; ++i){
            x = rayNode[i](0);
            y = rayNode[i](1);
            z = rayNode[i](2);
            segmentDense = segmentWet = 0.0;
            ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
            rayNodeWet.push_back(segmentDense / (segmentWet + 0.0001));
        }
    }

    void WeighRayValue(const VectorVec3d &rayNode, const Volume<unsigned short> &locOrigImg, VecXd &rayNodeWet)
    {
        typedef double spheredist;
        int nxss = int(rayNode.size());
        rayNodeWet.setZero(nxss);

        spheredist x, y, z;
        spheredist segmentDense, segmentWet;//, w1;
        for (int i = 0; i < nxss; ++i){
            x = rayNode[i](0);
            y = rayNode[i](1);
            z = rayNode[i](2);
            segmentDense = segmentWet = 0.0;
            ContourUtil::CalculateSphereOneNode(locOrigImg, 0.05, x, y, z, segmentDense, segmentWet);
            rayNodeWet(i) = segmentDense / (segmentWet + 0.0001);
        }
    }

    bool IsPointInRange(const Vec3d& pt, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax)
    {
        Vec3i tmpPt(int(pt(0)), int(pt(1)), int(pt(2)));
        if (tmpPt(0) < xMin || tmpPt(0) > xMax || tmpPt(1) < yMin || tmpPt(1) > yMax ||
            tmpPt(2) < zMin || tmpPt(2) > zMax) {
            return false;
        }
        return true;
    }

    bool IsPointInRange(const Vec5d& pt, int xMin, int xMax, int yMin, int yMax, int zMin, int zMax)
    {
        Vec3i tmpPt(int(pt(0)), int(pt(1)), int(pt(2)));
        if (tmpPt(0) < xMin || tmpPt(0) > xMax || tmpPt(1) < yMin || tmpPt(1) > yMax ||
            tmpPt(2) < zMin || tmpPt(2) > zMax) {
            return false;
        }
        return true;
    }

    template<>
    bool CompareVector(std::vector<double>& lhs, std::vector<double> &rhs)
    {
        if (lhs.size() != rhs.size()) return false;
        std::sort(lhs.begin(), lhs.end()); std::sort(rhs.begin(), rhs.end());
        for (size_t i = 0; i < lhs.size(); ++i) {
            if (NGUtility::IsEqualDouble(lhs[i], rhs[i])) return false;
        }
        return true;
    }

    template<>
    bool IsInVector(const std::vector<double>& lhs, const double& arg)
    {
        if (lhs.empty()) return false;
        for (size_t i = 0; i < lhs.size(); ++i) {
            if (NGUtility::IsEqualDouble(lhs[i], arg)) return true;
        }
        return false;
    }

    

}