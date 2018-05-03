//time : 2013-11-20
//this file include Eigen and std::vector to provide base data structure
//author:zhouhang,WLNO BMP

#ifndef BASETYPES_H
#define BASETYPES_H
#include <Eigen/Core>
#include <vector>
#include <stack>

#define NG_SMART_POINTER_TYPEDEF(type, name) typedef std::shared_ptr<type> name
#define NG_SMART_POINTER_DEFINE(type, var) std::shared_ptr<type> var
#define NG_SMART_POINTER_NEW_DEFAULT(type, var) var = std::shared_ptr<type>(new type)
#define NG_SMART_POINTER_NEW(type, var, para) var = std::shared_ptr<type>(new type(para))
#define NG_SMART_POINTER_DEFINE_NEW_DEFAULT(type, var) std::shared_ptr<type> var = std::shared_ptr<type>(new type)
#define NG_SMART_POINTER_DEFINE_NEW(type, var, para) std::shared_ptr<type> var = std::shared_ptr<type>(new type(para))
#define NG_CREATE_DYNAMIC_CAST(type, var1, var2)  std::shared_ptr<type> var1 = std::dynamic_pointer_cast<type>(var2)
#define NG_DYNAMIC_CAST(type, var1, var2)  var1 = std::dynamic_pointer_cast<type>(var2)
//macro:output error file and line
#define NG_ERROR_MESSAGE(userInfo) { printf(" %s %d:", __FILE__, __LINE__); printf(userInfo); printf("\n");}

typedef long long int Int64;

typedef unsigned char NGCHAR;
enum class DATATYPE{ IMAGE1, IMAGE8, IMAGE16, IMAGE32, SOMA, TREE };
enum class READMODE{ INVALID, TIFF3D, MOSTD, HDF5 };
enum class Record{ DeleteBranch, SmoothCurve };
enum class NEURITETYPE{ UNKNOWN, SOMA, AXON, DENDRITE, APICAL };

typedef Eigen::Vector2i Vec2i;
typedef Eigen::Vector3i Vec3i;
typedef Eigen::Vector4i Vec4i;
typedef Eigen::VectorXi VecXi;

typedef Eigen::Matrix2i Mat2i;
typedef Eigen::Matrix3i Mat3i;
typedef Eigen::Matrix4i Mat4i;
typedef Eigen::Matrix<int, 6, 1> Vec6i;

typedef Eigen::Vector2f Vec2f;
typedef Eigen::Vector3f Vec3f;
typedef Eigen::Vector4f Vec4f;
typedef Eigen::Matrix<float, 5, 1> Vec5f;

typedef Eigen::Matrix2f Mat2f;
typedef Eigen::Matrix3f Mat3f;
typedef Eigen::Matrix4f Mat4f;

typedef Eigen::Vector2d Vec2d;
typedef Eigen::Vector3d Vec3d;
typedef Eigen::Vector4d Vec4d;
typedef Eigen::Matrix<double, 5, 1> Vec5d;
typedef Eigen::Matrix<double, 6, 1> Vec6d;
typedef Eigen::Matrix<double, 7, 1> Vec7d;

typedef Eigen::Matrix2d Mat2d;
typedef Eigen::Matrix3d Mat3d;
typedef Eigen::Matrix4d Mat4d;

typedef Eigen::VectorXd VecXd;
typedef Eigen::MatrixXd MatXd;

typedef Eigen::VectorXi VecXi;
typedef Eigen::MatrixXi MatXi;

typedef std::vector<Vec2i, Eigen::aligned_allocator<Vec2i> > VectorVec2i;
typedef std::vector<Vec3i, Eigen::aligned_allocator<Vec3i> > VectorVec3i;
typedef std::vector<Vec4i, Eigen::aligned_allocator<Vec4i> > VectorVec4i;
typedef std::vector<Vec6i, Eigen::aligned_allocator<Vec6i> > VectorVec6i;

typedef std::vector<Vec3f, Eigen::aligned_allocator<Vec3f> > VectorVec3f;
typedef std::vector<Vec4f, Eigen::aligned_allocator<Vec4f> > VectorVec4f;

typedef std::vector<Vec2d, Eigen::aligned_allocator<Vec2d> > VectorVec2d;
typedef std::vector<Vec3d, Eigen::aligned_allocator<Vec3d> > VectorVec3d;
typedef VectorVec3d Line3d;
typedef std::vector<Vec4d, Eigen::aligned_allocator<Vec4d> > VectorVec4d;
typedef std::vector<Vec5d, Eigen::aligned_allocator<Vec5d> > VectorVec5d;
typedef VectorVec5d Line5d;
typedef std::vector<Vec6d, Eigen::aligned_allocator<Vec6d> > VectorVec6d;
typedef std::vector<Vec7d, Eigen::aligned_allocator<Vec7d> > VectorVec7d;
typedef std::vector<VecXd, Eigen::aligned_allocator<VecXd> > VectorVecXd;

typedef std::vector<Mat2i, Eigen::aligned_allocator<Mat2i> > VectorMat2i;

typedef std::vector<Mat2d, Eigen::aligned_allocator<Mat2d> > VectorMat2d;
typedef std::vector<Mat3d, Eigen::aligned_allocator<Mat3d> > VectorMat3d;
typedef std::vector<Mat4d, Eigen::aligned_allocator<Mat4d> > VectorMat4d;
typedef std::vector<Vec7d, Eigen::aligned_allocator<Vec7d> > VectorMat7d;
typedef std::vector<MatXd, Eigen::aligned_allocator<MatXd> > VectorMatXd;

typedef std::stack<Vec3d, Eigen::aligned_allocator<Vec3d> > StackVec3d;


typedef std::vector<VectorVec3d, Eigen::aligned_allocator<VectorVec3d> > CellVec3d;
typedef std::vector<CellVec3d, Eigen::aligned_allocator<CellVec3d> > XCellVec3d;
typedef std::vector<VectorVec4d, Eigen::aligned_allocator<VectorVec4d> > CellVec4d;
typedef std::vector<VectorVec7d, Eigen::aligned_allocator<VectorVec7d> > CellVec7d;

#endif // BASETYPES_H
