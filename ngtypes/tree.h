#ifndef TREE_H
#define TREE_H
#include <tree.hh>
#include <iostream>
#include <vector>
#include <tuple>
#include <memory>
#include "basetypes.h"
#include "ineurondataobject.h"

struct NGHashList{//x mod 1000
    NGHashList(int n = 1000) :modn(n){
        //v.resize(n);
        v = new VectorVec3i **[modn];
        for (int i = 0; i < modn; ++i) {
            v[i] = NULL;//new VectorVec3i[modn];
        }
        //for (auto &it : v) it.resize(n);
        sortFunc = [&](const Vec3i &lhs, const Vec3i &rhs)->bool{ \
            if (lhs(0) != rhs(0)) return lhs(0) < rhs(0); \
            else if (lhs(1) != rhs(1)) return lhs(1) < rhs(1); \
            else if (lhs(2) != rhs(2)) return lhs(2) < rhs(2); \
            else return false; };
        comp1 = [](const Vec3i &rhs, const Vec3i &val)->bool{
            return rhs(0) - val(0) < 0; };
        comp2 = [](const Vec3i &rhs, const Vec3i &val)->bool{
            return rhs(1) - val(1) < 0; };
        comp3 = [](const Vec3i &rhs, const Vec3i &val)->bool{
            return rhs(2) - val(2) < 0; };
    }
    ~NGHashList(){ Clear(); }
    std::function<bool(const Vec3i &lhs, const Vec3i &rhs)> sortFunc,
        comp1, comp2, comp3;
    //std::function<bool(const Vec4i &lhs, const Vec4i &rhs)> findFunc;
    int mod_hash(const int arg)const{ return arg%modn; }
    int modn;
    //std::vector<std::vector<VectorVec3i> > v;
    VectorVec3i ***v = NULL;
    void insert(const Vec3i& arg);
    void _insert(const Vec3i& arg);
    //void insert(const VectorVec3d& arg);
    void insert(const std::vector<Line5d> &arg);
    bool isExistInRange(const Vec3i &val, int thre = 4) const;
    void Clear();
    void Reset();
    bool IsEmpty(){ return v != NULL; }
};

class TreeCurve : public INeuronDataObject
{
public:
    TreeCurve(INeuronProcessObject *p = NULL);
    ~TreeCurve();
    bool IsEmpty()const;
    void Swap(std::vector<VectorVec5d>&);
    void SetCurve(std::vector<VectorVec5d>&);
    void Clear(){m_Curve.clear();}
    size_t size()const{return m_Curve.size();}
    const std::vector<VectorVec5d>& GetCurve() const{return m_Curve;}
    std::vector<VectorVec5d>& GetCurve(){return m_Curve;}
private:
    std::vector<VectorVec5d> m_Curve;
};
typedef std::shared_ptr<TreeCurve> TreeCurvePointer;
typedef std::shared_ptr<const TreeCurve> CTreeCurvePointer;

class TreeConnect : public INeuronDataObject
{
public:
    TreeConnect(INeuronProcessObject *p = NULL);
    ~TreeConnect();
    bool IsEmpty()const;
    void Swap(VectorMat2i&);
    void SetConnect(const VectorMat2i&);
    void Clear(){m_Connect.clear();}
    size_t size()const{return m_Connect.size();}
    const VectorMat2i& GetConnect() const{return m_Connect;}
private:
    VectorMat2i m_Connect;
};

typedef std::shared_ptr<TreeConnect> TreeConnectPointer;
typedef std::shared_ptr<const TreeConnect> CTreeConnectPointer;


struct LineID
{
    LineID(size_t arg = 0, NEURITETYPE t = NEURITETYPE::DENDRITE) :id(arg), type(t){}
    ~LineID(){}
    size_t id;
    NEURITETYPE type;
    bool operator==(const LineID& other){ return other.id == this->id; }
    friend std::ostream& operator <<(std::ostream& str, const LineID& arg){
        str << arg.id;
        return str;
    }
};

typedef std::vector < std::pair<Vec3d, Vec3d> > BOUNDINFOLIST;//pos dir

class NeuronTree;
typedef std::shared_ptr<NeuronTree> NeuronTreePointer;
typedef std::shared_ptr<const NeuronTree> CNeuronTreePointer;
class NeuronTree : public INeuronDataObject
{
public:
    NeuronTree(INeuronProcessObject *p = NULL) { m_ProcessObject = p; identifyName = "NeuronTree";  }//m_curInitPt.second = std::numeric_limits<size_t>::max();
    ~NeuronTree(){ Clear(); }
    virtual bool IsEmpty()const{ return m_curveList.empty(); }
    void Clear(){ m_curveList.clear(); m_Topo.clear(); }
    size_t size()const{ return m_curveList.size(); }
    int m_type =2;
    tree<LineID> m_Topo;
    Vec3d m_curInitPt=Vec3d(0.0,0.0,0.0);//pos and curve id
    BOUNDINFOLIST m_traceInitInfo;
    std::vector<Line5d>  m_curveList;//5th is traverse flag
    std::vector<std::pair<int, NeuronTreePointer> > m_DelBranch;
	std::vector<std::vector<std::pair<size_t, Line5d> >> m_ChangedBranch;

    NGHashList m_hashList;
    void BuildHash();
    //traverse flag
    //std::vector<std::pair<int,int> > traverseFlag_;
    //int curLineID, curBeg = std::numeric_limits<int>::max(), curEnd;
    //std::vector<std::vector<char> > traverseFlag_;
    //check flag
    VectorVec3d m_suspiciousPositions_;
    std::vector<int> m_suspiciousCheckState_;
} ;


class NeuronPopulation : public INeuronDataObject
{
public:
    NeuronPopulation(INeuronProcessObject *p = NULL){ m_ProcessObject = p; identifyName = "NeuronPopulation"; }
    ~NeuronPopulation(){}
    bool IsEmpty()const{ return m_pop.empty(); }
    void Clear(){ m_pop.clear(); }
    size_t size()const{ return m_pop.size(); }
    std::vector<NeuronTreePointer>  m_pop;
    std::vector<NeuronTreePointer> deletedPop;//just save once.
};

typedef std::shared_ptr<NeuronPopulation> NeuronPopulationPointer;
typedef std::shared_ptr<const NeuronPopulation> CNeuronPopulationPointer;

//extra function for tree library
namespace KPTREE{
    NEURITETYPE Num2NeuriteType(int num);
    int NeuriteType2Num(NEURITETYPE);
    size_t FindNearestID(const Vec3d& pt, const Line5d& line, double threv = 1.0);
	size_t FindNearestID2(const Vec3d& pt, const Line5d& line, double &mindist);
    bool SearchLineIDInTree(const NeuronTree& nTree, size_t id, tree<LineID>::iterator &resit);
    bool SearchLineIDInPopulation(const NeuronPopulation& treeList, size_t id, tree<LineID>::iterator &resit);
    void UpdateDeleteTreeIDList(NeuronTree& nTree, size_t deletedIDList);
    void UpdateDeleteTreeIDList(NeuronTree& nTree, std::vector<size_t>& deletedIDList);
    void UpdateDeleteSortedTreeIDList(NeuronTree& nTree, std::vector<size_t>& deletedIDList);
    bool SetTreeRoot(NeuronTree&nTree, const Vec3d& root, double threv = 0.5, bool = true, bool = true, bool = true);
    bool SetTreeRoot(std::vector<Line5d>& curve, tree<LineID>& topo, const Vec3d& root, double threv = 0.5, bool flag = true, bool = true, bool = true);
    bool MergeTwoLine(NeuronTree& nTree, size_t mergeID1, size_t megeID2);
    size_t FastSearchNearestLine(const Vec3d& pt, const std::vector<Line5d>& lines, const size_t = 3lu, const double = 10.0);
    size_t FastSearchNearestLine(const Vec3d& pt, const std::vector<Line5d>& lines, size_t &vertID, const size_t = 3lu, const double = 10.0);
    //void FindSymmericDiff(const std::vector<Line5d>&, const std::vector<Line5d>&, const NGHashList&, const NGHashList&, int rad, int threv, Line3d& res);
    void FindSymmericDiff(const NeuronTree &tree1, const NeuronTree &tree2, int rad, int threv, Line3d& res, std::vector<int> &);

	/*khtao*/
	void FindCurrentImagePoints(const std::vector<Line5d> &lines, std::vector<Line5d> &Currentlines, std::vector<VectorVec2i> &CurrentIDs,
		size_t x_min, size_t x_max, size_t y_min, size_t y_max, size_t z_min, size_t z_max);
	void FindCurrentImagePoints(const std::vector<Line5d> &lines, const std::vector<size_t> curBranchID, std::vector<Line5d> &Currentlines, std::vector<VectorVec2i> &CurrentIDs,
		size_t x_min, size_t x_max, size_t y_min, size_t y_max, size_t z_min, size_t z_max);

	void InsertCurrentLine(std::vector<Line5d> &lines, const CellVec3d &InsetLines, std::vector<VectorVec2i> &CurrentIDs);
	size_t FindNearestID(const Vec3d& pt, const Line5d& line, double &mindist, double threv);

    bool ConjunctChildCurve(NeuronTree& tree, size_t id);
    int RefreshTree(NeuronTree& tree, int xMax = -1, int yMax= -1, int zMax = -1, READMODE = READMODE::MOSTD);
    int RefreshTreeCastration(NeuronTree& tree, int xMax = -1, int yMax = -1, int zMax = -1, READMODE = READMODE::MOSTD);

   template<class T>
void print_tree_bracketed(const tree<T>& t, std::ostream& str=std::cout);

template<class T>
void print_subtree_bracketed(const tree<T>& t, typename tree<T>::iterator iRoot, 
                                      std::ostream& str=std::cout);



// Iterate over all roots (the head) and print each one on a new line
// by calling printSingleRoot.

template<class T>
void print_tree_bracketed(const tree<T>& t, std::ostream& str) 
    {
    int headCount = t.number_of_siblings(t.begin());
    int headNum = 0;
    for(typename tree<T>::sibling_iterator iRoots = t.begin(); iRoots != t.end(); ++iRoots, ++headNum) {
        print_subtree_bracketed(t,iRoots,str);
        if (headNum != headCount) {
            str << std::endl;
            }
        }
    }


// Print everything under this root in a flat, bracketed structure.

template<class T>
void print_subtree_bracketed(const tree<T>& t, typename tree<T>::iterator iRoot, std::ostream& str) 
    {
    if(t.empty()) return;
    if (t.number_of_children(iRoot) == 0) {
        str << *iRoot;	
        }
    else {
        // parent
        str << *iRoot;
        str << "(";
        // child1, ..., childn
        int siblingCount = t.number_of_siblings(t.begin(iRoot));
        int siblingNum;
        typename tree<T>::sibling_iterator iChildren;
        for (iChildren = t.begin(iRoot), siblingNum = 0; iChildren != t.end(iRoot); ++iChildren, ++siblingNum) {
            // recursively print child
            print_subtree_bracketed(t,iChildren,str);
            // comma after every child except the last one
            if (siblingNum != siblingCount ) {
                str << ", ";
                }
            }
        str << ")";
        }
    }
}
#endif // TREE_H
