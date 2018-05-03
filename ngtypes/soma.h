#ifndef SOMA_H
#define SOMA_H
#include "ineurondataobject.h"
#include <vector>
struct Cell{
    Cell(){}
    ~Cell(){}
    Cell(int i, double xx, double yy, double zz, double rr, double aa, double bb, double mrr,
         double ss, double vv):ID(i), x(xx), y(yy), z(zz), r(rr), a(aa), b(bb), mr(mrr), s(ss),
        v(vv){}
    int ID;
    double x, y, z, r, a, b, mr, s, v;//
};

class Soma : public INeuronDataObject
{
public:
    static Cell MakeCell(int i, double, double, double, double = 1.0, double = 0.0, double = 0.0, double = 0.0, double = 0.0, double = 0.0);
    static Cell MakeCell(double, double, double, double = 1.0, double = 0.0, double = 0.0, double = 0.0, double = 0.0, double = 0.0);
    Soma(INeuronProcessObject* p = NULL);
    bool IsEmpty()const;
    void push_back(const Cell&);
    void clear(){m_Source.clear();}
    const Cell& GetCell(size_t idx)const{return m_Source[idx];}
    Cell& GetCell(size_t idx){return m_Source[idx];}
    std::vector<Cell>& GetAllCell(){ return m_Source; }
    size_t size()const {return m_Source.size();}
    std::vector<Cell> m_Source; // why capsule?
private:
    Soma(const Soma&);//uncopyable
    Soma& operator=(const Soma&);
};

#endif // SOMA_H
