#include "soma.h"

Soma::Soma(INeuronProcessObject *p)
{
    m_ProcessObject = p;
    identifyName =  std::string("Soma");
}

bool Soma::IsEmpty() const
{
    return m_Source.empty();
}

void Soma::push_back(const Cell &obj)
{
    m_Source.push_back(obj);
}

Cell Soma::MakeCell(int i, double x, double y, double z, double r/*= 1.0*/,
    double a /*= 0.0*/, double b/*= 0.0*/, double mr/*= 0.0*/, double s/*= 0.0*/, double v/*= 0.0*/)
{
    return Cell(i, x, y, z, r, a, b, mr, s, v);
}

Cell Soma::MakeCell(double x, double y, double z, double r/*= 1.0*/,
    double a /*= 0.0*/, double b/*= 0.0*/, double mr/*= 0.0*/, double s/*= 0.0*/, double v/*= 0.0*/)
{
    return Cell(0, x, y, z, r, a, b, mr, s, v);
}
