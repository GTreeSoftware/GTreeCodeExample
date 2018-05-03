#ifndef SEMIAUTOCONNCET_H
#define SEMIAUTOCONNCET_H
#include "ngtypes\basetypes.h"
#include "ngtypes\volume.h"

struct pointTable
{
	double x, y, z, grayValue, costToStart, costToEnd;
	double positionToPrevious;//1 is equal level,sqrt(2) is intermediate level,sqrt(3) is vertex
	double getTotalCost(){ return (costToStart + costToEnd); };
	Vec3d positionVector;
	std::shared_ptr<pointTable> fatherPoint;

	pointTable()
	{
		x = 0; y = 0; z = 0;
		grayValue = 0;
		costToStart = 0;
		costToEnd = 0;
		positionToPrevious = 0;
	};
	pointTable(double arg1, double arg2, double arg3)
	{
		x = arg1; y = arg2; z = arg3;
		grayValue = 0;
		costToStart = 0;
		costToEnd = 0;
		positionToPrevious = 0;
	};

	pointTable(const pointTable &arg)
	{
		*this = arg;
	};

	pointTable operator+(pointTable &arg)const
	{
		pointTable sum(*this);
		//sum.x = x + arg.x;
		//sum.y = y + arg.y;
		//sum.z = z + arg.z;
		//sum.positionVector = positionVector;
		//sum.positionToPrevious = positionToPrevious;

		////sum.fatherPoint = std::make_shared<pointTable>(arg);//no copy function
		////sum.fatherPoint = arg;

		return sum += arg;
	};
	pointTable& operator+=(pointTable &arg)
	{
		x = x + arg.x;
		y = y + arg.y;
		z = z + arg.z;
		positionVector = positionVector;
		positionToPrevious = positionToPrevious;

		//sum.fatherPoint = std::make_shared<pointTable>(new pointTable(arg));//no copy function
		//sum.fatherPoint.get() = &arg;

		return *this;
	};

	pointTable& operator=(const pointTable &arg)
	{
		x = arg.x;
		y = arg.y;
		z = arg.z;
		grayValue = arg.grayValue;
		costToStart = arg.costToStart;
		costToEnd = arg.costToEnd;
		positionToPrevious = arg.positionToPrevious;
		positionVector = arg.positionVector;
		//deep copy
		//here recursion, need to control the end of the conditions
		if (arg.fatherPoint.get())
			fatherPoint = arg.fatherPoint;/////////////////////

		return *this;
	};

	bool operator==(const pointTable &arg)
	{
		return (x == arg.x && y == arg.y && z == arg.z);
	};

};

class SemiAutoConnect
{
public:
	SemiAutoConnect(Vec3d&, Vec3d&, const std::shared_ptr<const SVolume>&,double);
	~SemiAutoConnect();
	VectorVec3d pathTable, resampleLine;
	//std::vector<std::vector<std::shared_ptr<pointTable>>> allOpenTable;
	double  normThreshold = 0.07;
	bool errorFlag = false;
	void search();

private:
	double grayMax, grayMin,resampleThreshold;
	bool searchFlag = false;
	IVolume dataFlag;
	std::shared_ptr<const SVolume> sourceData;
	Vec3d sToendVec, vertexPoint;
	pointTable startPoint, endPoint;
	std::shared_ptr<pointTable> currentPoint;
	std::vector<std::shared_ptr<pointTable>>  allNeighborTable, currentOpenTable;
	std::vector<std::shared_ptr<pointTable>>::iterator currPointIte, endPointFinder;

	bool eveluationRange(pointTable&);
	double normaValue(double);
	double calcCostToStart(pointTable&);
	double calcCostToEnd(pointTable&);

	void getGrayRange();
	void getNeighborToOpenTable();
	void searchOpenTable();
	void ThreeDimdCurveResample(VectorVec3d&, double, VectorVec3d&);
};

#endif