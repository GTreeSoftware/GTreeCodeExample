#ifndef MOSTIMAGE_H
#define MOSTIMAGE_H
#include "imagecache.h"
#include <tiffio.h>
#include <QtGui>
#include<qdir.h>
#include <iostream>
#include <ctime>

#include <QDialog>
#include <fstream>
#include <sstream>
#define BLOCK_SIZE 512.0
#define DEPTH_SIZE 512
#define SIZE_TEST 1000
#pragma warning(disable:4996)
using  namespace std;

class MostImage 
{
	
public:
	int id;
	int img_type;           //图像类型，区分血管数据，胞体数据，荧光数据
	long sz0;            //三个方向的大小，像素
	long sz1;
	long sz2;
	int temp;
	int level_size;
	int block_num[10][3];
	int IOR_sz0;
	int IOR_sz1;
	int IOR_sz2;
	double rez_x, rez_y, rez_z; //the resolution of a image pixel along the 3 axes
	double rez_IOR_x, rez_IOR_y, rez_IOR_z; //the resolution of a image pixel along the 3 axes
	double origin_x, origin_y, origin_z;
	double origin_IOR_x, origin_IOR_y, origin_IOR_z;
	//unsigned char * data1d;  //图像数据
    unsigned short * data1d;  //图像数据
	//unsigned char * IORdata; //选取的感兴趣区域
    unsigned short * IORdata; //选取的感兴趣区域
	QString file_path;
    ImageCaches imageCaches_;

    MostImage():
        imageCaches_(10, 16)
	{
		id = 0;
		img_type = 0;
		sz0 = sz1 = sz2 = 0;
		IOR_sz0 = IOR_sz1 = IOR_sz2 = 0;
		rez_x = rez_y =rez_z = 1;
		rez_IOR_x = rez_IOR_y = rez_IOR_z = 0;
		origin_x = origin_y = origin_z =0;
		origin_IOR_x = origin_IOR_y = origin_IOR_z=0;
		temp =0;
		data1d=0;
		IORdata=0;
		file_path.clear();
		level_size=0;
		for (int i =0;i<10;i++)
		{
			block_num[i][0]=0;
			block_num[i][1]=0;
			block_num[i][2]=0;
		}

	}

	virtual ~MostImage()
	{
	}
		
	bool  select_IOR_file(int xbegin,int xend,int ybegin,int yend,int zbegin,int zend,int level);

	bool loadDataInfo(const char *path);      
	
	};


#endif
