#include<tiffio.h>
#include<iostream>
#include"most_data.h"
#include <sstream>
#include "QString"

using namespace std;
#define BLOCK_SIZE 512.0

int main(int argc, char *argv[])
{
			
	const char * open_path = "M:/mostd_test2/test.mostd";
	MostImage *mostimage = new MostImage();
	mostimage->loadDataInfo(open_path);

	mostimage->select_IOR_file(0,14225,0,10932,0,2416,5);
	return 0;

}
