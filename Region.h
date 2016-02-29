#ifndef __CS267_REGION_H__
#define __CS267_REGION_H__
#include <vector>

using namespace std;

class Region
{


public:
	double x;
	double y;
	int num;
	double dim;
	std::vector<int> neighbors;
	//vector<int> neighbors;
	std::vector<particle_t> particles;
	int width;
	int a;
	int b;

	//Region(double x, double y, int num, int totalNum, double regDim);
	Region(double x_in, double y_in, int num_in, int totalNum, double regDim);
	Region();



};

#endif
