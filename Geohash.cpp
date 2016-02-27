#include <math.h>
#include <float.h>
#include <cblas.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include <vector>
#include <iostream>

Region::Region(double x, double y, int num, int totalNum, double regDim){
	double ul.x = x;
	double ul.y = y;
	int num = num;
	double dim = regDim;
	std::vector<int> neighbors;
	int width = sqrt(totalNum);
	int a = num-width;
	int b = num+width;
	if ((num % width == 1) or (num % width == 0)){
		if (num == totalNum){
			neighbors.push_back(a-1);
			neighbors.push_back(a);
			neighbors.push_back(num-1);
		}else if (num == 1){
			neighbors.push_back(num+1);
			neighbors.push_back(b);
			neighbors.push_back(b+1);
		}else if (num == width){
			neighbors.push_back(num-1);
			neighbors.push_back(b-1);
			neighbors.push_back(b);
		}else if (num == totalNum - (width+1)){
			neighbors.push_back(a);
			neighbors.push_back(a+1);
			neighbors.push_back(num+1);
		}else if (num % width == 0){
			neighbors.push_back(a-1);
			neighbors.push_back(a);
			neighbors.push_back(num-1);
			neighbors.push_back(b-1);
			neighbors.push_back(b);
		}else if (num % width == 1){
			neighbors.push_back(a);
			neighbors.push_back(a+1);
			neighbors.push_back(num+1);
			neighbors.push_back(b);
			neighbors.push_back(b+1);
		}
	}else if (num < width){
		neighbors.push_back(num-1);
		neighbors.push_back(num+1);
		neighbors.push_back(b-1);
		neighbors.push_back(b);
		neighbors.push_back(b+1);
	}else if ((num < totalNum) and (num > (totalNum - width))){
		neighbors.push_back(a-1);
		neighbors.push_back(a);
		neighbors.push_back(a+1);
		neighbors.push_back(num-1);
		neighbors.push_back(num+1);
	}else{
		neighbors.push_back(a-1);
		neighbors.push_back(a);
		neighbors.push_back(a+1);
		neighbors.push_back(num-1);
		neighbors.push_back(num+1);
		neighbors.push_back(b-1);
		neighbors.push_back(b);
		neighbors.push_back(b+1);
	}
	

}


GeoRegion::GeoRegion(particle_t *inp_data, double size, int n)
{
	double dimension = size;
	double numRegions = 2*ceiling((log2(n)/log2(4)));//how to always get a multiple of 2
	double regionDim = size/sqrt(numRegions);
	Region regionsList[numRegions];
	double xdim;
	double ydim;
	for (int i = 1; i < (numRegions + 1; i++)){
		if (i%sqrt(numRegions) != 0){
			xdim = regionDim*((i%sqrt(numRegions))-1);
		}else{
			xdim = regionDim*(sqrt(numRegions)-1);
		}
		ydim = min(regionDim*floor(i/sqrt(numRegions)), regionDim*(sqrt(numRegions)-1));
		regionsList[(i-1)] = Region(xdim, ydim, i, numRegions, regionDim)
	}
	//check regions method
}
