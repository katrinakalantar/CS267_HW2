#include <math.h>
#include <float.h>
#include <cblas.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include <vector>
#include <iostream>
#include "Region.h"

Region::Region(double x_in, double y_in, int num_in, int totalNum, double regDim){
	//double ul.x = x;
	//double ul.y = y;
	double x = x_in;
	double y = y_in;
	int num = num_in;
	double dim = regDim;
	std::vector<int> neighbors;
	int width = sqrt(totalNum);
	int a = num-width;
	int b = num+width;
	std::vector<particle_t> particles;
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

Region::Region(){

}
