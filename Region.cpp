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
	double x1 = x_in;
	double y1 = y_in;
	int num1 = num_in;
	double dim1 = regDim;
	std::vector<int> neighbors1;
	int width1 = sqrt(totalNum);
	int a1 = num1-width1;
	int b1 = num1+width1;
	std::vector<particle_t> particles1;
	if ((num1 % width1 == 1) or (num1 % width1 == 0)){
		//printf("Region_a\n");
		if (num1 == totalNum){
			neighbors1.push_back(a1-1);
			neighbors1.push_back(a1);
			neighbors1.push_back(num1-1);
		}else if (num1 == 1){
			neighbors1.push_back(num1+1);
			neighbors1.push_back(b1);
			neighbors1.push_back(b1+1);
		}else if (num1 == width1){
			neighbors1.push_back(num1-1);
			neighbors1.push_back(b1-1);
			neighbors1.push_back(b1);
		}else if (num1 == totalNum - (width1)+1){
			neighbors1.push_back(a1);
			neighbors1.push_back(a1+1);
			neighbors1.push_back(num1+1);
		}else if (num1 % width1 == 0){
			neighbors1.push_back(a1-1);
			neighbors1.push_back(a1);
			neighbors1.push_back(num1-1);
			neighbors1.push_back(b1-1);
			neighbors1.push_back(b1);
		}else if (num1 % width1 == 1){
			neighbors1.push_back(a1);
			neighbors1.push_back(a1+1);
			neighbors1.push_back(num1+1);
			neighbors1.push_back(b1);
			neighbors1.push_back(b1+1);
		}
	}else if (num1 < width1){
		//printf("Region_b\n");
		neighbors1.push_back(num1-1);
		neighbors1.push_back(num1+1);
		neighbors1.push_back(b1-1);
		neighbors1.push_back(b1);
		neighbors1.push_back(b1+1);
	}else if ((num1 < totalNum) and (num1 > (totalNum - width1))){
		//printf("Region_c\n");
		neighbors1.push_back(a1-1);
		neighbors1.push_back(a1);
		neighbors1.push_back(a1+1);
		neighbors1.push_back(num1-1);
		neighbors1.push_back(num1+1);
	}else{
		//printf("Region_d\n");
		neighbors1.push_back(a1-1);
		neighbors1.push_back(a1);
		neighbors1.push_back(a1+1);
		neighbors1.push_back(num1-1);
		neighbors1.push_back(num1+1);
		neighbors1.push_back(b1-1);
		neighbors1.push_back(b1);
		neighbors1.push_back(b1+1);
	}
	//printf("num neighbors = %lu\n", neighbors1.size());
	init(x1,y1,num1,dim1,neighbors1,particles1,width1,a1,b1);

}

Region::Region(){

}

void Region::init(double inp_x, double inp_y, int inp_num, double inp_dim, std::vector<int> inp_neighbors, std::vector<particle_t> inp_particles, int inp_width, int inp_a, int inp_b){

	x = inp_x;
	y = inp_y;
	num = inp_num;
	dim = inp_dim;
	neighbors = inp_neighbors;
	particles = inp_particles;
	width = inp_width;
	a = inp_a;
	b = inp_b;

}

double Region::get_x(){
	return x;
}

double Region::get_y(){
	return y;
}

int Region::get_num(){
	return num;
}

double Region::get_dim(){
	return dim;
}
std::vector<int> Region::get_neighbors() {
	return neighbors;
}
std::vector<particle_t> Region::get_particles() {
	return particles;
}

void Region::add_particle(particle_t *inp_part) {
	particles.push_back(*inp_part);
}

int Region::get_width(){
	return width;
}

int Region::get_a(){
	return a;
}
int Region::get_b(){
	return b;
}
