#include <math.h>
#include <float.h>
//#include <home/cc/cs199/fa13/class/cs199-cly/Downloads/ATLAS/include/cblas.h>
//#include <cblas.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include <vector>
#include <iostream>
#include "GeoRegion.h"
#include "Region.h"
#include "common.h"
//#include <array>


GeoRegion::GeoRegion(particle_t *inp_data, int n)
{
	//printf("inside GeoRegion\n");
	double dimension1 = sqrt(.0005 * n);
	int sqrt_numRegions1 = (int) ceil(pow(n, .25));
	int numRegions1 = pow(sqrt_numRegions1, 2);
//	int numRegions1 = pow(ceil((log2(n)/log2(20))), 2);
//	if (numRegions1 > n){
//		numRegions1 = pow(floor(sqrt(n)),2);
//	}
	//int numRegions1 = 25; //PARAMETERIZED THIS
	//int sqrt_numRegions1 =5;
//	int sqrt_numRegions1 = (int) sqrt(numRegions1);
	double regionDim1 = dimension1/sqrt(numRegions1);
	double a = dimension1/5.0;
	//printf("regionDim = %f\n",regionDim1);
	//printf("a = %f\n",a);
	//printf("inside GeoRegion - 2\n");
	//Region regionsList[numRegions];
	std::vector<Region> regionsList1;
	std::vector<Region*> ptrregionsList1;
	//printf("inside GeoRegion - 3\n");
	//particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
	//Region *regionsList = (Region*) malloc( (numRegions) * sizeof(Region));
	double xdim1;
	double ydim1;
	for (int i = 1; i < (numRegions1 + 1); i++){
		//printf("region num= %d\n", i);
		//printf("inside GeoRegion - forLoop - \n");
		//printf("i=%d\n",i);
		if (i%sqrt_numRegions1 != 0){
			//printf("inside GeoRegion - forLoop - inside if\n");
			xdim1 = regionDim1*((i%sqrt_numRegions1)-1);
		}else{
			//printf("inside GeoRegion - forLoop - inside else\n");
			xdim1 = regionDim1*(sqrt_numRegions1-1);
		}
		//printf("inside GeoRegion - forLoop - final calculations %d\n",i);
		ydim1 = min(regionDim1*floor(i/sqrt_numRegions1), regionDim1*(sqrt_numRegions1-1));
		//regionsList[(i-1)] = Region(xdim, ydim, i, numRegions, regionDim);
		//printf("regionDim inside for loop %f\n",regionDim1);
		//printf("whee %d\n", i);
		//Region aRegion = Region(xdim1, ydim1, i, numRegions1, regionDim1);
		regionsList1.push_back(Region(xdim1, ydim1, i, numRegions1, regionDim1));
		//printf("whee %d\n", i);
	}
	//printf("regionDim last line of creation %f\n",regionDim1);

	init(dimension1, numRegions1, sqrt_numRegions1, regionDim1, regionsList1, xdim1, ydim1);
}

void GeoRegion::init(double inp_dimension, int inp_numRegions, int inp_sqrt_numRegions, double inp_regionDim, std::vector<Region> inp_regionsList, double inp_xdim, double inp_ydim){

	dimension = inp_dimension;
	numRegions = inp_numRegions;
	sqrt_numRegions = inp_sqrt_numRegions;
	regionDim = inp_regionDim;
	regionsList = inp_regionsList;
	xdim = inp_xdim;
	ydim = inp_ydim;

}
//do i need a constructor or init function?? how do i use classes in c >_<


double GeoRegion::get_dimension(){
	return dimension;
}

int GeoRegion::get_numRegions(){
	return numRegions;
}

int GeoRegion::get_sqrt_numRegions(){
	return sqrt_numRegions;
}

double GeoRegion::get_regionDim(){
	return regionDim;
}

double GeoRegion::get_xdim(){
	return xdim;
}

double GeoRegion::get_ydim(){
	return ydim;
}

std::vector<Region> GeoRegion::get_regionsList(){
	return regionsList;
}


void GeoRegion::update_particles(particle_t &p, int regionNum){
	//printf("update_particles - regionNum - %d\n",regionNum);
	//Region* a = &regionsList[regionNum-1]; //TRYING THIS TO SEE IF IT IS AN INDEX ERROR
	//printf("a.x = %f\n",a.x);
	//printf("region Num is %d\n", regionNum);
	regionsList[regionNum-1].add_particle(&p); //do i want a pointer to p (*p) or &p?
	//printf("pushed back particle p\n");
}
void GeoRegion::clear_particles(){
	//std::vector<particle_t> v; //empty vector for swapping
	for (int i = 0; i < (numRegions); i++){
		regionsList[i].Rclear_particles(); //clears vector and reallocates memory
		//printf("how many particles %lu\n", regionsList[i].get_particles().size());
	}
}

void GeoRegion::update_location(particle_t &p){//}, GeoRegion georeg){
	int pReg = get_region(p.x, p.y);//, georeg);
	//printf("pReg = %d\n",pReg);
	p.region = pReg;
	//printf("about to update_particles\n");
	update_particles(p, pReg); //is this how i call a function that is part of georeg?
	p.edge = check_edge(p.x, p.y);//, georeg);
	//printf("pEdge = %d\n",p.edge);
//	if (p.edge == 1){
//		printf("pEdge = %d\n",p.edge);
//	}
}

//void GeoRegion::update_location_fast(){
//	for (int regInd = 0; regInd < numRegions; regInd++ ){
//		curReg = regionsList[regInd];
//		regParts = curReg.get_particles();
//		for (partInd = 0; partInd < regParts.size(); partInd++){
//
//		}
//	}
//}

int GeoRegion::get_region(double x, double y){//}, GeoRegion georeg){
	//printf("regionDim = %f\n",regionDim);//georeg.regionDim);
	//printf("x= %f\n", x);
	//printf("y= %f\n", y);

	//int regNum = ceil(x/georeg.regionDim) + ((ceil(y/georeg.regionDim)-1)*sqrt(georeg.numRegions));
	int regNum = ceil(x/regionDim) + ((ceil(y/regionDim)-1)*sqrt(numRegions));

	//printf("refNum %d\n", regNum);
	return regNum;
}

int GeoRegion::check_edge(double x, double y){// GeoRegion georeg){
	int edge = 0;
	//double xDistFromEdge = min((georeg.regionDim - (x%georeg.regionDim)), (x%georeg.regionDim)); //check both left and right edge
	int colNum = floor(x/regionDim);//georeg.regionDim);
	int rowNum = floor(y/regionDim);//georeg.regionDim);
	double xDistFromEdge = min((x - colNum*regionDim), ((colNum + 1)*regionDim - x));
	//double xDistFromEdge = min((x - colNum*georeg.regionDim), ((colNum + 1)*georeg.regionDim - x));

	//double yDistFromEdge = min((georeg.regionDim - (y%georeg.regionDim)), (y%georeg.regionDim));
	double yDistFromEdge = min((y - rowNum*regionDim), ((rowNum + 1)*regionDim - y));
	//double yDistFromEdge = min((y - rowNum*georeg.regionDim), ((rowNum + 1)*georeg.regionDim - y));
	double cutoff = 0.01; //TOOK THIS FROM THE DEFINE CUTOFF=0.01 in common file
	
	if ((xDistFromEdge < cutoff) or (yDistFromEdge < cutoff)){
		edge = 1; //this means the particle is too near to the edge of the region
	}
	return edge;
}

//GeoRegion::~GeoRegion(){
//	delete dimension;
//	delete numRegions;
//	delete sqrt_numRegions;
//	delete regionsDim;
//	delete regionsList;
//	delete xdim;
//	delete ydim;
//}
