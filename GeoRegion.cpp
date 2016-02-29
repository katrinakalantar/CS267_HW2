#include <math.h>
#include <float.h>
#include <cblas.h>
#include <stdlib.h>
#include <stdio.h>
#include "common.h"
#include <vector>
#include <iostream>
#include "GeoRegion.h"
#include "Region.h"
#include "common.h"
//#include <array>


GeoRegion::GeoRegion(particle_t *inp_data, double size, int n)
{
	double dimension = size;
	//double numRegions = pow(ceil((log2(n)/log2(4))), 2);
	int numRegions = 25; //PARAMETERIZED THIS
	int sqrt_numRegions = 5;
	double regionDim = size/sqrt(numRegions);
	//Region regionsList[numRegions];
	std::vector<Region> regionsList;
	//particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
	//Region *regionsList = (Region*) malloc( (numRegions) * sizeof(Region));
	double xdim;
	double ydim;
	for (int i = 1; i < (numRegions + 1); i++){
		if (i%sqrt_numRegions != 0){
			xdim = regionDim*((i%sqrt_numRegions)-1);
		}else{
			xdim = regionDim*(sqrt_numRegions-1);
		}
		ydim = min(regionDim*floor(i/sqrt_numRegions), regionDim*(sqrt_numRegions-1));
		regionsList[(i-1)] = Region(xdim, ydim, i, numRegions, regionDim);
		//regionsList.push_back(Region(xdim, ydim, i, numRegions, regionDim));
	}
}
//do i need a constructor or init function?? how do i use classes in c >_<

void GeoRegion::update_particles(particle_t &p, int regionNum){
	Region a = regionsList[regionNum];
	regionsList[regionNum].particles.push_back(p); //do i want a pointer to p (*p) or &p?
}
void GeoRegion::clear_particles(){
	std::vector<particle_t> v; //empty vector for swapping
	for (int i = 1; i < (numRegions + 1); i++){
		regionsList[i].particles.swap(v); //clears vector and reallocates memory
	}
}

void GeoRegion::update_location(particle_t &p, GeoRegion georeg){
	int pReg = get_region(p.x, p.y, georeg);
	p.region = pReg;
	update_particles(p, pReg); //is this how i call a function that is part of georeg?
	p.edge = check_edge(p.x, p.y, georeg);
}

int GeoRegion::get_region(double x, double y, GeoRegion georeg){
	int regNum = ceil(x/georeg.regionDim) + ((ceil(y/georeg.regionDim)-1)*sqrt(georeg.numRegions));
	return regNum;
}

int GeoRegion::check_edge(double x, double y, GeoRegion georeg){
	int edge = 0;
	//double xDistFromEdge = min((georeg.regionDim - (x%georeg.regionDim)), (x%georeg.regionDim)); //check both left and right edge
	int colNum = floor(x/georeg.regionDim);
	int rowNum = floor(y/georeg.regionDim);
	double xDistFromEdge = min((x - colNum*georeg.regionDim), ((colNum + 1)*georeg.regionDim - x));

	//double yDistFromEdge = min((georeg.regionDim - (y%georeg.regionDim)), (y%georeg.regionDim));
	double yDistFromEdge = min((y - rowNum*georeg.regionDim), ((rowNum + 1)*georeg.regionDim - y));
	double cutoff = 0.01; //TOOK THIS FROM THE DEFINE CUTOFF=0.01 in common file
	
	if ((xDistFromEdge < cutoff) or (yDistFromEdge < cutoff)){
		edge = 1; //this means the particle is too near to the edge of the region
	}
	return edge;
}

