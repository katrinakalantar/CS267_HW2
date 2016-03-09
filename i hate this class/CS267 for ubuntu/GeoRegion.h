#ifndef __CS267_GEOREGION_H__
#define __CS267_GEOREGION_H__
#include "Region.h"
#include "common.h"

using namespace std;

class GeoRegion
{
	double dimension;
	int numRegions;
	double regionDim;
	int sqrt_numRegions;
	double xdim;
	double ydim;
	std::vector<Region> regionsList;
	void init(double inp_dimension, int inp_numRegions, int inp_sqrt_numRegions, double inp_regionDim, std::vector<Region> inp_regionsList, double inp_xdim, double inp_ydim);
//private:


public:

	GeoRegion(particle_t *inp_data, int n);
	void update_particles(particle_t &p, int regionNum);
	void clear_particles();
	void update_location(particle_t &p);//, GeoRegion georeg);
	int get_region(double x, double y);//, GeoRegion georeg);
	int check_edge(double x, double y);// GeoRegion georeg);
	double get_dimension();
	int get_numRegions();
	double get_regionDim();
	int get_sqrt_numRegions();
	double get_xdim();
	double get_ydim();
	std::vector<Region> get_regionsList();
};


#endif
