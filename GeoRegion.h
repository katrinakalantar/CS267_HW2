#ifndef __CS267_GEOREGION_H__
#define __CS267_GEOREGION_H__
#include "Region.h"
#include "common.h"

using namespace std;

class GeoRegion
{



public:
	double dimension;
	double numRegions;
	double regionDim;
	//Region regionsList;
	std::vector<Region> regionsList;
	double xdim;
	double ydim;
	
	GeoRegion(particle_t *inp_data, double size, int n);
	void update_particles(particle_t &p, int regionNum);
	void clear_particles();
	void update_location(particle_t &p, GeoRegion georeg);
	int get_region(double x, double y, GeoRegion georeg);
	int check_edge(double x, double y, GeoRegion georeg);
};


#endif
