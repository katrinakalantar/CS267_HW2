#ifndef __CS267_REGION_H__
#define __CS267_REGION_H__
#include <vector>

using namespace std;

class Region
{
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
	void init(double inp_x, double inp_y, int inp_num, double inp_dim, std::vector<int> inp_neighbors, std::vector<particle_t> inp_particles, int inp_width, int inp_a, int inp_b);

public:


	//Region(double x, double y, int num, int totalNum, double regDim);
	Region(double x_in, double y_in, int num_in, int totalNum, double regDim);
	Region();
	double get_x();
	double get_y();
	int get_num();
	double get_dim();
	std::vector<int> get_neighbors();
	std::vector<particle_t> get_particles();
	int get_width();
	int get_a();
	int get_b();
	void add_particle(particle_t *inp_part);
	void Rclear_particles();


};

#endif
