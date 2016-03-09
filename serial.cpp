#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "common.h"
#include "GeoRegion.h"
#include "Region.h"

//
//  benchmarking program
//
int main( int argc, char **argv )
{    

    //printf("inside main of serial.cpp\n");
    int navg,nabsavg=0;
    double davg,dmin, absmin=1.0, absavg=0.0;

    if( find_option( argc, argv, "-h" ) >= 0 )
    {
        printf( "Options:\n" );
        printf( "-h to see this help\n" );
        printf( "-n <int> to set the number of particles\n" );
        printf( "-o <filename> to specify the output file name\n" );
        printf( "-s <filename> to specify a summary file name\n" );
        printf( "-no turns off all correctness checks and particle output\n");
        return 0;
    }

    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;


    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    double size = set_size( n );
    init_particles( n, particles );
    GeoRegion aGeoRegion = GeoRegion(particles, n);

    for (int i = 0; i < n; i++){
    	aGeoRegion.update_location(particles[i]);
    }
    //
    //  simulate a number of time steps
    //
    double simulation_time = read_timer( );
	
    for( int step = 0; step < NSTEPS; step++ )
    {
	navg = 0;
        davg = 0.0;
	dmin = 1.0;
        //
        //  compute forces
        //
	std::vector<Region> regList = aGeoRegion.get_regionsList();
	for (int regInd = 0; regInd < regList.size(); regInd++){
		std::vector<particle_t> regPart;
		std::vector<int> neighborsInd;
		std::vector<particle_t> neighborParts;
        std::vector<particle_t> nhParts;
		std::vector<Region*> neighbors;
		Region* curReg = &regList[regInd];
		neighborsInd = curReg->get_neighbors();
		for (std::vector<int>::size_type nInd = 0; nInd < neighborsInd.size(); nInd++){
			neighbors.push_back(&regList[(neighborsInd[nInd]-1)]);
		}
		regPart = curReg->get_particles();
		for (std::vector<particle_t>::size_type partInd = 0; partInd < regPart.size(); partInd++){
			regPart[partInd].ax = regPart[partInd].ay = 0;
			for (std::vector<int>::size_type partBInd1 = 0; partBInd1 < partInd; partBInd1++){
				apply_force(regPart[partInd],regPart[partBInd1], &dmin, &davg, &navg);
			}
			for (std::vector<int>::size_type partBInd2 = (partInd+1); partBInd2 < regPart.size(); partBInd2++){
				apply_force(regPart[partInd],regPart[partBInd2], &dmin, &davg, &navg);
			}
			if (regPart[partInd].edge == 1){
				for (std::vector<Region*>::size_type neighInd = 0; neighInd < neighbors.size(); neighInd++){
					neighborParts = neighbors[neighInd]->get_particles();
					for (std::vector<particle_t>::size_type npInd = 0; npInd < neighborParts.size(); npInd++){
						apply_force(regPart[partInd],neighborParts[npInd], &dmin, &davg, &navg);
					}
					vector<particle_t>().swap(neighborParts);
				}
			}

		}
		vector<Region*>().swap(neighbors);
		vector<int>().swap(neighborsInd);
		vector<particle_t>().swap(regPart);
	}
 
        //
        //  move particles
        //
        aGeoRegion.clear_particles();

        for( int i = 0; i < n; i++ ){
        	move( particles[i] );
        	aGeoRegion.update_location(particles[i]);
        }

        if( find_option( argc, argv, "-no" ) == -1 )
        {
          //
          // Computing statistical data
          //
          if (navg) {
            absavg +=  davg/navg;
            nabsavg++;
          }
          if (dmin < absmin) absmin = dmin;
		
          //
          //  save if necessary
          //
          if( fsave && (step%SAVEFREQ) == 0 )
              save( fsave, n, particles );
        }
    }
    simulation_time = read_timer( ) - simulation_time;

    printf("printing simulation time in seconds\n");
    printf( "n = %d, simulation time = %g seconds", n, simulation_time);

    if( find_option( argc, argv, "-no" ) == -1 )
    {
      if (nabsavg) absavg /= nabsavg;
    // 
    //  -The minimum distance absmin between 2 particles during the run of the simulation
    //  -A Correct simulation will have particles stay at greater than 0.4 (of cutoff) with typical values between .7-.8
    //  -A simulation where particles don't interact correctly will be less than 0.4 (of cutoff) with typical values between .01-.05
    //
    //  -The average distance absavg is ~.95 when most particles are interacting correctly and ~.66 when no particles are interacting
    //
    printf( ", absmin = %lf, absavg = %lf", absmin, absavg);
    if (absmin < 0.4) printf ("\nThe minimum distance is below 0.4 meaning that some particle is not interacting");
    if (absavg < 0.8) printf ("\nThe average distance is below 0.8 meaning that most particles are not interacting");
    }
    printf("\n");     

    //
    // Printing summary data
    //
    if( fsum) 
        fprintf(fsum,"%d %g\n",n,simulation_time);
 
    //
    // Clearing space
    //
    if( fsum )
        fclose( fsum );    
    free( particles );
    if( fsave )
        fclose( fsave );
    
    return 0;
}
