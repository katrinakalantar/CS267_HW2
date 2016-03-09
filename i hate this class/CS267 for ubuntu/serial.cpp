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

    //printf("1: after help menu\n");
    
    int n = read_int( argc, argv, "-n", 1000 );

    char *savename = read_string( argc, argv, "-o", NULL );
    char *sumname = read_string( argc, argv, "-s", NULL );
    
    FILE *fsave = savename ? fopen( savename, "w" ) : NULL;
    FILE *fsum = sumname ? fopen ( sumname, "a" ) : NULL;

    //printf("2: assigned initial variables\n");

    particle_t *particles = (particle_t*) malloc( n * sizeof(particle_t) );
    //printf("3: created particle_t array\n");
    double size = set_size( n );
    //printf("size = %f\n",size);
    //printf("3.5\n");
    init_particles( n, particles );
    //printf("\n4\n");
    GeoRegion aGeoRegion = GeoRegion(particles, n);

    /*
    printf("regionDim in serial.cpp = %f\n",aGeoRegion.get_regionDim());
    printf("numRegions in serial.cpp = %d\n",aGeoRegion.get_numRegions());

    printf("sqrt_numRegions in serial.cpp = %d\n",aGeoRegion.get_sqrt_numRegions());
    printf("dimension in serial.cpp = %f\n",aGeoRegion.get_dimension());

    printf("5\n");
     */

    for (int i = 0; i < n; i++){
        //printf("particle num %d\n", i);
    	aGeoRegion.update_location(particles[i]);//, aGeoRegion);
    }
    //printf("6\n");
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
    int totalParts = 0;
    for (int rInd = 0; rInd < regList.size(); rInd++){
        Region* cReg = &regList[rInd];
        //printf("num region %d\n", rInd);
        std::vector<particle_t> partz = cReg->get_particles();
        //printf("num parts %lu\n", partz.size());
        totalParts += partz.size();
    }
    //printf("num particles total %d\n", totalParts);
    //printf("total num regions %lu\n", regList.size());

	//std::vector<Region> emptyReg;
	//std::vector<particle_t> emptyPart;
	//std::vector<int> emptyInt;
	for (int regInd = 0; regInd < regList.size(); regInd++){
		std::vector<particle_t> regPart;
		std::vector<int> neighborsInd;
		std::vector<particle_t> neighborParts;
        std::vector<particle_t> nhParts;
		std::vector<Region*> neighbors;
		//printf("region %d\n", regInd);
		Region* curReg = &regList[regInd];
		//printf("region assigned");
		neighborsInd = curReg->get_neighbors();
		//printf("neighbors assigned");
		//Region neighbors = (Region) malloc( neighborsInd.size() * sizeof(Region) );
		for (std::vector<int>::size_type nInd = 0; nInd < neighborsInd.size(); nInd++){
			//printf("neighbor %lu\n", nInd);
            //printf("neighborInd %d\n", neighborsInd[nInd]);
			neighbors.push_back(&regList[(neighborsInd[nInd]-1)]);
		}
		//printf("neighbors created");
		//printf("len neighbors %d\n", neighbors.size());
		regPart = curReg->get_particles();
        //printf("num reg parts %lu\n", regPart.size());
        //for (std::vector<Region*>::size_type nhInd = 0; nhInd < neighbors.size(); nhInd++){
                    //nhParts = neighbors[nhInd]->get_particles();
                    //printf("num neighbor parts %lu\n", nhParts.size());
                    //std::vector<particle_t>().swap(nhParts);
        //}
		//printf("num particles %d\n", regPart.size());
		for (std::vector<particle_t>::size_type partInd = 0; partInd < regPart.size(); partInd++){
			regPart[partInd].ax = regPart[partInd].ay = 0;
			for (std::vector<int>::size_type partBInd1 = 0; partBInd1 < partInd; partBInd1++){
				apply_force(regPart[partInd],regPart[partBInd1], &dmin, &davg, &navg);
			}
			for (std::vector<int>::size_type partBInd2 = (partInd+1); partBInd2 < regPart.size(); partBInd2++){
				apply_force(regPart[partInd],regPart[partBInd2], &dmin, &davg, &navg);
			}
			//printf("done with particles in region\n");
			if (regPart[partInd].edge == 1){
				//printf("edge\n");
				for (std::vector<Region*>::size_type neighInd = 0; neighInd < neighbors.size(); neighInd++){
					neighborParts = neighbors[neighInd]->get_particles();
                    //printf("num neighbor parts %lu\n", neighborParts.size());
					for (std::vector<particle_t>::size_type npInd = 0; npInd < neighborParts.size(); npInd++){
						apply_force(regPart[partInd],neighborParts[npInd], &dmin, &davg, &navg);
					}
					vector<particle_t>().swap(neighborParts);
				}
            
				//printf("done with neighbors\n");
			}

		}
		vector<Region*>().swap(neighbors);
		vector<int>().swap(neighborsInd);
		vector<particle_t>().swap(regPart);
		//neighborsInd.swap(emptyInt);
		//regPart.swap(emptyPart);
	}
//        for( int i = 0; i < n; i++ )
//        {
//            //printf("computing forces for i%d\n",i);
//            particles[i].ax = particles[i].ay = 0; //change following for loop to just use particles in same
//            //region or neighbors (if edge)
//            int reg = particles[i].region;
//            //printf("reg=%d\n",reg);
//            std::vector<int> relev_reg;
//            relev_reg.push_back(particles[i].region);
//            if((particles[i].edge == 1)){ // or (particles[1].edge == 0)){
//                std::vector<int> reg_neighbors = aGeoRegion.get_regionsList()[reg-1].get_neighbors();
//                relev_reg.insert(relev_reg.end(), reg_neighbors.begin(), reg_neighbors.end());
//                //printf("edge\n");
//            }else if (particles[i].edge != 0){
//                printf("aaaaaaaaaaah!!!!");
//            }
//
//            //printf("reg_neighbors size: %lu\n", relev_reg.size());
//
//            for(std::vector<int>::size_type k = 0; k != relev_reg.size(); k++) {
//                //printf("inside loop k\n");
//                //printf("k=%d\n",reg_neighbors[k]);
//                //std::vector<particle_t> parts = reg_neighbors[k].get_particles(); // aGeoRegion.get_regionsList()[reg_neighbors[k-1]].particles;
//                std::vector<particle_t> parts =  aGeoRegion.get_regionsList()[relev_reg[k]-1].get_particles();
//
//                for(std::vector<particle_t>::size_type j = 0; j != parts.size(); j++){
//                    if(particles[i].id!=parts[j].id){
//                        apply_force(particles[i],parts[j], &dmin, &davg, &navg);
//                    }
////                    else
////                        //printf("particles[i] == particles[j]\n");
//                }
//            }
//
//        }
 
        //
        //  move particles
        //
        //printf("just before moving particles\n");
        aGeoRegion.clear_particles();

        for( int i = 0; i < n; i++ ){
            double pxold = particles[i].x;
            double pyold = particles[i].y;
            //printf("p x %f\n", particles[i].x);
            //printf("p y %f\n", particles[i].y);
        	move( particles[i] );
            double pxnew = particles[i].x;
            double pynew = particles[i].y;
            if(pxnew == pxold){
                printf("shit, didn't move in x %d\n", i);
            }
            if(pynew == pyold){
                printf("shit, didn't move in y %d\n", i);
            }
            //printf("p new x %f\n", particles[i].x);
            //printf("p new y %f\n", particles[i].y);
        	aGeoRegion.update_location(particles[i]);
        }
        	//update_location(particles[i], aGeoRegion);
        //printf("finished moving particles\n");

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
    //free(aGeoRegion);
    if( fsave )
        fclose( fsave );
    
    return 0;
}
