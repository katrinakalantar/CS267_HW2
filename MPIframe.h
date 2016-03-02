/*
 *  quadtree.h
 *  Header file for a quadtree.
 *
 *  https://github.com/ninjin/barnes-hut-sne/blob/master/quadtree.h
 *
 *  Created by Laurens van der Maaten.
 *  Copyright 2012, Delft University of Technology. All rights reserved.
 *
 */

#ifndef QUADTREE_H
#define QUADTREE_H
#include "common.h"
#include <stdlib.h>
#include <unordered_map>
#include <vector>
#include <assert.h>
#include <iostream>
#include "frame.h"

#define density 0.0005

static inline double min(double x, double y) { return (x <= y ? x : y); }
static inline double max(double x, double y) { return (x <= y ? y : x); }
static inline double abs(double x) { return (x < .0 ? -x : x); }

class MPIFrame: public Frame {

public:
    MPIFrame(const int _rank, const int _block_stride,
             const int _n_x, const int _n_y,
             particle_t *particles, const int _n_particles) :
            Frame(_n_x, _n_y, particles, _n_particles),
            rank(_rank),
            block_stride(_block_stride),
            block_x(_rank / block_stride),
            block_y(_rank % block_stride){

        NW_arrival_buffer = new particle_t[_n_particles];
        NW_n_arrivals = 0;

        N_arrival_buffer = new particle_t[_n_particles];
        N_arrivals = 0;

        NE_arrival_buffer = new particle_t[_n_particles];
        NE_arrivals = 0;

        W_arrival_buffer = new particle_t[_n_particles];
        W_n_arrivals = 0;

        E_arrival_buffer = new particle_t[_n_particles];
        E_arrivals = 0;

        SW_arrival_buffer = new particle_t[_n_particles];
        SW_arrivals = 0;

        S_arrival_buffer = new particle_t[_n_particles];
        S_arrivals = 0;

        SE_arrival_buffer = new particle_t[_n_particles];
        SE_arrivals = 0;

        NW_departure_buffer = new particle_t[_n_particles];
        NW_n_departures = 0;

        N_departure_buffer = new particle_t[_n_particles];
        N_n_departures = 0;

        NE_departure_buffer = new particle_t[_n_particles];
        NE_n_departures = 0;

        W_departure_buffer = new particle_t[_n_particles];
        W_n_departures = 0;

        E_departure_buffer = new particle_t[_n_particles];
        E_n_departures = 0;

        SW_departure_buffer = new particle_t[_n_particles];
        SW_n_departures = 0;

        S_departure_buffer = new particle_t[_n_particles];
        S_n_departures = 0;

        SE_departure_buffer = new particle_t[_n_particles];
        SE_n_departures = 0;

    }

    inline void update_locations(particle_t* particles, const int n_particles, bool next_frame = true) {

        particle_t ****target_grid;
        int **target_size_grid;
        if (next_frame) {
            target_grid = next_part_grid;
            target_size_grid = next_size_grid;
        } else {
            target_grid = part_grid;
            target_size_grid = size_grid;
        }

        int x_idx, y_idx;
        for (int i = 0; i < n_particles; ++i) {
            get_idx(particles[i].x, particles[i].y, x_idx, y_idx);
            int &size_x_y = target_size_grid[x_idx][y_idx];
            if (size_x_y == -1) {
                target_grid[x_idx][y_idx] = new particle_t *[n_particles];
                size_x_y = 0;
            }
            target_grid[x_idx][y_idx][size_x_y] = particles + i;
            ++size_x_y;
        }

        if (next_frame) {
            for (int i = 0; i < n_x; ++i) {
                for (int j = 0; j < n_y; ++j) {
                    size_grid[i][j] = min(0, size_grid[i][j]);
                }
            }

            swap_part_grid = part_grid;
            swap_size_grid = size_grid;

            part_grid = target_grid;
            size_grid = target_size_grid;

            next_part_grid = swap_part_grid;
            next_size_grid = swap_size_grid;

        }
    }

protected:
    const int rank;
    const int block_stride;
    const int block_x;
    const int block_y;

    particle_t* NW_arrival_buffer;
    int NW_n_arrivals;

    particle_t* N_arrival_buffer;
    int N_arrivals;

    particle_t* NE_arrival_buffer;
    int NE_arrivals;

    particle_t* W_arrival_buffer;
    int W_n_arrivals;

    particle_t* E_arrival_buffer;
    int E_arrivals;

    particle_t* SW_arrival_buffer;
    int SW_arrivals;

    particle_t* S_arrival_buffer;
    int S_arrivals;

    particle_t* SE_arrival_buffer;
    int SE_arrivals;

    particle_t* NW_departure_buffer;
    int NW_n_departures;

    particle_t* N_departure_buffer;
    int N_n_departures;

    particle_t* NE_departure_buffer;
    int NE_n_departures;

    particle_t* W_departure_buffer;
    int W_n_departures;

    particle_t* E_departure_buffer;
    int E_n_departures;

    particle_t* SW_departure_buffer;
    int SW_n_departures;

    particle_t* S_departure_buffer;
    int S_n_departures;

    particle_t* SE_departure_buffer;
    int SE_n_departures;


};

#endif