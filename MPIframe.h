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
#include <mpi.h>
#include "frame.h"

static inline double min(double x, double y) { return (x <= y ? x : y); }
static inline double max(double x, double y) { return (x <= y ? y : x); }
static inline double abs(double x) { return (x < .0 ? -x : x); }

class MPIFrame{

public:
    MPIFrame(const int _block_stride,
             const int _block_x, const int _block_y,
             const int _n_block_x, const int _n_block_y,
             const int _n_x, const int _n_y,
             particle_t *particles, const int _n_particles) :
            rank(_block_x * _block_stride + block_y),
            block_stride(_block_stride),
            block_x(_block_x),
            block_y(_block_y),
            offset_x(_block_x * _n_x * (sqrt( density * _n_particles) / ((double) _n_x))),
            offset_y(_block_y * _n_y * (sqrt( density * _n_particles) / ((double) _n_y))){

        part_grid = new particle_t***[n_x];
        next_part_grid = new particle_t***[n_x];

        size_grid = new int*[n_x];
        next_size_grid = new int*[n_x];

        int i, j;

        for(i = 0; i < n_x; ++i){
            part_grid[i] = new particle_t**[n_y];
            next_part_grid[i] = new particle_t**[n_y];

            size_grid[i] = new int[n_y];
            next_size_grid[i] = new int[n_y];

            for(j = 0; j < n_y; ++j){
                size_grid[i][j] = -1;
                next_size_grid[i][j] = -1;
            }

        }

        update_locations(particles, n_particles, false);

        MPI_Datatype PARTICLE;
        MPI_Type_contiguous( 6, MPI_DOUBLE, &PARTICLE );
        MPI_Type_commit( &PARTICLE );

        /*
        * Init Send buffers
        */
        NWs_buffer = new particle_t[_n_particles];
        NWs_n = 0;

        Ns_buffer = new particle_t[_n_particles];
        Ns_n = 0;

        NEs_buffer = new particle_t[_n_particles];
        NEs_n = 0;

        Ws_buffer = new particle_t[_n_particles];
        Ws_n = 0;

        Es_buffer = new particle_t[_n_particles];
        Es_n = 0;

        SWs_buffer = new particle_t[_n_particles];
        SWs_n = 0;

        Ss_buffer = new particle_t[_n_particles];
        Ss_n = 0;

        SEs_buffer = new particle_t[_n_particles];
        SEs_n = 0;

        /*
         * Receive buffers
         */
        NWr_buffer = new particle_t[_n_particles];
        NWr_n = 0;

        Nr_buffer = new particle_t[_n_particles];
        Nr_n = 0;

        NEr_buffer = new particle_t[_n_particles];
        NEr_n = 0;

        Wr_buffer = new particle_t[_n_particles];
        Wr_n = 0;

        Er_buffer = new particle_t[_n_particles];
        Er_n = 0;

        SWr_buffer = new particle_t[_n_particles];
        SWr_n = 0;

        Sr_buffer = new particle_t[_n_particles];
        Sr_n = 0;

        SEr_buffer = new particle_t[_n_particles];
        SEr_n = 0;

    }

    /*
     * ------------------------------------------------
     *
     *              COMMUNICATION METHODS
     *
     * ------------------------------------------------
     *

    /*
     * Send/receive information about particles in NW corner
     */
    inline void comm_locations(int target, int source, bool should_send, bool should_recv, int &msg_idx,
                               int n_to_send, int & n_to_recv, particle_t* send_buffer, particle_t* recv_buffer){
        int msg_tag = ++msg_idx;
        // Send info about size to NW neighbor
        if (should_send) {
            MPI_Send(&(n_to_send), 1, MPI_INT, target, msg_tag, MPI_COMM_WORLD);
        }
        // Receive info from the SE neighbor
        if (should_recv) {
            MPI_Recv(&(n_to_recv), 1, MPI_INT, source, msg_tag, MPI_COMM_WORLD);
        }
        msg_tag = ++msg_idx;
        // Send particles
        if (should_send) {
            MPI_Send(send_buffer, n_to_send, PARTICLE, target, msg_tag, MPI_COMM_WORLD);
        }
        // Receive particles
        if (should_recv) {
            MPI_Recv(recv_buffer, n_to_recv, PARTICLE, source, msg_tag, MPI_COMM_WORLD);
        }
    }

    inline void comm_all_locations(){

        int msg_idx = 0;

        /*
         * Send to NW, receive from SE
         */
        comm_locations(
                (block_x - 1) * block_stride + block_y + 1,
                (block_x + 1) * block_stride + block_y - 1,
                block_x > 0 && block_y < n_block_y - 1,
                block_x < n_block_x - 1 && block_y > 0,
                msg_idx,
                NWs_n,
                SEr_n,
                NWs_buffer,
                SEr_buffer
        );

        /*
         * Send to N, receive from S
         */
        comm_locations(
                block_x * block_stride + block_y + 1,
                block_x * block_stride + block_y - 1,
                block_y < n_block_y - 1,
                block_y > 0,
                msg_idx,
                Ns_n,
                Sr_n,
                Ns_buffer,
                Sr_buffer
        );

        /*
         * Send to NE, receive from SW
         */
        comm_locations(
                (block_x + 1) * block_stride + block_y + 1,
                (block_x - 1) * block_stride + block_y - 1,
                block_x < n_block_x - 1 && block_y < n_block_y - 1,
                block_x > 0 && block_y > 0,
                msg_idx,
                NEs_n,
                SWr_n,
                NEs_buffer,
                SWr_buffer
        );

        /*
         * Send to W, receive from E
         */
        comm_locations(
                (block_x - 1) * block_stride + block_y,
                (block_x + 1) * block_stride + block_y,
                block_x > 0,
                block_x < n_block_x - 1,
                msg_idx,
                Ws_n,
                Er_n,
                Ws_buffer,
                Er_buffer
        );

        /*
         * Send to E, receive from W
         */
        comm_locations(
                (block_x + 1) * block_stride + block_y,
                (block_x - 1) * block_stride + block_y,
                block_x < n_block_x - 1,
                block_x > 0,
                msg_idx,
                Es_n,
                Wr_n,
                Es_buffer,
                Wr_buffer
        );

        /*
         * Send to SW, receive from NE
         */
        comm_locations(
                (block_x - 1) * block_stride + block_y - 1,
                (block_x + 1) * block_stride + block_y + 1,
                block_x > 0 && block_y > 0,
                block_x < n_block_x - 1 && block_y < n_block_y - 1,
                msg_idx,
                SWs_n,
                NEr_n,
                SWs_buffer,
                NEr_buffer
        );

        /*
         * Send to S, receive from N
         */
        comm_locations(
                block_x * block_stride + block_y - 1,
                block_x * block_stride + block_y + 1,
                block_x > 0,
                block_y < n_block_y - 1,
                msg_idx,
                Ss_n,
                Nr_n,
                Ss_buffer,
                Nr_buffer
        );

        /*
         * Send to SE, receive from NW
         */
        comm_locations(
                (block_x + 1) * block_stride + block_y - 1,
                (block_x - 1) * block_stride + block_y + 1,
                block_x < n_block_x - 1 && block_y > 0,
                block_x > 0 && block_y < n_block_y - 1,
                msg_idx,
                SEs_n,
                NWr_n,
                SEs_buffer,
                NWr_buffer
        );

    }


    /**
     * ----------------------------------------------------------
     *
     *                  SIMULATION METHODS
     *
     * ----------------------------------------------------------
     */

    inline void apply_forces(particle_t& part,
                             double *dmin,
                             double *davg,
                             int *navg){

        int x_idx, y_idx, size_x_y;
        particle_t** ptr;

        /*
         * Interaction with particles in the same cell
         */
        get_idx(part.x, part.y, x_idx, y_idx);
        size_x_y = size_grid[x_idx][y_idx];
        ptr = part_grid[x_idx][y_idx];
        for(int i = 0; i < size_x_y; ++i){
            apply_force(part, **(ptr++), dmin, davg, navg);
        }

        // South west neighboring cell
        if(x_idx > 0 && y_idx > 0){
            size_x_y = size_grid[x_idx - 1][y_idx - 1];
            ptr = part_grid[x_idx - 1][y_idx - 1];
            for(int i = 0; i < size_x_y; ++i){
                apply_force(part, **(ptr++), dmin, davg, navg);
            }
        }else if(block_x > 0 && block_y > 0){
            for(int i = 0; i < SW_n; ++i){
                apply_force(part, SW_buffer[i], dmin, davg, navg);
            }
        }

        // South neighboring cell
        if(y_idx > 0){
            size_x_y = size_grid[x_idx][y_idx - 1];
            ptr = part_grid[x_idx][y_idx - 1];
            for(int i = 0; i < size_x_y; ++i){
                apply_force(part, **(ptr++), dmin, davg, navg);
            }
        }else if(block_y > 0){
            for(int i = 0; i < S_n; ++i){
                apply_force(part, S_buffer[i], dmin, davg, navg);
            }
        }

        // Shout east neighboring cell
        if(x_idx < n_x - 1 && y_idx > 0){
            size_x_y = size_grid[x_idx + 1][y_idx - 1];
            ptr = part_grid[x_idx + 1][y_idx - 1];
            for(int i = 0; i < size_x_y; ++i){
                apply_force(part, **(ptr++), dmin, davg, navg);
            }
        }else if(block_x < n_block_x - 1 && block_y > 0){
            for(int i = 0; i < SE_n; ++i){
                apply_force(part, SE_buffer[i], dmin, davg, navg);
            }
        }

        // West neighboring cell
        if(x_idx > 0){
            size_x_y = size_grid[x_idx - 1][y_idx];
            ptr = part_grid[x_idx - 1][y_idx];
            for(int i = 0; i < size_x_y; ++i){
                apply_force(part, **(ptr++), dmin, davg, navg);
            }
        }else if(block_x > 0){
            for(int i = 0; i < W_n; ++i){
                apply_force(part, W_buffer[i], dmin, davg, navg);
            }
        }

        // East neighboring cell
        if(x_idx < n_x - 1){
            size_x_y = size_grid[x_idx + 1][y_idx];
            ptr = part_grid[x_idx + 1][y_idx];
            for(int i = 0; i < size_x_y; ++i){
                apply_force(part, **(ptr++), dmin, davg, navg);
            }
        }else if(block_x < n_block_x - 1){
            for(int i = 0; i < E_n; ++i){
                apply_force(part, E_buffer[i], dmin, davg, navg);
            }
        }

        // North west neighboring cell
        if(x_idx > 0 && y_idx < n_y - 1){
            size_x_y = size_grid[x_idx - 1][y_idx + 1];
            ptr = part_grid[x_idx - 1][y_idx + 1];
            for(int i = 0; i < size_x_y; ++i){
                apply_force(part, **(ptr++), dmin, davg, navg);
            }
        }else if(block_x > 0 && block_y < n_block_y - 1){
            for(int i = 0; i < NW_n; ++i){
                apply_force(part, NW_buffer[i], dmin, davg, navg);
            }
        }

        // North neighboring cell
        if(y_idx < n_y - 1){
            size_x_y = size_grid[x_idx][y_idx + 1];
            ptr = part_grid[x_idx][y_idx + 1];
            for(int i = 0; i < size_x_y; ++i){
                apply_force(part, **(ptr++), dmin, davg, navg);
            }
        }else if(block_y < n_block_y - 1){
            for(int i = 0; i < N_n; ++i){
                apply_force(part, N_buffer[i], dmin, davg, navg);
            }
        }

        // North east neighboring cell
        if(x_idx < n_x - 1 && y_idx < n_y - 1){
            size_x_y = size_grid[x_idx + 1][y_idx + 1];
            ptr = part_grid[x_idx + 1][y_idx + 1];
            for(int i = 0; i < size_x_y; ++i){
                apply_force(part, **(ptr++), dmin, davg, navg);
            }
        }else if(block_x < n_block_x - 1 && block_y < n_block_y - 1){
            for(int i = 0; i < NE_buffer; ++i){
                apply_force(part, NE_buffer[i], dmin, davg, navg);
            }
        }

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
    const double delta_x;
    const double delta_y;
    const double size;
    const int n_x;
    const int n_y;
    const int n_particles;

    particle_t**** part_grid;
    int** size_grid;

    particle_t**** next_part_grid;
    int** next_size_grid;

    particle_t**** swap_part_grid;
    int** swap_size_grid;

    const int rank;
    const int block_stride;
    const int block_x;
    const int block_y;
    const int n_block_x;
    const int n_block_y;
    const double x_offset;
    const double y_offset;

    const std::vector<particle_t> particles;
    const std::vector<int> to_remove;

    /*
     * -------------------------------------------
     *
     *              COMMUNICATION BUFFERS
     *
     * -------------------------------------------
     */

    /*
     * Send buffers
     */
    particle_t* NWs_buffer;
    int NWs_n;

    particle_t* Ns_buffer;
    int Ns_n;

    particle_t* NEs_buffer;
    int NEs_n;

    particle_t* Ws_buffer;
    int Ws_n;

    particle_t* Es_buffer;
    int Es_n;

    particle_t* SWs_buffer;
    int SWs_n;

    particle_t* Ss_buffer;
    int Ss_n;

    particle_t* SEs_buffer;
    int SEs_n;

    /*
     * Receive buffers
     */
    particle_t* NWr_buffer;
    int NWr_n;

    particle_t* Nr_buffer;
    int Nr_n;

    particle_t* NEr_buffer;
    int NEr_n;

    particle_t* Wr_buffer;
    int Wr_n;

    particle_t* Er_buffer;
    int Er_n;

    particle_t* SWr_buffer;
    int SWr_n;

    particle_t* Sr_buffer;
    int Sr_n;

    particle_t* SEr_buffer;
    int SEr_n;

    /**
     * ---------------------------------------
     *
     *              PROTECTED METHODS
     *
     * ---------------------------------------
     */

    /**
    * Look up the indices of that location.
    * IN: x, y. OUT: x_idx, y_idx
    */
    inline void get_idx(const double &x, const double &y, int &x_idx, int &y_idx){
        x_idx = (int) (x / delta_x);
        y_idx = (int) (y / delta_y);
    }

    /**
     * Add the content of a receive buffer to the list of particles
     */
    inline void add_content_to_mem(particle_t* buff, const int n_p){
        for(int i = 0; i < n_p; ++i){
            particles.push_back(buff);
        }
    }

};

#endif