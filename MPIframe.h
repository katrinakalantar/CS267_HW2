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

        mem = new particle_t[n_particles];
        n_particles = 0;

        next_mem = new particle_t[n_particles];
        next_n_particles = 0;

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
        * Init send buffers for locations
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
         * Init receive buffers for locations
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

        /*
         * Init send buffers for particles
         */
        pNWs_buffer = new particle_t[_n_particles];
        pNWs_n = 0;

        pNs_buffer = new particle_t[_n_particles];
        pNs_n = 0;

        pNEs_buffer = new particle_t[_n_particles];
        pNEs_n = 0;

        pWs_buffer = new particle_t[_n_particles];
        pWs_n = 0;

        pEs_buffer = new particle_t[_n_particles];
        pEs_n = 0;

        pSWs_buffer = new particle_t[_n_particles];
        pSWs_n = 0;

        pSs_buffer = new particle_t[_n_particles];
        pSs_n = 0;

        pSEs_buffer = new particle_t[_n_particles];
        pSEs_n = 0;

        /*
         * Init receive buffers for particles
         */
        pNWr_buffer = new particle_t[_n_particles];
        pNWr_n = 0;

        pNr_buffer = new particle_t[_n_particles];
        pNr_n = 0;

        pNEr_buffer = new particle_t[_n_particles];
        pNEr_n = 0;

        pWr_buffer = new particle_t[_n_particles];
        pWr_n = 0;

        pEr_buffer = new particle_t[_n_particles];
        pEr_n = 0;

        pSWr_buffer = new particle_t[_n_particles];
        pSWr_n = 0;

        pSr_buffer = new particle_t[_n_particles];
        pSr_n = 0;

        pSEr_buffer = new particle_t[_n_particles];
        pSEr_n = 0;

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
    inline void comm_buffer(int target, int source, bool should_send, bool should_recv, int &msg_idx,
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

    /**
     * Communicate particle locations so that
     * forces can be computed accurately for particles
     * near the border.
     */
    inline void comm_all_locations(){

        int msg_idx = 0;

        /*
         * Send to NW, receive from SE
         */
        comm_buffer(
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
        comm_buffer(
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
        comm_buffer(
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
        comm_buffer(
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
        comm_buffer(
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
        comm_buffer(
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
        comm_buffer(
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
        comm_buffer(
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
     * Communicate particles so that they are handled
     * by the processor magaging the region they belong to.
     */
    inline void comm_all_particles(){

        int msg_idx = 0;

        /*
         * Send to NW, receive from SE
         */
        comm_buffer(
                (block_x - 1) * block_stride + block_y + 1,
                (block_x + 1) * block_stride + block_y - 1,
                block_x > 0 && block_y < n_block_y - 1,
                block_x < n_block_x - 1 && block_y > 0,
                msg_idx,
                pNWs_n,
                pSEr_n,
                pNWs_buffer,
                pSEr_buffer
        );

        /*
         * Send to N, receive from S
         */
        comm_buffer(
                block_x * block_stride + block_y + 1,
                block_x * block_stride + block_y - 1,
                block_y < n_block_y - 1,
                block_y > 0,
                msg_idx,
                pNs_n,
                pSr_n,
                pNs_buffer,
                pSr_buffer
        );

        /*
         * Send to NE, receive from SW
         */
        comm_buffer(
                (block_x + 1) * block_stride + block_y + 1,
                (block_x - 1) * block_stride + block_y - 1,
                block_x < n_block_x - 1 && block_y < n_block_y - 1,
                block_x > 0 && block_y > 0,
                msg_idx,
                pNEs_n,
                pSWr_n,
                pNEs_buffer,
                pSWr_buffer
        );

        /*
         * Send to W, receive from E
         */
        comm_buffer(
                (block_x - 1) * block_stride + block_y,
                (block_x + 1) * block_stride + block_y,
                block_x > 0,
                block_x < n_block_x - 1,
                msg_idx,
                pWs_n,
                pEr_n,
                pWs_buffer,
                pEr_buffer
        );

        /*
         * Send to E, receive from W
         */
        comm_buffer(
                (block_x + 1) * block_stride + block_y,
                (block_x - 1) * block_stride + block_y,
                block_x < n_block_x - 1,
                block_x > 0,
                msg_idx,
                pEs_n,
                pWr_n,
                pEs_buffer,
                pWr_buffer
        );

        /*
         * Send to SW, receive from NE
         */
        comm_buffer(
                (block_x - 1) * block_stride + block_y - 1,
                (block_x + 1) * block_stride + block_y + 1,
                block_x > 0 && block_y > 0,
                block_x < n_block_x - 1 && block_y < n_block_y - 1,
                msg_idx,
                pSWs_n,
                pNEr_n,
                pSWs_buffer,
                pNEr_buffer
        );

        /*
         * Send to S, receive from N
         */
        comm_buffer(
                block_x * block_stride + block_y - 1,
                block_x * block_stride + block_y + 1,
                block_x > 0,
                block_y < n_block_y - 1,
                msg_idx,
                pSs_n,
                pNr_n,
                pSs_buffer,
                pNr_buffer
        );

        /*
         * Send to SE, receive from NW
         */
        comm_buffer(
                (block_x + 1) * block_stride + block_y - 1,
                (block_x - 1) * block_stride + block_y + 1,
                block_x < n_block_x - 1 && block_y > 0,
                block_x > 0 && block_y < n_block_y - 1,
                msg_idx,
                pSEs_n,
                pNWr_n,
                pSEs_buffer,
                pNWr_buffer
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

        comm_all_locations();

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
            for(int i = 0; i < SWr_n; ++i){
                apply_force(part, SWr_buffer[i], dmin, davg, navg);
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
            for(int i = 0; i < Sr_n; ++i){
                apply_force(part, Sr_buffer[i], dmin, davg, navg);
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
            for(int i = 0; i < SEr_n; ++i){
                apply_force(part, SEr_buffer[i], dmin, davg, navg);
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
            for(int i = 0; i < Wr_n; ++i){
                apply_force(part, Wr_buffer[i], dmin, davg, navg);
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
            for(int i = 0; i < Er_n; ++i){
                apply_force(part, Er_buffer[i], dmin, davg, navg);
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
            for(int i = 0; i < NWr_n; ++i){
                apply_force(part, NWr_buffer[i], dmin, davg, navg);
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
            for(int i = 0; i < Nr_n; ++i){
                apply_force(part, Nr_buffer[i], dmin, davg, navg);
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
            for(int i = 0; i < NEr_buffer; ++i){
                apply_force(part, NEr_buffer[i], dmin, davg, navg);
            }
        }

    }

    /**
     * Update the locations of the particles and
     * figure out which need to be sent to neighbors
     * for managing and which are in a bordering zone.
     */
    inline void update_locations(particle_t* particles, const int n_particles, bool next_frame = true) {

        clear_loc_buffers();

        particle_t ****target_grid;
        int **target_size_grid;
        particle_t* target_mem;
        int *target_n_particles;
        if (next_frame) {
            target_grid = next_part_grid;
            target_size_grid = next_size_grid;
            target_mem = next_mem;
            target_n_particles = &next_n_particles;
        } else {
            target_grid = part_grid;
            target_size_grid = size_grid;
            target_mem = mem;
            target_n_particles = &n_particles;
        }

        int x_idx, y_idx;
        for (int i = 0; i < n_particles; ++i) {

            const particle_t &part = particles[i];

            get_idx(part.x, part.y, x_idx, y_idx);

            if(y_idx < 0 && x_idx < 0){
                // This particle will move to SW neighbor
                pSWs_buffer[pSWs_n++] = part;
            }else if(y_idx == 0 && x_idx == 0){
                // This particle is now in SW corner
                SWs_buffer[SWs_n++] = part;
            }

            if(y_idx < 0){
                // This particle will move to S neighbor
                pSs_buffer[pSs_n++] = part;
            }else if(y_idx == 0){
                // This particle is now on S border (watch out for D. Trump...)
                Ss_buffer[Ss_n++] = part;
            }

            if(x_idx >= n_x && y_idx < 0){
                // This particle will move to SE neighbor
                pSEs_buffer[pSEs_n++] = part;
            }else if(x_idx == n_x - 1 && y_idx == 0){
                // This particle is now in SE corner
                SEs_buffer[SEs_n++] = part;
            }

            if(x_idx < 0){
                // This particle will move to W neighbor
                pWs_buffer[pWs_n++] = part;
            }else if(x_idx == 0){
                // This particle is now on W border
                Ws_buffer[Ws_n++] = part;
            }

            if(x_idx >= n_x){
                // This particle will move to E neighbor
                pEs_buffer[pEs_n++] = part;
            }else if(x_idx == n_x - 1){
                // This particle is now on E border
                Es_buffer[Es_n++] = part;
            }

            if(x_idx < 0 && y_idx >= n_y){
                // This particle will move to NW neighbor
                pNWs_buffer[pNWs_n++] = part;
            }else if(x_idx == 0 && y_idx == n_y){
                // This particle is now on NW corner
                NWs_buffer[NWs_n++] = part;
            }

            if(y_idx >= n_y){
                // This particle will move to N neighbor
                pNs_buffer[pNs_n++] = part;
            }else if(y_idx == n_y - 1){
                // This particle is now on North border
                Ns_buffer[Ns_n++] = part;
            }

            if(x_idx >= n_x && y_idx >= n_y){
                // This particle will move to NE neighbor
                pNEs_buffer[pNEs_n++] = part;
            }else if(x_idx == n_x - 1 && y_idx == n_y - 1){
                // This particle is now in NE corner
                NEs_buffer[NEs_n++] = part;
            }

            if(x_idx >= 0 && y_idx >= 0 && x_idx < n_x  && y_idx < n_x) {
                int &size_x_y = target_size_grid[x_idx][y_idx];
                if (size_x_y == -1) {
                    target_grid[x_idx][y_idx] = new particle_t *[n_particles];
                    size_x_y = 0;
                }
                target_grid[x_idx][y_idx][size_x_y] = particles + i;
                ++size_x_y;
            }

        }

        if(next_frame) {
            /**
             * Get particles from neighbors that have moved
             * within the region this processor manages
             * and add them to the memory and location grid
             */
            /**
             * Add particles from NW neighbor
             */
            add_particles_from_buffer(
                    pNWr_buffer,
                    pNWr_n,
                    target_grid,
                    target_size_grid,
                    target_mem,
                    *target_n_particles
            );

            /**
             * Add particles from N neighbor
             */
            add_particles_from_buffer(
                    pNr_buffer,
                    pNr_n,
                    target_grid,
                    target_size_grid,
                    target_mem,
                    *target_n_particles
            );

            /**
             * Add particles from NE neighbor
             */
            add_particles_from_buffer(
                    pNEr_buffer,
                    pNEr_n,
                    target_grid,
                    target_size_grid,
                    target_mem,
                    *target_n_particles
            );

            /**
             * Add particles form E neighbor
             */
            add_particles_from_buffer(
                    pEr_buffer,
                    pEr_n,
                    target_grid,
                    target_size_grid,
                    target_mem,
                    *target_n_particles
            );

            /**
             * Add particles from W neighbor
             */
            add_particles_from_buffer(
                    pWr_buffer,
                    pWr_n,
                    target_grid,
                    target_size_grid,
                    target_mem,
                    *target_n_particles
            );

            /**
             * Add particles from SW neighbor
             */
            add_particles_from_buffer(
                    pSWr_buffer,
                    pSWr_n,
                    target_grid,
                    target_size_grid,
                    target_mem,
                    *target_n_particles
            );

            /**
             * Add particles from S neighbor
             */
            add_particles_from_buffer(
                    pSr_buffer,
                    pSr_n,
                    target_grid,
                    target_size_grid,
                    target_mem,
                    *target_n_particles
            );

            /**
             * Add particles from SE neighbor
             */
            add_particles_from_buffer(
                    pSEr_buffer,
                    pSEr_n,
                    target_grid,
                    target_size_grid,
                    target_mem,
                    *target_n_particles
            );

            clear_p_buffers();
        }

        /**
         * Give the location of particles in bordering zones
         * to neighbors.
         */
        comm_all_locations();

        /**
         * Swap frames
         */

        if (next_frame) {

            particle_t**** swap_part_grid;
            int** swap_size_grid;
            particle_t* swap_mem;
            int swap_n_particles;

            for (int i = 0; i < n_x; ++i) {
                for (int j = 0; j < n_y; ++j) {
                    size_grid[i][j] = min(0, size_grid[i][j]);
                }
            }

            swap_part_grid = part_grid;
            swap_size_grid = size_grid;
            swap_mem = mem;
            swap_n_particles = n_particles;

            part_grid = target_grid;
            size_grid = target_size_grid;
            mem = target_mem;
            n_particles = target_n_particles;

            next_part_grid = swap_part_grid;
            next_size_grid = swap_size_grid;
            next_mem = swap_mem;
            next_n_particles = swap_n_particles;

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

    particle_t* mem;
    int n_particles;

    particle_t* next_mem;
    int next_n_particles;

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
     * Send buffers for locations
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
     * Receive buffers for locations
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

    /*
     * Send buffers for particles
     */
    particle_t* pNWs_buffer;
    int pNWs_n;

    particle_t* pNs_buffer;
    int pNs_n;

    particle_t* pNEs_buffer;
    int pNEs_n;

    particle_t* pWs_buffer;
    int pWs_n;

    particle_t* pEs_buffer;
    int pEs_n;

    particle_t* pSWs_buffer;
    int pSWs_n;

    particle_t* pSs_buffer;
    int pSs_n;

    particle_t* pSEs_buffer;
    int pSEs_n;

    /*
     * Receive buffers for particles
     */
    particle_t* pNWr_buffer;
    int pNWr_n;

    particle_t* pNr_buffer;
    int pNr_n;

    particle_t* pNEr_buffer;
    int pNEr_n;

    particle_t* pWr_buffer;
    int pWr_n;

    particle_t* pEr_buffer;
    int pEr_n;

    particle_t* pSWr_buffer;
    int pSWr_n;

    particle_t* pSr_buffer;
    int pSr_n;

    particle_t* pSEr_buffer;
    int pSEr_n;


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
        x_idx = (int) ((x - x_offset) / delta_x);
        y_idx = (int) ((y - y_offset) / delta_y);
    }

    /**
     * Clear location buffers
     */
    inline void clear_loc_buffers(){
        NWs_n   = 0;
        Ns_n    = 0;
        NEs_n   = 0;
        Ws_n    = 0;
        Es_n    = 0;
        SWs_n   = 0;
        Ss_n    = 0;
        SEs_n   = 0;

        NWr_n   = 0;
        Nr_n    = 0;
        NEr_n   = 0;
        Wr_n    = 0;
        Er_n    = 0;
        SWr_n   = 0;
        Sr_n    = 0;
        SEr_n   = 0;
    }

    /**
     * Clear particle buffers
     */
    inline void clear_p_buffers(){
        pNWs_n  = 0;
        pNs_n   = 0;
        pNEs_n  = 0;
        pWs_n   = 0;
        pEs_n   = 0;
        pSWs_n  = 0;
        pSs_n   = 0;
        pSEs_n  = 0;

        pNWr_n  = 0;
        pNr_n   = 0;
        pNEr_n  = 0;
        pWr_n   = 0;
        pEr_n   = 0;
        pSWr_n  = 0;
        pSr_n   = 0;
        pSEr_n  = 0;
    }

    /**
     * Add the content of a receive buffer to the list of particles
     */
    inline void add_particles_from_buffer(
            const particle_t* buffer,
            const int n_buffer,
            particle_t ****target_grid,
            int **target_size_grid,
            particle_t* target_mem,
            int &target_n_particles){

        int x_idx, y_idx;

        for(int i = 0; i < n_buffer; ++i) {

            const particle_t &part = particles[i];

            get_idx(part.x, part.y, x_idx, y_idx);

            target_mem[target_n_particles++] = part;

            assert(x_idx >= 0);
            assert(y_idx >= 0);
            assert(x_idx < n_x);
            assert(y_idx < n_y);

            if(y_idx == 0 && x_idx == 0){
                // This particle is now in SW corner
                SWs_buffer[SWs_n++] = part;
            }

            if(y_idx == 0){
                // This particle is now on S border (watch out for D. Trump...)
                Ss_buffer[Ss_n++] = part;
            }

            if(x_idx == n_x - 1 && y_idx == 0){
                // This particle is now in SE corner
                SEs_buffer[SEs_n++] = part;
            }

            if(x_idx == 0){
                // This particle is now on W border
                Ws_buffer[Ws_n++] = part;
            }

            if(x_idx == n_x - 1){
                // This particle is now on E border
                Es_buffer[Es_n++] = part;
            }

            if(x_idx == 0 && y_idx == n_y){
                // This particle is now on NW corner
                NWs_buffer[NWs_n++] = part;
            }

            if(y_idx == n_y - 1){
                // This particle is now on North border
                Ns_buffer[Ns_n++] = part;
            }

            if(x_idx == n_x - 1 && y_idx == n_y - 1){
                // This particle is now in NE corner
                NEs_buffer[NEs_n++] = part;
            }

            if (x_idx >= 0 && y_idx >= 0 && x_idx < n_x && y_idx < n_x) {
                int &size_x_y = target_size_grid[x_idx][y_idx];
                if (size_x_y == -1) {
                    target_grid[x_idx][y_idx] = new particle_t *[n_particles];
                    size_x_y = 0;
                }
                target_grid[x_idx][y_idx][size_x_y] = particles + i;
                ++size_x_y;
            }

        }

    }

};

#endif