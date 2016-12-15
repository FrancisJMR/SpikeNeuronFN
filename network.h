/** Spiking Neuron Neural Simulator
* Written by Francis Jeanson 2012
*/

#include <stdio.h> 
#include <math.h>
#include <time.h>
#include <stdlib.h>

#include <gsl/gsl_randist.h>
#include <gsl/gsl_cdf.h>

/* Code Utils */
#define TRUE 1
#define FALSE 0

/* Simulation Defs */
#define TIME_SCALE_FACTOR 1
#define SIMULTIME 1.0   // Simulation Duration
#define STEPSIZE 0.0001 // 0.00001
#define tau 0.4         // 0.04 /* integration constant */
#define TOTAL_STEPS (int)(SIMULTIME/STEPSIZE)

/* Network Model */
#define NET_SIZE 24
#define SHEET_WIDTH 5
#define SHEET_HEIGHT 5

#define INHIB_PERCENT 0.0               // percentage of inhibitory nodes
#define DEFAULT_THRESHOLD 12            // Number of necessary synchronous spikes to reach threshold on target
#define WEIGHT (0.12/DEFAULT_THRESHOLD) // Default connection weight // 0.12 minimum weight for excitation by one afferent

#define CONNECT_RANDOM TRUE
#define CONNECTION_PROB 1.0              
#define MIN_DELAY (int)(0.01/STEPSIZE)   // in slots: seconds/resolution -- with tau 0.4 STEPSIZE 0.0001: relative refractory 92ts,  absolute refractory 88ts
#define MAX_DELAY (int)(0.03/STEPSIZE)   // in slots: seconds/resolution

#define CONNECT_TOPOGRAPHIC FALSE
#define CONNECT_EXCIT_PROB 0.05         // probability that an excitatory connection is made
#define CONNECT_INHIB_PROB 0.125        // probability that an inhibitory connection is made
#define DELAY_EXCIT_DECAY_CONST 140     // Should be proportional to SHEET_WIDTH and SHEET_HEIGHT
#define DELAY_INHIB_DECAY_CONST DELAY_EXCIT_DECAY_CONST//(3*DELAY_EXCIT_DECAY_CONST);
#define MIN_TOPO_DELAY (int)(0.0005/STEPSIZE)
#define MAX_TOPO_DELAY (int)(0.01/STEPSIZE)

#define CONNECT_SMALL_WORLD FALSE   // connect network using preferential attachement
#define SMALL_WORLD_CONSTANT 0.02   // a smaller constant means smaller std

#define STDP_WINDOW (int)(0.03/STEPSIZE)
#define STDP_WEIGHT_FACTOR 0.03
#define MAX_WEIGHT 0.08

typedef int boolean;

typedef struct {
    float v;
    float r;
    float I;
    float xpos;
    float ypos;
    float axon[MAX_DELAY];
    int is_excit;
    int dpush; // the circular index to push output voltage on this axon
    float v_hist[STDP_WINDOW];
    int I_duration; //used to determine if input is long enough
} Neuron;

typedef struct {
    float weight;
    boolean are_connected;
    int delay;
} Connection;

const gsl_rng * r; // libgsl random var

void save_nodeconnections(int node_id, Connection Connections[][NET_SIZE], Neuron *Neurons, char *file_name);
double std_int(int *value_list, int length);
double std_float(float *value_list, int length);
int on_sheet_area(int i, int x, int y, int rad);
int get_node_x(int i);
int get_node_y(int i);
void generate_network(Connection Connections[][NET_SIZE], Neuron *Neurons, int mindelay, int maxdelay, float connect_prob, float weight, char *file_name);
int load_network_file(char *file_name, Connection Connections[][NET_SIZE], Neuron *Neurons);