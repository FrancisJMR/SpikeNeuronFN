/** Spiking Neuron Neural Simulator
* Written by Francis Jeanson 2012
*/
#include "network.h"

/* Neuron Model */
#define vi -1.5  /* initial value of v in Volt */
#define recovery -0.375
#define input 0.3 // lower bound: 0.3 upper bound: 1.3 (saturation)

#define fn_p 0.08
#define fn_a 0.7
#define fn_b 0.8

#define SPIKE_DETECT_THRESHOLD 1.0
#define SPIKE_DETECT_WIDTH 0.005 // time width of a spike. Used to evaluate spike times.

FILE *output_file;

float min(float a, float b);
float max(float a, float b);
void printout(float t, float v);
void saveout(float t, float v);
void save_spiketimes(float st, int n_id, char *spiketimes_filename);
void init_neurons(Neuron *Neurons, float inhib_perc);
void read_input_file(FILE *fp, Neuron *Neurons, int cellstart, int cellstop);
void set_poisson_input(int t, Neuron *Neurons, int nid, float spike_interval);
float fn_v(float t, float v, float r, float i);
float fn_r(float t, float v, float r);
void propagate_signals(Connection Connections[][NET_SIZE], Neuron *NSource, Neuron *NTarget);
void integrate_neurons(int discrete_ts, Neuron *Neurons, FILE *fp, int cellstart, int cellstop, int stim_on);
void perform_stdp(int discrete_ts, Connection Connections[][NET_SIZE], Neuron *NSource, Neuron *NTarget);
int store_neurone_activity(int discrete_ts, Neuron *Neurons, char *target_spiketimes_filename);
void detect_pattern(int discrete_ts, Neuron *Neurons, float detect_offset, int cellstart, int cellstop, int *pattern_counts);