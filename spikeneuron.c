/** Spiking Neuron Neural Simulator
* Written by Francis Jeanson 2012
*/

#include "spikeneuron.h"


float min(float a, float b) {
    return  a < b ? a : b;
}

float max(float a, float b) {
    return  a < b ? b : a;
}

void printout(float t, float v) {
    printf("%f\t%f\n",t,v); 
} 

void saveout(float t, float v) {
    fprintf(output_file, "%f\t%f\n", t/TIME_SCALE_FACTOR, v);
}

void save_spiketimes(float st, int n_id, char *spiketimes_filename) {
    FILE *spiketimes_file;
    spiketimes_file = fopen(spiketimes_filename, "a");
    
    fprintf(spiketimes_file, "%i\t%f\n", n_id, st/TIME_SCALE_FACTOR);
    
    fclose(spiketimes_file);
}

void init_neurons(Neuron *Neurons, float inhib_perc) {
    int i,j;
    for (i=0; i<NET_SIZE; i++) {
        // initialize the 'terms' that will differentiate over time.
        Neurons[i].v = vi; // init membrane voltage
        Neurons[i].r = recovery; // recovery factor
        Neurons[i].I = 0.0;
        Neurons[i].I_duration = 0;
        Neurons[i].xpos = 0.0;
        Neurons[i].ypos = 0.0;
        for (j=0; j<MAX_DELAY; j++) {
            Neurons[i].axon[j] = 0.0;
        }
        Neurons[i].dpush = MAX_DELAY-1;
        Neurons[i].is_excit = drand48() < inhib_perc ? FALSE : TRUE;
        for (j=0; j<STDP_WINDOW; j++) {
            Neurons[i].v_hist[j] = 0.0;
        }
    }
}

void read_input_file(FILE *fp, Neuron *Neurons, int cellstart, int cellstop){
    int input_network_size = 64;
    float input_vals[input_network_size];
    int i, n;
    int line_length = 15*input_network_size;
    char line[line_length];
    char field[15];
    
    if (fgets(line, line_length, fp) != NULL) {
        const char *c = line;
        i = 0;
        while ( sscanf(c, "%9[^\t]%n", field, &n) == 1 ) {
            input_vals[i] = atof(field);
            c+=n;
            if (*c != '\t') break;
            ++c; // skip delimiter
            i++;
        }
        
        // Add file value to neurone input
        for (i=0; i<NET_SIZE; i++) {
            
            // if cellstop larger then read-in those cells only
            if (cellstart < cellstop) {
                
                if (input_vals[i] > 0 && cellstart <= i && i < cellstop) {
                    Neurons[i].I += input_vals[i];
                }
                
            } else if (input_vals[i] > 0) {
                Neurons[i].I += input_vals[i];
            }
            
            // limit the new input
            Neurons[i].I = min(Neurons[i].I, 1.3);
        }
        //printf("net vals %g %g %g %g %g\n", input_vals[0], input_vals[1], input_vals[2], input_vals[3], input_vals[4]);
    }
}


void set_poisson_input(int t, Neuron *Neurons, int nid, float spike_interval){

    int si = spike_interval/STEPSIZE; // convert the mean period
    int ps = round(-log(drand48()) * si); // compute poisson probability
    // NOTE: ps is actually some poisson distributed time within the interval si
    if (ps == si/2 && Neurons[nid].I_duration < t) {
        Neurons[nid].I = input;
        Neurons[nid].I_duration = t+(SPIKE_DETECT_WIDTH/STEPSIZE);
    }
}

/** FitHug-Nagumo model
 * Vs(j) = Vs(j) + tau * ( Vs(j) - (Vs(j).^3)/3 - Rs(j) + Is(j) );
 * Rs(j) = Rs(j) + tau * p*(1.0*Vs(j) + a - b * Rs(j));
 */
float fn_v(float t, float v, float r, float i) {
    return v + tau * (v - pow(v, 3)/3 - r + i);
}

float fn_r(float t, float v, float r) {
    return r + tau * fn_p * (v + fn_a - fn_b * r);
}

/** Do all the network signal propagation and neuron model integration
 * Neurons: the 1D list of neurons
 * Connections: the 2D list of connections from one neuron to another.
 */
void propagate_signals(Connection Connections[][NET_SIZE], Neuron *NSource, Neuron *NTarget) {
    // propagate connection firing values
    int i, j, dpull;
    
    // transmit activity from pre-synaptic to post-synaptic cells when spike arrives
    for (j=0; j<NET_SIZE; j++) {
        for (i=0; i<NET_SIZE; i++) {
            if (Connections[i][j].are_connected) {
                // dpull: the circular index to pull voltage from this axon according to delay
                dpull = (NSource[i].dpush+Connections[i][j].delay)%MAX_DELAY;
                if (NSource[i].is_excit) NTarget[j].I += Connections[i][j].weight * max(NSource[i].axon[dpull], 0.0);
                else NTarget[j].I -= Connections[i][j].weight * max(NSource[i].axon[dpull], 0.0);
            }
            // limit total upper bound input to cell
            NTarget[j].I = min(NTarget[j].I, 1.3);
        }
    }
}

//void integrate_neurons(float ts, Neuron *Neurons, int input_patt[PATT_TIME][NET_SIZE]) {
void integrate_neurons(int discrete_ts, Neuron *Neurons, FILE *fp, int cellstart, int cellstop, int stim_on) {
    
    int i;
    
    // store history before integration
    float ts = discrete_ts*STEPSIZE;
    int spike_count = 0;
    int hist_index = discrete_ts%STDP_WINDOW;
    
    for (i=0; i<NET_SIZE; i++) {
        Neurons[i].v_hist[hist_index] = Neurons[i].v;
    }
    
    // get FN spikes from extenal file
    //if (fp != NULL) read_input_file(fp, Neurons, cellstart, cellstop);
    if (stim_on && cellstart < cellstop && cellstop <= NET_SIZE) {
        for (i=cellstart; i<cellstop; i++){
            set_poisson_input(discrete_ts, Neurons, i, 0.005);
        }
    }
    
    // integrate the voltage and conductance of each neurone
    for (i=0; i<NET_SIZE; i++) {
        Neurons[i].v = fn_v(ts, Neurons[i].v, Neurons[i].r, Neurons[i].I);
        Neurons[i].r = fn_r(ts, Neurons[i].v, Neurons[i].r);
	}    
    
    // store new voltage on axon and update dpush index
    for (i=0; i<NET_SIZE; i++) {
            Neurons[i].axon[Neurons[i].dpush] = Neurons[i].v;
            Neurons[i].dpush--;
            if (Neurons[i].dpush < 0) Neurons[i].dpush = MAX_DELAY-1;
    }
    
    // reset input to neurones
    for (i=0; i<NET_SIZE; i++) {
        if (Neurons[i].I_duration < discrete_ts){
            Neurons[i].I = 0.0;
        }
    }
        
}

void perform_stdp(int discrete_ts, Connection Connections[][NET_SIZE], Neuron *NSource, Neuron *NTarget) {
    int tmp;
    
    int i,j;
    float ts = discrete_ts*STEPSIZE;
    int k, spike_found;
    int spiketime_diff; // the discrete pre or post synaptic spiketime
    float new_weight;
    int hist_offset, hist_offset_prev;
    int hist_index     = discrete_ts%STDP_WINDOW;
    int hist_prev_step = (discrete_ts-1)%STDP_WINDOW;
    for (i=0; i<NET_SIZE; i++) {
        
        // only apply STDP if current cell just spiked & enough time has passed to look back
        if (NTarget[i].v >= SPIKE_DETECT_THRESHOLD && 
            NTarget[i].v_hist[hist_prev_step] < SPIKE_DETECT_THRESHOLD && 
            discrete_ts > STDP_WINDOW) {
            
            //printf("plastic!\n");
            
            // check all connected cells to see if spiked within STDP_WINDOW
            for (j=0; j<NET_SIZE; j++) {
                
                // consider only pre-synaptic case and post-synaptic cases
                if (Connections[i][j].are_connected || Connections[j][i].are_connected) {
                    
                    // scan other neurone history for last spike
                    k=0;
                    spike_found = FALSE;
                    while(k<STDP_WINDOW-1 && spike_found == FALSE){
                        // determine position to fetch activation history on circular history list
                        if (hist_index - k < 0) hist_offset = STDP_WINDOW + hist_index - k;
                        else hist_offset = hist_index - k;
                        if (hist_index - k - 1 < 0) hist_offset_prev = STDP_WINDOW + hist_index - k - 1;
                        else hist_offset_prev = hist_index - k - 1;
                        
                        // determine if other neurone spiked or not
                        if (NSource[j].v_hist[hist_offset] >= SPIKE_DETECT_THRESHOLD && NSource[j].v_hist[hist_offset_prev] < SPIKE_DETECT_THRESHOLD) {
                            spike_found = TRUE;
                            spiketime_diff = k;
                        }
                        k++;
                    }
                                         
                    if (spike_found) {
                        // if j is target and spiked just before afferent i then depress
                        if (Connections[i][j].are_connected) {
                            new_weight = max(Connections[i][j].weight - (STDP_WEIGHT_FACTOR * STDP_WINDOW * 1.0)/spiketime_diff, 0.0);
                            //printf("dW %f\n", new_weight);
                            Connections[i][j].weight = new_weight;
                        }
                        // if j is afferent and spiked just before target i then reinforce
                        if (Connections[j][i].are_connected) {
                            new_weight = min(Connections[j][i].weight + (STDP_WEIGHT_FACTOR * STDP_WINDOW * 1.0)/spiketime_diff, MAX_WEIGHT);
                            //printf("dW %f\n", new_weight);
                            Connections[j][i].weight = new_weight;
                        }
						
                    }
                }
            }
        }
    }
}

int store_neurone_activity(int discrete_ts, Neuron *Neurons, char *target_spiketimes_filename) {
    
    // Store spike times
    int i;
    float ts = discrete_ts*STEPSIZE;
    int spike_count = 0;
    int hist_index = discrete_ts%STDP_WINDOW;
    int prev_step     = (discrete_ts-1)%STDP_WINDOW;
    for (i=0; i<NET_SIZE; i++) {
        if (Neurons[i].v >= SPIKE_DETECT_THRESHOLD && Neurons[i].v_hist[prev_step] < SPIKE_DETECT_THRESHOLD) {
            //printf("prev v %f, current v %f\n", Neurons[i].v_hist[prev_step], Neurons[i].v);
            save_spiketimes(ts, i, target_spiketimes_filename);
            spike_count++;
        }
    }
    return spike_count;
}

/** detect_pattern
 Add a tick in pattern_counts each time a cell spikes.
 Start counting only after detect_offset time has passed.
 */
void detect_pattern(int discrete_ts,  Neuron *Neurons, float detect_offset, int cellstart, int cellstop, int *pattern_counts) {

    int i;
    float ts = discrete_ts*STEPSIZE;
    int spike_count = 0;
    int hist_index = discrete_ts%STDP_WINDOW;
    int prev_step     = (discrete_ts-1)%STDP_WINDOW;
        
    for (i=0; i<NET_SIZE; i++) {
        if (i >= cellstart && i<cellstop && Neurons[i].v >= SPIKE_DETECT_THRESHOLD && Neurons[i].v_hist[prev_step] < SPIKE_DETECT_THRESHOLD) {
            if (ts >= detect_offset) pattern_counts[i-cellstart] += 1;
        }
    }
}