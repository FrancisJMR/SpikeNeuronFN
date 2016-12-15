/** Spiking Neuron Neural Simulator
* Written by Francis Jeanson 2012
*/

#include "network.h"

void save_nodeconnections(int node_id, Connection Connections[][NET_SIZE], Neuron *Neurons, char *file_name) {
    FILE *nodeconnections_file;
    nodeconnections_file = fopen(file_name, "a");
    int i,c;
    char connection_str[50];
    for (i=0; i<NET_SIZE; i++) {
        // c_net[source][target]
        if (Connections[i][node_id].are_connected && nodeconnections_file != NULL) {
            if (Neurons[i].is_excit) c = 1; else c = 0;
            sprintf(connection_str, "%d%s%d%s%d%s%d%s%g%s", node_id," ", i, " ", c ," ",Connections[i][node_id].delay," ", Connections[i][node_id].weight,"\n");
            //sprintf(connection_str, "%d%s%d%s%d%s%d%s", node_id," ", i, " ", c ," ",Connections[i][node_id].delay, "\n");
            fputs ( connection_str, nodeconnections_file);
        }
    }
    
    fclose(nodeconnections_file);
}

/** Compute the standard deviation of a list of ints.
 */
double std_int(int *value_list, int length) {
    
    int i;
    double mean = 0.0;
    double sumdiffsq = 0.0;
    
    for (i=0; i<length; i++) {
        mean += value_list[i]*1.0;
    }
    mean /= length;
    
    
    for (i=0; i<length; i++) {
        sumdiffsq += pow( value_list[i] - mean, 2);
    }
    
    return sqrt(sumdiffsq / length);
    
}

double std_float(float *value_list, int length) {
    
    int i;
    double mean = 0.0;
    double sumdiffsq = 0.0;
    
    for (i=0; i<length; i++) {
        mean += value_list[i]*1.0;
    }
    mean /= length;
    
    
    for (i=0; i<length; i++) {
        sumdiffsq += pow( value_list[i] - mean, 2);
    }
    
    return sqrt(sumdiffsq / length);
    
}

int on_sheet_area(int i, int x, int y, int rad) {
    
    int j,k;    
    for (j=0; j<SHEET_WIDTH; j++) {
        for (k=0; k<SHEET_HEIGHT; k++) {
            if (get_node_x(i) >= x - rad 
                && get_node_x(i) <= x + rad 
                && get_node_y(i) >= y - rad
                && get_node_y(i) <= y + rad)
                return TRUE;
        }
    }
    
    return FALSE;
}

int get_node_x(int i) {
    return i%SHEET_WIDTH;
}

int get_node_y(int i) {
    return i/SHEET_WIDTH;
}

/** void generate_network()
 * Generate a delay network with uniform or preferential attachment
 * for a small world topology.
 */
void generate_network(Connection Connections[][NET_SIZE], Neuron *Neurons, int mindelay, int maxdelay, float connect_prob, float weight, char *file_name) {
    
    int i,j,k,s,t;
    int afferents_arr[NET_SIZE];
    double target_score_arr[NET_SIZE], max;
    
    // init afferents_array and Connections
    for (i=0; i<NET_SIZE; i++) {
        
        afferents_arr[i] = 0;
        
        for (j=0; j<NET_SIZE; j++) {
            Connections[i][j].are_connected = FALSE;
            Connections[i][j].weight = weight;
        }
        
        // init the axon potential
        for (k=0; k<MAX_DELAY; k++) {
            Neurons[i].axon[k] = 0.0;
        }
        
    }
        
    if (CONNECT_SMALL_WORLD) {
        printf("Connecting small world...\n");
        int connections_bin = connect_prob * NET_SIZE * NET_SIZE;
        printf("Connections to make: %i\n", connections_bin);
        for (i=0; i<connections_bin; i++) {
            
            // pick random source
            s = (rand() / ((double)RAND_MAX + 1)) * NET_SIZE;
            t = (rand() / ((double)RAND_MAX + 1)) * NET_SIZE;
            
            if (Connections[t][s].are_connected) {
                i--;
                continue;
            }
            
            // use montecarlo to determine target cell selection probability
            // target with highest score on this round gets selected
            for (j=0; j<NET_SIZE; j++) {
                target_score_arr[j] = drand48() * (connect_prob * afferents_arr[j]/NET_SIZE + 1.0);//randdouble( CONNECTION_PROB + afferents_arr[j] * SMALL_WORLD_CONSTANT / (NET_SIZE - (CONNECTION_PROB * NET_SIZE)) );
                //else target_score_arr[j] = randdouble( CONNECTION_PROB );
            }
            
            // get max target_score_arr value
            max = 0.0;
            for (j=0; j<NET_SIZE; j++) {
                if (target_score_arr[j] > max) {
                    max = target_score_arr[j];
                    t = j;
                }
            }
            
            // set connection with probability
            if (t == s) { 
                Connections[t][s].are_connected = FALSE; // prevent connection-to-self
                i--; // don't use up the connection from the bin
            } else if (Connections[t][s].are_connected) {
                i--; // don't re-set and existing connection
            } else {
                Connections[t][s].are_connected = TRUE;
                
                afferents_arr[t]++;
                //printf("Connections made: %i\n", i);
            }        
        }
    } else if (CONNECT_TOPOGRAPHIC) {
        
        printf("Connecting randomly with open boundaries...\n");
        int mx = SHEET_WIDTH/2;
        int my = SHEET_HEIGHT/2;
        
        float px, py, p_exc, p_inh;
        
        for (s=0; s<NET_SIZE; s++) {
            
            // determine if the connection should be set or not based on topological rules 
            // corner cells connect 4x less than center cells            
            // if source is away from center then reduce the probability of making a connection proportionally
            px = 1 - (abs(mx - get_node_x(s))*0.5 / mx);
            py = 1 - (abs(my - get_node_y(s))*0.5 / my);
            
            p_exc = px*py*CONNECT_EXCIT_PROB;
            p_inh = px*py*CONNECT_INHIB_PROB;
            
            //printf("neuron %i, excit prob %f, inhib prob %f\n", s, p_exc, p_inh);
            
            for (t=0; t<NET_SIZE; t++) {
                
                if (s != t) {
                    // if source node is excitatory
                    if (Neurons[s].is_excit) {
                        if (drand48() < p_exc) {
                            Connections[s][t].are_connected = TRUE;
                            afferents_arr[t]++;
                        }
                    // if source is inhibitory
                    } else {
                        if (drand48() < p_inh) {
                            Connections[s][t].are_connected = TRUE;
                            afferents_arr[t]++;
                        }
                    }
                }
            }
        }
        
    } else {
        
        // connect with uniform random distribution
        printf("Connecting randomly...\n");
        for (t=0; t<NET_SIZE; t++) {
            for (s=0; s<NET_SIZE; s++) {
                if (drand48() <= connect_prob && s != t) {
                    Connections[s][t].are_connected = TRUE;
                    afferents_arr[t]++;
                }
            }
        }
        
    }
    
    // assign delays    
    int linear_delay;
    int max_exc_delay = 0;
    int max_inh_delay = 0;
    int decay_const;
    int del_missing = 0;
    double delay;
    //int random_delay;
    for (i=0; i<NET_SIZE; i++) {
        
        for (j=0; j<NET_SIZE; j++) {
            
            // compute topological delay
            if (Neurons[i].is_excit) decay_const = DELAY_EXCIT_DECAY_CONST;
            else decay_const = DELAY_INHIB_DECAY_CONST;
            
            delay = decay_const * 
            sqrt( pow( (i/SHEET_WIDTH)/(double)SHEET_WIDTH - (j/SHEET_WIDTH)/(double)SHEET_WIDTH, 2 ) + 
                 pow( (i%SHEET_HEIGHT)/(double)SHEET_HEIGHT - (j%SHEET_HEIGHT)/(double)SHEET_HEIGHT, 2 ) 
                 );
            
            linear_delay = (int) delay + MIN_TOPO_DELAY;
            
            // override topological
            //if (CONNECT_RANDOM) linear_delay = MIN_DELAY+(drand48()*((MAX_DELAY*STEPSIZE)-(MIN_DELAY*STEPSIZE)))/STEPSIZE; // convert delay in s to number of slots as per stepsize
            if (CONNECT_RANDOM) linear_delay = mindelay+(drand48()*((maxdelay*STEPSIZE)-(mindelay*STEPSIZE)))/STEPSIZE; // convert delay in s to number of slots as per stepsize
            //random_delay = MIN_DELAY+(drand48()*((MAX_DELAY*STEPSIZE)-(MIN_DELAY*STEPSIZE)))/STEPSIZE; // convert delay in s to number of slots as per stepsize
            
            if (Connections[i][j].are_connected) {
                
                if (CONNECT_TOPOGRAPHIC && linear_delay > MAX_TOPO_DELAY) {
                    Connections[i][j].are_connected = FALSE;
                    afferents_arr[j]--;
                    del_missing++;
                }

                Connections[i][j].delay = linear_delay;
                //printf("Connection delay from %i to %i: %i (rand %i)\n", i, j, linear_delay, random_delay);
            }            
            
            if (Neurons[i].is_excit && linear_delay > max_exc_delay) max_exc_delay = linear_delay;
            if (!Neurons[i].is_excit && linear_delay > max_inh_delay) max_inh_delay = linear_delay;
        }
    }

    // get network stats
    int total_connections = 0;
    for (t=0; t<NET_SIZE; t++) {
        for (s=0; s<NET_SIZE; s++) {
            if (Connections[t][s].are_connected) {
                total_connections++;
            }
        }        
    }
    
    double aff_mean;
    aff_mean = total_connections*1.0/NET_SIZE;
    
    // compute Delay SD
    int delay_list[total_connections];
    double delay_mean = 0.0;
    i=0;
    for (t=0; t<NET_SIZE; t++) {
        for (s=0; s<NET_SIZE; s++) {
            if (Connections[t][s].are_connected) {
                delay_list[i] = Connections[t][s].delay;
                i++;
                delay_mean += Connections[t][s].delay;
            }
        }
    }    
    
    delay_mean /= total_connections*1.0;
    
    //save network info for target cell t
    for (t=0; t<NET_SIZE; t++) {
        save_nodeconnections( t, Connections, Neurons, file_name );
    }
    printf("Total afferents: %i\nMissing due to Delay Cutoff: %i\nMean afferents per node: %g, Connections SD: %g\nDelay mean: %g (%gs), Delay SD: %g\n", total_connections, del_missing, aff_mean, std_int(afferents_arr, NET_SIZE), delay_mean, delay_mean*STEPSIZE/TIME_SCALE_FACTOR, std_int(delay_list, total_connections));
    printf("Default_Threshold %g, Max EXC delay in network: %i (%f),  Max INH delay in network: %i (%f)\n", DEFAULT_THRESHOLD*1.0, max_exc_delay, max_exc_delay*STEPSIZE/TIME_SCALE_FACTOR, max_inh_delay, max_inh_delay*STEPSIZE/TIME_SCALE_FACTOR);
    
}

int load_network_file(char *file_name, Connection Connections[][NET_SIZE], Neuron *Neurons) {
    
    FILE *fp;
    
    fp = fopen(file_name, "r");
    
    if (fp == NULL) {
        printf("Can't open network file...\n");
        return FALSE;
    }
    printf("Loading network file...\n");
    
    
    // RESET NETWORK
    int i,j;
    for (i=0; i<NET_SIZE; i++){
        for (j=0; j<NET_SIZE; j++) {
            Connections[i][j].are_connected = FALSE;
            Connections[i][j].weight = WEIGHT;
            Connections[i][j].delay = MAX_DELAY;
        }
    }
    
    int net_vals[5];
    float weight;
    int n;
    char line[100];
    char field[10];
    
    while (fgets(line, 100, fp) != NULL) {
//        printf("line: %s", line);
        const char *c = line;
        i = 0;
        while ( sscanf(c, "%9[^ ]%n", field, &n) == 1 ) {
            if (i<4) net_vals[i] = atoi(field);
            else weight = atof(field);
            c+=n;
            if (*c != ' ') break;
            ++c; // skip delimiter
            i++;
        }
        
        Connections[net_vals[1]][net_vals[0]].are_connected = TRUE;
        Connections[net_vals[1]][net_vals[0]].delay         = net_vals[3];
        Connections[net_vals[1]][net_vals[0]].weight        = weight;
        Neurons[net_vals[1]].is_excit = net_vals[2];
        
        //printf("net vals %i %i %i %i %g\n", net_vals[0], net_vals[1],net_vals[2],net_vals[3], weight);
    }
    fclose(fp);
    
    // PRINT NETWORK STATS
    int total_connections   = 0;
    float total_weight      = 0.0;
    int total_delay         = 0;
    int total_excitatory    = 0;
    int max_delay = 0;
    int min_delay = MAX_DELAY;
    
    for (i=0; i<NET_SIZE; i++){
        for (j=0; j<NET_SIZE; j++) {
            if (Connections[i][j].are_connected) {
                total_connections++;
                total_weight += Connections[i][j].weight;
                total_delay += Connections[i][j].delay;
                
                if (Connections[i][j].delay>max_delay) max_delay = Connections[i][j].delay;
                if (Connections[i][j].delay<min_delay) min_delay = Connections[i][j].delay;
            }
        }
        
        if (Neurons[i].is_excit) total_excitatory++;
    }    
    
    printf("Total connections: %i\n", total_connections);
    printf("Mean afferents per node: %g, Delay mean: %g (%g)\n", total_connections*1.0/NET_SIZE, total_delay*1.0/total_connections, (total_delay*1.0/total_connections)*STEPSIZE/TIME_SCALE_FACTOR);
    printf("Default_Threshold %g, Max delay: %i (%f),  Min delay in network: %i (%f)\n", DEFAULT_THRESHOLD*1.0, max_delay, max_delay*STEPSIZE/TIME_SCALE_FACTOR, min_delay, min_delay*STEPSIZE/TIME_SCALE_FACTOR);
    
    return TRUE;
}
