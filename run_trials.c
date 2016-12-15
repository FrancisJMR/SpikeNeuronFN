/** Spiking Neuron Neural Simulator
* Written by Francis Jeanson 2012
*/

#include "spikeneuron.h"

// STATISTIC MAIN: Run N times: generate network, apply STDP, test recall
int main(int argc, char *argv[]) {
    
    int i, j, t, c, p;
    float ts;
    int over_excit = 0;
    
    // init random seeds
    srand(time(NULL));
    srand48(time(NULL));
    
    char spiketimes_output_filename[100];
    sprintf(spiketimes_output_filename, "out_spiketimes_.txt");
    
    char network_filename[100];
    sprintf(network_filename, "out_nodeconnCs_gen.txt");
    
    // delete old files
    FILE *fp;    
    FILE *fp_null   = NULL;
    
	int sample_size     = atoi(argv[1]); // (INT) the number of trials to run
	int hold_factor     = atoi(argv[2]); // (INT) the sleep hold factor
    int patterns_number = atoi(argv[3]); // (INT) the number of patterns to test for
    int patterns_size   = atoi(argv[4]); // (INT) the number of cells per pattern to assign
    int recall_size     = atoi(argv[5]); // (INT) the number of cells to stimulate for recall testing
    int tot_learn_time  = atoi(argv[6]); // (INT) the total time of simulation for training
    
    float patterns_recall[patterns_number][sample_size]; // The final score for each pattern
    
    for (p=0; p<patterns_number; p++) {
        for (i=0; i<sample_size; i++){
            patterns_recall[p][i] = 0.0;
        }
    }
    
    int self_sustained[sample_size];
    for (i=0; i<sample_size; i++){
        self_sustained[i] = 0;
    }
    
    int hold_steps     = hold_factor*STDP_WINDOW+MAX_DELAY;
    int pat_learn_time = (tot_learn_time/patterns_number) - hold_steps;
    if (pat_learn_time <= 0) {
        printf("Not enough time to train each pattern!\n");
        exit(1);
    }

    if (patterns_number*patterns_size > NET_SIZE) {
        printf("Too many patterns, or size of patterns too big!\n");
        exit(1);
    }

    // Initialize all neurons of the network.
    Neuron Ns[NET_SIZE]; // input neurons
    Connection Cs[NET_SIZE][NET_SIZE]; // input network
        
    // Main statistic loop
    for (i=0; i<sample_size; i++) {
        printf("\nTest %i\n", i);
		
        // reset neurones and axons
        init_neurons(Ns, INHIB_PERCENT);
        
        // GENERATE THE NETWORK
        generate_network(Cs,  Ns,  MIN_DELAY, MAX_DELAY, CONNECTION_PROB, WEIGHT, network_filename);
        printf("Base Weight %g\n", WEIGHT);
        //init learn
        int pat_learn_counter = 0;
		int stimulation_on = 1;
        
        int cell_start = 0;
        int cell_stop  = cell_start+patterns_size;
        //printf("cell start stop %i %i\n", cell_start, cell_stop);
    
        // reset the output file
        fp = fopen(spiketimes_output_filename, "w");
        fclose(fp);
        
        // train patterns
        for (t=0; t<tot_learn_time; t++){

            ts = t*STEPSIZE;
            
            if (t>0 && t%(pat_learn_time+hold_steps) == 0) {
                cell_start += patterns_size;
                cell_stop  += patterns_size;
                pat_learn_counter = 0;
                stimulation_on = 1;
                //printf("Learn time %i  cell start stop %i %i\n", t, cell_start, cell_stop);
            }
            
            // hold pattern
            if (pat_learn_counter == pat_learn_time) {
                stimulation_on = 0;
            }
            
            propagate_signals(Cs, Ns, Ns);
            integrate_neurons(t, Ns, fp_null, cell_start, cell_stop, stimulation_on);
            perform_stdp(t, Cs, Ns, Ns);
            //store_neurone_activity(t, Ns, spiketimes_output_filename);
            
            pat_learn_counter++;
        }
        
        printf("Pattern learn time %i, hold steps %i\n", pat_learn_time, hold_steps);
        
        
        //determine if sustained
        for (j=0; j<NET_SIZE; j++) {
            for (c=0; c<STDP_WINDOW; c++){
                if (Ns[j].v_hist[c] > SPIKE_DETECT_THRESHOLD) self_sustained[i] = 1;
            }
        }
        
        // Save learned node connections NB: This saves the network from the last trial only!
        fp = fopen("out_nodeconnCs_stdp_test.txt", "w");
        fclose(fp);
        for (c=0;c<NET_SIZE;c++){
            save_nodeconnections( c, Cs, Ns, "out_nodeconnCs_stdp_test.txt");
        }
        //*/
        
        printf("\nLearning Complete...\n\n");

        // BEGIN RECALL
        cell_start      = 0;
        cell_stop       = cell_start+patterns_size;
        int recall_start = cell_stop-recall_size;
        int recall_stop = cell_start+recall_size;
        if (recall_stop > cell_stop) printf("Warning: Recall Size Bigger than Pattern Size.\n");
		int pattern_counts[NET_SIZE]; // spike counter
		
        int discrete_ts;
        // test each pattern
        for (p=0; p<patterns_number; p++){
            
            //printf("Recall Start , Stop %i %i\n", recall_start, cell_stop);
            
            for (j=0; j<NET_SIZE; j++) {
                pattern_counts[j] = 0;
            }
            
            init_neurons(Ns, INHIB_PERCENT);
            
            for (t=0; t<10000; t++) {
         
                propagate_signals(Cs, Ns, Ns);
                integrate_neurons(t, Ns,  fp_null, recall_start, cell_stop, 1);
                detect_pattern(t, Ns, 0.0, 0, NET_SIZE, pattern_counts);
                //store_neurone_activity(t, Ns, spiketimes_output_filename);
            }
            
            for (j=0; j<NET_SIZE; j++) {
				// pattern good if cell in the subset spiked more than 20 times
                if (j >= cell_start && j < cell_stop && pattern_counts[j] > 20) patterns_recall[p][i] += 1.0; 
                if ( (j < cell_start || j >= cell_stop) && pattern_counts[j] > 20) over_excit += 1;
            }
            
            cell_start  += patterns_size; // always start from the first cell in the pattern
            cell_stop   += patterns_size;
            recall_start+= patterns_size;
            recall_stop += patterns_size;
        }
        printf("\n");
                
    }
    
    for (p=0; p<patterns_number; p++) {
        float avg = 0.0;
        printf("Pattern sample recall rate [");
        for (i=0; i<sample_size; i++){
            avg += patterns_recall[p][i];
            printf("%g, ", patterns_recall[p][i]/(patterns_size*1.0));
        }
        printf("]\n");
        printf("Pattern %i recall: %g \n", p, avg/(sample_size*1.0*patterns_size));
    }
    printf("Over Excits %i\n", over_excit);
    
    printf("Self sustained trials: ");
    int tot=0;
    for (i=0; i<sample_size; i++) {
        printf("%i ", self_sustained[i]);
        tot += self_sustained[i];
    }
    printf("\n");
    printf("Self Sustained total %i/%i\n", tot, sample_size);
    
    return 0;
}