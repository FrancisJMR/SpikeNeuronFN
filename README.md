Spikeneuron Fitzhug-Nagumo Model
version: 0.1
by Francis Jeanson 2012

Core Files: network.c/.h spikeneuron.c/.h

Features
- FitzHugh-Nagumo Neurone Model
- Discrete Transmission Delays
- Learning via STDP
- Option for Random, 2D Topological, or Small-World Connectivity

Most parameters are adjustable in network.h, neurone specific parameters adjustable in spikeneuron.h

Additional Files: run_trial.c

This file contains the main function to run trial experiments for multi-pattern learning and recall. 

Compile with gcc using:
gcc spikeneuron.c network.c run_trials.c -o run_trials

Run example:
$ ./run_trials 10 2 3 8 4 10000
This:
- runs 10 trials
- Uses a silence period factor of 2 (2*300+300 time steps)
- Generate 3 consecutive patterns
- Each pattern has 8 cells
- 4 neurones will be used to test recall of each pattern
- Learning period will have a total of 10000 time steps

Example output:
Pattern sample recall rate [0.75, 0.75, 0.875, 0.875, 0.875, 1, 0.625, 0.625, 0.875, 1, ]
Pattern 0 recall: 0.825 
Pattern sample recall rate [0.875, 0.875, 0.75, 0.75, 0.875, 1, 1, 0.75, 1, 1, ]
Pattern 1 recall: 0.8875 
Pattern sample recall rate [0.875, 1, 0.875, 0.625, 0.875, 0.875, 0.5, 0.75, 0.5, 0.875, ]
Pattern 2 recall: 0.775 
Over Excits 0
Self sustained trials: 0 0 0 0 0 0 0 0 0 0 
Self Sustained total 0/10

Pattern sample recall rate lists the average (among the 3 patterns) recall for each trial.
The total average is displayed below.
Over Excits displays the number of cells that do not belong to the pattern that fired during recall.
Self sustained trials lists a 1 if least one of the three patterns self sustained for each trial.


