#!/bin/bash

# For standard I3322
for i in {1..20}
do
    ./run --i3322 -l 3 -d -m $i >> data/i3322.dat
done

# For randomized I3322
#for j in {1..30}
#do
    #for i in {5..80}
    #do
        #./run -S $j --i3322 -l 3 -d -m $i -i 1000 | tee -a data/r3322.dat
    #done
#done

