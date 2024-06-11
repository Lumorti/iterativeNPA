#!/bin/bash

# generic loop
> out.log
for i in {1..10000}
do
    echo "Running iteration $i"
    echo "Running iteration $i" >> out.log
    ./run -S $i --r3322 -t 1 -i 10000 >> out.log
done
