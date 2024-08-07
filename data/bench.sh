#!/bin/bash

#for size in 10 20 30 40 50 60 70 80 90 100 110 120 130 140
#do
    #for iter in {1..100}
    #do
        #outputOptim=$(/usr/bin/time -v ../run -S $iter --rxx22 $size -c 4 -D -s o &> /dev/stdout)
        #outputMosek=$(/usr/bin/time -v ../run -S $iter --rxx22 $size -c 4 -D -s m &> /dev/stdout)
        #outputSCS=$(/usr/bin/time -v ../run -S $iter --rxx22 $size -c 4 -D -s s &> /dev/stdout)
        #memoryOptim=$(echo "$outputOptim" | grep "Maximum resident set size" | awk '{print $6}')
        #memoryMosek=$(echo "$outputMosek" | grep "Maximum resident set size" | awk '{print $6}')
        #memorySCS=$(echo "$outputSCS" | grep "Maximum resident set size" | awk '{print $6}')
        #timeOptim=$(echo "$outputOptim" | grep "Time to solve" | awk '{print $4}')
        #timeMosek=$(echo "$outputMosek" | grep "Time to solve" | awk '{print $4}')
        #timeSCS=$(echo "$outputSCS" | grep "Time to solve" | awk '{print $4}')
        #timeOptim=$(echo "$timeOptim" | tr -dc '0-9')
        #timeMosek=$(echo "$timeMosek" | tr -dc '0-9')
        #timeSCS=$(echo "$timeSCS" | tr -dc '0-9')
        #resultOptim=$(echo "$outputOptim" | grep "Result" | awk '{print $2}')
        #resultMosek=$(echo "$outputMosek" | grep "Result" | awk '{print $2}')
        #resultSCS=$(echo "$outputSCS" | grep "Result" | awk '{print $2}')
        #linOptim=$(echo "$outputOptim" | grep "Final linear error" | awk '{print $4}')
        #linMosek=$(echo "$outputMosek" | grep "Final linear error" | awk '{print $4}')
        #linSCS=$(echo "$outputSCS" | grep "Final linear error" | awk '{print $4}')
        #eigOptim=$(echo "$outputOptim" | grep "Final min eig" | awk '{print $4}')
        #eigMosek=$(echo "$outputMosek" | grep "Final min eig" | awk '{print $4}')
        #eigSCS=$(echo "$outputSCS" | grep "Final min eig" | awk '{print $4}')
        #timeGen=$(echo "$outputOptim" | grep "Time to generate" | awk '{print $4}')
        #timeGen=$(echo "$timeGen" | tr -dc '0-9')
        #echo "$iter $size $timeGen | $memoryOptim $timeOptim $resultOptim $linOptim $eigOptim | $memoryMosek $timeMosek $resultMosek $linMosek $eigMosek | $memorySCS $timeSCS $resultSCS $linSCS $eigSCS"
        #echo "$iter $size $timeGen | $memoryOptim $timeOptim $resultOptim $linOptim $eigOptim | $memoryMosek $timeMosek $resultMosek $linMosek $eigMosek | $memorySCS $timeSCS $resultSCS $linSCS $eigSCS" >> rxx.dat
    #done
#done

for size in 150 160 170 180 190 200 210
do
    for iter in {1..100}
    do
        outputOptim=$(/usr/bin/time -v ../run -S $iter --rxx22 $size -c 4 -D -s o &> /dev/stdout)
        memoryOptim=$(echo "$outputOptim" | grep "Maximum resident set size" | awk '{print $6}')
        timeOptim=$(echo "$outputOptim" | grep "Time to solve" | awk '{print $4}')
        timeOptim=$(echo "$timeOptim" | tr -dc '0-9')
        resultOptim=$(echo "$outputOptim" | grep "Result" | awk '{print $2}')
        linOptim=$(echo "$outputOptim" | grep "Final linear error" | awk '{print $4}')
        eigOptim=$(echo "$outputOptim" | grep "Final min eig" | awk '{print $4}')
        timeGen=$(echo "$outputOptim" | grep "Time to generate" | awk '{print $4}')
        timeGen=$(echo "$timeGen" | tr -dc '0-9')
        echo "$iter $size $timeGen | $memoryOptim $timeOptim $resultOptim $linOptim $eigOptim | 0 0 0 0 0 | 0 0 0 0 0"
        echo "$iter $size $timeGen | $memoryOptim $timeOptim $resultOptim $linOptim $eigOptim | 0 0 0 0 0 | 0 0 0 0 0" >> rxx.dat
    done
done

git add rxx*.dat -f
git commit -m "automatic commit of data"
git push
