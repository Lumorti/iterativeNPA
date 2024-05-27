
# Repeat many times
> data.dat
for i in {1..1000}
do

    # Run the program
    ./run --i3322 -l 2 -v 3 -t 3 -S $i > output.log

    # Extract the values from "Result: 0.34" using awk
    result=$(awk '/Result: / {print $2}' output.log)

    # Extract the values from "Num zero eig: 11"
    zeroEig=$(awk '/Num zero eig: / {print $4}' output.log)

    # Extract the values from "Num pos eig: 11"
    posEig=$(awk '/Num pos eig: / {print $4}' output.log)

    # Write to file
    echo "$i $result $zeroEig $posEig"
    echo "$result $zeroEig $posEig" >> temp.dat

done


