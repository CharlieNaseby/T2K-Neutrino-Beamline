#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <argument>"
    exit 1
fi
argument="$1"
counter=1989
for (( i=0; i<10; i++ )); do
# Loop to run the program with the argument and an incrementing counter
    echo "Running: $argument $counter"
    rebdsimOptics "${argument}/${argument}_${counter}.root" "${argument}/${argument}_${counter}_optics.root" &
    ((counter++))
done

wait
pattern="${argument}/${argument}*_optics.root"
files=($pattern)
rebdsimCombine "${argument}_optics_combined.root" "${files[@]}"

