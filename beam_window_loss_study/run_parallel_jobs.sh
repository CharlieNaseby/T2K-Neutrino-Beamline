#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <argument>"
    exit 1
fi


argument="$1"

mkdir -p "${argument}"

counter=1989
for (( i=0; i<10; i++ )); do
# Loop to run the program with the argument and an incrementing counter
    echo "Running: $argument $counter"
    bdsim --file=target_beamline.gmad --outfile="${argument}/${argument}_${counter}" --batch --ngenerate=50000 --seed="$counter" &
    ((counter++))
done

wait

bash create_blm_histos.sh ${argument}


