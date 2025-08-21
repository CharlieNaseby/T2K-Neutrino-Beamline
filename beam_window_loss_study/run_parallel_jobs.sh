#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <argument>"
    exit 1
fi


argument="$1"
counter=1989
for (( i=0; i<4; i++ )); do
# Loop to run the program with the argument and an incrementing counter
    echo "Running: $argument $counter"
    bdsim --file=target_beamline.gmad --outfile="${argument}_${counter}" --batch --ngenerate=10000 --seed="$counter" &
    ((counter++))
done

wait

bash create_blm_histos.sh ${argument}


