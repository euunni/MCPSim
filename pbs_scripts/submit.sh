#!/bin/bash

NUM_JOBS=40
INTERVAL=7200 # 2h

while true; do
    echo "Submitting $NUM_JOBS jobs..."
    for ((i=0; i<$NUM_JOBS; i++)); do
        qsub script_multiRun.pbs
        sleep 1
    done
    echo "Waiting for 2 hours before submitting the next batch..."
    sleep $INTERVAL
done
