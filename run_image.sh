#!/bin/bash

SCHEDULE=$1
PATTERN=$2
BUCKETS=$3
CORES=$4

docker run -v $(pwd)/nonbdna_pipeline:/nonbdna_pipeline \
    nonbdna-workflow:latest \
    $SCHEDULE $PATTERN $BUCKETS $CORES
