#!/bin/bash

cat * | grep passed | grep IR | awk -F ' ' '{ print $7 }' > processed_succesfully_shuffled_IR.txt
