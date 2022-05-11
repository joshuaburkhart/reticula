#!/bin/bash
# go interactive before doing anything on the hpc
srun -I30 -p sandbox -N 1 -c 1 --mem=6G -t 0-01:00:00 --pty /bin/bash


