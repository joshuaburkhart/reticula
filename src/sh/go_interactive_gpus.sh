#!/bin/bash
# go interactive before doing anything on the hpc
# srun -I30 -p kill-shared --gpus-per-node=8 -N 1 -c 8 --mem=32G -t 0-01:00:00 --pty /bin/bash
srun -I30 -p gpu --gres=gpu:8 -N 1 -c 8 --mem=32G -t 0-06:00:00 --pty /bin/bash


