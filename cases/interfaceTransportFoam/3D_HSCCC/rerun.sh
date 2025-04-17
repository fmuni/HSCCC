#!/bin/bash

rm -r slurm* process*
decomposePar
sbatch runjob
