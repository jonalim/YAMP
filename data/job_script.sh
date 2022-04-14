#!/bin/bash
#$ -cwd -V
#$ -P jravel-lab
#$ -q threaded.q 
#$ -pe thread 7 
#$ -m ea
#$ -l mem_free=52G

/local/projects-t3/MSL/pipelines/bin/nextflow run /local/projects-t3/MSL/pipelines/packages/YAMP/YAMP.nf -profile base,conda -with-dag dag.svg --reads /local/scratch/jolim/temp/mg/input_1_median  -resume
