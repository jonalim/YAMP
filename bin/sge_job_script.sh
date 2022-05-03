#!/bin/bash
#$ -cwd -V
#$ -P jravel-lab
#$ -q all.q  
#$ -m ea
#$ -l mem_free=4G

/local/projects-t3/MSL/pipelines/bin/nextflow run /local/projects-t3/MSL/pipelines/packages/YAMP/YAMP.nf -profile base,sge,conda -with-dag dag.svg --reads /local/projects-t3/XDENT/covid_samples
