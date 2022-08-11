#!/bin/bash
#$ -cwd -V
#$ -P jravel-lab
#$ -q all.q  
#$ -m ea
#$ -l mem_free=4G

#cp -r /autofs/scratch/jolim/persist/mg_mt/mg/test_data/test_data2/* ./
/local/projects-t3/MSL/pipelines/bin/nextflow run /local/projects-t3/MSL/pipelines/packages/YAMP/YAMP.nf -profile base,sge,conda -with-dag dag.svg --indir /local/projects-t3/<project>
