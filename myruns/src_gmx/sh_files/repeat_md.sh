#!/bin/bash 

# submit the run_md job (run_md.sh)
ID=$(sbatch run_md.sh | awk '{print $4}')
echo "$ID"
# submit the MD run job in a loop (run_md.sh)
#for ((i=1;i<=3;i+=1))
#do 
#ID=$(sbatch -w "ended("$ID")" run_md.sh | awk '{print $2}' | tr "<" " " | tr ">" " " | awk '{print $1}')
#done
