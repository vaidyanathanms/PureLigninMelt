#!/bin/bash 

# submit the preprocessing job (run_preprocess.sh)
ID=$(bsub run_preprocess.sh | awk '{print $2}' | tr "<" " " | tr ">" " " | awk '{print $1}')

# submit the MD run job in a loop (run_md.sh)
<<<<<<< HEAD
for ((i=1;i<=3;i+=1))
=======
for ((i=1;i<=2;i+=1))
>>>>>>> origin/master
do 
ID=$(bsub -w "ended("$ID")" run_md.sh | awk '{print $2}' | tr "<" " " | tr ">" " " | awk '{print $1}')
done
