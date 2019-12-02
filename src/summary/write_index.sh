#!/bin/bash

### Variables
index_file="$1"
index_new="$2"
n="$#"
files="${@:3:n}"
cp $index_file $index_new

echo ########### ------- Files to include in index.md ------------#########
echo 
echo $files
echo $index_file


### Add html files to script
for sample in $files
do
    link="+ [${sample}](batch_effect_${sample}.html)"
    link_sim="+ [${sample}](simulation_${sample}.html)"
    link_qc="+ [${sample}](countSimQC_${sample}.html)"
    echo $link
    echo $link_sim
    echo $link_qc
    sed -i "/^## Batch effects/a $link" $index_new
    sed -i "/^## Simulation/a $link_sim" $index_new
    sed -i "/^## CountSimQC/a $link_qc" $index_new
done


echo ###########-------- finished index file ----------------###########
echo "added $n files to index.md"
