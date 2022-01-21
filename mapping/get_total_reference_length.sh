#!/bin/bash -l

## script can be run with an srun call
## srun -n 4 -p interactive ./get_total_reference_length.sh
## this script outputs a file with all the reference lengths in a directory with subdirectories holding reference genomes
## note: these subdirectories should already have *.dict files (since they have been used for mapping) 

for d in ./*/ ; do
    cd $d 
    name=$(echo $d | sed 's/\.//1' | sed 's/\///g')
    path=$(awk 'NR == 2 {print $5}' *.dict | sed 's/UR:file://g')
    num_contigs=$(wc -l <*.dict)
    if (( ${num_contigs} > 1 )) ; then
        len_total=0
        for i in $(seq 2 ${num_contigs})
            do
            len=$(awk -v i="$i" 'NR == i {print $3}' *.dict | sed 's/LN://g')
            len_total=$((len_total + len))
        done
        awk -v name="$name" -v path="$path" -v len_total="$len_total" 'NR == 1 {print name, len_total, path}' *.dict >> ../reference_lengths
    else
        awk -v name="$name" -v path="$path" 'NR == 1 {print name, $2, path}' *.dict >> ../reference_lengths
    fi
    cd ..
done
