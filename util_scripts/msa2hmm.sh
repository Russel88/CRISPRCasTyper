#!/bin/bash

# Turn MSAs into HMMs

while IFS= read -r i
do
    k1=$(echo $i | awk '{print $1}')
    k2=$(echo $i | awk '{print $2}')
    echo $i
    hmmbuild --amino data/Profiles/${k2}.hmm data/trunc_Variants/${k1}.FASTA
done < data/Naming.txt
