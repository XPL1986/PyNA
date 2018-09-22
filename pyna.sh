#!/usr/bin/env bash

for bam in "$@"
    do
        if [[ ${bam} = *"DNA"* ]]; then
            ./dna.sh ${bam}
            ./filter.sh ${bam}
        elif [[ ${bam} = *"RNA"* ]]; then
            ./rna.sh ${bam}
            ./filter.sh ${bam}
        else
            echo "${bam} is an invalid input."
        fi
    done