#!/bin/bash

DIRS=$(echo star/batch_2/*_pass2)
for d in ${DIRS[@]}; do
  #echo $d >> metadata/batch_2_mappingrates.txt
  grep -H "Uniquely mapped reads % |" $d/Log.final.out >> metadata/batch_2_mappingrates.txt
done
