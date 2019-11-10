#!/usr/bin/bash
#sed -i -e 's/\r$//' test1.sh 
cd /project/hackathon/hackers11/shared/GettingStarted/TestSet/eggNOG_aligns/proteins
for file in `ls *.fasta`
do
echo ${file}
output_name=${file:0:4}
echo ${output_name}
./splitAlign.py ${file} -o ${output_name}
done

##this directory deposit the split align files
mkdir /project/hackathon/hackers11/shared/GettingStarted/TestSet/eggNOG_aligns/split_aligns

mv *half* /project/hackathon/hackers11/shared/GettingStarted/TestSet/eggNOG_aligns/split_aligns