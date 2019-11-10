#!/bin/bash
#compare.sh

##Use Case
# bash compare.sh -t 511145 -f cydA_trimmed_ENOG4105C4M.fasta -d 0.9

usage() {
  echo "-t --taxid"
  echo "-f --fasta file"
  echo "-d --depth"
  exit 1
}

OPTIND=1
while getopts :t:f:d: opt
do
  case ${opt} in
    t) tid=${OPTARG};;
    f) fasta=${OPTARG};;
    d) depth=${OPTARG};;
  esac
done

shift $((${OPTIND} -1))

pname=$(echo $fasta | awk -F'[_.]' '{print $1}')
filename=$(echo ${pname}_${depth})

depth=$(bc <<< "${depth}")
Ecoli=$(grep ${tid} ${fasta} -A1);

while read coline
do
  coline_name=$(echo ${coline} | tr -d ">")
  read coline_seq
  while read line
  do
    read line_seq
    i=0
    match=0
    total=0
    while (( i++ < ${#line_seq} ))
    do
      EP=$(expr substr "${coline_seq}" $i 1)
      XP=$(expr substr "${line_seq}" $i 1)
      if [[ ${EP} != ${XP} ]]
      then
        total=$((total + 1))
      elif [[ ${EP} != "-" && ${XP} != "-" ]]
      then
        total=$((total + 1))
        match=$((match + 1))
      fi
    done
    Adepth=$(bc <<< "scale = 5; ${match}/${total}")
    if (( $(echo "${Adepth} > ${depth}" | bc) ))
    then
      echo ${line} >> ${filename}_${coline_name}.txt
      echo ${line_seq} >> ${filename}_${coline_name}.txt
    fi
  done < "${fasta}"
done <<< "${Ecoli}"
