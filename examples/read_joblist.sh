#!/bin/bash

joblist=joblist #0

#declare -a myarray
declare -a myarray0
declare -a myarray1
declare -a myarray2

# Load file into array.
let i=0
let i0=0
let i1=0
let i1=0
while IFS=$'\n' read -r line_data; do
    # Parse “${line_data}” to produce content 
    # that will be stored in the array.
    # (Assume content is stored in a variable 
    # named 'array_element'.)
    # ...

#sbatch slurm_batch_pathnames10.sl
#Submitted batch job 113835
#    if [ $((i%2)) -eq 0 ] 
#    then 
##    	echo $i is "even" 
#    	echo line_data=$line_data
#    	array_element0=${i#*pathnames}
#        echo array_element0=$array_element0
#    	myarray0[i]="${array_element0}"
#    	((++i0))
#    fi
    echo line_data read=$line_data
    #array_element1=${line_data#*pathnames}
    array_element1=${line_data#*pathnames}
    echo array_element1=$array_element1    
    array_element1=${array_element1%*.sl}
    echo array_element1=$array_element1    

    myarray1[i]="${array_element1}" # Populate array.
#    myarray[i]="${line_data}"
    ((++i))
done < ${joblist}1$joblistgroup
#done < ${joblist}1


let i=0
while IFS=$'\n' read -r line_data; do
	echo line_data read=$line_data
	array_element2bis=${line_data##job}
	array_element2=${line_data:20:26}
        
	echo array_element2bis=$array_element2
	echo array_element2=$array_element2
    	myarray2[i]="${array_element2}" # Populate array.
    	((++i))
done < ${joblist}2$joblistgroup
#done < ${joblist}2





let i=0
while (( ${#myarray1[@]} > i )); do
    printf "${myarray1[i++]}\n"
done

let i=0
while (( ${#myarray2[@]} > i )); do
    printf "${myarray2[i++]}\n"
done

let i=0
while (( ${#myarray2[@]} > i )); do
    echo 'test' "${myarray1[$i]}" 'is job' "${myarray2[$i]}"
    ((++i))
done

let i=0
while (( ${#myarray2[@]} > i )); do
    
    echo 'test' "${myarray1[$i]}" 'is job' "${myarray2[$i]}"
    ((++i))
done
#exit

tot_successes="$(cat output-*.out  | grep CONGRATULATIONS  | wc -l)"
tot_runs="$(ls -1 output-*.out | wc -l)"
echo $tot_successes successful out of $tot_runs

let isuc=0
let ifail=0
let i=0
while (( ${#myarray2[@]} > i )); do

#    stringin < output-*myarray1[$i].out   
#    issuccess="$(cat output-*myarray1[$i].out  | grep CONGRATULATIONS  | wc -l)"

#echo search output-*${myarray1[$i]}.out
#outputlist=$(ls -1 output-*.out | grep ${myarray1[$i]})
#echo outputlist = $outputlist

nummer=${myarray2[$i]}
echo nummer=$nummer
capture=$(find *"$nummer"*out)
capture2=`find *"$nummer"*out`

echo capture=$capture
echo capture2=$capture2

issuccess="$(cat "$capture"  | grep CONGRATULATIONS  | wc -l)"
issuccess2=`cat "$capture"  | grep CONGRATULATIONS  | wc -l`
echo $issuccess
echo $issuccess2



#[ $1 -gt 100 ]
#if [[grep -Fxq CONGRATULATIONS "$capture"]]

echo issuccess=$issuccess
echo [ $issuccess -eq 1 ]

  if [ $issuccess == 1 ]
   then
    # code if found
    ((++isuc))
    echo 'test' "${myarray1[$i]}" 'is job' "${myarray2[$i]}" success=$issuccess OK
  else
     ((++ifail))
    echo 'test' "${myarray1[$i]}" 'is job' "${myarray2[$i]}" success=$issuccess not found
    # code if not found
  fi


    ((++i))
done

echo i=$i
echo isuc=$isuc
echo ifail=$ifail
#output-prod-cn08-113896.out



