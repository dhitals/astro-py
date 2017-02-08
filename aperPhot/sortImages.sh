#!/bin/bash

# make a list of all images and their RAs
fitskey RA n*/*fits > allfiles.tmp

declare -a obj ra

awk '{print $1}' objects.lis > obj.tmp
awk '{print $2}' objects.lis > ra.tmp

# read in the object names
let count=0
while read LINE; do
    obj[$count]=$LINE
    ((count++))
done < "obj.tmp"

# read in the object RAs
let count=0
while read LINE; do
    ra[$count]=$LINE
    ((count++))
done < "ra.tmp"

echo ${ra[*]}
# for each object, make a list of files and move them
size=${#obj[@]}
for i in $( seq $size ); do

    echo $i ${obj[$i]} ${ra[$i]}
    grep ${ra[$i]} allfiles.tmp > 1.tmp       # grep all files with the RA into a file
    sed 's/\(\.fits\).*/\1/' 1.tmp > 2.tmp    # keep only the filenames

    mkdir ${obj[$i]}
    while read LINE; do
	mv $LINE ${obj[$i]}/ # move all files to the new dir	
    done < "2.tmp"

done

rm *.tmp
