#!/bin/bash
#
# dokfwhm.sh -- given a directory full of FITS images, measure the FWHM
#
# usage: dokfwhm.sh path-name
# e.g. : dokfwhm.sh 'images/PG*.fits'


function join { local IFS="$1"; shift; badImagesDir="$*"; }

outfile='kfwhm.dat'
## backup $outfile if it already exists
if [ -a $outfile ]; then
    mv $outfile $outfile.bak
    echo "    Existing  $outfile moved to $outfile.bak"
fi
echo "frame       <FWHM>   skylevel  medianFWHM  NSTARS" > $outfile

## create a badImagesDir if it does not exist

IFS='/' read -ra array <<< "$1"
unset array[${#array[@]}-1]
join / "${array[@]}" badimages

if [ ! -d $badImagesDir ]; then
   mkdir -p $badImagesDir
fi

for eachImage in $1; do
    echo $eachImage
    
    # split the path+image name by '/' --imName will the last element
    IFS='/.' read -ra array <<< "$eachImage"
    obj=${array[${#array[@]}-2]}  # get obj

    # do the kfwhm
    kfwhm -p kfwhm.par -v $eachImage > $obj.fwhm
    tail -1 $obj.fwhm >> $outfile

    # read the FWHM and move to baddir if < 0.0
    out="$(tail -1 $obj.fwhm)"
    IFS=' ' read -ra array <<< "${out}"
    fwhm="$(bc -l <<< ${array[1]})"

    crit=3
    if [[ $(echo "$fwhm < 0" | bc) -eq 1 ]]; then
	mv $eachImage $badImagesDir
    fi
done
	    
