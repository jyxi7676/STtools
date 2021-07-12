#!/bin/bash
miseq=
hiseq=
hdmilength=20
outdir=
miseq_pos="./spatialcoordinates.txt"
whitelists="./whitelist.txt"

progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
Firstly, we need to process our data with bash script extractCoord.sh, which can be found under script folder in this repository.
-m    <miseq>         : Path of read file from 1st-Seq.  Required.
-h    <hiseq>         : Path of read1 from 2nd-Seq.  Required.
-l    <hdmilength>    : An integer indicating the length of the HDMIs; For now, it can only take 20 or 30. In default, we assume if MiSeq is used for 1st-Seq, then hdmilength=20; if HiSeq is used for 1st-Seq, then hdmilength=30.
-o    <outdir>        : Path of output files. Required.
#-p    <miseq_pos>     : Five columns representing 1st-Seq HDMIs, lane, tile, X, Y.  Required.
#-w    <whitelists>    : This is the whitelists of HDMIs used for STARsolo alignment. If MiSeq is used for 1st-Seq, then whitelists are the reverse complementary of HDMIs in bottom tiles from 1st-Seq ; if HiSeq is used for 1st-Seq,  whitelists are the reverse complementary of HDMIs in all tiles in lane 2 from 1st-Seq.  Required.
EOF
}

function error_exit() {
    echo "ERROR: $1
    " >&2
    usage
    exit 1
}

function check_set() {
    value=$1
    name=$2
    flag=$3

    if [[ -z "$value" ]]
    then error_exit "$name has not been specified.  $flag flag is required"
    fi
}

set -e
# Fail if any of the commands in a pipeline fails
set -o pipefail

while getopts ":m:h:l:o:" options; do
  #echo $options
  case ${options} in
    m ) miseq=$OPTARG;;
    h ) hiseq=$OPTARG;;
    l ) hdmilength=$OPTARG;;
    o ) outdir=$OPTARG;;
    #p ) miseq_pos=$OPTARG;;
   # w ) whitelists=$OPTARG;;
    h ) usage
          exit 1;;
    \? ) usage
         exit 1;;
    * ) usage
          exit 1;;
  esac
done

shift $(($OPTIND - 1))

#echo "number: "$#
echo 'here'
check_set "$miseq" "Path of read file from 1st-Seq" "-m"
check_set "$hiseq" "Path of read1 from 2nd-Seq" "-h"
check_set "$outdir" "Path to output files" "-o"
#check_set "$whitelists" "Whitelists of barcodes from extractCoord.sh"  "-w"
#check_set "$miseq_pos" "Spatialcoordinates representing 1st-Seq HDMIs, lane, tile, X, Y." "-p"

#####Nueed to double check the checking commands
if (( $# != 0 ))
then error_exit "Incorrect number of arguments"
fi

#if hdmilength is not provided, provide a default value of 20
if [[ -z "$hdmilength" ]]
then 
  hdmilength=20
  echo "the length of the HDMIs is assigned to 20."
fi

#check if file1 and file2 exist
[ -f "$miseq" ] && echo "$miseq exists." || echo "$miseq does not exist!"
[ -f "$hiseq" ] && echo "$hiseq exists." || echo "$hiseq does not exist!"
[ -d "$outdir" ] && echo "$outdir exists." || echo "$outdir does not exist!"

#check if whitelists exsit
#[ -f "$whitelists" ] && echo "$whitelists exists." || echo "$whitelists does not exist!"

#check if whitelists exsit
#[ -f "$miseq_pos" ] && echo "$miseq_pos exists." || echo "$miseq_pos does not exist!"



#extract coordinates
echo "print here"
cd $outdir
echo " Extract HDMI sequences, position, and STARsolo Whitelists"
zcat $miseq | sed -n '1~4s/:/ /gp' | cut -d ' ' -f 4-7 >  ./pos-MiSeq-temp.txt

zcat $miseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-20  > ./HDMIs-MiSeq-temp.txt
cat ./HDMIs-MiSeq-temp.txt | rev | tr ACGTN TGCAN > ./HDMIs-MiSeq-temp-rev.txt
paste  ./HDMIs-MiSeq-temp-rev.txt ./pos-MiSeq-temp.txt | column -s $'\t' -t > ./MiSeq-temp-revHDMIs-pos.txt
#remove duplicated HDMIs
awk '!seen[$1]++' ./MiSeq-temp-revHDMIs-pos.txt  > $miseq_pos
##bottom tiles for MiSeq
if [  "$hdmilength" -eq 30 ]
then
  zcat $hiseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-30 > ./HDMI_SeqScope_2nd.txt
  zcat $miseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 3-32  > ./HDMIs-MiSeq-temp.txt
  cat ./HDMIs-MiSeq-temp.txt | rev | tr ACGTN TGCAN > ./HDMIs-MiSeq-temp-rev.txt
  paste  ./HDMIs-MiSeq-temp-rev.txt ./pos-MiSeq-temp.txt | column -s $'\t' -t > ./MiSeq-temp-revHDMIs-pos.txt
  #remove duplicated HDMIs
  awk '!seen[$1]++' ./MiSeq-temp-revHDMIs-pos.txt  > $miseq_pos
  cat $miseq_pos | awk '{ print $1 }' > $whitelists
else 
  zcat $hiseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-20 > ./HDMI_SeqScope_2nd.txt
  zcat $miseq | perl -lane 'print $_ if ( $. % 4 == 2 )'  | cut -c 1-20  > ./HDMIs-MiSeq-temp.txt
  cat ./HDMIs-MiSeq-temp.txt | rev | tr ACGTN TGCAN > ./HDMIs-MiSeq-temp-rev.txt
  paste  ./HDMIs-MiSeq-temp-rev.txt ./pos-MiSeq-temp.txt | column -s $'\t' -t > ./MiSeq-temp-revHDMIs-pos.txt
  #remove duplicated HDMIs
  awk '!seen[$1]++' ./MiSeq-temp-revHDMIs-pos.txt  > $miseq_pos
  cat $miseq_pos | awk '{ if ($3 > 2100) { print $1 } }' > $whitelists

fi
#rm ./pos-MiSeq-temp.txt
#rm ./HDMIs-MiSeq-temp.txt
#rm ./MiSeq-temp-revHDMIs-pos.txt
