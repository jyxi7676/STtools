#!/bin/bash
file1=
file2=
hdmilength=20
whitelists=
outprefix=
starpath=
seqtkpath=
geneIndex=

progname=`basename $0`

function usage () {
    cat >&2 <<EOF
USAGE: $progname [options] <unmapped-queryname-sorted.bam>
This step is to preprocess the fastq files and to align the data to reference genome.The bash script takes in several user defined parameters and outputs STARsolo summary statistics, and DGE in the current directory. Note: Here we assume you already have the reference genome that is needed for STARsolo alignment. If not please refer to https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/lessons/03_alignment.html.
-a    <file1>         : Read1 from 2nd-Seq.  Required.
-b    <file2>         : Read2 from 2nd-Seq.  Required.
-l    <hdmilength>    : An integer indicating the length of the HDMIs; For now, it can only take 20 or 30. In default, we assume if MiSeq is used for 1st-Seq, then hdmilength=20; if HiSeq is used for 1st-Seq, then hdmilength=30.
-w    <whitelists>    : Whitelists of barcodes from extractCoord.sh.  Required.
-o    <outprefix>     : Prefix for STARsolo output.  Required.
-t    <starpath>      : Full path for STAR software.  Required.
-q    <seqtkpath>     : Full path for seqtk tool.  Required.
-g    <geneIndex>     : Reference genome directory.  Required.
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

while getopts ":a:b:l:w:o:t:q:g:" options; do
  #echo $options
  case ${options} in
    a ) file1=$OPTARG;;
    b ) file2=$OPTARG;;
    l ) hdmilength=$OPTARG;;
    w ) whitelists=$OPTARG;;
    o ) outprefix=$OPTARG;;
    t ) starpath=$OPTARG;;
    q ) seqtkpath=$OPTARG;;
    g ) geneIndex=$OPTARG;;
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

check_set "$file1" "Read1 from 2nd-Seq" "-a"
check_set "$file2" "Read2 from 2nd-Seq" "-b"
check_set "$whitelists" "Whitelists of barcodes from extractCoord.sh"  "-w"
check_set "$starpath" "Path for STAR software" "-t"
check_set "$seqtkpath" "Path for seqtk tool" "-q"
check_set "$geneIndex" "Reference genome directory" "-g"

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

#if outprefix is not provided, then make a temp directory
if [[ -z "$outprefix" ]]
then outprefix=`mktemp -d`
     echo "Using temporary directory $outprefix"
fi

#check if file1 and file2 exist
[ -f "$file1" ] && echo "$file1 exists." || echo "$file1 does not exist!"
[ -f "$file2" ] && echo "$file2 exists." || echo "$file2 does not exist!"

#check if whitelists exsit
[ -f "$whitelists" ] && echo "$whitelists exists." || echo "$whitelists does not exist!"

#check if reference genome directory exist. 
[ -d "$geneIndex" ] && echo "$geneIndex exists." || echo "$geneIndex does not exist!"

#check if STAR exist or not
star_executable=$starpath$"/STAR"
if [[ ! ( -x $star_executable && -f $star_executable ) ]]
then error_exit "STAR executable $star_executable passed via -s does not exist or is not executable"
#elif which STAR > /dev/null
#then echo > /dev/null
else echo "$star_executable exists"
fi

#check if seqtk tool exist or not
seqtk_executable=$seqtkpath$"/seqtk"
if [[ ! ( -x $seqtk_executable && -f $seqtk_executable ) ]]
then error_exit "seqtk executable $seqtk_executable passed via -s does not exist or is not executable"
#elif which STAR > /dev/null
#then echo > /dev/null
else echo "$seqtk_executable exists"
fi

file1_basename="$(basename -s .fastq.gz $file1)"
#echo $f1
file1_dir="$(dirname "${file1}")"
file1_final_postfix="_modified.fastq"
file1_final="$file1_dir/$file1_basename$file1_final_postfix"
file1_final_gz="$file1_final.gz"


#cd $seqtkpath
echo 'trimming R1'
$seqtk_executable trimfq -q 0 -l $hdmilength $file1 > ./file1_trim.fastq
pigz -p 8 -f ./file1_trim.fastq
#gzip ./file1_trim.fastq
echo $file1_final

if [  "$hdmilength" -eq 30 ]
then 
  paste <(zcat ./file1_trim.fastq.gz) <(zcat $file2) | perl -lane 'if ( $. % 4 == 1 ) { print "$F[0] $F[1]"; } elsif ( $. % 4 == 3 ) { print "+"; } else { print substr($F[0],0,30).substr($F[1],0,9).substr($F[0],50); }' > $file1_final
else 
  paste <(zcat ./file1_trim.fastq.gz) <(zcat $file2) | perl -lane 'if ( $. % 4 == 1 ) { print "$F[0] $F[1]"; } elsif ( $. % 4 == 3 ) { print "+"; } else { print substr($F[0],0,20).substr($F[1],0,9).substr($F[0],50); }' > $file1_final
fi
#gzip $file1_final
pigz -p 8 -f $file1_final
rm ./file1_trim.fastq.gz
echo 'finish trimming'
#echo $file1_final
rndstart=`expr 1 + $hdmilength`

echo 'Start alignment'
$star_executable    --genomeDir  $geneIndex \
                    --readFilesIn  $file2 $file1_final_gz  \
                    --outSAMtype BAM SortedByCoordinate  \
                    --readFilesCommand zcat \
                    --runDirPerm All_RWX \
                    --outFileNamePrefix $outprefix  \
                    --soloType CB_UMI_Simple \
                    --soloCBstart 1 --soloCBlen $hdmilength \
                    --soloUMIstart $rndstart --soloUMIlen 9 \
                    --soloCBwhitelist $whitelists \
                    --runThreadN 6 \
                    --soloBarcodeReadLength 0 \
                    --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
                    --outFilterScoreMinOverLread 0 \
                    --outFilterMatchNminOverLread 0 \
                    --clip3pAdapterSeq AAAAAAAAAA \
                    --clip3pAdapterMMp 0.1 \
                    --soloFeatures Gene GeneFull SJ Velocyto \
                    --limitOutSJcollapsed 1000000 \
                    --soloCellFilter None
