READS_FILE='chr22.fastq'

if [ "$#" -eq 7 ]; then
    FREQ_MINIMIZERS=$1
    K=$2
    THRESHOLD=$3
    SAMPLES=$4
    LIMIT=$5
    THREADS=$6
    PREPROCESS=$7

elif [ "$#" -eq 8 ]; then
    READS_FILE=$1
    FREQ_MINIMIZERS=$2
    K=$3
    THRESHOLD=$4
    SAMPLES=$5
    LIMIT=$6
    THREADS=$7
    PREPROCESS=$8
else
    echo 'Usage: ./canu.sh [reads file] freq_minimizers k threshold samples limit threads preprocess'
    echo 'Example: ./canu.sh chr22.fastq 8 20 85 -3 100 8 5000'
    exit
fi

javac *.java
readlistfile=`java HashContainment $READS_FILE $FREQ_MINIMIZERS $K $THRESHOLD seed=$SAMPLES limit=$LIMIT threads=24 preprocess=5000`
echo 'Extracting reads from: ' $readlistfile
java ExtractReads $readlistfile $READS_FILE
newreadsfile=$readlistfile'.fastq'
#rm -r 'canu_'$newreadsfile
#canu-1.8/*/bin/canu -d 'canu_'$newreadsfile -p chr22 genomeSize=35m -useGrid=false -stopOnLowCoverage=1 -pacbio-corrected $newreadsfile
