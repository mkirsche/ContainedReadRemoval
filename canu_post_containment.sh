READS_FILE='chr22.fastq'

if [ "$#" -eq 4 ]; then
    FREQ_MINIMIZERS=$1
    K=$2
    THRESHOLD=$3
    SAMPLES=$4
fi

if [ "$#" -eq 5 ]; then
    READS_FILE=$1
    FREQ_MINIMIZERS=$2
    K=$3
    THRESHOLD=$4
    SAMPLES=$4
fi

if [ "$#" -eq 6 ]; then
    READS_FILE=$1
    FREQ_MINIMIZERS=$2
    K=$3
    THRESHOLD=$4
    LIMIT=$5
fi

javac *.java
readlistfile=`java HashContainment $READS_FILE $FREQ_MINIMIZERS $K $THRESHOLD seed=$SAMPLES --fnOnly`
java ExtractReads $readlistfile $READS_FILE
newreadsfile=$readlistfile'.fastq'
rm -r 'canu_'$newreadsfile
canu-1.8/*/bin/canu -d 'canu_'$newreadsfile -p chr22 genomeSize=35m -useGrid=false -stopOnLowCoverage=1 -pacbio-corrected $newreadsfile
