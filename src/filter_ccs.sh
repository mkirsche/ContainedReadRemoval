if [ "$#" -eq 8 ]; then
    READS_FILE=$1
    FREQ_MINIMIZERS=$2
    K=$3
    THRESHOLD=$4
    SAMPLES=$5
    LIMIT=$6
    THREADS=$7
    PREPROCESS=$8
else
    echo 'Usage: $0 reads file freq_minimizers k threshold samples limit threads preprocess'
    echo 'Example: $0 chr22.fastq 8 20 85 -3 100 8 5000'
    exit
fi

javac *.java
readlistfile=`java -Xmx64G HashContainment $READS_FILE $FREQ_MINIMIZERS $K $THRESHOLD seed=$SAMPLES limit=$LIMIT threads=$THREADS preprocess=$PREPROCESS`
java ExtractReads $readlistfile $READS_FILE
echo $readlistfile'.fastq'
