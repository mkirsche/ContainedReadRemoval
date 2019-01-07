javac *.java

WTDIR='wtdbg2'

if [ "$#" -eq 6 ]; then
    READS_FILE=$1
    W1=$2
    K1=$3
    W2=$4
    K2=$5
    CT=$6
else
    echo 'Usage: ./error.sh <reads_file> <w1> <k1> <w2> <k2> <ct>'
    echo 'Example: ./error.sh chr22.fastq 5 8 50 14 10'
    exit
fi

java PB_FilterContainedReads $READS_FILE w1=$W1 k1=$K1 w2=$W2 k2=$K2 ct=$CT
ofn=`java PB_FilterContainedReads $READS_FILE w1=$W1 k1=$K1 w2=$W2 k2=$K2 ct=$CT fnonly`
java ExtractReads $ofn $READS_FILE

NEW_READS_FILE=$ofn'.fastq'
$WTDIR/wtdbg2 -t 16 -i $NEW_READS_FILE -L 5000 -fo $WTDIR'/'$NEW_READS_FILE
gunzip $WTDIR'/'$NEW_READS_FILE'.ctg.lay.gz'; 
$WTDIR/wtpoa-cns -t 16 -i $WTDIR'/'$NEW_READS_FILE'.ctg.lay' -fo $WTDIR'/'$NEW_READS_FILE'.ctg.lay.fa'



