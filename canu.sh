READS_FILE='chr22.fastq'
FREQ_MINIMIZERS=$1
K=$2
THRESHOLD=$3
SAMPLES=$4

javac *.java
readlistfile=`java HashContainment $READS_FILE $FREQ_MINIMIZERS $K $THRESHOLD seed=$SAMPLES`
java ExtractReads $readlistfile $READS_FILE
newreadsfile=$readlistfile'.fastq'
rm -r 'canu_'$newreadsfile
canu-1.8/*/bin/canu -d 'canu_'$newreadsfile -p chr22 genomeSize=35m -useGrid=false -stopOnLowCoverage=1 -pacbio-corrected $newreadsfile
