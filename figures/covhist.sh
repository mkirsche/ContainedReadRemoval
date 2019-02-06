BINDIR=`dirname $(readlink -f "$0")`
MINIMAP='/scratch/groups/mschatz1/mkirsche/github/minimap2/minimap2'
READ_FILE=$1
REF_FILE=$2
PAF_FILE=$READ_FILE'.paf'
COV_FILE=$READ_FILE'.cov.txt'

$MINIMAP $REF_FILE $READ_FILE > $PAF_FILE
javac $BINDIR/*.java
java -cp $BINDIR CoverageDistribution $PAF_FILE > $COV_FILE
python coverage_histogram $COV_FILE $COV_FILE'.png'
