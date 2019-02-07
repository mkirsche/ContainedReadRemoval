BINDIR=`dirname $(readlink -f "$0")`
MINIMAP='/scratch/groups/mschatz1/mkirsche/github/minimap2/minimap2'

while getopts r:g:p: opt; do
    case $opt in 
        p) PAF_FILE=${OPTARG} ;;
        r) READ_FILE=${OPTARG} ;;
        g) REF_FILE=${OPTARG} ;;
    esac
done

if [ ! -z $PAF_FILE ]
then
    echo 'Using existing paf file'
else
    PAF_FILE=$READ_FILE'.paf'
    $MINIMAP $REF_FILE $READ_FILE > $PAF_FILE
fi

COV_FILE=$PAF_FILE'.cov.txt'

javac $BINDIR/*.java
java -cp $BINDIR CoverageDistribution $PAF_FILE > $COV_FILE
python $BINDIR/coverage_histogram.py $COV_FILE $PAF_FILE'.png'
