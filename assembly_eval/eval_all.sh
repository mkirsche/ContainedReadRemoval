assemblydir=$1
reffile=$2
buscolineage=$3

BINDIR=`dirname $(readlink -f "$0")`

mode='normal'

while getopts d:r:l:m: option
do
    case "${option}"
        in
        d) assemblydir=${OPTARG};;
        r) reffile=${OPTARG};;
        l) buscolineage=${OPTARG};;
        m) mode=${OPTARG}
    esac
done

echo $mode

ls $assemblydir | parallel --gnu --jobs 4 $BINDIR/eval.sh $assemblydir'/'{} $reffile {}'_stats' -l $buscolineage -m $mode

