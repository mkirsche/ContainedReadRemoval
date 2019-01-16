assemblydir=$1
reffile=$2
buscolineage=$3

BINDIR=`dirname $(readlink -f "$0")`

ls $assemblydir | parallel --gnu --jobs 4 $BINDIR/eval.sh $assemblydir'/'{} $reffile {}'_stats' -l $buscolineage

