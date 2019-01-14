assemblydir=$1
reffile=$2
buscolineage=$3

BINDIR=`dirname $(readlink -f "$0")`

for i in `ls $assemblydir`
do
    echo 'Evaluating ' $i
    $BINDIR/eval.sh $assemblydir'/'$i $reffile $i'_stats' -l $buscolineage
done

