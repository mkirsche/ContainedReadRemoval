assemblydir=$1
reffile=$2
buscolineage=$3
for i in `ls $assemblydir`
do
    echo 'Evaluating ' $i
    ./eval.sh $assemblydir/$i $reffile $i_stats -l $buscolineage
done

