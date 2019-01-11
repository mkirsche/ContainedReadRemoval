# Performs assembly evaluation on a single assembly

quastfile=/home-3/mkirsche@jhu.edu/build/quast-5.0.1/quast-lg.py
buscofile=/home-3/mkirsche@jhu.edu/build/busco/scripts/run_BUSCO.py
buscolineage=/home-3/mkirsche@jhu.edu/build/busco/embryophyta_odb9

BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`

assembly=$1
ref=$2
outdir=$3

mkdir $outdir

if [ ! -f "$BINDIR/AssemblyStats.class" ]
then
    javac $BINDIR/AssemblyStats.java
fi

java -cp $BINDIR AssemblyStats $WORKINGDIR/$assembly

python $buscofile -i $WORKINGDIR/$assembly -l $buscolineage -o $WORKINGDIR/$outdir/busco -m genome

$quastfile -o $WORKINGDIR/$outdir/quast -r $ref $WORKINGDIR/$assembly

