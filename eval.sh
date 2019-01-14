# Performs assembly evaluation on a single assembly

quastfile=/home-3/mkirsche@jhu.edu/build/quast-5.0.1/quast-lg.py
buscofile=/home-3/mkirsche@jhu.edu/build/busco/scripts/run_BUSCO.py
buscolineage=/home-3/mkirsche@jhu.edu/build/busco/embryophyta_odb9

BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`

assembly=$1
ref=$2
outdir=$3

while getopts q:b:l option
do
    case "${option}"
        in
        q) quastfile=${OPTARG};;
        b) buscofile=${OPTARG};;
        l) buscolineage=${OPTARG};;
    esac
done

mkdir $outdir

if [ ! -f "$BINDIR/AssemblyStats.class" ]
then
    javac $BINDIR/AssemblyStats.java
fi

java -cp $BINDIR AssemblyStats $WORKINGDIR/$assembly | tee $outdir/assemblystats.txt

cd $outdir
python $buscofile -i $WORKINGDIR/$assembly -l $buscolineage -o busco -m genome
cd $WORKINGDIR

$quastfile -o $WORKINGDIR/$outdir/quast -r $ref $WORKINGDIR/$assembly

