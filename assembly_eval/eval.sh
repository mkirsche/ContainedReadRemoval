# Performs assembly evaluation on a single assembly

quastfile=/home-3/mkirsche@jhu.edu/build/quast-5.0.1/quast-lg.py
buscofile=/home-3/mkirsche@jhu.edu/build/busco/scripts/run_BUSCO.py
buscolineage=/home-3/mkirsche@jhu.edu/build/busco/embryophyta_odb9

BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`

mode='normal'

while getopts a:r:o:q:b:l:m: option
do
    case "${option}"
        in
        a) assembly=${OPTARG};;
        r) ref=${OPTARG};;
        o) outdir=${OPTARG};;
        q) quastfile=${OPTARG};;
        b) buscofile=${OPTARG};;
        l) buscolineage=${OPTARG};;
        m) mode=${OPTARG}
    esac
done

echo 'Evaluating: '$assembly
echo '  Reference genome: '$ref
echo '  Output directory: '$outdir
echo '  Busco lineage: '$buscolineage
echo '  Mode: '$mode

if [ -d $outdir ]
then
    rm -r $outdir
fi

mkdir $outdir

if [ ! -f "$BINDIR/AssemblyStats.class" ]
then
    javac $BINDIR/AssemblyStats.java
fi

java -cp $BINDIR AssemblyStats $WORKINGDIR/$assembly | tee $outdir/assemblystats.txt

if [ "$mode" = "fast" ]; then
    exit 0;
fi

cd $outdir
python $buscofile -f --blast_single_core -i $WORKINGDIR/$assembly -l $buscolineage -o busco -m genome
cd $WORKINGDIR

$quastfile -o $WORKINGDIR/$outdir/quast -r $ref $WORKINGDIR/$assembly

