BINDIR=`dirname $(readlink -f "$0")`
WORKINGDIR=`pwd`
assembler='wtdbg2'
wtdir='/scratch/groups/mschatz1/mkirsche/hashing/wtdbg2'
canufile='$BINDIR/../canu-1.8/*/bin/canu'
usage() 
{ 
    echo "Usage: $0 -r <readfile> -o <outfile> [-c -t <readtype> -l <genome length>]" 1>&2;
    echo "  -c is an option for using the Canu assembler" 1>&2;
    echo "  If using Canu, the following parameters are also needed: " 1>&2;
    echo "    -t <readtype> having value pacbio-raw, nanopore-raw, or pacbio-corrected" 1>&2;
    echo "    -l <genome length>" 1>&2
    exit 1; 

}

while getopts c:o:r:t:l option
do
    case "${option}"
        in
        c) assembler='canu';;
        o) outfile=${OPTARG};;
        r) readfile=${OPTARG};;
        t) readtype=${OPTARG};;
        l) length=${OPTARG};;
    esac
done

if [ -z "${outfile}" ] || [ -z "${readfile}" ]; then
    usage
fi

basenamereadfile=${readfile##*/}

if [ "$assembler" = "canu" ]; then
    if [ -z "${length}" ] || [ -z "${readtype}" ]; then
        usage
    fi
    OUTDIR='canu_'$outfile
    if [ -d $OUTDIR ]; then
        rm -r $OUTDIR
    fi
    $canufile -d $OUTDIR -p $outfile genomeSize=$length -useGrid=false -stopOnLowCoverage=1 -$type $readfile
else
    OUTDIR=$WORKINGDIR'/wtdbg2_assemblies'
    mkdir -p $OUTDIR
    $wtdir/wtdbg2 -t 16 -i $readfile -L 5000 --rescue-low-cov-edges -fo $OUTDIR'/'$basenamereadfile
    gunzip -f $OUTDIR'/'$basenamereadfile'.ctg.lay.gz'; 
    $wtdir/wtpoa-cns -t 16 -i $OUTDIR'/'$basenamereadfile'.ctg.lay' -fo $OUTDIR'/'$basenamereadfile'.ctg.lay.fa'
fi


