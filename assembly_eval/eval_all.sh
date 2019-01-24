assemblydir=$1
reffile=$2
buscolineage=$3

BINDIR=`dirname $(readlink -f "$0")`

mode='normal'

while getopts d:r:l:m:f: option
do
    case "${option}"
        in
        d) assemblydir=${OPTARG};;
        r) reffile=${OPTARG};;
        l) buscolineage=${OPTARG};;
        m) mode=${OPTARG};;
        f) fileassemblylist=${OPTARG}
    esac
done

echo $mode

if [ -z "${fileassemblylist}" ]
then
    ls $assemblydir | parallel --gnu --jobs 4 $BINDIR/eval.sh -a $assemblydir'/'{} -r $reffile -o {}'_stats' -l $buscolineage -m $mode
else
    cat $fileassemblylist | parallel --gnu --jobs 4 $BINDIR/eval.sh -a $assemblydir'/'{} -r $reffile -o {}'_stats' -l $buscolineage -m $mode
fi

