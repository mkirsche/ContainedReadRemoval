usage() 
{ 
    echo "Usage: $0 -r <readfile> -t <readtype> -o <outfile>" 1>&2;
    echo "  -t <readtype> has value pacbio, nanopore, or ccs" 1>&2;
    exit 1; 

}

while getopts t:o:r:a option
do
    case "${option}"
        in
        o) outfile=${OPTARG};;
        r) readfile=${OPTARG};;
        t) readtype=${OPTARG};;
        a) assembler=${OPTARG};;
    esac
done

if [ -z "${outfile}" ] || [ -z "${readfile}" ] || [ -z "${readtype}" ]; then
    usage
fi

if [ "$readtype" = "ccs" ]; then
    newreadsfile=`./filter_ccs.sh 8 25 95 3 100 24 5000`
elif [ "$readtype" = "pacbio" ]; then
    # TODO: Find better parameters for Raw Pacbio reads
    newreadsfile=`./filter_pb_np.sh $readfile 5 12 50 18 25`
elif [ "$readtype" = "nanopore" ]; then
    newreadsfile=`./filter_pb_np.sh $readfile 5 50 130 19 24`
else
    usage  
fi
