BINDIR=`dirname $(readlink -f "$0")`
javac -cp "$BINDIR:$BINDIR/commons-math3-3.6.1.jar" SimulateReadsFromGenome.java 

java -cp "$BINDIR:$BINDIR/commons-math3-3.6.1.jar" SimulateReadsFromGenome coverage=30 maxLength=5000000

java -cp ../src PB_FilterContainedReads simulatedreads.fa nt=2 ct=.022 ofn=uncontained.txt k1=15
