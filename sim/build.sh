BINDIR=`dirname $(readlink -f "$0")`
javac -cp "$BINDIR:$BINDIR/commons-math3-3.6.1.jar" SimulateReadsFromGenome.java 

java -cp "$BINDIR:$BINDIR/commons-math3-3.6.1.jar" SimulateReadsFromGenome
