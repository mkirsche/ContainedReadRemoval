BINDIR=`dirname $(readlink -f "$0")`
javac -cp "$BINDIR:$BINDIR/commons-math3-3.6.1.jar" SimulateReadsFromGenome.java 

java -cp "$BINDIR:$BINDIR/commons-math3-3.6.1.jar" SimulateReadsFromGenome coverage=30 maxLength=5000000 error=0.12 sample=ERsample.fastq

#java -cp "$BINDIR:$BINDIR/commons-math3-3.6.1.jar" SimulateReadsFromGenome coverage=30 maxLength=5000000 error=0.05 sample=ERsample.fastq

#java -cp "$BINDIR:$BINDIR/commons-math3-3.6.1.jar" SimulateReadsFromGenome coverage=30 maxLength=5000000 error=0.12

#javac ../src/*.java
#java -cp ../src PB_FilterContainedReads simulatedreads.fa nt=2 ct=.015 et=0.05 w1=5 ofn=uncontained.txt k1=15 rt=.1 lf=15000 method=rectangle logging

#java -cp ../src PB_FilterContainedReads simulatedreads.fa nt=2 ct=.21 ct2=.21 w1=5 ofn=uncontained.txt k1=15 rt=.4
