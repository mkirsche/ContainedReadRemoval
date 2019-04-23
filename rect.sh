len=$1
ct=$2
name=$3
javac src/*.java; time java -cp src PB_FilterContainedReads ../ERR2173373.fastq nt=24 ct=$ct w1=5 ofn=/home-3/mkirsche@jhu.edu/hashing/CCS/results/$name.txt dfn=/home-3/mkirsche@jhu.edu/hashing/CCS/results/$name.debug.txt k1=15 rt=.6 lf=$len method=rectangle  2>&1 | tee /home-3/mkirsche@jhu.edu/hashing/CCS/results/$name.log.txt

java -cp /home-3/mkirsche@jhu.edu/hashing/CCS/src ExtractReads /home-3/mkirsche@jhu.edu/hashing/CCS/results/$name.txt /home-3/mkirsche@jhu.edu/hashing/ERR2173373.fastq

/home-3/mkirsche@jhu.edu/hashing/CCS/src/assemble.sh -r /home-3/mkirsche@jhu.edu/hashing/CCS/results/$name.txt.fastq -o fuzzy_$name

java -cp /home-3/mkirsche@jhu.edu/hashing/CCS/assembly_eval AssemblyStats wtdbg2_assemblies/$name.txt.fastq.ctg.lay.fa | tee fuzzy_$name.stats

