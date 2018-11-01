javac *.java
for i in `ls -d canu_*/`; do java GetCanuResults $i; echo ''; done
