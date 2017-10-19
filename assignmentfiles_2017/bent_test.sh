javac -cp 'contest.jar' *.java
jar cmf MainClass.txt submission.jar *.class
echo "Bent Cigar: "
java -jar testrun.jar -submission=player41 -evaluation=BentCigarFunction -seed=1
