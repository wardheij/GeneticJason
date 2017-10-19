javac -cp 'contest.jar' ec/*.java
jar cmf MainClass.txt submission.jar ec/*.class
echo "Bent Cigar: "
java -jar testrun.jar -submission=ec.player41 -evaluation=BentCigarFunction -seed=1
