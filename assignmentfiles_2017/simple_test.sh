javac -cp 'math.jar;contest.jar' player41.java
jar cmf MainClass.txt submission.jar player41.class
echo "Sphere: "
java -jar testrun.jar -submission=player41 -evaluation=SphereEvaluation -seed=1
echo ""
echo "Bent Cigar: "
java -jar testrun.jar -submission=player41 -evaluation=BentCigarFunction -seed=1
echo ""
echo "Schaffers: "
java -jar testrun.jar -submission=player41 -evaluation=SchaffersEvaluation -seed=1
echo ""
echo "Katsuura: "
java -jar testrun.jar -submission=player41 -evaluation=KatsuuraEvaluation -seed=1
