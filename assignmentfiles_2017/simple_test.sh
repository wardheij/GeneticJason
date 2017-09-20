javac -cp contest.jar player0.java
jar cmf MainClass.txt submission.jar player0.class
echo "Sphere: "
java -jar testrun.jar -submission=player0 -evaluation=SphereEvaluation -seed=1
echo ""
echo "Bent Cigar: "
java -jar testrun.jar -submission=player0 -evaluation=BentCigarFunction -seed=1
echo ""
echo "Schaffers: "
java -jar testrun.jar -submission=player0 -evaluation=SchaffersEvaluation -seed=1
echo ""
echo "Katsuura: "
java -jar testrun.jar -submission=player0 -evaluation=KatsuuraEvaluation -seed=1
