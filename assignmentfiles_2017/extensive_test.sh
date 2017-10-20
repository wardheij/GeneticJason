javac -cp 'contest.jar' *.java
jar cmf MainClass.txt submission.jar *.class

for (( i = 0; i < 10; i++ )); do
  echo "Bent Cigar: "
  java -jar testrun.jar -submission=player41 -evaluation=BentCigarFunction -seed=$i
  echo ""
  echo "Schaffers: "
  java -jar testrun.jar -submission=player41 -evaluation=SchaffersEvaluation -seed=$i
  echo ""
  echo "Katsuura: "
  java -jar testrun.jar -submission=player41 -evaluation=KatsuuraEvaluation -seed=$i
  echo ""
done
