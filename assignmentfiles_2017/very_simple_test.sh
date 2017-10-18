javac -cp 'contest.jar' ec/*.java
jar cmf MainClass.txt submission.jar ec/*.class
echo "Sphere: "
java -jar testrun.jar -submission=ec.player41 -evaluation=SphereEvaluation -seed=1
