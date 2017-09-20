javac -cp contest.jar player0.java
jar cmf MainClass.txt submission.jar player0.class
java -jar testrun.jar -submission=player0 -evaluation=BentCigarFunction -seed=1
