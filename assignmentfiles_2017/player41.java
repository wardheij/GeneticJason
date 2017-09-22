import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import static java.lang.System.out;

public class player41 implements ContestSubmission
{
	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;

	public player41()
	{
		rnd_ = new Random();
	}

	public void setSeed(long seed)
	{
		// Set seed of algortihms random process
		rnd_.setSeed(seed);
	}

	public void setEvaluation(ContestEvaluation evaluation)
	{
		// Set evaluation problem used in the run
		evaluation_ = evaluation;

		// Get evaluation properties
		Properties props = evaluation.getProperties();
    // Get evaluation limit
    evaluations_limit_ = Integer.parseInt(props.getProperty("Evaluations"));
		// Property keys depend on specific evaluation
		// E.g. double param = Double.parseDouble(props.getProperty("property_name"));
    boolean isMultimodal = Boolean.parseBoolean(props.getProperty("Multimodal"));
    boolean hasStructure = Boolean.parseBoolean(props.getProperty("Regular"));
    boolean isSeparable = Boolean.parseBoolean(props.getProperty("Separable"));

		// Do sth with property values, e.g. specify relevant settings of your algorithm
    if(isMultimodal){
      out.println("Multimodal");
    }
		if(hasStructure){
			out.println("Regular");
		}
		if(isSeparable){
			out.println("Separable");
		}
  }

	public void run()
	{
		// Run your algorithm here

		hillClimber();

    // int evals = 0;
    // // init population
    // // calculate fitness
    // while(evals<evaluations_limit_){
    //     // Select parents
    //     // Apply crossover / mutation operators
    //     double child[] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
    //     // Check fitness of unknown fuction
    //     Double fitness = (double) evaluation_.evaluate(child);
    //     evals++;
    //     // Select survivors
    // }

	}

	public void hillClimber()
	{
		int evals = 0;
		double currBest[] = randomStart(10);
		double currBestFitness = (double) evaluation_.evaluate(currBest);
		evals++;

		// Test randoms for 1% of evals
		for(int i = 1; i < evaluations_limit_ * 0.01; i++, evals++)
		{
			double child[]= randomStart(10);

			double fitness = (double) evaluation_.evaluate(child);

			if (fitness >= currBestFitness) {
				currBestFitness = fitness;
				currBest = child;
			}
		}

		// Hillclimb
		for(; evals < evaluations_limit_; evals++)
		{
		    double child[];
				do
				{
				 	child = randomArray(10);
					child = sumArray(currBest, child);
				} while(!verify(child));

		    double fitness = (double) evaluation_.evaluate(child);

		    if (fitness >= currBestFitness) {
		    	currBestFitness = fitness;
					currBest = child;
		    }
		}

	}

	public double[] randomArray(int n){
			double arr[] = new double[n];

			for (int i = 0; i < n; i++) {
				// arr[i] = rnd_.nextDouble() - 0.5;
				arr[i] = rnd_.nextGaussian()*0.1;
			}

			return arr;
	}

	public double[] sumArray(double[] a, double[] b)
	{
		for (int i = 0; i < a.length; i++) {
			b[i] += a[i];
		}

		return b;
	}

	public double[] randomStart(int n)
	{
			double arr[] = new double[n];

			for (int i = 0; i < n; i++) {
				arr[i] = (rnd_.nextDouble() - 0.5) * 10;
			}

			return arr;
	}

	public boolean verify(double[] a)
	{
		for (int i = 0; i < a.length; i++) {
			if(a[i] < -5.0 || a[i] > 5.0)
			{
				return false;
			}
		}
		return true;
	}

}
