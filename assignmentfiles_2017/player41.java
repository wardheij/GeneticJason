import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;

import java.util.Random;
import java.util.Properties;
import java.util.ArrayList;
import static java.lang.System.out;

public class player41 implements ContestSubmission
{
	Random rnd_;
	ContestEvaluation evaluation_;
    private int evaluations_limit_;

	int doThis;

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

		doThis = 0;
		// Do sth with property values, e.g. specify relevant settings of your algorithm
    if(isMultimodal){
      out.println("Multimodal");
			doThis = 0;
    }
		if(hasStructure){
			out.println("Regular");
			doThis = 0;
		}
		if(isSeparable){
			out.println("Separable");
		}


  }

	public void run()
	{
		// Run your algorithm here
		if(doThis == 0)
		{
			//hillClimber();

			// NOTE: This is basically a hillclimber with 100 startingpoints..
			// The best one gets taken and gets 10 children etc.
			plantPropagation(100,1,10);
		}
		else if (doThis == 1)
		{
			plantPropagation(100,100,5);
		}
		else if (doThis == 2)
		{
			// NOTE: can't get this to work currently. Imports.
			fireworks();
		}

		// NOTE: THINGS WE NEED ASAP:
		// Gradient ascent: https://en.wikipedia.org/wiki/Gradient_descent

		// NOTE: other interesting (and quite easily implemented) options include:
		// CMA-ES: https://en.wikipedia.org/wiki/CMA-ES
		// ACO: https://en.wikipedia.org/wiki/Ant_colony_optimization_algorithms
		// PSO: https://en.wikipedia.org/wiki/Particle_swarm_optimization
		// FA: https://en.wikipedia.org/wiki/Firefly_algorithm
		// EDA: https://en.wikipedia.org/wiki/Estimation_of_distribution_algorithm
	}

	public void fireworks()
	{
		double mins[] = {-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0,-5.0};
		double maxs[] = {5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0,5.0};
		FA fire = new FA(10,10,1,1,1,10,maxs,mins,"",evaluation_);
		fire.FAframework();
	}

	public void hillClimber()
	{
		int evals = 0;
		double currBest[] = randomStart(10);
		double currBestFitness = (double) evaluation_.evaluate(currBest);
		evals++;

		// Test randoms for 1% of evals
		for(int i = 0; i < evaluations_limit_ * 0.01; i++, evals++)
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
				 	child = randomArray(10, currBestFitness);
					child = sumArray(currBest, child);
				} while(!verify(child));

		    double fitness = (double) evaluation_.evaluate(child);

		    if (fitness >= currBestFitness) {
		    	currBestFitness = fitness;
					currBest = child;
		    }
		}

	}

	public double[] randomArray(int n, double fitness){
			double arr[] = new double[n];

			if (fitness <= 1){
				fitness = 1;
			} else {
				// fitness *= fitness;
				// fitness = 1 << (int)fitness;
				// fitness = factorial((int)fitness);
				// fitness = gammaFact(fitness);
				fitness = Math.pow(1.5, fitness);
			}

			for (int i = 0; i < n; i++) {
				// arr[i] = (rnd_.nextDouble() - 0.5)*2 / fitness;

				arr[i] = rnd_.nextGaussian() / fitness;
			}

			return arr;
	}

	public double[] sumArray(double[] a, double[] b)
	{
		for (int i = 0; i < a.length; i++) {
			b[i] += a[i];
			// b[i] = (double)Math.round(b[i] * 10000d) / 10000d;
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

	public int factorial(int n)
	{
		int total = 1;
		for (int i = 1; i < n; i++) {
			total *= i;
		}
		return total;
	}

	public double gammaFact(double n)
	{
		return Math.sqrt(2 * Math.PI * n)* Math.pow(n / Math.E, n);
	}

	public ArrayList<double[]> createpopulation(int popSize, int dimension)
	{
		ArrayList<double[]> population = new ArrayList<double[]>();
		for(int i = 0; i < popSize; i++)
		{
			population.add(randomStart(dimension));
		}
		return population;
	}

	public double[] getFitnessPopulation(ArrayList<double[]> population, int popSize)
	{
		double fitness[] = new double[popSize];
		for(int i = 0; i < popSize; i++)
		{
			fitness[i] = (double) evaluation_.evaluate(population.get(i));
			if(fitness[i] < 0.0)
			{
				fitness[i] = 0.0;
			}
		}
		return fitness;
	}

	public void sortPopulation(ArrayList<double[]> population, double[] fitness, int popSize)
	{
		double tmp[];
		double temp;

		// Sinksort
		for(int i = 0; i < popSize; i++)
		{
			for(int j = 1; j < popSize - i; j++)
			{
				if(fitness[j - 1] <= fitness[j])
				{
					temp = fitness[j - 1];
					fitness[j - 1] = fitness[j];
					fitness[j] = temp;

					tmp = population.get(j - 1);
					population.set(j - 1, population.get(j));
					population.set(j, tmp);
				}
			}
		}
	}

	public double[] getD(int n, double fitness){
			double arr[] = new double[n];

			for (int i = 0; i < n; i++) {
				// arr[i] = 2 * (1 - fitness) * (rnd_.nextDouble() - 0.5);
				// arr[i] = rnd_.nextGaussian() / Math.pow(1.5, fitness);

			}

			return arr;
	}

	public ArrayList<double[]> createRunners(double[] origin, double fitness, int maxRunners, int dimensions)
	{
		ArrayList<double[]> children = new ArrayList<double[]>();

		//NOTE: unsure about this part in the paper... can't get the actual thing to
		// so I just did this..
		// double correctedFitness = (1. / 2.) * (Math.tanh(4. * (fitness / 10.) - 2.) + 1.);
		double correctedFitness = fitness / 10.;
		int numberOfRunners = (int)Math.ceil(maxRunners * correctedFitness * rnd_.nextDouble());

		if (numberOfRunners == 0) {
			numberOfRunners = 1;
		}

		for(int i = 0; i < numberOfRunners; i++)
		{
			double child[] = new double[dimensions];

			do
			{
				// child = sumArray(origin, getD(10, fitness));
				child = sumArray(origin, randomArray(10, fitness));

			} while(!verify(child));

			children.add(child);
		}

		return children;
	}

	public double[] copyArray(double[] in)
	{
		double out[] = new double[in.length];

		for(int i = 0; i < in.length; i++)
		{
			out[i] = in[i];
		}

		return out;
	}


	public void plantPropagation(int a, int b, int c)
	{
		int dimensions = 10;
		// SPPA parameters
		int startPopSize = a;
		int popSelection = b;
		int maxRunners = c;
		int generations = evaluations_limit_;

		ArrayList<double[]> population = createpopulation(startPopSize, dimensions);
		double fitness[];
		ArrayList<double[]> newPopulation = new ArrayList<double[]>();

		for(int i = 0; i < generations; i++)
		{
			//NOTE: we do have the fitness for some of these already, TODO; fix this
			fitness = getFitnessPopulation(population, population.size());

			// Sort population by fitness descending
			sortPopulation(population, fitness, population.size());
			// out.println("Generation: " + i + "   Best fitness: " + fitness[0] + "   Current population: " + population.size());
			// out.print("Best: [");
			// for (int j = 0; j < dimensions - 1; j++)
			// {
			// 	out.print(population.get(0)[j] + ", ");
			// }
			// out.print(population.get(0)[dimensions - 1] + "]\n");

			// We're done here.
			if(fitness[0] == 10.0)
			{
				break;
			}

			// This population sucks major ass. Lets rebuild society...
			if(fitness[0] <= 0.00001)
			{
				population = createpopulation(startPopSize, dimensions);
			}

			newPopulation.clear();

			// limit the amount of new stuff...
			// NOTE: We could add probabilities here so that lower scoring population also has a chance.
			for(int j = 0; j < Math.min(population.size(), popSelection); j++)
			{
				ArrayList<double[]> runners = createRunners(population.get(j), fitness[j], maxRunners, dimensions);

				// New population
				// NOTE: normally you just add the runners only; weird results if you do that
				newPopulation.add(population.get(j));
				newPopulation.addAll(runners);
			}

			// Create new copy.
			population.clear();
			for (int j = 0; j < newPopulation.size(); j++)
			{
				population.add(copyArray(newPopulation.get(j)));
			}

		}
	}
}
