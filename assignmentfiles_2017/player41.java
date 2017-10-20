// package ec;

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
		// SPHERE
		if (hasStructure && isSeparable && !isMultimodal) {
			doThis = 0;
		}
		// BENT CIGAR
		else if (!hasStructure && !isSeparable && !isMultimodal) {
			doThis = 1;
		}
		// SCHAFFERS
		else if (hasStructure && !isSeparable && isMultimodal) {
			doThis = 2;
		}
		// KATSUURA
		else if (!hasStructure && !isSeparable && isMultimodal) {
			doThis = 3;
		}

  }

	public void run()
	{
		// COMPETITION:
		// SPHERE
		if (doThis == 0) {
			hillClimber(randomStart(10), true);
		}
		// BENT CIGAR
		else if (doThis == 1) {
			hillClimber(randomStart(10), true);
		}
		// SCHAFFERS
		else if (doThis == 2) {
			plantPropagation(100,20,10,true);
		}
		// KATSUURA
		else if (doThis == 3) {
			plantPropagation(500,10,10,true);
		}

	}

	public void hillClimber(double[] start, boolean randomStart)
	{
		int evals = 0;
		double currBest[] = start;
		double currBestFitness = (double) evaluation_.evaluate(currBest);
		evals++;

		// Test randoms for 5% of evals
		// System.out.println("Start Random search");
		if(randomStart) {
			for(int i = 0; i < evaluations_limit_ * 0.05; i++, evals++)
			{
				double child[]= randomStart(10);

				double fitness = (double) evaluation_.evaluate(child);

				if (fitness >= currBestFitness) {
					currBestFitness = fitness;
					currBest = child;
				}
			}
		}
		// System.out.println("Evals done: " + evals + "\t currBest: " + currBestFitness);

		// Hillclimb
		// System.out.println("Start Hillclimber");
		for(int i = 0; evals < evaluations_limit_ *2 ; i++, evals++)
		{
		    double child[];
				do
				{
				 	child = randomArray(10, currBestFitness);
					child = sumArray(currBest, child);
				} while(!verify(child));

		    double fitness = (double) evaluation_.evaluate(child);

		    if (fitness > currBestFitness) {
		    	currBestFitness = fitness;
					currBest = child;

					// System.out.println("Evals done: " + evals + "\t Iterations done:" + i + "\t currBest: " + currBestFitness);
		    }

		}

	}

	public double[] randomArray(int n, double fitness){
			double arr[] = new double[n];
			double change;
			double rand;


			if (fitness <= 9.999 && doThis == 3){
				change = 9.999;
			} else if (fitness <= 1.) {
				change = 1.;
			} else {
				// fitness *= fitness;
				// fitness = 1 << (int)fitness;
				// fitness = factorial((int)fitness);
				// fitness = gammaFact(fitness);
				change = Math.pow(1.6, fitness);
			}

			for (int i = 0; i < n; i++) {
				rand = rnd_.nextGaussian();

				if (1. / change < (10. - fitness)) {
					arr[i] = rand / change;
				} else {
					arr[i] = rand * (10. - fitness) / 2.16;
				}

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

	public ArrayList<double[]> createRunners(double[] origin, double fitness, int maxRunners, int dimensions)
	{
		ArrayList<double[]> children = new ArrayList<double[]>();

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

	public void gradientAscent(double[] givenState)
	{
		double maxIterations = evaluations_limit_ / 11;

		// double[] oldState = randomStart(10);
		double[] oldState = givenState;
		double oldFitness = (double) evaluation_.evaluate(oldState);

		int evals = 1;

		double[] newState = sumArray(oldState, randomArray(10, oldFitness));
		double newFitness;

		double[] oldGradient= new double[10];
		double[] newGradient = new double[10];

		// System.out.println("Start Gradient Ascent");
		for (int i = 0; i < maxIterations; i++, evals += 11) {
			newFitness = (double) evaluation_.evaluate(newState);

			// upon convergence, break
			if (newFitness == 10.0) {
				break;
			}

			double dy = newFitness - oldFitness;

			// Prevent devision by 0. NOTE: dit kan ook anders
			if (dy != 0.0)
			{
				// calculate gradient
				newGradient = calculateGradients(newState, newFitness);

				// Shift back the new values as old
				oldGradient = copyArray(newGradient);
				oldState = copyArray(newState);

				// if (newFitness != oldFitness) {
				// 	System.out.println("Evals done: " + evals + "\t Iterations done:" + i + "\t currBest: " + newFitness);
				// }

				oldFitness = newFitness;

				newState = getNewState(newState, newGradient, newFitness);

				// The new state must be correct.
				if(!verify(newState)){
					newState = randomStart(10);
				}
			} else {
				newState = randomStart(10);
			}

		}
	}

	private double[] calculateGradients(double[] state, double fitness) {
		double[] gradient = new double[10];
		double change;

		change = 0.001;

		// System.out.println("Begin \t change: " + change);

		for (int i = 0; i < state.length; i++) {

			state[i] += change;
			gradient[i] = ((double) evaluation_.evaluate(state) - fitness) / change;
			// System.out.println(gradient[i]);
			state[i] -= change;
		}

		return gradient;
	}


	private double[] getNewState(double[] oldState, double[] gradient, double fitness)
	{
		double[] newState = new double[10]; // copy vector

		double change = (10 - fitness) / 100.;

		for (int i = 0; i < oldState.length; i++) {

			double x = oldState[i]; // remember value

			// adjust value by gradient
			if (gradient[i] != 0.0) {
				newState[i] = oldState[i] + change * gradient[i];
			}
		}
		return newState;
	}

	public double[] concat(double[] a, double[] b) {
	   int aLen = a.length;
	   int bLen = b.length;
	   double[] c= new double[aLen+bLen];
	   System.arraycopy(a, 0, c, 0, aLen);
	   System.arraycopy(b, 0, c, aLen, bLen);
	   return c;
	}


	public void plantPropagation(int a, int b, int c, boolean optimise)
	{
		int dimensions = 10;
		// SPPA parameters
		int startPopSize = a;
		int popSelection = b;
		int maxRunners = c;
		int generations = evaluations_limit_;
		int currEvals = 0;

		ArrayList<double[]> population = createpopulation(startPopSize, dimensions);
		double fitness[] = new double[10];
		double parentFitness[] = new double[10];
		ArrayList<double[]> newPopulation = new ArrayList<double[]>();
		ArrayList<double[]> parentPopulation = new ArrayList<double[]>();

		// System.out.println("Start PPA \t optimise == " + optimise);
		for(int i = 0; i < generations; i++)
		{
			fitness = getFitnessPopulation(population, population.size());
			currEvals += population.size();

			fitness = concat(fitness, parentFitness);
			population.addAll(parentPopulation);

			parentFitness = new double[0];

			// Sort population by fitness descending
			sortPopulation(population, fitness, population.size());

			// System.out.println("Evals done: " + currEvals + "\t Iterations done:" + i + "\t currBest: " + fitness[0]);

			// We're done here.
			if((fitness[0] >= 9.9 || currEvals > evaluations_limit_ * 0.8) && optimise )
			{
				break;
			}

			// This population sucks major ass. Lets rebuild society...
			if(fitness[0] <= 0.00001)
			{
				population = createpopulation(startPopSize, dimensions);
			}

			newPopulation.clear();
			parentPopulation.clear();

			// limit the amount of new stuff...
			for(int j = 0; j < Math.min(population.size(), popSelection); j++)
			{
				ArrayList<double[]> runners = createRunners(population.get(j), fitness[j], maxRunners, dimensions);

				// New population
				parentPopulation.add(population.get(j));
				double[] temp = {fitness[j]};
				parentFitness = concat(parentFitness, temp);

				newPopulation.addAll(runners);
			}

			// Create new copy.
			population.clear();
			for (int j = 0; j < newPopulation.size(); j++)
			{
				population.add(copyArray(newPopulation.get(j)));
			}

		}

		// System.out.println("Fitness after PPA: " + fitness[0]);

		if (optimise) {
			// System.out.println("Start Optimising");
			hillClimber(population.get(0),false);
			// gradientAscent(population.get(0));
		}
	}



}
