package fr.inria.optimization.cmaes.examples;

import fr.inria.optimization.cmaes.CMAEvolutionStrategy;
import fr.inria.optimization.cmaes.fitness.IObjectiveFunction;

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
	private static final MAX_DIM = 10;
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
        else if (doThis == 3)
        {
            CMA_ES.optimze();
        }

		// NOTE: THINGS WE NEED ASAP:
		// Gradient ascent: https://en.wikipedia.org/wiki/Gradient_descent
		// CMA-ES: https://en.wikipedia.org/wiki/CMA-ES <<<<<<< WARD <<<<<<<

		// NOTE: other interesting (and quite easily implemented) options include:
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
			// TODO: random restart on !verify(child)


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

	public double[] gradientAscent(double[] oldState, double[] currentState, int maxIterations)
	{
		double[] state = currentState;
		double alpha = 0.01; // learning rate
		double oldGradient;
		double gradient;
		// TODO: init random oldState?
		// NOTE: Ja. Ik denk dat het voor de eerste iteratie gewoon het makkelijkst
		// is om een random positie te kiezen oid. Hebben we het even over.

		for (int i = 0; i < maxIterations, i++) {

			// upon convergence, break
			if (gradient - oldGradient < 0) {
				break;
			}
			// calculate gradient
			gradient = calculateGradient(oldState, state);

			// NOTE: dit hieronder werkt niet lekker in Java. Gebruik hier de copyArray() die ik gemaakt heb
			oldState = state;
			state = getNewState(state, gradient, alpha);
			oldGradient = gradient;
		}
	}

	// NOTE: Idealiter wil je dit zo min mogelijk doen.
	// Nu gebeuren er 4 evaluaties per dim per iteratie, terwijl er 1 nodig is
	// (er is immers maar 1 nieuwe state)
	// Daarnaast berekend dit niet de gradient, maar is dit de delta y.
	// (verandering over de score)
	// De gradient krijg je door de verandering in de dimensies te delen door
	// de verandering in score. dx/dy (oftewel de numerieke afgeleide)
	private double calculateGradient(double[] oldState, double[] newState)
	{
		return (double) evaluation_.evaluate(newState) - evaluation_.evaluate(oldState);
	}

	// NOTE: Dit is wel goed geloof ik.
	// Gradient moet een array zijn van gradients; gradient[i]
	// We moeten even nadenken over alpha. In principe is wat die if nu doet
	// een soort van indicator geven dat alpha te groot is.
	private double[] getNewState(double[] oldState, double gradient, double alpha)
	{
		double[] newState = oldState.clone(); // copy vector
		for (int i = 0; i < oldState.size(); i++)
		{
			double x = oldState[i]; // remember value

			// adjust value by gradient
			oldState[i] = oldState[i] + (alpha * gradient);

			// if the adjustment doesn't improve the gradient, do not adjust
			if (calculateGradient(oldState, newState) < 0) {
				oldState[i] = x;
			}
		}
		return newState;
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

    /*
     * @see CMAEvolutionStrategy
     *
     * @author Nikolaus Hansen, released into public domain.
     * got this from https://github.com/okanasik/cma-es
     */

    public class CMA_ES {
        public static void optimize(String[] args) {
            // IObjectiveFunction fitfun = new Rosenbrock();

            // new a CMA-ES and set some initial values
            CMAEvolutionStrategy cma = new CMAEvolutionStrategy();
            cma.readProperties(); // read options, see file CMAEvolutionStrategy.properties
            cma.setDimension(10); // overwrite some loaded properties
            cma.setInitialX(0.05); // in each dimension, also setTypicalX can be used
            cma.setInitialStandardDeviation(0.2); // also a mandatory setting
            cma.options.stopFitness = 1e-14;       // optional setting

            // initialize cma and get fitness array to fill in later
            double[] fitness = cma.init();  // new double[cma.parameters.getPopulationSize()];

            // initial output to files
            cma.writeToDefaultFilesHeaders(0); // 0 == overwrites old files

            // iteration loop
            while(cma.stopConditions.getNumber() == 0) {

                // --- core iteration step ---
                double[][] pop = cma.samplePopulation(); // get a new population of solutions
                for(int i = 0; i < pop.length; ++i) {    // for each candidate solution i
                    // a simple way to handle constraints that define a convex feasible domain
                    // (like box constraints, i.e. variable boundaries) via "blind re-sampling"
                                                           // assumes that the feasible domain is convex, the optimum is
                    while (!fitfun.isFeasible(pop[i]))     //   not located on (or very close to) the domain boundary,
                        pop[i] = cma.resampleSingle(i);    //   initialX is feasible and initialStandardDeviations are
                                                           //   sufficiently small to prevent quasi-infinite looping here
                    // compute fitness/objective value
                    // fitness[i] = fitfun.valueOf(pop[i]); // fitfun.valueOf() is to be minimized
                        fitness[i] = getFitnessPopulation(population[i], population.size());
                }
                cma.updateDistribution(fitness);         // pass fitness array to update search distribution
                // --- end core iteration step ---

                // output to files and console
                cma.writeToDefaultFiles();
                int outmod = 150;
                if (cma.getCountIter() % (15*outmod) == 1)
                    cma.printlnAnnotation(); // might write file as well
                if (cma.getCountIter() % outmod == 1)
                    cma.println();
            }
            // evaluate mean value as it is the best estimator for the optimum
            cma.setFitnessOfMeanX(fitfun.valueOf(cma.getMeanX())); // updates the best ever solution

            // final output
            cma.writeToDefaultFiles(1);
            cma.println();
            cma.println("Terminated due to");
            for (String s : cma.stopConditions.getMessages())
                cma.println("  " + s);
            cma.println("best function value " + cma.getBestFunctionValue()
                    + " at evaluation " + cma.getBestEvaluationNumber());

            // we might return cma.getBestSolution() or cma.getBestX()

        } // main
    } // class

}
