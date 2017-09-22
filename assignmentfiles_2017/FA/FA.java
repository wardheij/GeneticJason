/**
 * 
 */
package old;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.util.Random;

/**
 * @author ºúÎ¡Î¡
 *
 */
public class FA {
	//parameters of FA
	private int numLocations;
	private int numMaxSparks;
	private double numBoundA;
	private double numBoundB;
	private double numMaxAmplitude;
	private int numGaussianSparks;
	private double [] maxbound;
	private double [] minbound;
	private int dimension;
	private String infopath; 
	
	//objects of FA
	private Spark [] fireworks;
	private Spark [][] sparks;
	private Spark [] gaussiansparks;
	
	//information of FA
	private int numGenerations;
	private double optimumvalue;
	private int numFunctionEvaluations;
	
	//minimum value of double
	private double eps;
	
	//fitness function
	Function func;
	
	//constructor
	public FA(int n,int m,double a,double b,double am,int mg,double [] maxb,double [] minb,String info,Function fun) {
		numLocations = n;
		numMaxSparks = m;
		numBoundA = a;
		numBoundB = b;
		numMaxAmplitude = am;
		numGaussianSparks = mg;
		maxbound = maxb;
		minbound = minb;
		dimension = maxbound.length;
		infopath = info;
		
		eps = 1e-38;
		
		func = fun;
	}
	
	//fireworks algorithm framework
	public double FAframework() {
		//select n initial locations
		selectinitiallocations();
		while(stopcriteria() == false) {
			//set off n fireworks
			setoff();
			//select n locations
			selectlocations();
		}
		return optimumvalue;
	}
	//select n initial locations
	private void selectinitiallocations() {
		numGenerations = 0;
		numFunctionEvaluations = 0;
		fireworks = new Spark [numLocations];
		//random position
		double [] randpos = new double [dimension];
		//set random position
		int i,j;
		//set random positions for all fireworks
		for(i = 0;i < numLocations;i ++) {
			fireworks[i] = new Spark();
			//generate a random position
			for(j = 0;j < dimension;j ++) {
				randpos[j] = maxbound[j] - Math.random() * (maxbound[j] - minbound[j]) * 0.25;
			}
			fireworks[i].setposition(randpos);
		}
	}
	//set off n fireworks
	private void setoff() {
		numGenerations ++;
		//get max(worst) and min(best) value
		double maxvalue = fireworks[0].getvalue(func);
		double minvalue = fireworks[0].getvalue(func);
		int i;
		for(i = 1;i < numLocations;i ++) {
			if(fireworks[i].getvalue(func) > maxvalue) {
				maxvalue = fireworks[i].getvalue(func);
			}
			if(fireworks[i].getvalue(func) < minvalue) {
				minvalue = fireworks[i].getvalue(func);
			}
		}
		double summaxdiff = 0.0;
		double summindiff = 0.0;
		for(i = 0;i < numLocations;i ++) {
			summaxdiff += maxvalue - fireworks[i].getvalue(func);
			summindiff += fireworks[i].getvalue(func) - minvalue;
		}
		
		//get number of sparks for all fireworks
		int [] numSparks = new int [numLocations];
		double tmpcoef;
		for(i = 0;i < numLocations;i ++) {
			tmpcoef = (maxvalue - fireworks[i].getvalue(func) + eps) / (summaxdiff + eps);
			if(tmpcoef < numBoundA) {
				tmpcoef = numBoundA;
			}
			if(tmpcoef > numBoundB) {
				tmpcoef = numBoundB;
			}
			numSparks[i] = (int) (numMaxSparks * tmpcoef);
		}
		
		//get amplitude of explosion for all fireworks
		double ampExplosition [] = new double [numLocations];
		for(i = 0;i < numLocations;i++) {
			ampExplosition[i] =(fireworks[i].getvalue(func) - minvalue + eps) / (summindiff + eps) * numMaxAmplitude;
		}
		
		//generate sparks for all fireworks
		sparks = new Spark [numLocations] [];
		//temporary position
		double [] tmppos = new double [dimension];
		double [] fireworkpos;
		for(i = 0;i < numLocations;i++) {
			sparks[i] = new Spark [numSparks[i]];
			fireworkpos = fireworks[i].getposition();
			//get all sparks' position
			int k;
			for(k = 0 ;k < numSparks[i];k ++) {
				sparks[i][k] = new Spark();
				//select z directions
				Random rand=new Random();
				boolean [] randflag = new boolean [dimension];
				int j;
				for(j = 0;j < dimension;j ++) {
					randflag[j] = false;
				}
				int numExplosionDirections = (int) (dimension * Math.random());
				int randomcount = 0;
				int tmprand;
				while(randomcount < numExplosionDirections) {
					tmprand = rand.nextInt(dimension);
					if(randflag[tmprand] == false) {
						randflag[tmprand] = true;
						randomcount ++;
					}
				}
				//explode
				double displacement = ampExplosition[i] * (Math.random() - 0.5) * 2;
				for(j = 0 ;j < dimension;j ++) {
					if(randflag[j] == true) {
						tmppos[j] = fireworkpos[j] + displacement;
						//out of bound
						if(tmppos[j] < minbound[j] || tmppos[j] > maxbound[j]) {
							double abspos = Math.abs(tmppos[j]);
							while(abspos >= 0) {
								abspos -= (maxbound[j] - minbound[j]);
							}
							abspos += (maxbound[j] - minbound[j]);
							tmppos[j] = minbound[j] + abspos;
						}
					}
					else {
						tmppos[j] = fireworkpos[j];
					}
				}
				//set position of the spark
				sparks[i][k].setposition(tmppos);
			}
		}
		
		//gaussian explode
		gaussiansparks = new Spark [numGaussianSparks];
		int k;
		for(k = 0;k < numGaussianSparks;k ++) {
			gaussiansparks[k] = new Spark();
			//randomly select a firework
			Random rand=new Random();
			i = Math.abs(rand.nextInt()) % numLocations;
			fireworkpos = fireworks[i].getposition();
			//select z directions
			boolean [] randflag = new boolean [dimension];
			int j;
			for(j = 0;j < dimension;j ++) {
				randflag[j] = false;
			}
			int numExplosionDirections = (int) (dimension * Math.random());
			int randomcount = 0;
			int tmprand;
			while(randomcount < numExplosionDirections) {
				tmprand = Math.abs(rand.nextInt()) % dimension;
				if(randflag[tmprand] == false) {
					randflag[tmprand] = true;
					randomcount ++;
				}
			}
			//explode
			double gaussiancoef = 1.0 + rand.nextGaussian();
			for(j = 0 ;j < dimension;j ++) {
				if(randflag[j] == true) {
					tmppos[j] = fireworkpos[j] * gaussiancoef;
					//out of bound
					if(tmppos[j] < minbound[j] || tmppos[j] > maxbound[j]) {
						double abspos = Math.abs(tmppos[j]);
						while(abspos >= 0) {
							abspos -= (maxbound[j] - minbound[j]);
						}
						abspos += (maxbound[j] - minbound[j]);
						tmppos[j] = minbound[j] + abspos;
					}
				}
				else {
					tmppos[j] = fireworkpos[j];
				}
			}
			//set position of the spark
			gaussiansparks[k].setposition(tmppos);
		}
	}
	//select n locations
	private void selectlocations() {
		//select the best location
		Spark bestspark = fireworks[0];
		int i,j,k;
		for(i = 1;i < numLocations;i ++) {
			if(fireworks[i].getvalue(func) < bestspark.getvalue(func)) {
				bestspark = fireworks[i];
			}
		}
		for(i = 0;i < numLocations;i ++) {
			for(j = 0;j < sparks[i].length;j ++) {
				if(sparks[i][j].getvalue(func) < bestspark.getvalue(func)) {
					bestspark = sparks[i][j];
				}
			}
		}
		for(i = 0;i < numGaussianSparks;i ++) {
			if(gaussiansparks[i].getvalue(func) < bestspark.getvalue(func)) {
				bestspark = gaussiansparks[i];
			}
		}
		optimumvalue = bestspark.getvalue(func);
		//output the best value
		PrintStream info = null;
		try {
			info = new PrintStream(new FileOutputStream(infopath,true));
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		info.println("best value = " + (1 - optimumvalue));
		info.println("best position = " + bestspark.getposition().toString());
		info.println("---------------------------------------------------------");
		info.close();
		//select the rest n-1 locations
		//count the number of fireworks and sparks
		int numFireworksSparks = numLocations + numGaussianSparks;
		for(i = 0;i < numLocations;i ++) {
			for(j = 0;j < sparks[i].length;j ++) {
				numFireworksSparks ++;
			}
		}
		//calculate the number of function evaluations
		numFunctionEvaluations += numFireworksSparks;
		//put all the fireworks and sparks in an array
		double [][] fireworkspos = new double [numFireworksSparks][];
		int idx = 0;
		for(i = 0;i < numLocations;i ++) {
			fireworkspos[idx] = fireworks[i].getposition();
			idx ++;
		}
		for(i = 0;i < numLocations;i ++) {
			for(j = 0;j < sparks[i].length;j ++) {
				fireworkspos[idx] = sparks[i][j].getposition();
				idx ++;
			}
		}
		for(i = 0;i < numGaussianSparks;i ++) {
			fireworkspos[idx] = gaussiansparks[i].getposition();
			idx ++;
		}
		//calculate the selection probability of each location
		double [] selectionprobability = new double [numFireworksSparks];
		double sumprob = 0;
		for(i = 0;i < numFireworksSparks;i ++) {
			selectionprobability[i] = 0;
			for(j = 0;j < numFireworksSparks;j ++) {
				double tmpdis = 0;
				for(k = 0;k < dimension;k ++) {
					tmpdis += (fireworkspos[i][k] - fireworkspos[j][k]) * (fireworkspos[i][k] - fireworkspos[j][k]);
				}
				selectionprobability[i] += Math.sqrt(tmpdis);
			}
			sumprob += selectionprobability[i];
		}
		double [] cumulativeprobability = new double [numFireworksSparks];
		for(i = 0;i < numFireworksSparks;i ++) {
			if(sumprob < eps) {
				selectionprobability[i] = 1.0 / numFireworksSparks;
			}
			else {
				selectionprobability[i] /= sumprob;
			}
			if(i == 0) {
				cumulativeprobability[i] =selectionprobability[i];
			}
			else {
				cumulativeprobability[i] = cumulativeprobability[i - 1] + selectionprobability[i];
			}
		}
		//select n-1 locations according to the selection probability
		int [] nextlocations = new int [numLocations - 1];
		for(k = 0;k < numLocations - 1;k ++) {
			double randpointer = Math.random();
			for(i = 0;i < numFireworksSparks;i ++) {
				if(randpointer <= cumulativeprobability[i]) {
					break;
				}
			}
			nextlocations[k] = i;
		}
		//set next generations
		Spark[] nextfireworks = new Spark [numLocations];
		nextfireworks[numLocations - 1] = bestspark;
		boolean breakflag;
		for(k = 0;k < numLocations - 1;k ++) {
			idx = 0;
			breakflag = false;
			for(i = 0;i < numLocations;i ++) {
				if(idx == nextlocations[k]) {
					nextfireworks[k] = fireworks[i];
					breakflag = true;
					break;
				}
				idx ++;
			}
			if(breakflag == true) {
				continue;
			}
			for(i = 0;i < numLocations;i ++) {
				for(j = 0;j < sparks[i].length;j ++) {
					if(idx == nextlocations[k]) {
						nextfireworks[k] = sparks[i][j];
						breakflag = true;
						break;
					}
					idx ++;
				}
				if(breakflag == true) {
					break;
				}
			}
			if(breakflag == true) {
				continue;
			}
			for(i = 0;i < numGaussianSparks;i ++) {
				if(idx == nextlocations[k]) {
					nextfireworks[k] = gaussiansparks[i];
					breakflag = true;
					break;
				}
				idx ++;
			}
		}
		fireworks = nextfireworks;
	}
	private boolean stopcriteria() {
		//if(numGenerations < 2000) {
		if(numFunctionEvaluations < 300000) {
			return false;
		}
		else {
			//System.out.println("numGenerations=" + numGenerations + ",numFunctionEvaluations=" + numFunctionEvaluations);
			return true;
		}
	}
}
