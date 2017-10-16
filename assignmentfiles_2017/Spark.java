/**
 *
 */
package EC;

import org.vu.contest.ContestSubmission;
import org.vu.contest.ContestEvaluation;
/**
 * @author ��ΡΡ
 *
 */

public class Spark {
	//position of the spark
	private double [] position;
	//if this spark has been evaluated
	private boolean bEvaluated;
	//the function value of this spark
	private double value;

	//constructor
	public Spark() {
		bEvaluated = false;
	}

	//set the position of the spark
	public void setposition(double [] pos) {
		position = new double [pos.length];
		int i;
		for (i = 0;i < pos.length ;i++) {
			position [i] = pos [i];
		}
		bEvaluated = false;
	}

	//get the position of the spark
	public double[] getposition() {
		return position;
	}

	//get the function value of the spark
	public double getvalue(ContestEvaluation func) {
		if(bEvaluated == false) {
			bEvaluated = true;
			value = (double) func.evaluate(position);
			return value;
		}
		else {
			return value;
		}
	}

}
