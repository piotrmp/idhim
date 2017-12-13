package io.github.piotrmp.idhim;

import java.util.HashMap;
import java.util.Map;

import org.ejml.simple.SimpleMatrix;

/**
 * Stores observed identity values in one of two ways
 */
public class ObservedIdentities {
	/**
	 * As a list of distributions for each moment
	 */
	public Map<Moment, SimpleMatrix[]> W;
	/**
	 * As a mean value for every moment
	 */
	public Map<Moment, Double> mean;

	/**
	 * Create empty object
	 */
	public ObservedIdentities() {
		W = new HashMap<Moment, SimpleMatrix[]>();
		mean = new HashMap<Moment, Double>();
	}

	/**
	 * Add identity variable distribution for a given moment
	 * 
	 * @param moment
	 *            moment associated with the added W
	 * @param x
	 *            identity variable values
	 * @param y
	 *            number of occurrences of respective values
	 */
	public void add(Moment moment, double[] x, double[] y) {
		double[][] tmp = new double[2][];
		tmp[0] = x;
		tmp[1] = y;
		double[][] x1 = { x };
		double[][] y1 = { y };
		SimpleMatrix xx = (new SimpleMatrix(x1)).transpose();
		SimpleMatrix yy = (new SimpleMatrix(y1)).transpose();
		SimpleMatrix[] xxyy = { xx, yy };
		W.put(moment, xxyy);
		// compute means
		mean.put(moment, xx.elementMult(yy).elementSum() / yy.elementSum());
	}

	/**
	 * Add identity variable values by providing a mean
	 * 
	 * @param moment
	 *            moment associated with the added W
	 * @param mean
	 *            mean value of this identity variable
	 */
	public void add(Moment moment, Double mean) {
		this.mean.put(moment, mean);
	}
}
