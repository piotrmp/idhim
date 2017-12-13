package io.github.piotrmp.idhim;

import java.util.List;

import org.ejml.simple.SimpleMatrix;

/**
 * Stores probability density of mass (or other quantity) values for a single particle type
 */
public class MassDensities {
	/**
	 * List of coordinates at which observations were made
	 */
	public List<Coordinates> x;
	/**
	 * Number of observations
	 */
	public SimpleMatrix y;

	/**
	 * Constructor converting double[][] to SimpleMatrix
	 * @param densX
	 * @param densY
	 */
	public MassDensities(List<Coordinates> densX, double[][] densY) {
		x = densX;
		y = (new SimpleMatrix(densY)).transpose();
	}

	/**
	 * Simple constructor
	 * @param x
	 * @param y
	 */
	public MassDensities(List<Coordinates> x, SimpleMatrix y) {
		this.x = x;
		this.y = y;
	}

}
