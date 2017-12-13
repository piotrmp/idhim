package io.github.piotrmp.idhim;

/**
 * Represents coordinates in n-dimensional space
 */
public class Coordinates implements Comparable<Coordinates> {

	/**
	 * Value of the coordinates
	 */
	private double[] coors;

	/**
	 * Creates new objects from given coordinates
	 * 
	 * @param coors
	 *            value of the coordinates
	 */
	public Coordinates(double[] coors) {
		if (coors == null)
			throw new IllegalArgumentException("Coordinates can't be empty.");
		else
			this.coors = coors;
	}

	/**
	 * Compares to other coordinates. Objects are equal if they have the same
	 * naumber of dimensions and all the coordinates are equal.
	 * 
	 * @param coordinates
	 *            to compare to
	 */
	public int compareTo(Coordinates o) {
		if (coors.length != o.coors.length)
			return coors.length - o.coors.length;
		for (int i = 0; i < coors.length; ++i)
			if (coors[i] < o.coors[i])
				return -1;
			else if (coors[i] > o.coors[i])
				return 1;
		return 0;
	}

	/**
	 * Return a String representation of coordinates: series of comma-separated
	 * numbers in square brackets
	 */
	public String toString() {
		StringBuilder sb = new StringBuilder("[");
		for (int i = 0; i < coors.length; ++i)
			sb.append((i == 0 ? "" : ",") + coors[i]);
		sb.append("]");
		return sb.toString();
	}

}
