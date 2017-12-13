package io.github.piotrmp.idhim;

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Represents a polynomial of k variables 
 */
public class Polynomial {
	/**
	 * Polynomial representation, mapping from variable product (e.g. x*y^2) to coefficient 
	 */
	public Map<Moment, Double> m;

	/**
	 * Create zero polynomial
	 */
	public Polynomial() {
		m = new HashMap<Moment, Double>();
	}

	/** Create zero polynomial of rank n
	 * @param n rank of created polynomial
	 */
	public Polynomial(int n) {
		this();
		m.put(new Moment(new int[n]), 1.0);
	}

	/** Initiate a polynomial with a constant 
	 * @param n rank of created polynomial
	 * @param val value of free constant
	 */
	public Polynomial(int n, double val) {
		this();
		m.put(new Moment(new int[n]), val);
	}

	/** Add another polynomial to this one: values at the same moments will be summed 
	 * @param other polynomial to be added
	 */
	public void add(Polynomial other) {
		for (Moment otherM : other.m.keySet())
			if (m.containsKey(otherM))
				m.put(otherM, m.get(otherM) + other.m.get(otherM));
			else
				m.put(otherM, other.m.get(otherM));
		return;
	}
	
	/** Multiply this polynomial by the other
	 * @param other polynomial to multiply with
	 */
	public void multiply(Polynomial other) {
		Map<Moment, Double> old = m;
		m = new HashMap<Moment, Double>();
		for (Moment m1 : old.keySet())
			for (Moment m2 : other.m.keySet()) {
				Moment newM = new Moment(m1);
				// add moments
				newM.add(m2);
				// multiply values
				double newV = old.get(m1) * other.m.get(m2);
				if (m.containsKey(newM))
					m.put(newM, m.get(newM) + newV);
				else
					m.put(newM, newV);
			}

	}

	/** Multiply with a constant value
	 * @param val
	 */
	public void multiply(double val) {
		for (Moment m1 : m.keySet())
			m.put(m1, m.get(m1) * val);

	}

	/**
	 * Create string representation
	 */
	public String toString() {
		StringBuilder sb = new StringBuilder();
		for (Moment m1 : m.keySet())
			sb.append("N"+m1 + "*" + m.get(m1) + "+");
		if (sb.length() > 0)
			sb.deleteCharAt(sb.length() - 1);
		return sb.toString();
	}

	
	/**
	 * Clean from moment with zero coefficient
	 */
	public void clean() {
		for (Moment moment : new HashSet<Moment>(m.keySet()))
			if (m.get(moment) == 0.0)
				m.remove(moment);
	}
}
