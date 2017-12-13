package io.github.piotrmp.idhim;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedList;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.SortedMap;
import java.util.TreeMap;

/**
 * Represents a moment of a series of random variables, e.g. [0,2,1,0]=>
 * <(x_0^0)*(x_1^2)*(x_2^1)*(x_3^0)> = <x_1*x_1*x_2> The sum of the values in
 * array is the moment rank, e.g. 0+2+1+0=3
 */
public class Moment implements Comparable<Moment> {
	/**
	 * Rank of each of variables in the moment, i.e. a power it is assigned
	 */
	public int[] n;

	/**
	 * Creates a new moment from a given rank array
	 * 
	 * @param n
	 *            array of ranks (powers in the moment)
	 */
	public Moment(int[] n) {
		this.n = n;
	}

	/**
	 * Creates empty moment, i.e. with all powers set to zero
	 * 
	 * @param k
	 *            length of the moment, i.e. number of random variables
	 */
	public Moment(int k) {
		this(new int[k]);
	}

	/**
	 * Copy from another moment
	 * 
	 * @param other
	 *            moment to copy from
	 */
	public Moment(Moment other) {
		this.n = Arrays.copyOf(other.n, other.n.length);
	}

	/**
	 * Checks if a moment equals another one. Two moments are equal if the have
	 * the same length and values in the arrays.
	 */
	@Override
	public boolean equals(Object obj) {
		if (this == obj)
			return true;
		if (obj == null)
			return false;
		if (getClass() != obj.getClass())
			return false;
		Moment other = (Moment) obj;
		return Arrays.equals(n, other.n);
	}

	/**
	 * Moment comparison is based on comparison of their rank arrays
	 */
	@Override
	public int compareTo(Moment arg0) {
		if (arg0 == null)
			throw new NullPointerException();
		if (arg0.n.length != n.length)
			return n.length - arg0.n.length;
		for (int i = 0; i < n.length; ++i)
			if (n[i] != arg0.n[i])
				return n[i] - arg0.n[i];
		return 0;
	}

	/**
	 * Hash code computation is based on the included rank array
	 */
	@Override
	public int hashCode() {
		return Arrays.hashCode(n);
	}

	/**
	 * A moment is represented in String as series of numbers in angular
	 * brackets, e.g. '<0210>'
	 */
	public String toString() {
		StringBuilder sb = new StringBuilder();
		if (n.length == 0)
			sb.append("<empty>");
		else {
			sb.append("<");
			for (int i = 0; i < n.length; ++i)
				sb.append(n[i]);
		}
		sb.append(">");
		return sb.toString();
	}

	/**
	 * An unambiguous conversion to string adding underscore between values
	 * 
	 * @return String representation of the moment
	 */
	public String toStringU() {
		StringBuilder sb = new StringBuilder();
		for (int i = 0; i < n.length; ++i) {
			if (sb.length() > 0)
				sb.append('_');
			sb.append(n[i]);
		}
		return sb.toString();
	}

	/**
	 * Generates every possible permutation of the non-zero elements of the
	 * moment
	 * 
	 * @return list of permutations, each represented as array of non-zero
	 *         integers
	 */
	public List<int[]> nonZeroPermutations() {
		List<Integer> nonzeros = new ArrayList<>();
		for (int i = 0; i < n.length; ++i)
			if (n[i] != 0)
				nonzeros.add(i);
		List<List<Integer>> permutations = permute(nonzeros);
		List<int[]> result = new ArrayList<int[]>(permutations.size());
		for (List<Integer> permutation : permutations) {
			int[] x = new int[permutation.size()];
			for (int i = 0; i < permutation.size(); ++i)
				x[i] = permutation.get(i);
			result.add(x);
		}
		return result;
	}

	/**
	 * Generates every possible moment of given rank and length
	 * 
	 * @param n
	 *            rank of the generated moments (sum of powers)
	 * @param k
	 *            length of the generated moments (number of variables involved)
	 * @return list of moments generated for given criteria
	 */
	public static List<Moment> all(int n, int k) {
		Set<Moment> result = new HashSet<>();
		if (n < 0)
			throw new IllegalArgumentException("Can't generate moments for rank less than 0.");
		else if (n == 0)
			// only one possible moment of rank 0
			result.add(new Moment(new int[k]));
		else {
			// generate moments of lower rank
			List<Moment> lower = all(n - 1, k);
			for (Moment lowerM : lower)
				// for each of them increment one of the values
				for (int i = 0; i < k; ++i) {
					Moment copy = new Moment(lowerM);
					copy.n[i]++;
					result.add(copy);
				}
		}
		return new ArrayList<Moment>(result);
	}

	/**
	 * Returns every possible permutation of a given list of items (including
	 * the original)
	 * 
	 * @param items
	 *            list of items to permute
	 * @return list of permuted lists of items
	 */
	private static List<List<Integer>> permute(List<Integer> items) {
		List<List<Integer>> result = new ArrayList<List<Integer>>();
		if (items.size() == 1)
			// only way to permute single item
			result.add(items);
		else
			// for every possible division of the list..
			for (int i = 0; i < items.size(); ++i) {
				List<Integer> copy = new ArrayList<Integer>(items);
				Integer last = copy.remove(i);
				// permute the first part
				for (List<Integer> shorter : permute(copy)) {
					shorter.add(last);
					// leave the rest unchanged
					result.add(shorter);
				}
			}
		return result;
	}

	/**
	 * Generate all possible sequences of k (length of n[]) moments, such that
	 * they all sum to this moment
	 * 
	 * @return list of possible distribution of n[] over k moments
	 */
	public List<Moment[]> getDistributions() {
		@SuppressWarnings("unchecked")
		// generate every possible distribution of n[i] over a single moment
		List<Moment>[] columnCombinations = new List[n.length];
		int[] columnCombinationCounts = new int[n.length];
		for (int i = 0; i < n.length; ++i) {
			columnCombinations[i] = all(n[i], n.length);
			columnCombinationCounts[i] = columnCombinations[i].size();
		}
		// for every possible choice of moments in each column..
		List<int[]> choices = generateChoices(columnCombinationCounts, 0);
		List<Moment[]> result = new ArrayList<Moment[]>();
		for (int[] choice : choices) {
			// generate the result
			Moment[] thischoice = new Moment[n.length];
			for (int i = 0; i < choice.length; ++i)
				thischoice[i] = columnCombinations[i].get(choice[i]);
			result.add(thischoice);
		}
		return (result);
	}

	/**
	 * Generate all possible choices where at each step you can choose from a
	 * given number of options (starting from a given index)
	 * 
	 * @param options
	 *            array of numbers denoting number of options at every step
	 * @param i
	 *            index of the first step that should be taken into account
	 *            while generating choices
	 * @return list of choices as arrays of indices for the chosen options at
	 *         each step, e.g.
	 *         [1,3,2]=>[0,0,0],[0,0,1],[0,1,0],[0,1,1],[0,2,0],[0,2,1]
	 */
	private static List<int[]> generateChoices(int[] options, int i) {
		List<int[]> results = new ArrayList<int[]>();
		if (i == options.length) {
			// steps left to choose from, generate a single default option
			results.add(new int[0]);
			return results;
		}
		// generate choices for the next steps
		List<int[]> furtherResults = generateChoices(options, i + 1);
		// for every possible choice at this step..
		for (int j = 0; j < options[i]; ++j)
			// for every possible choice at further steps..
			for (int[] further : furtherResults) {
				// combine the two
				int[] thisResult = new int[options.length - i];
				thisResult[0] = j;
				for (int k = 0; k < further.length; ++k)
					thisResult[k + 1] = further[k];
				results.add(thisResult);
			}
		return results;
	}

	/**
	 * Add a value of another moment to this one
	 * 
	 * @param other
	 *            the other moment
	 */
	public void add(Moment other) {
		for (int i = 0; i < n.length; ++i)
			n[i] += other.n[i];
	}

	/**
	 * Get all possible distributions of this moment's n[] into an unlimited
	 * number of other (possible repeating) moments, e.g.
	 * [210]=>{2x100+010,200+010,110+100,210}
	 * 
	 * @return a list of possible distributions, each expressed as map from a
	 *         moment constituent to how many times it occurs for this
	 *         combination
	 */
	public List<Map<Moment, Integer>> getDistributionsInf() {
		// list of possible configurations
		List<Map<Moment, Integer>> configurations = new ArrayList<Map<Moment, Integer>>();
		// initial configuration: one empty model
		Map<Moment, Integer> init = new HashMap<Moment, Integer>();
		init.put(new Moment(0), 1);
		configurations.add(init);
		int k = n.length;
		// build configurations going through moment elements
		for (int i = 0; i < k; ++i) {
			List<Map<Moment, Integer>> newConfigurations = new ArrayList<Map<Moment, Integer>>();
			// for each configuration up to now...
			for (Map<Moment, Integer> conf : configurations) {
				// moments present in this configuration
				Moment[] moments = conf.keySet().toArray(new Moment[conf.size()]);
				// generate variations representing possible distribution of
				// current element between moments of this configuration
				List<int[]> variations = getVariations(n[i], conf.size());
				for (int[] variation : variations) {
					@SuppressWarnings("unchecked")
					List<int[][]>[] combinations = new List[variation.length];
					// for each moment set
					for (int j = 0; j < moments.length; ++j)
						// generate possible distributions of the value added
						// (variation[j])
						// there into these several (conf.get(moments[j]))
						// identical moments
						combinations[j] = getCombinations(variation[j],
								moments[j].sum() > 0 ? conf.get(moments[j]) : variation[j]);
					// number of combinations for each moment set is a number of possible choices there
					int[] possibleChoices = new int[combinations.length];
					for (int j = 0; j < combinations.length; ++j)
						possibleChoices[j] = combinations[j].size();
					// for each choice from a given sequance
					for (int[] choice : generateChoices(possibleChoices, 0)) {
						// generate new configuration
						SortedMap<Moment, Integer> newconf = new TreeMap<Moment, Integer>();
						for (int j = 0; j < moments.length; ++j) {
							int[][] combination = combinations[j].get(choice[j]);
							for (int l = 0; l < combination.length; ++l) {
								int[] ma = null;
								ma = Arrays.copyOf(moments[j].n, moments[j].n.length + 1);
								ma[ma.length - 1] = combination[l][0];
								Moment m = new Moment(ma);
								newconf.put(m, combination[l][1]);
							}
						}
						// add new empty moment so that there always is one
						newconf.put(new Moment(new int[i + 1]), 1);
						newConfigurations.add(newconf);
					}
				}
			}
			configurations = newConfigurations;
		}
		List<Map<Moment, Integer>> result = new ArrayList<Map<Moment, Integer>>();
		for (Map<Moment, Integer> conf : configurations) {
			// remove empty moments
			conf.remove(new Moment(k));
			// if configuration nonempty, add it to results
			if (!conf.isEmpty())
				result.add(conf);
		}
		return result;
	}

	/**
	 * Generate all possible distributions of n elements into k boxes
	 * 
	 * @param n
	 *            number of items to be distributed
	 * @param k
	 *            number of boxes
	 * @return a list of possible distributions: vectors of length k, each
	 *         summing up to n
	 */
	private static List<int[]> getVariations(int n, int k) {
		List<int[]> results = new LinkedList<int[]>();
		if (n == 0)
			// if no elements only one possibility
			results.add(new int[k]);
		else if (k == 1) {
			// if one box only one possibility
			int[] res = { n };
			results.add(res);
		} else
			// for every possible number of items in last box...
			for (int i = 0; i <= n; ++i) {
				// generate variation using the remaining boxes and items
				for (int[] lower : getVariations(n - i, k - 1)) {
					int[] upper = Arrays.copyOf(lower, lower.length + 1);
					// and put i items in the last box
					upper[k - 1] = i;
					results.add(upper);
				}
			}
		return results;
	}

	/**
	 * Generate all possible distributions of n elements into k
	 * indistinguishable boxes
	 * 
	 * @param n
	 *            number of elements to distribute
	 * @param k
	 *            number of boxes
	 * @return a list of possible distribution, each as two-dimensional array,
	 *         where each row has two columns, containing number of items in box
	 *         and number of boxes
	 */
	private static List<int[][]> getCombinations(int n, int k) {
		// naive and inefficient approach
		List<int[][]> results = new LinkedList<int[][]>();
		if (n == 0) {
			// if no elements to distribute we just have k empty boxes
			int[][] result = new int[1][];
			int[] zerozero = { 0, k };
			result[0] = zerozero;
			results.add(result);
		} else {
			List<int[]> preresults = new LinkedList<int[]>();
			// get all variations for the same problem
			List<int[]> variations = getVariations(n, k);
			ext: for (int[] variation : variations) {
				int[] sorted = Arrays.copyOf(variation, variation.length);
				// sort each variation
				Arrays.sort(sorted);
				for (int[] result : preresults)
					if (Arrays.equals(sorted, result))
						// if this combination already has been added from
						// different variation, ignore it
						continue ext;
				preresults.add(sorted);
			}
			// convert the combinations to the desired format
			for (int[] preresult : preresults) {
				Map<Integer, Integer> map = new HashMap<Integer, Integer>();
				for (int i : preresult)
					if (map.containsKey(i))
						map.put(i, map.get(i) + 1);
					else
						map.put(i, 1);
				int[][] array = new int[map.size()][];
				int counter = 0;
				for (Integer i : map.keySet()) {
					array[counter] = new int[2];
					array[counter][0] = i;
					array[counter][1] = map.get(i);
					counter++;
				}
				results.add(array);
			}
		}
		return results;
	}

	/**
	 * Sum of moment elements - its rank
	 * 
	 * @return moment sum
	 */
	public int sum() {
		int r = 0;
		for (int i : n)
			r += i;
		return r;
	}

}
