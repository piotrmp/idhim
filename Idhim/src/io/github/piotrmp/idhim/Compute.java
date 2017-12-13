package io.github.piotrmp.idhim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.ejml.factory.SingularMatrixException;
import org.ejml.simple.SimpleMatrix;

/**
 * Reads the data and computes higher-order moments (see paper for explanation
 * of calculation)
 */
public class Compute {
	/**
	 * List of particle types
	 */
	private static List<String> types = new ArrayList<String>();
	/**
	 * List of measurement settings
	 */
	private static List<String> bins = new ArrayList<String>();
	/**
	 * Mass (or other quantity) densities for each particle type
	 */
	private static Map<String, MassDensities> rho = new HashMap<String, MassDensities>();
	/**
	 * Identity variable values for particle types at mass points
	 */
	private static Map<String, MassDensities> w = new HashMap<String, MassDensities>();
	/**
	 * User-provided observed W values
	 */
	private static ObservedIdentities W = new ObservedIdentities();
	/**
	 * Unnormalised multiplicity values
	 */
	private static Map<Moment, Double> UN = new HashMap<Moment, Double>();
	/**
	 * Event multiplicity values
	 */
	private static Map<Moment, Double> N = new HashMap<Moment, Double>();
	/**
	 * Number of particle types
	 */
	private static int k;
	/**
	 * Maximum rank for which moments will be computed
	 */
	private static int nMax;
	/**
	 * Rank of currently available moments
	 */
	private static int n;
	/**
	 * Computed new multiplicity moments
	 */
	private static List<Moment> newN;
	/**
	 * Formulas linking moments of W with moments of N
	 */
	private static Map<Moment, Polynomial> formulas = new HashMap<Moment, Polynomial>();
	/**
	 * Cache structure to avoid re-calculating u() for the same moment
	 */
	private static Map<Moment, double[]> uCache = new HashMap<Moment, double[]>();
	/**
	 * b vector for linear problem
	 */
	private static SimpleMatrix b;
	/**
	 * A matrix for linear problem
	 */
	private static SimpleMatrix A;
	/**
	 * Correction for malformed input files
	 */
	private static int correction = 1;
	/**
	 * Provide more information to standard output
	 */
	private static boolean verbose = false;
	/**
	 * Whether mass distributions are to be treated as normalised
	 */
	private static boolean normalisedMasses = false;
	/**
	 * File with particle types
	 */
	private static File typesF;
	/**
	 * File with bin names
	 */
	private static File binsF;
	/**
	 * Directory with mass densities
	 */
	private static File massDir;
	/**
	 * File with means of W
	 */
	private static File meansF;
	/**
	 * File for output
	 */
	private static File outF;

	/**
	 * Main function
	 * 
	 * @param args
	 *            command-line arguments to parse
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {
		boolean go = prepareCLI(args);
		if (go)
			run();
	}

	/**
	 * Parse command line arguments
	 * 
	 * @param args
	 *            arguments
	 * @return if a parse is successful and the program should proceed
	 */
	private static boolean prepareCLI(String[] args) {
		// define command-line options
		Options options = new Options();
		Option help = new Option("h", "help", false, "print this message");
		Option verboseO = new Option("v", "verbose", false, "be extra verbose");
		Option typesFileO = new Option("t", "types", true,
				"file with newline-separated <k> particle types (defaults to './types.tsv')");
		typesFileO.setArgName("file");
		Option binsFileO = new Option("b", "bins", true,
				"file with newline-separated bin names (possibly as multiple tab-separated columns, defaults to './bins.tsv')");
		binsFileO.setArgName("file");
		Option rhoDirO = new Option("r", "rhos", true,
				"directory containing files with mass densities named 'rho_<type>_<bin>.tsv' (defaults to './rho')");
		rhoDirO.setArgName("dir");
		Option meansFileO = new Option("W", "meanW", true,
				"file with mean sums of identity variables with moment provided as tab-separated <k> columns, followed by a single moment value (defaults to './meanW.tsv')");
		meansFileO.setArgName("file");
		Option outFileO = new Option("o", "out", true,
				"path to TSV file to which the computed moment values will be written (none if undefined)");
		outFileO.setArgName("file");
		options.addOption(help);
		options.addOption(verboseO);
		options.addOption(typesFileO);
		options.addOption(binsFileO);
		options.addOption(rhoDirO);
		options.addOption(meansFileO);
		options.addOption(outFileO);
		// parse input
		CommandLineParser parser = new DefaultParser();
		CommandLine line = null;
		try {
			line = parser.parse(options, args);
		} catch (ParseException e) {
			System.out.println("Unable to parse the command-line arguments");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(" ", options);
			return false;
		}
		if (line.hasOption("h") || (!line.hasOption("v") && !line.hasOption("o"))) {
			// doesn't make sense to run if no output will be displayed
			if (!line.hasOption("h"))
				System.out.println("No output selected: choose either -v or -o to see the results.");
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(" ", options);
			return false;
		}
		verbose = line.hasOption("v");
		// default files: to be overridden by user-specified
		String typesFilePath = "./types.tsv";
		String binsFilePath = "./bins.tsv";
		String rhoDirPath = "./rho";
		String meansFilePath = "./meanW.tsv";
		String outFilePath = null;
		if (line.hasOption("t"))
			typesFilePath = line.getOptionValue("t");
		if (line.hasOption("b"))
			binsFilePath = line.getOptionValue("b");
		if (line.hasOption("r"))
			rhoDirPath = line.getOptionValue("r");
		if (line.hasOption("W"))
			meansFilePath = line.getOptionValue("W");
		if (line.hasOption("o"))
			outFilePath = line.getOptionValue("o");
		typesF = new File(typesFilePath);
		binsF = new File(binsFilePath);
		massDir = new File(rhoDirPath);
		meansF = new File(meansFilePath);
		if (outFilePath != null)
			outF = new File(outFilePath);
		return true;
	}

	/**
	 * Run the computations
	 * 
	 * @throws IOException
	 */
	public static void run() throws IOException {
		int nStart=1;
		System.out.println("Reading basic data...");
		readBasic();
		System.out.println("Reading mass densities...");
		readMasses();
		System.out.println("Computing expected identity values...");
		computeExpected();
		System.out.println("Reading observed identity values...");
		readObservationsMeans();
		if (normalisedMasses){
			N=UN;
			nStart=2;
		}
		BufferedWriter bw = null;
		if (outF != null)
			bw = new BufferedWriter(new FileWriter(outF));
		// add higher moments iteratively to reduce errors
		for (n = nStart; n <= nMax; ++n) {
			System.out.println("n=" + n + ", Generating moments ...");
			newN = Moment.all(n, k);
			generateFormulas();
			System.out.println("n=" + n + ", Preparing the linear problem...");
			prepareMatrices();
			System.out.println("n=" + n + ", Solving the linear problem...");
			solve(bw);
		}
		if (bw != null)
			bw.close();
	}

	/**
	 * Read basic data from files
	 * 
	 * @throws IOException
	 */
	private static void readBasic() throws IOException {
		if (!typesF.exists() || typesF.isDirectory())
			throw new RuntimeException("Unable to locate file: " + typesF.getName());
		try {
			BufferedReader reader = new BufferedReader(new FileReader(typesF));
			String line = null;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				if (!line.isEmpty())
					types.add(line);
			}
			reader.close();
		} catch (FileNotFoundException e) {
			throw new IOException("Unable to read file: " + typesF.getName(), e);
		}
		if (types.size() == 0)
			throw new RuntimeException("Unable to find types in the file.");
		k = types.size();
		if (verbose) {
			System.out.print("Found types:");
			for (String type : types)
				System.out.print(" " + type);
			System.out.println();
		}
		if (!binsF.exists() || binsF.isDirectory())
			throw new RuntimeException("Unable to locate file: " + binsF.getName());
		try {
			BufferedReader reader = new BufferedReader(new FileReader(binsF));
			String line = null;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				String[] parts = line.split("\\s+");
				StringBuilder sb = new StringBuilder();
				for (String part : parts) {
					if (sb.length() != 0)
						sb.append("_");
					sb.append(part);
				}
				bins.add(sb.toString());
			}
			reader.close();
		} catch (FileNotFoundException e) {
			throw new IOException("Unable to read file: " + binsF.getName(), e);
		}
		if (bins.size() == 0)
			throw new RuntimeException("Unable to find bins in the file.");
		if (verbose) {
			System.out.print("Found bins:");
			for (String bin : bins)
				System.out.print(" " + bin);
			System.out.println();
		}
	}

	/**
	 * Read mass densities from files
	 * 
	 * @throws IOException
	 */
	private static void readMasses() throws IOException {
		if (!massDir.exists() || !massDir.isDirectory())
			throw new RuntimeException("Unable to locate directory with mass densities: " + massDir.getName());
		for (String bin : bins) {
			List<Coordinates> densX = null;
			double[][] densY = null;
			int j = 0;
			for (String type : types) {
				String filename = "rho_" + type + "_" + bin + ".tsv";
				File file = new File(massDir, filename);
				if (!file.exists() || file.isDirectory())
					throw new RuntimeException("Unable to locate file with mass densities: " + filename);
				if (densX == null) {
					// First read X coordinates
					densX = new ArrayList<Coordinates>();
					try {
						BufferedReader reader = new BufferedReader(new FileReader(file));
						String line = null;
						while ((line = reader.readLine()) != null) {
							line = line.trim();
							String[] parts = line.split("\\s+");
							if (parts.length < 2)
								throw new RuntimeException("Unable to parse mass density file: " + file);
							// coordinates can be multi-dimensional (not
							// necessarily mass only)
							double[] xx = new double[parts.length - 1];
							for (int i = 0; i < xx.length; ++i)
								xx[i] = Double.parseDouble(parts[i]);
							densX.add(new Coordinates(xx));
						}
						reader.close();
					} catch (IOException e) {
						throw new IOException("Unable to read file: " + filename, e);
					}
					densY = new double[k][];
					for (int i = 0; i < k; ++i)
						densY[i] = new double[densX.size()];
				}
				try {
					// read Y values
					BufferedReader reader = new BufferedReader(new FileReader(file));
					String line = null;
					int i = 0;
					while ((line = reader.readLine()) != null) {
						line = line.trim();
						String[] parts = line.split("\\s+");
						if (parts.length < 2)
							throw new RuntimeException("Unable to parse mass density file: " + file);
						double[] xx = new double[parts.length - 1];
						for (int k = 0; k < xx.length; ++k)
							xx[k] = Double.parseDouble(parts[k]);
						Coordinates x = new Coordinates(xx);
						double y = Double.parseDouble(parts[parts.length - 1]) / correction;
						// if coordinates are different this time, we can't
						// proceed: files for different particle types have to
						// be aligned to compute w
						if (x.compareTo(densX.get(i)) != 0)
							throw new RuntimeException("Unexpected x coordinates in file " + file + ", expected "
									+ densX.get(i) + ", got " + x);
						densY[j][i] = y;
						++i;
					}
					reader.close();
				} catch (IOException e) {
					throw new IOException("Unable to read file: " + filename, e);
				}
				++j;
			}
			MassDensities densities = new MassDensities(densX, densY);
			rho.put(bin, densities);
		}
		// compute existing multiplicity moments by integrating (summing) the
		// rhos
		for (int i = 0; i < k; ++i) {
			double sum = 0;
			for (String bin : bins)
				sum += rho.get(bin).y.extractVector(false, i).elementSum();
			int[] mom = new int[k];
			mom[i] = 1;
			Moment moment = new Moment(mom);
			UN.put(moment, sum);
			if (verbose)
				System.out.println("A" + moment + "=" + sum);
		}
	}

	/**
	 * Read means of observed W
	 * 
	 * @throws IOException
	 */
	private static void readObservationsMeans() throws IOException {
		if (!meansF.exists() || meansF.isDirectory())
			throw new RuntimeException("Unable to locate file: " + meansF.getName());
		Map<String, Double> means = new HashMap<String, Double>();
		BufferedReader reader = null;
		// read all means
		try {
			reader = new BufferedReader(new FileReader(meansF));
			String line = null;
			while ((line = reader.readLine()) != null) {
				line = line.trim();
				String[] parts = line.split("\\s+");
				if (parts.length < 2) {
					reader.close();
					throw new RuntimeException("Unexpected format of file: " + meansF.getName());
				}
				StringBuilder sb = new StringBuilder();
				for (int i = 0; i < parts.length - 1; ++i) {
					if (sb.length() != 0)
						sb.append("_");
					sb.append(parts[i]);
				}
				Double value = Double.parseDouble(parts[parts.length - 1]);
				means.put(sb.toString(), value);
			}
			reader.close();
		} catch (FileNotFoundException e) {
			reader.close();
			throw new IOException("Unable to read file: " + meansF.getName(), e);
		} catch (NumberFormatException e) {
			reader.close();
			throw new RuntimeException("Unexpected format of file: " + meansF.getName());
		}
		if (means.size() == 0)
			throw new RuntimeException("Unable to find mean values in file: " + meansF.getName());
		int i = 0;
		// check to what rank moments are available
		all: while (true) {
			i++;
			System.out.println("Reading moments of rank " + i + " ...");
			List<Moment> moments = Moment.all(i, k);
			for (Moment moment : moments) {
				String name = moment.toStringU();
				if (!means.containsKey(name)) {
					if (i > 1)
						System.out.println("Not found for " + name + ", will compute moments up to rank " + (i - 1));
					else
						System.out.println("Not found for " + name + ", will not compute anything.");
					break all;
				}
				W.add(moment, means.get(name));
				if (verbose)
					System.out.println("W" + moment + "=" + W.mean.get(moment));
			}
		}
		nMax = i - 1;
	}

	/**
	 * Read full W distributions from files in a given directory (not necssary
	 * if we have means)
	 * 
	 * @param directory
	 *            where the files with W distributions reside
	 * @throws IOException
	 */
	@SuppressWarnings("unused")
	private static void readObservationsDistributions(File directory) throws IOException {
		File obsDir = new File(directory, "W");
		if (!obsDir.exists() || !obsDir.isDirectory())
			throw new RuntimeException("Unable to locate directory with mass densities: W");
		int i = 0;
		all: while (true) {
			i++;
			System.out.println("Reading moments of rank " + i + " ...");
			List<Moment> moments = Moment.all(i, k);
			for (Moment moment : moments) {
				int counter = 0;
				// check all possible names
				List<String> names = possibleNames(moment);
				File file = null;
				for (String filename : names) {
					File file1 = new File(obsDir, filename);
					if (file1.exists() && !file1.isDirectory()) {
						file = file1;
						break;
					}
				}
				if (file == null) {
					System.out
							.println("Not found for " + names.get(0) + ", will compute moments up to rank " + (i - 1));
					break all;
				}
				List<Double> xx = new ArrayList<Double>();
				List<Double> yy = new ArrayList<Double>();
				try {
					BufferedReader reader = new BufferedReader(new FileReader(file));
					String line = null;
					while ((line = reader.readLine()) != null) {
						line = line.trim();
						String[] parts = line.split("\\s");
						if (parts.length < 2)
							throw new RuntimeException("Unable to parse observation file: " + file);
						xx.add(Double.parseDouble(parts[parts.length - 2]));
						yy.add(Double.parseDouble(parts[parts.length - 1]));
					}
					reader.close();
				} catch (IOException e) {
					throw new IOException("Unable to read file: " + file, e);
				}
				double[] x = new double[xx.size()];
				double[] y = new double[xx.size()];
				for (int j = 0; j < xx.size(); j++) {
					// correction to avoid negative Ws
					x[j] = xx.get(j) >= 0 ? xx.get(j) : 0.0;
					y[j] = yy.get(j);
					counter += y[j];
				}
				W.add(moment, x, y);
				if (verbose)
					System.out.println("W" + moment + "=" + W.mean.get(moment) + " (" + counter + " samples)");
			}
		}
		nMax = i - 1;
	}

	/**
	 * Generate all possible file names for a given moment: different ordering
	 * is possible
	 * 
	 * @param moment
	 * @return possible names, e.g. (x,y,z)=[0,1,2] =>
	 *         ["W_y_z2.txt","W_z2_y.txt"]
	 */
	private static List<String> possibleNames(Moment moment) {
		List<String> result = new ArrayList<String>();
		for (int[] labeling : moment.nonZeroPermutations()) {
			StringBuilder sb = new StringBuilder();
			for (int i : labeling)
				sb.append("W_" + types.get(i) + (moment.n[i] > 1 ? moment.n[i] : ""));
			String filename = sb.toString() + ".txt";
			result.add(filename);
		}
		return result;
	}

	/**
	 * Compute expected value of identity variables from raw densities
	 */
	private static void computeExpected() {
		for (String bin : bins) {
			SimpleMatrix rhos = rho.get(bin).y;
			SimpleMatrix div = new SimpleMatrix(rhos.numRows(), rhos.numCols());
			for (int i = 0; i < rhos.numRows(); ++i) {
				double sum = rhos.extractVector(true, i).elementSum() + Double.MIN_VALUE;
				for (int j = 0; j < rhos.numCols(); ++j)
					div.set(i, j, sum);
			}
			SimpleMatrix ws = rhos.elementDiv(div);
			w.put(bin, new MassDensities(rho.get(bin).x, ws));
		}
	}

	/**
	 * Generate polynomials representing dependence between sought multiplicity
	 * distribution moments and known W moments
	 */
	private static void generateFormulas() {
		// for each moment of higher order..
		for (Moment moment : newN) {
			// generate a formula defining W using Ns of lower rank
			Polynomial formula = generateFormula(moment);
			if (verbose)
				System.out.println("W" + moment + " = " + formula);
			formulas.put(moment, formula);
		}
	}

	/**
	 * Generate a formula defining a given moment of W with a polynomial defined
	 * over moments of N (see paper for explanation)
	 * 
	 * @param moment
	 *            Moment of W to be computed
	 * @return A polynomial - coefficients for each N of lower rank
	 */
	private static Polynomial generateFormula(Moment moment) {
		List<Moment[]> distributionColumns = moment.getDistributions();
		List<Moment[]> distributionRows = new ArrayList<Moment[]>();
		for (Moment[] distributionColumn : distributionColumns)
			distributionRows.add(transpose(distributionColumn));
		Polynomial result = new Polynomial();
		for (Moment[] bigDistribution : distributionRows) {
			Polynomial bigFormula = new Polynomial(k);
			for (int i = 0; i < bigDistribution.length; ++i) {
				Polynomial localFormula = generateLocalFormula(bigDistribution[i], i);
				bigFormula.multiply(localFormula);
			}
			int multiplier = 1;
			for (int i = 0; i < k; ++i) {
				List<Integer> counts = new ArrayList<Integer>();
				for (Moment momentI : bigDistribution)
					counts.add(momentI.n[i]);
				multiplier *= newton(counts);
			}
			bigFormula.multiply(multiplier);
			result.add(bigFormula);
		}
		result.clean();
		return (result);
	}

	/**
	 * Transpose a matrix defined by a list of moments
	 * 
	 * @param distributionColumn
	 *            List of moments to be treated as columns
	 * @return Array of moments corresponding to rows
	 */
	private static Moment[] transpose(Moment[] distributionColumn) {
		Moment[] result = new Moment[k];
		for (int i = 0; i < k; ++i) {
			Moment resultI = new Moment(k);
			for (int j = 0; j < k; ++j)
				resultI.n[j] = distributionColumn[j].n[i];
			result[i] = resultI;
		}
		return result;
	}

	/**
	 * Generate a sub-formula contributing to a definition of a given moment of
	 * W (see paper for explanation)
	 * 
	 * @param moment
	 *            Moment of W to be defined
	 * @param i
	 * @return part of a polynomial
	 */
	private static Polynomial generateLocalFormula(Moment moment, int i) {
		if (moment.sum() == 0)
			return new Polynomial(moment.n.length);
		List<Map<Moment, Integer>> configurations = moment.getDistributionsInf();
		Polynomial result = new Polynomial();
		for (Map<Moment, Integer> configuration : configurations) {
			double val = 1;
			List<Integer> bigCounts = new ArrayList<Integer>(configuration.size());
			for (Moment momentL : configuration.keySet()) {
				for (int j = 0; j < configuration.get(momentL); ++j)
					val *= u(momentL, i);
				bigCounts.add(configuration.get(momentL));
			}
			Polynomial a = newtonOpen(i, bigCounts);
			for (int jj = 0; jj < k; ++jj) {
				List<Integer> smallCounts = new ArrayList<Integer>();
				for (Moment momentL : configuration.keySet())
					for (int j = 0; j < configuration.get(momentL); ++j)
						if (momentL.n[jj] != 0)
							smallCounts.add(momentL.n[jj]);
				if (!smallCounts.isEmpty())
					a.multiply(newton(smallCounts));
			}
			a.multiply(val);
			result.add(a);
		}
		result.clean();
		return result;
	}

	/**
	 * Compute function u for a given moment of w and a chosen rho (see paper
	 * for explanation)
	 * 
	 * @param moment
	 * @param i
	 * @return
	 */
	private static double u(Moment moment, int i) {
		if (moment.sum() == 0)
			return 1;
		// if u() has already been computed for these parameters, avoid doing it
		// again
		if (uCache.containsKey(moment) && !Double.isNaN(uCache.get(moment)[i]))
			return uCache.get(moment)[i];
		double sum = 0;
		for (int j = 0; j < bins.size(); ++j) {
			SimpleMatrix product = rho.get(bins.get(j)).y.extractVector(false, i);
			for (int kk = 0; kk < k; ++kk)
				for (int kkk = 0; kkk < moment.n[kk]; ++kkk)
					product = product.elementMult(w.get(bins.get(j)).y.extractVector(false, kk));
			sum += product.elementSum();
		}
		Moment Ni = new Moment(k);
		Ni.n[i] = 1;
		double result;
		if (!UN.containsKey(Ni))
			throw new RuntimeException("Moment not found: " + Ni);
		result = sum / UN.get(Ni);
		if (!uCache.containsKey(moment)) {
			double[] array = new double[k];
			for (int j = 0; j < k; ++j)
				array[j] = Double.NaN;
			uCache.put(moment, array);
		}
		uCache.get(moment)[i] = result;
		return result;
	}

	/**
	 * Computes a multinomial coefficient
	 * \binom{Ni}{a1,a2,a3,..,ak,Ni-sum(a1..ak)} as a polynomial with variable
	 * Ni
	 * 
	 * @param i
	 *            points to variable Ni to be used in the polynomial
	 * @param counts
	 *            parameters of a coefficient
	 * @return Multinomial coefficient as a polynomial for Ni
	 */
	private static Polynomial newtonOpen(int i, List<Integer> counts) {
		int sum = 0;
		for (Integer j : counts)
			sum += j;
		Polynomial result = new Polynomial(k);
		for (int j = 0; j < sum; ++j) {
			Polynomial Ni = new Polynomial();
			int[] nn = new int[k];
			nn[i] = 1;
			Ni.m.put(new Moment(nn), 1.0);
			Ni.add(new Polynomial(k, -j));
			result.multiply(Ni);
		}
		int factorials = 1;
		for (Integer j : counts)
			factorials *= factorial(j);
		result.multiply(1.0 / factorials);
		return result;
	}

	/**
	 * Compute a value of multinomial coefficient
	 * \binom{a1+a2+a3+..+ak}{a1,a2,a3,..,ak}
	 * 
	 * @param smallCounts
	 *            parameters of a coefficient
	 * @return value of a multinomial coefficient
	 */
	private static int newton(List<Integer> smallCounts) {
		int sum = 0;
		for (Integer i : smallCounts)
			sum += i;
		int result = factorial(sum);
		for (Integer i : smallCounts)
			result /= factorial(i);
		return result;
	}

	/**
	 * Factorial calculation
	 * 
	 * @param n
	 *            Factorial argument
	 * @return factorial value
	 */
	private static int factorial(int n) {
		int result = 1;
		for (int i = 1; i <= n; ++i)
			result *= i;
		return result;
	}

	/**
	 * Prepare a representation of the task as a linear problem
	 */
	private static void prepareMatrices() {
		double[][] bb = new double[newN.size()][];
		double[][] AA = new double[newN.size()][];
		for (int i = 0; i < newN.size(); ++i) {
			AA[i] = new double[newN.size()];
			// any ordering of moments is good, if used consistently
			Moment moment = newN.get(i);
			double[] val = { W.mean.get(moment) };
			for (Moment other : formulas.get(moment).m.keySet())
				if (other.sum() < n)
					val[0] -= N.get(other) * formulas.get(moment).m.get(other);
				else
					AA[i][newN.indexOf(other)] = formulas.get(moment).m.get(other);
			bb[i] = val;
		}
		b = new SimpleMatrix(bb);
		A = new SimpleMatrix(AA);
	}

	/**
	 * Solve a linear problem and write results as TSV
	 * 
	 * @param bw
	 *            Writer to write results to
	 * @throws IOException
	 */
	private static void solve(BufferedWriter bw) throws IOException {
		SimpleMatrix newNval;
		try {
			// simple!
			newNval = A.solve(b);
		} catch (SingularMatrixException e) {
			throw new RuntimeException("Singular matrix");
		}
		System.out.println("Solved!");
		for (int i = 0; i < newN.size(); ++i) {
			if (verbose)
				System.out.println("N" + newN.get(i) + "=" + newNval.get(i));
			if (bw != null) {
				for (int j = 0; j < newN.get(i).n.length; ++j)
					bw.write("" + newN.get(i).n[j] + "\t");
				bw.write("" + newNval.get(i) + "\n");
			}
			N.put(newN.get(i), newNval.get(i));
		}
		newN = null;
		newNval = null;
	}

}
