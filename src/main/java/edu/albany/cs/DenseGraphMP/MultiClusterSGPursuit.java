package edu.albany.cs.DenseGraphMP;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.ObjectOutputStream;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;
import org.apache.commons.math3.stat.descriptive.rank.Median;

import com.google.common.collect.Sets;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.base.Matrix;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.MultiDataGraph;
import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.PCAScore;
import edu.albany.cs.scoreFuncs.Stat;
import edu.albany.cs.tailApprox.TailApprox;

/**
 * Algorithm 1 : S2Graph-Mp Aglorithm .
 *
 * @author
 */
public class MultiClusterSGPursuit {

	/**
	 * graph info
	 */
	private final String fileName;
	private final int n;
	private final int p;
	private final MultiDataGraph graph;
	/**
	 * the total sparsity of x
	 */
	private final int k;
	/**
	 * the maximum number of |y|<s
	 */
	private final int s;
	/**
	 * the maximum number of connected components formed by F
	 */
	private final int g;
	/**
	 * bound on the total weight w(F) of edges in the forest F
	 */
	private final double B;
	/**
	 * number of iterations
	 */
	private final int iterations;
	/**
	 * Function to be optimized.
	 */
	private final Function function;
	private double lambda;
	/**
	 * results
	 */
	public double[] xi;
	public double[] yi;
	public Set<Integer> psiX;
	public Set<Integer> psiY;
	public double funcValue = -1.0D;
	public double runTime;

	public static int verboseLevel = 0;



	/*
	 * k: sparsity of nodes s: sparsity of attributes g: number of connected
	 * components B: budget iterations: The maximum number of iterations.
	 */
	
	public MultiClusterSGPursuit() {
		this.graph = null;
		this.n = 0;
		this.p = 0;
		this.k = 0; // sparsity of x
		this.s = 0; // sparsity of y
		this.g = 0;
		this.B = 0;
		this.fileName = null;
		this.iterations = 0;
		this.function = null;
		this.lambda = 0;
	}
	public MultiClusterSGPursuit(String fileName,
			MultiDataGraph graph, int k, int s, int g, double B,
			int iterations, Function func, double lambda) {
		this.graph = graph;
		this.n = graph.numOfNodes;
		this.p = graph.numOfFeas;
		this.k = k; // sparsity of x
		this.s = s; // sparsity of y
		this.g = g;
		this.B = B;
		this.fileName = fileName;
		this.iterations = iterations;
		this.function = func;
		this.lambda = lambda;

		if (run()) {
			// System.out.println("S2GraphMPDense successfully finished.");
		} else {
			System.out.println("S2GraphMPDense is failed.");
		}
	}

	private int[] ranknodes(ArrayList<Integer> nodes, int k) {
		double[] degrees = new double[nodes.size()];
		int[] rank = new int[Math.min(nodes.size(), k)];
		Arrays.fill(degrees, 0);
		for (int i = 0; i < nodes.size(); i++) {
			for (int j : nodes) {
				if (graph.AdjMatrix[nodes.get(i)][j] > 0) {
					degrees[i] += 1;
				}
			}
		}
		ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(degrees);
		Integer[] indexes = arrayIndexComparator.getIndices();
		Arrays.sort(indexes, arrayIndexComparator);
		for (int i = 0; i < rank.length; i++) {
			rank[i] = nodes.get(indexes[i]);
		}
		return rank;
	}

	/*
	 * This function identifies good initial values for x and y. k: sparsity of
	 * x (nodes) s: sparsity of y (attributes) OUTPUT: return the initial values
	 * of x and y.
	 */
	private List<double[]> CalcInitilvalues(int k, int s) {
		double[] x0 = new double[n];
		double[] y0 = new double[p];
		Arrays.fill(x0, 0);
		Arrays.fill(y0, 0);
		double threshold = 0.11;
		ArrayList<ArrayList<Integer>> res = new ArrayList<ArrayList<Integer>>();
		for (int j = 0; j < p; j++) {
			for (int i = 0; i < n; i++) {
				ArrayList<Integer> nns = graph.arrayListAdj.get(i);
				Collections.sort(nns);
				double[] vals = new double[nns.size()];
				for (int k1 = 0; k1 < nns.size(); k1++) {
					vals[k1] = graph.W[nns.get(k1)][j];
				}
				if (nns.size() > 6) {
					x0[i] = Math.abs(graph.W[i][j]
							- StatUtils.percentile(vals, 50));
				} else {
					x0[i] = Math.abs(graph.W[i][j] - StatUtils.mean(vals));
				}
			}
			ArrayList<Integer> S = new ArrayList<Integer>();
			for (int i = 0; i < n; i++) {
				if (x0[i] < threshold) {
					S.add(i);
				}
			}
			// if(S.size() > 0){
			// S = Utils.getLargestConnectedComponant(graph.arrayListAdj, S);
			// }else{
			// S = new ArrayList<Integer>();
			// }
			res.add(S);
		}
		Integer[] indexes = Utils.sortArrayList1(res);
		for (int i = 0; i < s; i++) {
			y0[indexes[i]] = 1;
		}
		Arrays.fill(x0, 0);
		int[] rank = ranknodes(res.get(indexes[0]), k);
		for (int i = 0; i < k; i++) {
			x0[rank[i]] = 1;
		}

		// String line = "features: ";
		// double[] xx = new double[n];
		// Arrays.fill(xx, 0);
		// for(int i = 0; i < s; i++){
		// line += indexes[i] + ",";
		// y0[indexes[i]] = 1;
		// for(int j:res.get(indexes[i])){
		// xx[j] += 1;
		// }
		// }
		// // System.out.println(line);
		// Arrays.fill(x0, 0);
		// line = "nodes: ";
		// int cnt = 0;
		// for(int i=0; i<n; i++){
		// if(xx[i] == s){
		// line += i + ",";
		// x0[i] = 1;
		// cnt += 1;
		// }
		// }
		// if(cnt <5){
		// for(int i:res.get(indexes[0])){
		// line += i + ",";
		// x0[i] = 1;
		// }
		// }
		// System.out.println(line);
		List<double[]> list = new ArrayList<>();
		list.add(x0);
		list.add(y0);
		return list;
	}

	private List<double[]> CalcInitilvalues1(int k, int s) {
		double threshold = -0.7;
		double[] x0 = new double[n];
		double[] y0 = new double[p];
		Arrays.fill(x0, 0);
		Arrays.fill(y0, 0);
		ArrayList<int[]> res = new ArrayList<int[]>();
		double[] scores = new double[p];
		for (int j = 0; j < p; j++) {
			ArrayList<Integer> S = new ArrayList<Integer>();
			for (int i = 0; i < n; i++) {
				ArrayList<Integer> nns = graph.arrayListAdj.get(i);
				// Collections.sort(nns);
				double[] dists = new double[nns.size()];
				for (int k1 = 0; k1 < nns.size(); k1++) {
					dists[k1] = Math.abs(graph.W[nns.get(k1)][j]
							- graph.W[i][j]);
				}
				if (nns.size() > 1) {
					x0[i] = StatUtils.percentile(dists, 50) * -1;
				} else {
					x0[i] = StatUtils.mean(dists) * -1;
				}
				if (x0[i] > threshold)
					S.add(i);
			}
			ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(x0);
			Integer[] indexes = arrayIndexComparator.getIndices();
			Arrays.sort(indexes, arrayIndexComparator);
			S = new ArrayList<Integer>();
			for (int i = 0; i < k * 5; i++) {
				S.add(indexes[i]);
			}
			if (S.size() > 0) {
				int[] rank = ranknodes(S, k);
				double[] values = new double[rank.length];
				for (int i = 0; i < rank.length; i++) {
					values[i] = graph.W[rank[i]][j];
				}
				double score = StatUtils.populationVariance(values) * -1;
				Arrays.fill(x0, 0);
				for (int i : rank) {
					x0[i] = 1;
				}
				Arrays.fill(y0, 0);
				y0[j] = 1;
				funcValue = function.getFuncValue(x0, y0, lambda);
				if (funcValue >= 0)
					scores[j] = -1000;
				else
					scores[j] = score;
				res.add(rank);
			} else {
				int[] rank = new int[] {};
				scores[j] = -1000;
				res.add(rank);
			}
		}
		ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(scores);
		Integer[] indexes = arrayIndexComparator.getIndices();
		Arrays.sort(indexes, arrayIndexComparator);
		Arrays.fill(y0, 0);
		for (int i = 0; i < s; i++) {
			if (scores[indexes[i]] > -1000)
				y0[indexes[i]] = 1;
		}
		if (StatUtils.sum(y0) == 0)
			y0[indexes[0]] = 1;
		Arrays.fill(x0, 0);
		for (int i : res.get(indexes[0])) {
			x0[i] = 1;
		}
		if (StatUtils.sum(x0) == 0) {
			System.out
					.println("NOTICE: no good initialization can be idetnfied!!!!!!!!!!!!!!!!!!!!!!!!!");
			x0[0] = 1;
			x0[1] = 1;
			x0[2] = 1;
		}
		List<double[]> list = new ArrayList<>();
		list.add(x0);
		list.add(y0);
		return list;
	}

	private List<double[]> CalcInitilvalues2(int k, int s) {
		double[] x0 = new double[n];
		double[] y0 = new double[p];
		Arrays.fill(x0, 0);
		Arrays.fill(y0, 0);
		ArrayList<ArrayList<int[]>> res = new ArrayList<ArrayList<int[]>>();
		ArrayList<double[]> scores = new ArrayList<double[]>();
		int[] trials = new int[] { 2, 3, 4, 5, 6 };
		for (int i = 0; i < trials.length; i++) {
			res.add(new ArrayList<int[]>());
			scores.add(new double[p]);
		}
		for (int j = 0; j < p; j++) {
			ArrayList<Integer> S = new ArrayList<Integer>();
			for (int i = 0; i < n; i++) {
				ArrayList<Integer> nns = graph.arrayListAdj.get(i);
				double[] dists = new double[nns.size()];
				for (int k1 = 0; k1 < nns.size(); k1++) {
					dists[k1] = Math.abs(graph.W[nns.get(k1)][j]
							- graph.W[i][j]);
				}
				if (nns.size() > 1) {
					x0[i] = StatUtils.percentile(dists, 50) * -1;
				} else {
					x0[i] = StatUtils.mean(dists) * -1;
				}
			}
			ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(x0);
			Integer[] indexes = arrayIndexComparator.getIndices();
			Arrays.sort(indexes, arrayIndexComparator);

			int ii = 0;
			for (int r : trials) {
				if (k * r > n) {
					continue;
				}
				S = new ArrayList<Integer>();
				for (int i = 0; i < k * r; i++) {
					S.add(indexes[i]);
				}
				if (S.size() > 0) {
					int[] rank = ranknodes(S, k);
					double fval = getFun(rank, new int[] { j }, 0);
					scores.get(ii)[j] = fval * -1;
					res.get(ii).add(rank);
					// double[] values = new double[rank.length];
					// for(int i=0; i<rank.length; i++){
					// values[i] = graph.W[rank[i]][j];
					// }
					// double score = 0;
					// if(rank.length > 0)
					// score = StatUtils.populationVariance(values) * -1;
					// else
					// score = -1000;
					// scores.get(ii)[j] = score;
					// res.get(ii).add(rank);
				} else {
					int[] rank = new int[] {};
					scores.get(ii)[j] = -1000;
					res.get(ii).add(rank);
				}
				ii += 1;
			}
		}
		double fValue = -1;
		ArrayList<Integer> Y = new ArrayList<Integer>();
		ArrayList<Integer> X = new ArrayList<Integer>();
		for (int ii = 0; ii < trials.length; ii++) {
			ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(
					scores.get(ii));
			Integer[] indexes = arrayIndexComparator.getIndices();
			Arrays.sort(indexes, arrayIndexComparator);
			ArrayList<Integer> Y1 = new ArrayList<Integer>();
			for (int i = 0; i < s; i++) {
				if (scores.get(ii)[indexes[i]] > -1000)
					Y1.add(indexes[i]);
			}

			if (res.get(ii).size() < 1)
				continue;
			int[] X1 = res.get(ii).get(indexes[0]);
			double fValue1 = getFun(X1, Y1, 0);
			if (fValue == -1 || fValue > fValue1) {
				Y.clear();
				for (int i : Y1)
					Y.add(i);
				X.clear();
				for (int i : X1)
					X.add(i);
				fValue = fValue1;
			}
			for (int j = 0; j < s; j++) {
				X1 = res.get(ii).get(indexes[j]);
				Y1.clear();
				Y1.add(indexes[j]);
				fValue1 = getFun(X1, Y1, 0);
				// fValue1 = scores.get(ii)[j];
				if (fValue == -1 || fValue > fValue1) {
					Y.clear();
					for (int i : Y1)
						Y.add(i);
					X.clear();
					for (int i : X1)
						X.add(i);
					fValue = fValue1;
				}
			}
			for (int j = 0; j < s * 2; j++) {
				if (j >= p) {
					continue;
				}
				if (Y.contains(indexes[j]) == false) {
					ArrayList<Integer> YY = new ArrayList<Integer>();
					for (int i : Y)
						YY.add(i);
					YY.add(indexes[j]);
					fValue1 = getFun(X, YY, 0);
					if (fValue == -1 || fValue > fValue1) {
						Y.clear();
						for (int i : YY)
							Y.add(i);
						fValue = fValue1;
					}
				}
			}
		}
		Arrays.fill(y0, 0);
		for (int i : Y)
			y0[i] = 1;
		Arrays.fill(x0, 0);
		for (int i : X)
			x0[i] = 1;
		if (StatUtils.sum(x0) == 0) {
			System.out
					.println("NOTICE: no good initialization of x can be idetnfied!!!!!!!!!!!!!!!!!!!!!!!!!");
			x0[0] = 1;
			x0[1] = 1;
			x0[2] = 1;
		}
		if (StatUtils.sum(y0) == 0) {
			System.out
					.println("NOTICE: no good initialization of y can be idetnfied!!!!!!!!!!!!!!!!!!!!!!!!!");
			y0[0] = 1;
		}
		List<double[]> list = new ArrayList<>();
		list.add(x0);
		list.add(y0);
		return list;
	}

	public double getFun(int[] X, ArrayList<Integer> Y, double lambda) {
		double[] x = new double[n];
		double[] y = new double[p];
		Arrays.fill(x, 0);
		for (int i : X) {
			x[i] = 1;
		}
		Arrays.fill(y, 0);
		for (int i : Y) {
			y[i] = 1;
		}
		return function.getFuncValue(x, y, lambda);
	}

	public double getFun(ArrayList<Integer> X, ArrayList<Integer> Y,
			double lambda) {
		double[] x = new double[n];
		double[] y = new double[p];
		Arrays.fill(x, 0);
		for (int i : X) {
			x[i] = 1;
		}
		Arrays.fill(y, 0);
		for (int i : Y) {
			y[i] = 1;
		}
		return function.getFuncValue(x, y, lambda);
	}

	public double getFun(int[] X, int[] Y, double lambda) {
		double[] x = new double[n];
		double[] y = new double[p];
		Arrays.fill(x, 0);
		for (int i : X) {
			x[i] = 1;
		}
		Arrays.fill(y, 0);
		for (int i : Y) {
			y[i] = 1;
		}
		return function.getFuncValue(x, y, lambda);
	}

	public double sum(double[] x) {
		double sum = 0;
		for (int i = 0; i < x.length; i++) {
			sum += x[i];
		}
		return sum;
	}

	public static void showstat(double[] x, String title, double lowerbound,
			double upperbound) {
		String line = "";
		for (int j = 0; j < x.length; j++) {
			if (x[j] != 0 && x[j] > lowerbound && x[j] < upperbound) {
				line += ", (" + j + ", " + x[j] + ")";
			}
		}
		System.out.println(title + line);
	}

	public static void showstat(double[] x, String title) {
		String line = "";
		for (int j = 0; j < x.length; j++) {
			if (x[j] != 0) {
				// line += ", (" + j + ", " + x[j] + ")";
				line += "," + x[j];
			}
		}
		System.out.println(title + line);
	}

	public void showstat(int[] x, String title) {
		String line = "";
		for (int j = 0; j < x.length; j++) {
			if (x[j] != 0) {
				line += ", (" + j + ", " + x[j] + ")";
			}
		}
		System.out.println(title + line);
	}

	private void initilization1(double[] xi, double[] yi, MultiDataGraph graph) {
		int cnt = 0;
		for (int j : graph.trueSubGraphNodes) {
			xi[j] = 1;
			cnt += 1;
			if (cnt > 5)
				break;
		}
		cnt = 0;
		for (int j : graph.trueFeas) {
			yi[j] = 1;
			cnt += 1;
			if (cnt > 2)
				break;
		}
	}

	private void StoreTestCaseData(double[] xi, double[] yi,
			Set<Integer> omegaX, Set<Integer> omegaY, double lambda,
			String fileName, String outputfilename) {
		try {
			FileOutputStream fos = new FileOutputStream(outputfilename);
			ObjectOutputStream oos = new ObjectOutputStream(fos);
			oos.writeObject(xi); // write MenuArray to ObjectOutputStream
			oos.writeObject(yi);
			oos.writeObject(Set2IntArray(omegaX));
			oos.writeObject(Set2IntArray(omegaY));
			oos.writeObject(lambda);
			oos.writeObject(fileName);
			oos.close();
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		Path p = Paths.get(fileName);
		String file = p.getFileName().toString();
		try {
			copyFile(new File(fileName), new File("testcases/getArgMinFxy/"
					+ file));
		} catch (Exception ex) {
			ex.printStackTrace();
		}
		// try {
		// FileInputStream istream = new FileInputStream(outputfilename);
		// ObjectInputStream ois = new ObjectInputStream(istream);
		// xi = (double []) ois.readObject();
		// yi = (double []) ois.readObject();
		// omegaX = IntArray2Set((int[]) ois.readObject());
		// omegaY = IntArray2Set((int[]) ois.readObject());
		// lambda = (double) ois.readObject();
		// fileName = (String) ois.readObject();
		// } catch(Exception ex) {
		// ex.printStackTrace();
		// }

	}

	public double[] median_WTx(double[] x) {
		int p = graph.W[0].length;
		double sum = StatUtils.sum(x);
		double[] medians = new double[p];
		ArrayList<Integer> S = new ArrayList<Integer>();
		for (int i = 0; i < x.length; i++) {
			if (x[i] != 0) {
				S.add(i);
			}
		}
		for (int j = 0; j < p; j++) {
			double[] vals = new double[S.size()];
			for (int k = 0; k < S.size(); k++) {
				vals[k] = x[S.get(k)] * graph.W[S.get(k)][j];
			}
			medians[j] = StatUtils.percentile(vals, 50) * sum;
		}
		return medians;
	}

	public double[] squaretimevector(double[][] M, double[] v) {
		double[] out = new double[M.length];
		for (int i = 0; i < M.length; i++) {
			out[i] = 0;
			for (int j = 0; j < v.length; j++) {
				out[i] += M[i][j] * M[i][j] * v[j];
			}
		}
		return out;
	}

	// double[] Ix = new double[n];
	// for(int i:graph.trueNodesSet){
	// xi[i] = 1;
	// }
	//
	// for(int i:graph.trueFeas){
	// yi[i] = 1;
	// }
	// Arrays.fill(Ix, 1.0D);
	// double[][] W_transpose = Matrix.transpose(graph.W);
	// double xT1 = Matrix.dot(xi, Ix);
	// double[] WTx;
	// WTx=Matrix.VecMultiplyMat(xi,graph.W);
	// double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
	// showstat(WTx, "WTx");
	// showstat(WTx_1Tx, "WTx_1Tx");
	// WTx = median_WTx(xi);
	// WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
	// showstat(WTx, "WTx");
	// showstat(WTx_1Tx, "WTx_1Tx");
	// double[] comp1 = Matrix.VarMultiplyVec(1.0/0.1,
	// squaretimevector(Matrix.subtract(W_transpose,Matrix.VecOutterPro(WTx_1Tx,
	// Ix)), xi));
	// double[] comp2 = Matrix.VarMultiplyVec(1.0/1,
	// squaretimevector(W_transpose, xi));
	// double[] gradient = Matrix.VecSubstract(comp1, comp2);
	// showstat(gradient, "gradient");
	// initilization1(xi, yi, graph);

	private boolean run() {
		long startTime = System.nanoTime();
		List<double[]> xy0 = CalcInitilvalues2(20, 3);
		runTime = (System.nanoTime() - startTime) / 1e9;
		// System.out.println("CalcInitilvalues2 run time:" + runTime);
		double[] xi = new double[n];
		double[] yi = new double[p];
		System.arraycopy(xy0.get(0), 0, xi, 0, xi.length);
		System.arraycopy(xy0.get(1), 0, yi, 0, yi.length);
		List<Double> fValues = new ArrayList<>();
		double oldFunValue = -1;
		for (int i = 0; i < this.iterations; i++) {
			// System.out.println("iteration: " + i);
			fValues.add(function.getFuncValue(xi, yi, lambda));
			// if(verboseLevel>0){
			// System.out.println(sum(xi) + ", " + sum(yi));
			// }
			double[] gradientFx = function.getGradientX(xi, yi, lambda);
			double[] gradientFy = function.getGradientY(xi, yi, lambda);
			if (verboseLevel > 0) {
				// System.out.println(StatUtils.min(gradientFx) + ", " +
				// StatUtils.max(gradientFx));
				// System.out.println(StatUtils.min(gradientFy) + ", " +
				// StatUtils.max(gradientFy));
				showstat(xi, "x:");
				showstat(yi, "y:");
				showstat(gradientFx, "gradientFx: ", -1000, 0);
				showstat(gradientFy, "gradientFy: ", -1000, 0);
			}
			// if(true) return true;

			gradientFx = Stat.normalizeGradient(xi, gradientFx, n);
			gradientFy = Stat.normalizeGradient(yi, gradientFy, p);
			if (verboseLevel > 0) {
				System.out.println(StatUtils.min(gradientFx) + ", "
						+ StatUtils.max(gradientFx));
				System.out.println(StatUtils.min(gradientFy) + ", "
						+ StatUtils.max(gradientFy));
				for (int j = 0; j < gradientFy.length; j++) {
					if (gradientFy[j] != 0) {
						System.out.println(j + ", " + gradientFy[j]);
					}
				}
			}
			HeadApprox head = new HeadApprox(graph.edges, graph.edgeCosts,
					gradientFx, k, g, B, graph.trueSubGraphNodes, false);
			Set<Integer> gammaX = new HashSet<>(head.bestForest.nodesInF);
			Set<Integer> gammaY = identifyDirection(gradientFy, 2 * s);
			Set<Integer> omegaX = Sets.union(gammaX, Stat.supp(xi));
			// System.out.println(omegaX.size());
			Set<Integer> omegaY = Sets.union(gammaY, Stat.supp(yi));
			// StoreTestCaseData(xi, yi, omegaX, omegaY, lambda, fileName,
			// "testcases/getArgMinFxy/test-case2.dat");
			List<double[]> b = function.getArgMinFxy(xi, yi, omegaX, omegaY,
					lambda);
			/**
			 * (b_x,b_y)= argMax f(x,y) s.t. supp(x) \in omegaX, supp(y) \in
			 * omegaY
			 */
			double[] bx = b.get(0);
			double[] by = b.get(1);
			TailApprox tail = new TailApprox(graph.edges, graph.edgeCosts, bx,
					k, g, B, graph.trueSubGraphNodes, false);
			psiX = new HashSet<>(tail.bestForest.nodesInF);
			psiY = identifyDirection(by, s);
			xi = projectionOnVector(bx, psiX);
			yi = projectionOnVector(by, psiY);
			funcValue = function.getFuncValue(xi, yi, lambda);
			if (oldFunValue == -1) {
				oldFunValue = funcValue;
			} else {
				if (Math.abs(oldFunValue - funcValue) < 0.001) {
					break;
				} else {
					oldFunValue = funcValue;
				}
			}
		}

		// for (int i = 0; i < this.iterations; i++) {
		// fValues.add(function.getFuncValue(xi, yi,lambda));
		// double[] gradientFy = function.getGradientY(xi, yi,lambda);
		// gradientFy = Stat.normalizeGradient(yi, gradientFy, p);
		// Set<Integer> gammaY = identifyDirection(gradientFy, 2 * s);
		// Set<Integer> omegaY = Sets.union(gammaY, Stat.supp(yi));
		// double[] by = function.getArgMinFy(xi, omegaY,lambda); /** (b_x,b_y)=
		// argMax f(x,y) s.t. supp(x) \in omegaX, supp(y) \in omegaY*/
		// psiY = identifyDirection(by, s);
		// yi = projectionOnVector(by, psiY);
		// funcValue = function.getFuncValue(xi, yi,lambda);
		// }
		runTime = (System.nanoTime() - startTime) / 1e9;
		return true;
	}

	private static void copyFile(File sourceFile, File destFile)
			throws IOException {
		if (!sourceFile.exists()) {
			return;
		}
		if (!destFile.exists()) {
			destFile.createNewFile();
		}
		FileChannel source = null;
		FileChannel destination = null;
		source = new FileInputStream(sourceFile).getChannel();
		destination = new FileOutputStream(destFile).getChannel();
		if (destination != null && source != null) {
			destination.transferFrom(source, 0, source.size());
		}
		if (source != null) {
			source.close();
		}
		if (destination != null) {
			destination.close();
		}

	}

	public int[] Set2IntArray(Set<Integer> x) {
		int[] data = new int[x.size()];
		int i = 0;
		for (int j : x) {
			data[i] = j;
			i += 1;
		}
		return data;
	}

	public Set<Integer> IntArray2Set(int[] x) {
		Set<Integer> mySet = new HashSet<Integer>();
		for (int i : x) {
			mySet.add(i);
		}
		return mySet;
	}



	private Set<Integer> identifyDirection(double[] gradientFy, int s) {
		Set<Integer> GammaY = new HashSet<>();
		double[] absGradientFy = new double[gradientFy.length];
		for (int i = 0; i < gradientFy.length; i++) {
			absGradientFy[i] = Math.abs(gradientFy[i]);
		}
		ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(absGradientFy);
		Integer[] indexes = arrayIndexComparator.getIndices();
		Arrays.sort(indexes, arrayIndexComparator);
		for (int i = 0; i < s; i++) {
			if (i < indexes.length && gradientFy[indexes[i]] != 0) {
				GammaY.add(indexes[i]);
			}
		}
		return GammaY;
	}

	private double[] projectionOnVector(double[] vector, Set<Integer> set) {
		if (vector == null || vector.length == 0) {
			return null;
		}
		double[] projectionVector = new double[vector.length];
		for (int i = 0; i < vector.length; i++) {
			if (set.contains(i)) {
				projectionVector[i] = vector[i];
			} else {
				projectionVector[i] = 0.0D;
			}
		}
		return projectionVector;
	}

	/**
	 * *************************************************************************
	 * *****************************************************************
	 */
	public static void MultiClusterDetection(String fileName,
			String resultFileName, int k, int s, double lambda, int clusterNum) {
		long startTimeAll = System.nanoTime();
		APDMInputFormat apdm = null;
		MultiDataGraph mgdGraph = null;
		int g = 1;
		int iterations = 10;
		double B = k - g + 0.0D;

		ArrayList<Set<Integer>> clusters = new ArrayList<Set<Integer>>();
		NormalDistribution norm =null;
		Random r=new Random();
		for (int c = 0; c < clusterNum; c++) {
			
			apdm = new APDMInputFormat(fileName);

			mgdGraph = new MultiDataGraph(apdm);
			PCAScore func = new PCAScore(mgdGraph.W, mgdGraph.AdjMatrix);

			MultiClusterSGPursuit s2GraphMPDense = new MultiClusterSGPursuit(
					fileName, mgdGraph, k, s, g, B, iterations, func, lambda);
			clusters.add(s2GraphMPDense.psiX);
			

			for (int i : s2GraphMPDense.psiX) {
				// System.out.println("Old:" + i + " "
				// + ArrayUtils.toString(mgdGraph.W[i]));
				for (int j : s2GraphMPDense.psiY) {
					double mu = r.nextInt(99) + 10;
					norm = new NormalDistribution(mu, 5);
					// System.out.print(mgdGraph.W[i][j] + " ");
					mgdGraph.W[i][j] = norm.sample();
					// System.out.print(mgdGraph.W[i][j] + "\n");
				}
				// System.out.println("New:" +
				// ArrayUtils.toString(mgdGraph.W[i]));
				// Utils.stop();
			}

		}
		/** Calculate Normalized cut for top 1 result ***/
		double nCut = 0.0D;
		int N = apdm.data.numNodes;
		double[] I = Matrix.identityVector(N);
		double tempNCut = 100000.0D;
		double[][] adjMatrix = new double[N][N];
		for (int i = 0; i < N; i++) {
			for (int j : apdm.data.graphAdj[i]) {
				adjMatrix[i][j] = 1.0;
				adjMatrix[j][i] = 1.0;
			}

		}
		for (Set<Integer> cluster : clusters) {
			double[] a1 = indicateVector(cluster, N); // top
																// 1
																// results
																// indicatr
		// vector
		double[] a1C = Matrix.VecSubstract(I, a1);
		double aW1_a = Matrix.dot(Matrix.VecMultiplyMat(a1, adjMatrix), a1C);
		double aw1 = Matrix.dot(Matrix.VecMultiplyMat(a1, adjMatrix), I);
		nCut = aW1_a / aw1;
			System.out.println("X:-" + cluster.toString());
		System.out.println("Normalized Cut=:" + nCut);
		}
		double runTime = (System.nanoTime() - startTimeAll) / 1e9;
		System.out.println("Total running time: " + runTime);

		/*
		 * double time = (System.nanoTime() - startTimeAll) / 1e9; try { File f
		 * = new File(resultFileName); if(f.exists()){ f.delete(); }
		 * FileOutputStream fos = new FileOutputStream(resultFileName);
		 * ObjectOutputStream oos = new ObjectOutputStream(fos);
		 * oos.writeObject(fileName); // write MenuArray to ObjectOutputStream
		 * oos.writeObject(stat); oos.close(); } catch(Exception ex) {
		 * System.out.println("!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! error");
		 * ex.printStackTrace(); }
		 */
	}

	private static void MultiClusterDetectionOn3Measure(String fileName,
			String resultFileName, int k, int s, double lambda, int clusterNum,
			String dataset) {
		long startTimeAll = System.nanoTime();
		APDMInputFormat apdm = new APDMInputFormat(fileName);
		MultiDataGraph mgdGraph = new MultiDataGraph(apdm);
		int g = 1;
		int iterations = 5;
		double B = k - g + 0.0D;

		ArrayList<SubSpace> clusterRes = new ArrayList<SubSpace>();
		NormalDistribution norm = null;

		Random r = new Random();
		System.out.print("Cluster(" + clusterNum + "): ");
		for (int c = 0; c < clusterNum; c++) {
			System.out.print(c + " ");

		
			PCAScore func = new PCAScore(mgdGraph.W, mgdGraph.AdjMatrix);

			MultiClusterSGPursuit s2GraphMPDense = new MultiClusterSGPursuit(
					fileName, mgdGraph, k, s, g, B, iterations, func, lambda);

			clusterRes.add(new SubSpace(s2GraphMPDense.psiX,
					s2GraphMPDense.psiY));
			norm = new NormalDistribution(0, 1);
			for (int i : s2GraphMPDense.psiX) {
				for (int j : s2GraphMPDense.psiY) {
					mgdGraph.W[i][j] = norm.sample();
				}
			}
		}
		System.out.println();
		int N = apdm.data.numNodes;
		double[][] adjMatrix = new double[N][N];
		for (int i = 0; i < N; i++) {
			for (int j : apdm.data.graphAdj[i]) {
				adjMatrix[i][j] = 1.0;
				adjMatrix[j][i] = 1.0;
			}

		}
		/********************************************************
		 * Calculate : 1. Node density: x^TAx/1^Tx 2. Attribute coherence: Sum
		 * of Euclideant 3. Cluster size
		 ********************************************************/
		// laod original atribute matrix W
		mgdGraph = new MultiDataGraph(apdm);
		String outResult = dataset + " " + k + " " + s + " ";
		String scores = "";
		List<Integer> c_num = (List<Integer>) Arrays.asList(5, 10, 15, 20);
		for (int t : c_num) {
			double avg_density = 0.0D;
			double avg_clusterSize = 0.0D;
			double avg_euclidean = 0.0D;

			double temp_density = 0.0D;
			double temp_clusterSize = 0.0D;
			double temp_euclidean_dist = 0.0D;
			int TOPK = t;
			int top = TOPK;

			for (SubSpace cluster : clusterRes) {
				if (top < 1) {
					break;
				}
				top--;
				Set<Integer> topRes = new HashSet<Integer>();
				for (int i : cluster.x) {
					topRes.add(i);
				}

				double[] xx = indicateVector(topRes, N); // top 1 results
															// indicator
															// vector

				double density = (1.0 * Matrix.dot(
						Matrix.VecMultiplyMat(xx, adjMatrix), xx))
						/ StatUtils.sum(xx);

				if (cluster.x.size() < 2) {
					TOPK--;
					continue;
				}
				temp_density += density;

				temp_clusterSize += cluster.x.size();

				// System.out.println(cluster.x.size() + " " +
				// cluster.x.toString());
				Set<ArrayList<Integer>> pairs = new HashSet<ArrayList<Integer>>();
				for (int i = 0; i < cluster.x.size() - 1; i++) {
					for (int j = i + 1; j < cluster.x.size(); j++) {
						ArrayList<Integer> temp = new ArrayList<Integer>();
						temp.add(cluster.x.get(i));
						temp.add(cluster.x.get(j));
						pairs.add(temp);
					}
				}
				double tempSum = 0.0D;
				for (ArrayList<Integer> pair : pairs) {
					double dist = EuclideanDist(mgdGraph.W[pair.get(0)],
							mgdGraph.W[pair.get(1)], cluster.y);
					tempSum += dist;
				}
				// System.out.println(tempSum + " / " + pairs.size() + " = " +
				// tempSum
				// / pairs.size());
				temp_euclidean_dist += 1.0 * (tempSum / pairs.size());

			}// for subspace

			// 1. average density
			avg_density = temp_density / TOPK;
			// 2. average cluster size
			avg_clusterSize = temp_clusterSize / TOPK;
			// 3. average euclidean distance
			avg_euclidean = temp_euclidean_dist / TOPK;
			scores += round(avg_density, 2) + " " + round(avg_clusterSize, 2)
					+ " " + round(avg_euclidean, 2);
			log("Top" + t + "/" + TOPK + " " + round(avg_density, 2) + " "
					+ round(avg_clusterSize, 2) + " " + round(avg_euclidean, 2));

		}
		log(outResult + scores + "\n", dataset);
	}

	private static double EuclideanDist(double[] ds1, double[] ds2,
			ArrayList<Integer> y) {
		// TODO Auto-generated method stub
		double eDist = 0.0D;
		for (int i : y) {

			eDist += Math.pow(Math.abs(ds1[i] - ds2[i]), 2);
		}
		return Math.sqrt(eDist);
	}

	public static double[] indicateVector(Set<Integer> cluster, int n) {
		double[] x = new double[n];
		for (int i : cluster) {
			x[i] = 1.0;
		}

		return x;
	}

	public static void log(String line) {
		FileWriter fstream;
		System.out.println(line);
		try {
			fstream = new FileWriter("data/DenseGraph/Log/log_DFB.txt", true);
			fstream.write(line + "\r\n");
			fstream.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void log(String line, String dataset) {
		FileWriter fstream;
		// System.out.println(line);
		try {
			fstream = new FileWriter(
					"data/DenseGraph/Log/log_SG-Pursuit_output-" + dataset
							+ ".txt", true);
			fstream.write(line);
			fstream.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void experiments_RealdataSets() {

		double lambda = 5.00D;
		String root = "data/DenseGraph/Realdata/apdm/";
		String outroot = "outputs/GAMerResult/";
		// "dblp", "genes", "imdb","arxivSmall", "patents", "arxivLarge", "dfb"
		for (String dataset : new String[] { "dblp" }) {
			dataset += "_APDM.txt";
			long startTimeAll = System.nanoTime();
			System.out.println(dataset);
			MultiClusterDetectionOn3Measure(root + dataset, "", 5, 3, lambda,
					5, dataset);

			System.out.println("Total running time: "
					+ (System.nanoTime() - startTimeAll) / 1e9);
		}
	}

	/*
	 *  DBLP         (774,20,1757)  
		DFB          (100,5,1106)          
		GENES        (2900,115,8264)    
		IMDB         (862,21,4388)        
		Patents      (100000,5,188631)                               
		ArxivSmall   (856,30,2660)                          
		ArxivLArge   (11989,300,119258)

	 * **/

	public static void experiments_RealdataSets3Measure() {

		ExecutorService pool = Executors.newFixedThreadPool(1);

		Map<String, Integer> attnumber = new HashMap<String, Integer>() {
			{
				put("dblp", 20);
				put("dfb", 5);
				put("genes", 115);
				put("imdb", 21);
				put("arxivLarge", 300);
			}
		};

		String root = "data/DenseGraph/Realdata/apdm/";
		String outroot = "outputs/GAMerResult/";
		int maxCluster = 20;

		// "dblp", "genes", "imdb", "dfb","arxivSmall", "patents", "arxivLarge"
		for (String dataset : new String[] { "dfb", "dblp", "genes", "imdb",
				"arxivLarge" }) {
			for (double lambda : new double[] { 5.0 }) {
				for (int s : new int[] { 3, 5, 10, 15, 20 }) {
					for (int k : new int[] { 5, 10, 15, 30 }) {
						if (attnumber.get(dataset) < s) {
							System.out.println(attnumber.get(dataset)
									+ " ************ "
									+ s);
							log(dataset + " " + k + " " + s + " -2 -2 -2 \n",
									dataset); // Empty results
							continue;
						}
						
						log("k=" + k + " s=" + s + " lambda=" + lambda);


						String inputFileString = root + dataset
								+ "_norm_APDM.txt";



						long startTimeAll = System.nanoTime();

							pool.execute(new Thread() {
								public void run() {
									MultiClusterDetectionOn3Measure(
											inputFileString, "", k, s, lambda,
										maxCluster, dataset);
								}
							});

						log("Total running time: "
								+ (System.nanoTime() - startTimeAll) / 1e9);
					}// k
				}// s
			}// lambda
		}// dataset

	}
	public static double mad(double[] autoCorrelationValues) {
		double[] tempTable = new double[autoCorrelationValues.length];
		Median m = new Median();
		double medianValue = m.evaluate(autoCorrelationValues);
		for (int i = 0; i < autoCorrelationValues.length; i++) {
			tempTable[i] = Math.abs(autoCorrelationValues[i] - medianValue);
		}
		return m.evaluate(tempTable);
	}

	public static void validateAPDMFiles() {
		String inRoot = "data/DenseGraph/DenseSubgraph_APDM/";
		File folder = new File(inRoot);
		File[] listOfFiles = folder.listFiles();
		for (int i = 0; i < listOfFiles.length; i++) {
			if (listOfFiles[i].isFile()) {
				// System.out.println("File " + listOfFiles[i].getName());
				APDMInputFormat apdm = null;
				MultiDataGraph mgdGraph = null;
				String fileName = listOfFiles[i].getName();
				try {
					apdm = new APDMInputFormat(inRoot + fileName);
					mgdGraph = new MultiDataGraph(apdm);
				} catch (Exception ex) {
					System.out
							.println("!!!!!!!!!!!!!! apdm file reading error..."
									+ fileName);
					File f = new File(inRoot + fileName);
					f.delete();
					System.out.println("!!!!!!!!!!!!!! apdm file deleted..."
							+ fileName);
					// ex.printStackTrace();
					// return;
				}

			} else if (listOfFiles[i].isDirectory()) {
				// System.out.println("Directory " + listOfFiles[i].getName());
			}
		}
	}

	public static double round(double value, int places) {
		if (places < 0)
			throw new IllegalArgumentException();

		long factor = (long) Math.pow(10, places);
		value = value * factor;
		long tmp = Math.round(value);
		return (double) tmp / factor;
	}

	public static class SubSpace {
		public ArrayList<Integer> x = new ArrayList<Integer>();
		public ArrayList<Integer> y = new ArrayList<Integer>();

		public SubSpace(Set<Integer> x_star, Set<Integer> y_star) {
			for (Integer xi : x_star) {
				x.add(xi);
			}
			for (Integer yi : y_star) {
				y.add(yi);
			}
		}
	}

	public static void main(String[] args) {

		experiments_RealdataSets3Measure();

	}
}
