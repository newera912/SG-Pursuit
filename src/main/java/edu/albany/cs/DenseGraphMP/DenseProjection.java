package edu.albany.cs.DenseGraphMP;

import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.Matrix;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.scoreFuncs.DensityProjectScore;
import edu.albany.cs.scoreFuncs.Stat;

import org.apache.commons.math3.distribution.BinomialDistribution;

import java.util.*;

/**
 * Algorithm 1 : S2Graph-Mp Aglorithm .
 *
 * @author
 */
public class DenseProjection {

	/**
	 * graph info
	 */
	private final int n;
	private final int p;
	private final Graph graph;
	/**
	 * the total sparsity of x
	 */
	private final int k;
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
	private final int t;
	/**
	 * Function to be optimized.
	 */
	// private final FunctionX func;

	/** Weight verctor */
	private double[] lambdas;
	private double[][] A;
	private double[] b;
	private double gamma;

	/** Possible subsets */
	private ArrayList<Set<Integer>> S;

	/**
	 * results
	 */
	public double[] x;
	public Set<Integer> bestS;
	public double funcValue = -1.0D;
	public double runTime;

	private int verboseLevel = 1;

	// public DenseGraphMP(Graph graph, int k, int s, int g, Function func) {
	// this(graph, k, s, g, k - g + 0.0D, func);
	// }
	//
	// public DenseGraphMP(Graph graph, int k, int s, int g, double B, Function
	// func) {
	// this(graph, k, s, g, B, 5, func);
	// }

	public DenseProjection(Graph graph, int k, int g, double B, int t,
			double[] lambdas, double[] b, double gamma) {
		this.graph = graph;
		this.n = graph.numOfNodes;
		this.p = graph.numOfFeas;
		this.k = k;
		this.g = g;
		this.B = B;
		this.t = t;
		this.lambdas = lambdas;
		this.b = b;
		this.gamma = gamma;
		this.bestS = new HashSet<Integer>();
		this.S = new ArrayList<Set<Integer>>();		
		//A = graph.AdjMatrix;
		A=new double[n][n];
		double p_in=0.35;
		double p_out=0.05;			
		for(int p=0;p<1000;p++){
			for(int j=p+1;j<1000;j++){
				if((p>=400 && p<500)){				
					
						if((j>=400 && j<500)){
							double succ=new BinomialDistribution(1,p_in).sample();
							A[p][j]=succ;
							A[j][p]=succ;
						}else{
							double succ=new BinomialDistribution(1,p_out).sample();
							A[p][j]=succ;
							A[j][p]=succ;
						}
					
				}else{
					double succ=new BinomialDistribution(1,p_out).sample();
						A[p][j]=succ;
						A[j][p]=succ;
				}
			}
		}
		ConnectedComponents cc = new ConnectedComponents(A);
	    System.out.println("is connected: " + cc.checkConnectivity());

		//Utils.stop();
		// run the algorithm
		if (run()) {
			System.out.println("S2GraphMPDense successfully finished.");
		} else {
			System.out.println("S2GraphMPDense is failed.");
		}
	}

	private boolean run() {

		long startTime = System.nanoTime();

		List<Double> fValues = new ArrayList<>();
		for (double lambda : lambdas) {
			long lamStartTime = System.nanoTime();
			x = Stat.getRandomVector(n, k);
			if (verboseLevel > 0) {
				System.out.println("\n\nk=" + k + " B=" + B
						+ "------------lambda: " + lambda + "------------");
			}

			/** line5 argmax x^T*b.b +lambda_i*(x^TAx/1^Tx) */
			DensityProjectScore dpScore = new DensityProjectScore(b, lambda, A);
			//DensityScore dpScore = new DensityScore(A);
			// ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, double[]
			// c, int s, int g, double B, int t,int[] trueSubGraph,
			// DensityProjectScore func
			GraphMP graphMP = new GraphMP(graph, x, k, g, B, t, dpScore, false);
			/** line6 */
			if (graphMP.resultNodes_supportX != null) {
				S.add(Stat.Array2Set(graphMP.resultNodes_supportX));
			}
			System.out.println("Lambda=" + lambda + " " + (System.nanoTime() - lamStartTime) / 1e9 + " sec ...");
			funcValue = dpScore.getFuncValue(graphMP.x);
		}// lambdas

		/** line 9: j = */
		bestS = argmin();
		System.out.println("BestS: " + bestS.toString());
		runTime = (System.nanoTime() - startTime) / 1e9;
		return true;
	}

	public Set<Integer> argmin() {
		Set<Integer> S_j = new HashSet<Integer>();
		double best_bs = Double.NEGATIVE_INFINITY;
		for (Set<Integer> S_i : S) {
			System.out.println("Density: " + density(S_i) + " size:"+ S_i.size() + " S=" + S_i.toString());
			double temp_bs = 0.0D;
			double density_Si = density(S_i);
			if (density_Si > gamma) {
				/** ||b_s||^2_2 */
				for (int i : S_i) {
					temp_bs += b[i] * b[i];
				}
				System.out.println(">>> temp_bs:" + density_Si + " best_bs:"
						+ best_bs);
				// if(temp_bs>best_bs){
				// System.out.println("temp_bs "+temp_bs+" best_bs"+best_bs);
				// best_bs=temp_bs;
				// S_j=S_i;
				// }
				if (density_Si > best_bs) {
					//System.out.println("temp_bs " + density_Si + " best_bs"+ best_bs);
					best_bs = density_Si;
					S_j = S_i;
				}
			}

		}

		return S_j;
	}

	public double density(Set<Integer> S_i) {
		double density = 0.0;
		double[] xS = Stat.getIndicateVec(S_i, n);
		//System.out.println("xAx" + Matrix.dot(Matrix.VecMultiplyMat(xS, A), xS)+ " S_i size:" + S_i.size());
		density = Matrix.dot(Matrix.VecMultiplyMat(xS, A), xS) / S_i.size();

		return density;

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
			GammaY.add(indexes[i]);
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

}
