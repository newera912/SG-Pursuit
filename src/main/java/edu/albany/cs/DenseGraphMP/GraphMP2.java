package edu.albany.cs.DenseGraphMP;

import edu.albany.cs.DenseGraphMP.FunctionX;
import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.MultiDataGraph;
import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.scoreFuncs.DensityProjectScore;
import edu.albany.cs.scoreFuncs.DensityScore;
import edu.albany.cs.scoreFuncs.Stat;
import edu.albany.cs.tailApprox.TailApprox;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Algorithm 1 : Graph-Mp Aglorithm in our IJCAI paper.
 *
 * @author Baojian bzhou6@albany.edu
 */
public class GraphMP2 {

	/** 1Xn dimension, input data */
	private final double[] c;
	private final int numOfNodes;
	/** graph info */
	private final HashSet<Integer> nodes;
	private final ArrayList<Integer[]> edges;
	private final ArrayList<Double> edgeCosts;
	/** the total sparsity of S */
	private final int s;
	/** the maximum number of connected components formed by F */
	private final int g;
	/** bound on the total weight w(F) of edges in the forest F */
	private final double B;
	/** number of iterations */
	private final int t;
	
	//private final DensityProjectScore function;
	private final DensityProjectScore function;
//	/** project vector b , adjacency Matix A */
//	private final double[] b;
//	private final double[] A;
	/** results */
	public double[] x;
	public int[] resultNodes_supportX;
	public int[] resultNodes_Tail = null;
	public double funcValue = -1.0D;
	public double runTime;
	public Graph graph;
	private int verboseLevel = 1;

	

	public GraphMP2(Graph graph, double[] c, int s, int g, double B, int t, DensityProjectScore func) {
		this.edges = graph.edges;
		this.nodes = new HashSet<Integer>(graph.nodes);
		this.numOfNodes = nodes.size();		
		this.edgeCosts = graph.edgeCosts;		
		this.c = c;
		this.s = s;
		this.g = g;
		this.B = B;
		this.t = t;
		
		this.function = func;
		/** run Graph-MP algorithm. */
		x = run();
	}
	
	private double[] run() {

		long startTime = System.nanoTime();
		double[] x = Stat.getRandomVector(numOfNodes, s);
	
		System.out.println("initial X: "+ArrayUtils.toString(Stat.supp(x)));
		ArrayList<Double> fValues = new ArrayList<>();
		for (int i = 0; i < this.t; i++) { // t iterations
			if (verboseLevel > 0) {
				System.out.println("s="+s+" B="+B+" GraphMP------------iteration: " + i + "------------");
				System.out.println("FunVal: "+function.getFuncValue(x));
			}
			fValues.add(function.getFuncValue(x));
			double[] gradientF = function.getGradient(x);
			System.out.println("Head Input gradientF="+Stat.printWithIndex(gradientF));
			gradientF = Stat.normalizeGradient(x, gradientF, numOfNodes);
			System.out.println("Norm headInput gradientF="+Stat.printWithIndex(gradientF));
			/** head approximation */
			HeadApprox head = new HeadApprox(edges, edgeCosts, gradientF, s, g, B,false);
			ArrayList<Integer> S = Utils.unionSets(head.bestForest.nodesInF, support(x));
			
			long argminStartTime = System.nanoTime();
			double[] b = function.getArgMinFx(S);
			System.out.println("argmin(S) S: "+ArrayUtils.toString(S));
			System.out.println("Tail Input b "+Stat.printWithIndex(b));
			System.out.println("Time for argmin= "+(System.nanoTime() - argminStartTime) / 1e9+" sec ...");
			/** tail approximation */
			TailApprox tail = new TailApprox(edges,edgeCosts, b, s, g, B,false);
			/** calculate x^{i+1} */
			for (int j = 0; j < b.length; j++) {
				x[j] = 0.0D;
			}
			for (int j : tail.bestForest.nodesInF) {
				x[j] = b[j];
			}
			if (verboseLevel > 0) {
				System.out.println("number of head nodes : " + head.bestForest.nodesInF.size()+" "+ArrayUtils.toString(head.bestForest.nodesInF));
				System.out.println("Tail Input supp(b)   : " + Stat.supp(b).size()+" "+Stat.supp(b));
				System.out.println("number of tail nodes : " + tail.bestForest.nodesInF.size()+" "+ArrayUtils.toString(tail.bestForest.nodesInF));
			}
			resultNodes_Tail = Utils.getIntArrayFromIntegerList(tail.bestForest.nodesInF);
			//System.out.println("X_i+1:" +Stat.supp(x));
			//Utils.stop();
		}
		resultNodes_supportX = getSupportNodes(x);
		funcValue = function.getFuncValue(x);
		runTime = (System.nanoTime() - startTime) / 1e9;
		System.out.println("Projection 1 iteration time: "+runTime+" sec...");
		return x;
	}




	private double[] initializeRandom() {
		double[] x0 = new double[c.length];
		Random rand = new Random();
		while(StatUtils.sum(x0)==0.0D){
		for (int i = 0; i < c.length; i++) {
			if (rand.nextDouble() < 0.5D) {
				x0[i] = 1.0D;
			} else {
				x0[i] = 0.0D;
			}
		}
		}
		return x0;
	}

	/**
	 * get a support of a vector
	 *
	 * @param x
	 *            array x
	 * @return a subset of nodes corresponding the index of vector x with
	 *         entries not equal to zero
	 */
	public ArrayList<Integer> support(double[] x) {
		if (x == null) {
			return null;
		}
		ArrayList<Integer> nodes = new ArrayList<>();
		for (int i = 0; i < x.length; i++) {
			if (x[i] != 0.0D) {
				nodes.add(i);
			}
		}
		return nodes;
	}

	/**
	 * get the nodes returned by algorithm
	 *
	 * @param x
	 *            array x
	 * @return the result nodes
	 */
	private int[] getSupportNodes(double[] x) {
		int[] nodes = null;
		for (int i = 0; i < x.length; i++) {
			if (x[i] != 0.0D) {
				/** get nonzero nodes */
				nodes = ArrayUtils.add(nodes, i);
			}
		}
		return nodes;
	}
	
	
	public static void UnitTestGraphMP(String inputFilePath) {

		System.out.println("\n------------------------------ test starts ------------------------------");
		System.out.println("testing file path: " + inputFilePath);
		/** step0: data file */
		APDMInputFormat apdm = new APDMInputFormat(inputFilePath);
		Graph graph=new MultiDataGraph(apdm);
		
		ArrayList<Integer[]> edges = apdm.data.intEdges;
		ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;
		Set<Integer> trueSet=new HashSet<Integer>();
		double[] b=new double[apdm.data.numNodes];
		double lambda=15;
		double mu=10.0D;
		int[] candidateS = new int[] { 35,40,45};
		
		
		 for(int q=0;q<b.length;q++){
	         	if(q>=400 && q<500){  
	         		trueSet.add(q);
	         		b[q]=new NormalDistribution(mu, 1.0D).sample();					
	         	}else{
	         		b[q]=new NormalDistribution(0, 1.0D).sample();
	         	}
	         }
		 
		/** step1: score function */
		DensityProjectScore func = new DensityProjectScore(b, lambda,graph.AdjMatrix);
		/** step2: optimization */
		
		double optimalVal = Double.MAX_VALUE;
		int g=1;
		PreRec bestPreRec = new PreRec();
		GraphMP2 bestGraphMP = null;
		for (int s : candidateS) {
			double B = s - g + 0.0D;
			int t = 3;
//			Graph graph, double[] c, int s, int g, double B, int t, DensityProjectScore func
			
			GraphMP2 graphMP = new GraphMP2(graph, b, s, g, B, t, func);
			double[] yx = graphMP.x;
			if (func.getFuncValue(yx) < optimalVal) {
				optimalVal = func.getFuncValue(yx);
				bestPreRec = new PreRec(graphMP.resultNodes_supportX, trueSet);
				bestGraphMP = graphMP;
				
			}
		}
		System.out.println("\n\nprecision : " + bestPreRec.pre + " ; recall : " + bestPreRec.rec);
		System.out.println("result subgraph is: " + Arrays.toString(bestGraphMP.resultNodes_Tail));
		System.out.println("true subgraph is: " + Arrays.toString(trueSet.toArray()));
		System.out.println("------------------------------ test ends --------------------------------\n");
	}
	
	
	public static void main(String[] args){
		String inputfile="data/DenseGraph/Dense_APDM/APDM_r_5_Cluster_10_in_0.35_out_0.05_case_0.txt";
		UnitTestGraphMP(inputfile);
	}
}
