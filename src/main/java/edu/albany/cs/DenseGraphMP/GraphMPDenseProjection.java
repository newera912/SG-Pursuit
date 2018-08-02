package edu.albany.cs.DenseGraphMP;
import java.util.Date;
import java.text.SimpleDateFormat;
import java.text.DateFormat;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.MultiDataGraph;
import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.scoreFuncs.DensityProjectScore;
import edu.albany.cs.scoreFuncs.Stat;
import edu.albany.cs.tailApprox.TailApprox;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.StatUtils;

import java.util.*;
import java.util.Map.Entry;

/**
 * Algorithm 1 : Graph-Mp Aglorithm in our IJCAI paper.
 *
 * @author Baojian bzhou6@albany.edu
 */
public class GraphMPDenseProjection {

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
	private final boolean debug;
	private final int[] ranks;
	
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

	
	//graph.edges, graph.edgeCosts, x, k, g, B, t, new int[0], dpScore
//	public GraphMP(Graph graph, double[] c, int s, int g, double B, int t, DensityProjectScore func) {
//		this.edges = graph.edges;
//		this.nodes = new HashSet<Integer>(graph.nodes);
//		this.numOfNodes = nodes.size();		
//		this.edgeCosts = graph.edgeCosts;		
//		this.c = c;
//		this.s = s;
//		this.g = g;
//		this.B = B;
//		this.t = t;
//		
//		//this.function = func;
//		this.function = null;
//		/** run Graph-MP algorithm. */
//		x = run();
//	}
	
	public GraphMPDenseProjection(Graph graph, double[] c, int s, int g, double B, int t, DensityProjectScore func, int[] ranks, boolean debug) {
		this.edges = graph.edges;
		this.nodes = new HashSet<Integer>(graph.nodes);
		this.numOfNodes = nodes.size();		
		this.edgeCosts = graph.edgeCosts;		
		this.c = c;
		this.s = s;
		this.g = g;
		this.B = B;
		this.t = t;
		this.debug = debug;
		this.ranks = ranks;
		
		this.function = func;
		/** run Graph-MP algorithm. */
		x = run();
	}
	public static void log(String line){
		FileWriter fstream;
		System.out.println(line);		
		try {
			fstream = new FileWriter("data/DenseGraph/Log/log.txt", true);
			fstream.write(line);
			fstream.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	private double[] run() {
		boolean debug = false;
		long startTime = System.nanoTime();
//		double[] x = Stat.getRandomVector(numOfNodes, s);
		double[] x = new double[numOfNodes];
		Arrays.fill(x, 0);
		for(int i=0; i< 2 * s; i++){
			x[ranks[i]] = 1; 
		}
	
		if (debug == true){
			System.out.println("initial X: "+ArrayUtils.toString(Stat.supp(x)));
		}
		ArrayList<Double> fValues = new ArrayList<>();
		for (int i = 0; i < this.t; i++) { // t iterations
			if (debug == true){
				if (verboseLevel > 0) {
					System.out.println("s="+s+" B="+B+" GraphMP------------iteration: " + i + "------------");
					System.out.println("FunVal: "+function.getFuncValue(x));
				}
			}
			fValues.add(function.getFuncValue(x));
			double[] gradientF = function.getGradient(x);
			if (debug == true){
				System.out.println("Head Input gradientF="+Stat.printWithIndex(gradientF));
			}
			gradientF = Stat.normalizeGradient(x, gradientF, numOfNodes);
			if (debug == true){
				System.out.println("Norm headInput gradientF="+Stat.printWithIndex(gradientF));
			}
			/** head approximation */
			HeadApprox head = new HeadApprox(edges, edgeCosts, gradientF, s, g, B, false);
			ArrayList<Integer> S = Utils.unionSets(head.bestForest.nodesInF, support(x));
			
			long argminStartTime = System.nanoTime();
			double[] b = function.getArgMinFx(S,false);

			if (debug == true){
				System.out.println("argmin(S) S: "+ArrayUtils.toString(S));
				System.out.println("Tail Input b "+Stat.printWithIndex(b));
				System.out.println("Time for argmin= "+(System.nanoTime() - argminStartTime) / 1e9+" sec ...");
			}
			/** tail approximation */
			TailApprox tail = new TailApprox(edges,edgeCosts, b, s, g, B, false);
			/** calculate x^{i+1} */
			for (int j = 0; j < b.length; j++) {
				x[j] = 0.0D;
			}
			for (int j : tail.bestForest.nodesInF) {
				x[j] = b[j];
			}
			if (verboseLevel > 0 && debug == true) {
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
		if (debug == true){
			System.out.println("Projection 1 iteration time: "+runTime+" sec...");
		}
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

	public static int[] getIndicesInOrder(double[] array) {
	    Map<Integer, Double> map = new HashMap<Integer, Double>(array.length);
	    for (int i = 0; i < array.length; i++)
	        map.put(i, array[i]);

	    List<Entry<Integer, Double>> l = 
	                           new ArrayList<Entry<Integer, Double>>(map.entrySet());

	    Collections.sort(l, new Comparator<Entry<?, Double>>() {
	            @Override
	            public int compare(Entry<?, Double> e1, Entry<?, Double> e2) {
	                return e2.getValue().compareTo(e1.getValue());
	            }
	        });

	    int[] result = new int[array.length];
	    for (int i = 0; i < result.length; i++)
	        result[i] = l.get(i).getKey();

	    return result;
	}	
	public static void UnitTestGraphMP1(String inputFilePath) {
		APDMInputFormat apdm = new APDMInputFormat(inputFilePath);
		Graph graph=new MultiDataGraph(apdm);
		double[] weights = new double[graph.AdjMatrix.length];
		for(int i = 0; i < graph.AdjMatrix.length; i++){
			for(int j = 0; j < graph.AdjMatrix[i].length; j++){
				weights[i] = weights[i] + graph.AdjMatrix[i][j];
			}
		}
		int[] ranks = getIndicesInOrder(weights);
//		for(int i=0; i < graph.numOfNodes; i++){
//			System.out.println(weights[ranks[i]]);
//		}
//		Arrays.sort(weights);
//		ArrayUtils.reverse(weights);
		
//		Set<Integer> trueSet=new HashSet<Integer>();
		double[] b = new double[apdm.data.numNodes];
		double lambda=0.0000;
		int[] candidateS = new int[] {35,40,45, 50, 55};
		for(int q=0;q<b.length;q++){
			b[q]=new NormalDistribution(Math.random(), 1.0D).sample();
		}
		for(int q : graph.trueSubGraphNodes){
			b[q]=new NormalDistribution(0.6, 0.01D).sample();
		}
		
		/** step1: score function */
		DensityProjectScore func = new DensityProjectScore(b, lambda,graph.AdjMatrix);
		/** step2: optimization */
		
		double optimalVal = Double.MAX_VALUE;
		int g=1;
		PreRec bestPreRec = new PreRec();
		GraphMPDenseProjection bestGraphMP = null;
		for (int s : candidateS) {
			double B = s - g + 0.0D;
			int t = 3;
			GraphMPDenseProjection graphMP = new GraphMPDenseProjection(graph, b, s, g, B, t, func, ranks, false);
			double[] yx = graphMP.x;
			PreRec prf=new PreRec(graphMP.resultNodes_supportX, graph.trueSubGraphNodes);
			log("\n Sparsity:" + s + " FValue:"+func.getFuncValue(yx)+" pre="+prf.pre+" rec="+prf.rec);
			if (func.getFuncValue(yx) < optimalVal) {
				optimalVal = func.getFuncValue(yx);
				bestPreRec = new PreRec(graphMP.resultNodes_supportX, graph.trueSubGraphNodes);
				bestGraphMP = graphMP;
			}
		}
		log("\n\nFValue:"+optimalVal+" Best precision : " + bestPreRec.pre + " ; Best recall : " + bestPreRec.rec);
		log("\nresult subgraph is: " + Arrays.toString(bestGraphMP.resultNodes_Tail));
		log("\ntrue subgraph is: " + Arrays.toString(graph.trueSubGraphNodes));
		log("\n------------------------------ test ends --------------------------------\n");
	}
	
	
//	public static void UnitTestGraphMP(String inputFilePath) {
//
//		System.out.println("\n------------------------------ Unit test 1 starts ------------------------------");
//		System.out.println("testing file path: " + inputFilePath);
//		/** step0: data file */
//		APDMInputFormat apdm = new APDMInputFormat(inputFilePath);
//		Graph graph=new MultiDataGraph(apdm);
//		
////		graph.
//		
//		ArrayList<Integer[]> edges = apdm.data.intEdges;
//		ArrayList<Double> edgeCosts = apdm.data.identityEdgeCosts;
//		Set<Integer> trueSet=new HashSet<Integer>();
//		double[] b=new double[apdm.data.numNodes];
//		double lambda=0.0000;
//		double mu=10.0D;
//		int[] candidateS = new int[] { 35,40,45};
//
//		 for(int q=0;q<b.length;q++){
//	         	if(q>=400 && q<500){  
//	         		trueSet.add(q);
//	         		b[q]=new NormalDistribution(mu, 0.01D).sample();					
//	         	}else{
//	         		b[q]=new NormalDistribution(0, 1.0D).sample();
//	         	}
//	         }
//		 
//		/** step1: score function */
//		DensityProjectScore func = new DensityProjectScore(b, lambda,graph.AdjMatrix);
//		/** step2: optimization */
//		
//		double optimalVal = Double.MAX_VALUE;
//		int g=1;
//		PreRec bestPreRec = new PreRec();
//		GraphMPDenseProjection bestGraphMP = null;
//		for (int s : candidateS) {
//			double B = s - g + 0.0D;
//			int t = 3;
////			Graph graph, double[] c, int s, int g, double B, int t, DensityProjectScore func
//			
//			GraphMPDenseProjection graphMP = new GraphMPDenseProjection(graph, b, s, g, B, t, func, false);
//			double[] yx = graphMP.x;
//			PreRec prf=new PreRec(graphMP.resultNodes_supportX, trueSet);
//			System.out.println("S="+s+" FValue:"+func.getFuncValue(yx)+" pre="+prf.pre+" rec="+prf.rec);
//			if (func.getFuncValue(yx) < optimalVal) {
//				optimalVal = func.getFuncValue(yx);
//				bestPreRec = new PreRec(graphMP.resultNodes_supportX, trueSet);
//				bestGraphMP = graphMP;
//				
//			}
//		}
//		System.out.println("\n\nFValue:"+optimalVal+" Best precision : " + bestPreRec.pre + " ; Best recall : " + bestPreRec.rec);
//		System.out.println("\n\nresult subgraph is: " + Arrays.toString(bestGraphMP.resultNodes_Tail));
//		System.out.println("\n\ntrue subgraph is: " + Arrays.toString(trueSet.toArray()));
//		System.out.println("\n\n------------------------------ test ends --------------------------------\n");
//	}
//	
	
	public static void main(String[] args){
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date = new Date();
		log("\n\n========================" + dateFormat.format(date) + "=============================\n\n");
		for(int i=0; i < 10; i++){
			String inputfile="data/DenseGraph/DenseSubgraph_APDM/APDM_Dense_subgraph_in_0.35_out_0.1_trueFeasNum_5_case_" + i + ".txt";
			log("\n\n\n\n" + inputfile);
			UnitTestGraphMP1(inputfile);
		}
	}
}
