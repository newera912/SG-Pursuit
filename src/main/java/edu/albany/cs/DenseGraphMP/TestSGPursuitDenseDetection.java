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
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Date;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.lang3.ArrayUtils;
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
public class TestSGPursuitDenseDetection {

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

//    public S2GraphMPDense(Graph graph, int k, int s, int g, Function func) {
//        this(graph, k, s, g, k - g + 0.0D, func);
//    }
//
//    public S2GraphMPDense(Graph graph, int k, int s, int g, double B, Function func) {
//        this(graph, k, s, g, B, 5, func);
//    }

    /*
     * k: sparsity of nodes
     * s: sparsity of attributes
     * g: number of connected components
     * B: budget
     * iterations: The maximum number of iterations. 
     */
    public TestSGPursuitDenseDetection(String fileName, MultiDataGraph graph, int k, int s, int g, double B, int iterations, Function func,double lambda) {
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
        this.lambda=lambda;
		System.out.println(n + " " + p);
        if (run()) {
//            System.out.println("S2GraphMPDense successfully finished.");
        } else {
            System.out.println("S2GraphMPDense is failed.");
        }
    }
    
    private int[] ranknodes(ArrayList<Integer> nodes, int k){
    	double[] degrees = new double[nodes.size()];
    	int[] rank = new int[Math.min(nodes.size(), k)];
    	Arrays.fill(degrees, 0);
    	for(int i=0; i<nodes.size(); i++){
    		for(int j:nodes){
    			if(graph.AdjMatrix[nodes.get(i)][j] > 0){
    				degrees[i] += 1;
    			}
    		}
    	}
        ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(degrees);
        Integer[] indexes = arrayIndexComparator.getIndices();
        Arrays.sort(indexes, arrayIndexComparator);
        for(int i=0; i<rank.length; i++){
        	rank[i] = nodes.get(indexes[i]);
        }
        return rank;
    }

    /*
     * This function identifies good initial values for x and y.  
     * k: sparsity of x (nodes)
     * s: sparsity of y (attributes)
     * OUTPUT: return the initial values of x and y. 
     */
    private List<double[]> CalcInitilvalues(int k, int s){
        double[] x0 = new double[n];
        double[] y0 = new double[p];
        Arrays.fill(x0, 0);
        Arrays.fill(y0, 0);
        double threshold = 0.11;
        ArrayList<ArrayList<Integer>> res = new ArrayList<ArrayList<Integer>>();
        for (int j = 0; j < p; j++){
        	for(int i = 0; i < n; i++){ 
	        	ArrayList<Integer> nns = graph.arrayListAdj.get(i);
	        	Collections.sort(nns);
	        	double[] vals = new double[nns.size()]; 
	        	for(int k1 = 0; k1 < nns.size(); k1++){
	        		vals[k1] = graph.W[nns.get(k1)][j];
	        	}
	        	if(nns.size() > 6){
		        	x0[i] = Math.abs(graph.W[i][j] - StatUtils.percentile(vals, 50));
	        	}else{
		        	x0[i] = Math.abs(graph.W[i][j] - StatUtils.mean(vals));
	        	}
        	}
            ArrayList<Integer> S = new ArrayList<Integer>();
            for(int i = 0; i < n; i++){
            	if(x0[i] < threshold){
                	S.add(i);
            	}
            }
//            if(S.size() > 0){
//                S = Utils.getLargestConnectedComponant(graph.arrayListAdj, S);
//            }else{
//            	S = new ArrayList<Integer>();
//            }
            res.add(S);
        } 
        Integer[] indexes = Utils.sortArrayList1(res);
        for(int i = 0; i < s; i++){
            y0[indexes[i]] = 1;
        }        
        Arrays.fill(x0, 0);
        int[] rank = ranknodes(res.get(indexes[0]), k);
        for(int i=0; i<k; i++){
        	x0[rank[i]] = 1;
        }
        
//        String line = "features: ";
//        double[] xx = new double[n];
//        Arrays.fill(xx, 0);
//        for(int i = 0; i < s; i++){
//        	line += indexes[i] + ",";
//            y0[indexes[i]] = 1;
//            for(int j:res.get(indexes[i])){
//            	xx[j] += 1;
//            }
//        }
////        System.out.println(line);
//        Arrays.fill(x0, 0);
//        line = "nodes: ";
//        int cnt = 0;
//        for(int i=0; i<n; i++){
//        	if(xx[i] == s){
//            	line += i + ",";
//        		x0[i] = 1;
//        		cnt += 1;
//        	}
//        }
//        if(cnt <5){
//            for(int i:res.get(indexes[0])){
//        	    line += i + ",";
//        	    x0[i] = 1;
//            }
//        }
//        System.out.println(line);
        List<double[]> list = new ArrayList<>();
        list.add(x0);
        list.add(y0);
        return list;
    }
    
    private List<double[]> CalcInitilvalues1(int k, int s){
    	double threshold = -0.7;
        double[] x0 = new double[n];
        double[] y0 = new double[p];
        Arrays.fill(x0, 0);
        Arrays.fill(y0, 0);
        ArrayList<int[]> res = new ArrayList<int[]>();
        double[] scores = new double[p];
        for (int j = 0; j < p; j++){
        	ArrayList<Integer> S = new ArrayList<Integer>();
        	for(int i = 0; i < n; i++){ 
	        	ArrayList<Integer> nns = graph.arrayListAdj.get(i);
//	        	Collections.sort(nns);
	        	double[] dists = new double[nns.size()]; 
	        	for(int k1 = 0; k1 < nns.size(); k1++){
	        		dists[k1] = Math.abs(graph.W[nns.get(k1)][j] - graph.W[i][j]);
	        	}
	        	if(nns.size() > 1){
		        	x0[i] = StatUtils.percentile(dists, 50) * -1;
	        	}else{
		        	x0[i] = StatUtils.mean(dists) * -1;
	        	}
	        	if(x0[i] > threshold) S.add(i);
        	}
            ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(x0);
            Integer[] indexes = arrayIndexComparator.getIndices();
            Arrays.sort(indexes, arrayIndexComparator);
            S = new ArrayList<Integer>();
            for(int i=0; i<k * 5; i++){
            	S.add(indexes[i]);
            }
        	if(S.size() > 0){
            	int[] rank = ranknodes(S, k);
            	double[] values = new double[rank.length];
            	for(int i=0; i<rank.length; i++){
            		values[i] = graph.W[rank[i]][j];
            	}
            	double score = StatUtils.populationVariance(values) * -1;
                Arrays.fill(x0, 0);
                for(int i:rank){
                	x0[i] = 1;
                }
                Arrays.fill(y0, 0);
                y0[j] = 1;
                funcValue = function.getFuncValue(x0, y0,lambda);
                if(funcValue >= 0) scores[j] = -1000;
                else scores[j] = score;
                res.add(rank);
        	}else{
        		int[] rank = new int[]{};
        		scores[j] = -1000;
        		res.add(rank);
        	}
        }
        ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(scores);
        Integer[] indexes = arrayIndexComparator.getIndices();
        Arrays.sort(indexes, arrayIndexComparator);
        Arrays.fill(y0, 0);
        for(int i = 0; i < s; i++){
        	if(scores[indexes[i]] > -1000) y0[indexes[i]] = 1;
        }        
        if(StatUtils.sum(y0) == 0) y0[indexes[0]] = 1;
        Arrays.fill(x0, 0);
        for(int i:res.get(indexes[0])){
        	x0[i] = 1;
        }
        if(StatUtils.sum(x0) == 0){
        	System.out.println("NOTICE: no good initialization can be idetnfied!!!!!!!!!!!!!!!!!!!!!!!!!");
        	x0[0] = 1;
        	x0[1] = 1;
        	x0[2] = 1;
        }
        List<double[]> list = new ArrayList<>();
        list.add(x0);
        list.add(y0);
        return list;
    }    

    private List<double[]> CalcInitilvalues2(int k, int s){
        double[] x0 = new double[n];
        double[] y0 = new double[p];
        Arrays.fill(x0, 0);
        Arrays.fill(y0, 0);
        ArrayList<ArrayList<int[]>> res = new ArrayList<ArrayList<int[]>>();
        ArrayList<double[]> scores = new ArrayList<double[]>();
        int[] trials = new int[]{2, 3, 4, 5, 6}; 
        for(int i=0; i<trials.length; i++){
        	res.add(new ArrayList<int[]>());
        	scores.add(new double[p]);
        }
        for (int j = 0; j < p; j++){
        	ArrayList<Integer> S = new ArrayList<Integer>();
        	for(int i = 0; i < n; i++){ 
	        	ArrayList<Integer> nns = graph.arrayListAdj.get(i);
	        	double[] dists = new double[nns.size()]; 
	        	for(int k1 = 0; k1 < nns.size(); k1++){
	        		dists[k1] = Math.abs(graph.W[nns.get(k1)][j] - graph.W[i][j]);
	        	}
	        	if(nns.size() > 1){
		        	x0[i] = StatUtils.percentile(dists, 50) * -1;
	        	}else{
		        	x0[i] = StatUtils.mean(dists) * -1;
	        	}
        	}
            ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(x0);
            Integer[] indexes = arrayIndexComparator.getIndices();
            Arrays.sort(indexes, arrayIndexComparator);

            int ii=0;
            for(int r:trials){
				if (k * r > n) {
					continue;
				}
                S = new ArrayList<Integer>();
				for (int i = 0; i < k * r; i++) {
                	S.add(indexes[i]);
                }
            	if(S.size() > 0){
                	int[] rank = ranknodes(S, k);
                	double fval = getFun(rank, new int[]{j}, 0);
                	scores.get(ii)[j] = fval * -1;
                    res.get(ii).add(rank);
//                	double[] values = new double[rank.length];
//                	for(int i=0; i<rank.length; i++){
//                		values[i] = graph.W[rank[i]][j];
//                	}
//                	double score = 0;
//                	if(rank.length > 0)
//                		score = StatUtils.populationVariance(values) * -1;
//                	else
//                		score = -1000;
//                    scores.get(ii)[j] = score;
//                    res.get(ii).add(rank);
            	}else{
            		int[] rank = new int[]{};
            		scores.get(ii)[j] = -1000;
            		res.get(ii).add(rank);
            	}
            	ii += 1;
            }
        }
        double fValue = -1;
        ArrayList<Integer> Y = new ArrayList<Integer>();
        ArrayList<Integer> X = new ArrayList<Integer>();
        for(int ii=0; ii<trials.length; ii++){
            ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(scores.get(ii));
            Integer[] indexes = arrayIndexComparator.getIndices();
            Arrays.sort(indexes, arrayIndexComparator);
            ArrayList<Integer> Y1 = new ArrayList<Integer>();
            for(int i = 0; i < s; i++){
            	if(scores.get(ii)[indexes[i]] > -1000) Y1.add(indexes[i]);
            }
			System.out.println(ii + "------------------");
            int[] X1 = res.get(ii).get(indexes[0]);
            double fValue1 = getFun(X1, Y1, 0);
            if(fValue == -1 || fValue > fValue1){
            	Y.clear();
            	for(int i:Y1) Y.add(i);
            	X.clear();
            	for(int i:X1) X.add(i);
            	fValue = fValue1;
            }
            for(int j=0; j<s; j++){
            	X1 = res.get(ii).get(indexes[j]);
            	Y1.clear();
            	Y1.add(indexes[j]);
                fValue1 = getFun(X1, Y1, 0);
//            	fValue1 = scores.get(ii)[j];
                if(fValue == -1 || fValue > fValue1){
                	Y.clear();
                	for(int i:Y1) Y.add(i); 
                	X.clear();
                	for(int i:X1) X.add(i); 
                	fValue = fValue1;
                }
            }
            for(int j=0; j<s*2; j++){
				if (j >= p) {
					continue;
				}
            	if(Y.contains(indexes[j]) == false){
            		ArrayList<Integer> YY = new ArrayList<Integer>();
            		for(int i:Y) YY.add(i);
            		YY.add(indexes[j]);
            		fValue1 = getFun(X, YY, 0);
                    if(fValue == -1 || fValue > fValue1){
                    	Y.clear();
                    	for(int i:YY) Y.add(i); 
                    	fValue = fValue1;
                    }
            	}
            }
        }
        Arrays.fill(y0, 0);
        for(int i:Y) y0[i] = 1; 
        Arrays.fill(x0, 0);
        for(int i:X) x0[i] = 1; 
        if(StatUtils.sum(x0) == 0){
        	System.out.println("NOTICE: no good initialization of x can be idetnfied!!!!!!!!!!!!!!!!!!!!!!!!!");
        	x0[0] = 1;
        	x0[1] = 1;
        	x0[2] = 1;
        }
        if(StatUtils.sum(y0) == 0){
        	System.out.println("NOTICE: no good initialization of y can be idetnfied!!!!!!!!!!!!!!!!!!!!!!!!!");
        	y0[0] = 1;
        }
        List<double[]> list = new ArrayList<>();
        list.add(x0);
        list.add(y0);
        return list;
    }        
    
    
    
    public double getFun(int[] X, ArrayList<Integer> Y, double lambda){
    	double[] x = new double[n];
    	double[] y = new double[p];
    	Arrays.fill(x, 0);
    	for(int i:X){
    		x[i] = 1;
    	}
    	Arrays.fill(y, 0);
    	for(int i:Y){
    		y[i] = 1;
    	}
    	return function.getFuncValue(x, y,lambda);
    }
    
    public double getFun(ArrayList<Integer> X, ArrayList<Integer> Y, double lambda){
    	double[] x = new double[n];
    	double[] y = new double[p];
    	Arrays.fill(x, 0);
    	for(int i:X){
    		x[i] = 1;
    	}
    	Arrays.fill(y, 0);
    	for(int i:Y){
    		y[i] = 1;
    	}
    	return function.getFuncValue(x, y,lambda);
    }
    
   public double getFun(int[] X, int[] Y, double lambda){
    	double[] x = new double[n];
    	double[] y = new double[p];
    	Arrays.fill(x, 0);
    	for(int i:X){
    		x[i] = 1;
    	}
    	Arrays.fill(y, 0);
    	for(int i:Y){
    		y[i] = 1;
    	}
    	return function.getFuncValue(x, y,lambda);
    }

    public double sum(double[] x){
    	double sum = 0;
    	for(int i=0; i< x.length; i++){
    		sum += x[i];
    	}
    	return sum;
    }

    public static void showstat(double[] x, String title, double lowerbound, double upperbound){
        String line = "";
        for(int j = 0; j < x.length; j++){
        	if(x[j] != 0 && x[j] > lowerbound && x[j] < upperbound){
        		line += ", (" + j + ", " + x[j] + ")";
        	}
        }
        System.out.println(title + line);
    }
    
    public static void showstat(double[] x, String title){
        String line = "";
        for(int j = 0; j < x.length; j++){
        	if(x[j] != 0){
//        		line += ", (" + j + ", " + x[j] + ")";
        		line += "," + x[j];
        	}
        }
        System.out.println(title + line);
    }
    
    public void showstat(int[] x, String title){
        String line = "";
        for(int j = 0; j < x.length; j++){
        	if(x[j] != 0){
        		line += ", (" + j + ", " + x[j] + ")";
        	}
        }
        System.out.println(title + line);
    }

    private void initilization1(double[] xi, double[] yi, MultiDataGraph graph){
	    int cnt = 0;
	    for(int j:graph.trueSubGraphNodes){
	    	xi[j] = 1;
	    	cnt += 1;
	    	if(cnt > 5)
	    		break; 
	    }
	    cnt = 0;
	    for(int j:graph.trueFeas){
	    	yi[j] = 1;
	    	cnt += 1;
	    	if(cnt > 2)
	    		break; 
	    }
    }
    
    private void StoreTestCaseData(double[] xi, double[] yi, Set<Integer> omegaX, Set<Integer> omegaY, double lambda, String fileName, String outputfilename){
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
        } catch(Exception ex) {
            ex.printStackTrace();
        }
        Path p = Paths.get(fileName);
        String file = p.getFileName().toString();
        try{
        	copyFile(new File(fileName), new File("testcases/getArgMinFxy/" + file));
        } catch(Exception ex) {
            ex.printStackTrace();
        }                        	  
//        try {
//            FileInputStream istream = new FileInputStream(outputfilename);
//            ObjectInputStream ois = new ObjectInputStream(istream);
//            xi = (double []) ois.readObject(); 
//            yi = (double []) ois.readObject(); 
//            omegaX = IntArray2Set((int[]) ois.readObject());
//            omegaY = IntArray2Set((int[]) ois.readObject());
//            lambda = (double) ois.readObject(); 
//            fileName = (String) ois.readObject();
//        } catch(Exception ex) {
//            ex.printStackTrace();
//        }
        
    }

    
    public double[] median_WTx(double[] x){
    	int p = graph.W[0].length;
    	double sum = StatUtils.sum(x);
    	double[] medians = new double[p];
    	ArrayList<Integer> S = new ArrayList<Integer>();
    	for(int i = 0; i < x.length; i++){
    		if(x[i] != 0){
    			S.add(i);
    		}
    	}
    	for(int j=0; j<p; j++){
    		double[] vals = new double[S.size()];
    		for(int k=0; k<S.size(); k++){
    			vals[k] = x[S.get(k)] * graph.W[S.get(k)][j];
    		}
    		medians[j] = StatUtils.percentile(vals, 50) * sum;
    	}
    	return medians;
    }
    
    public double[] squaretimevector(double[][] M, double[] v){
    	double[] out = new double[M.length];
    	for(int i = 0; i < M.length; i++){
    		out[i] = 0;
    		for(int j = 0; j < v.length; j++){
    			out[i] += M[i][j] * M[i][j] * v[j];
    		}
    	}
    	return out;
    }    

    
//  double[] Ix = new double[n];
//	for(int i:graph.trueNodesSet){
//		xi[i] = 1;
//	}
//
//	for(int i:graph.trueFeas){
//		yi[i] = 1;  
//	}
//  Arrays.fill(Ix, 1.0D); 
//  double[][] W_transpose = Matrix.transpose(graph.W);
//  double xT1 = Matrix.dot(xi, Ix);
//  double[] WTx;
//  WTx=Matrix.VecMultiplyMat(xi,graph.W);
//  double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
//  showstat(WTx, "WTx");
//  showstat(WTx_1Tx, "WTx_1Tx"); 
//  WTx = median_WTx(xi);
//  WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
//  showstat(WTx, "WTx");
//  showstat(WTx_1Tx, "WTx_1Tx");
//  double[] comp1 = Matrix.VarMultiplyVec(1.0/0.1, squaretimevector(Matrix.subtract(W_transpose,Matrix.VecOutterPro(WTx_1Tx, Ix)), xi));
//  double[] comp2 = Matrix.VarMultiplyVec(1.0/1, squaretimevector(W_transpose, xi));
//  double[] gradient = Matrix.VecSubstract(comp1, comp2);
//  showstat(gradient, "gradient");
//  initilization1(xi, yi, graph);

    
    private boolean run() {
        long startTime = System.nanoTime();
		List<double[]> xy0 = CalcInitilvalues2(20, 3);
        runTime = (System.nanoTime() - startTime) / 1e9;
        System.out.println("CalcInitilvalues2 run time:" + runTime);
        double[] xi = new double[n];
        double[] yi = new double[p];
        System.arraycopy(xy0.get(0), 0, xi, 0, xi.length);
        System.arraycopy(xy0.get(1), 0, yi, 0, yi.length);
        List<Double> fValues = new ArrayList<>();
        double oldFunValue = -1;
        for (int i = 0; i < this.iterations; i++) { 
//        	System.out.println("iteration: " + i);
            fValues.add(function.getFuncValue(xi, yi,lambda));
//            if(verboseLevel>0){
//                System.out.println(sum(xi) + ", " + sum(yi));
//            }
            double[] gradientFx = function.getGradientX(xi, yi,lambda);
            double[] gradientFy = function.getGradientY(xi, yi,lambda);
            if(verboseLevel>0){
//	            System.out.println(StatUtils.min(gradientFx) + ", " + StatUtils.max(gradientFx));
//	            System.out.println(StatUtils.min(gradientFy) + ", " + StatUtils.max(gradientFy));
	            showstat(xi, "x:");
	            showstat(yi, "y:");
	            showstat(gradientFx, "gradientFx: ", -1000, 0);
	            showstat(gradientFy, "gradientFy: ", -1000, 0);
            }
//            if(true) return true;

            gradientFx = Stat.normalizeGradient(xi, gradientFx, n);
            gradientFy = Stat.normalizeGradient(yi, gradientFy, p);
            if(verboseLevel>0){
	            System.out.println(StatUtils.min(gradientFx) + ", " + StatUtils.max(gradientFx));
	            System.out.println(StatUtils.min(gradientFy) + ", " + StatUtils.max(gradientFy));
	            for(int j = 0; j < gradientFy.length; j++){
	            	if(gradientFy[j] != 0){
	            		System.out.println(j + ", " + gradientFy[j]);
	            	}
	            }
            }
            HeadApprox head = new HeadApprox(graph.edges, graph.edgeCosts, gradientFx, k, g, B, graph.trueSubGraphNodes, false);
            Set<Integer> gammaX = new HashSet<>(head.bestForest.nodesInF);
            Set<Integer> gammaY = identifyDirection(gradientFy, 2 * s);
            Set<Integer> omegaX = Sets.union(gammaX, Stat.supp(xi));
//            System.out.println(omegaX.size());
            Set<Integer> omegaY = Sets.union(gammaY, Stat.supp(yi));
//            StoreTestCaseData(xi, yi, omegaX, omegaY, lambda, fileName, "testcases/getArgMinFxy/test-case2.dat");
            List<double[]> b = function.getArgMinFxy(xi, yi, omegaX, omegaY,lambda);   /** (b_x,b_y)= argMax f(x,y) s.t. supp(x) \in omegaX, supp(y) \in omegaY*/
            double[] bx = b.get(0);
            double[] by = b.get(1);
            TailApprox tail = new TailApprox(graph.edges, graph.edgeCosts, bx, k, g, B, graph.trueSubGraphNodes,false);
            psiX = new HashSet<>(tail.bestForest.nodesInF);
            psiY = identifyDirection(by, s);
            xi = projectionOnVector(bx, psiX);
            yi = projectionOnVector(by, psiY);
            funcValue = function.getFuncValue(xi, yi,lambda);
            if(oldFunValue == -1){
            	oldFunValue = funcValue;
            }else{
            	if(Math.abs(oldFunValue - funcValue) < 0.001){
            		break;
            	}else{
            		oldFunValue = funcValue;
            	}
            }
        }

//        for (int i = 0; i < this.iterations; i++) { 
//            fValues.add(function.getFuncValue(xi, yi,lambda));
//            double[] gradientFy = function.getGradientY(xi, yi,lambda);
//            gradientFy = Stat.normalizeGradient(yi, gradientFy, p);
//            Set<Integer> gammaY = identifyDirection(gradientFy, 2 * s);
//            Set<Integer> omegaY = Sets.union(gammaY, Stat.supp(yi));
//            double[] by = function.getArgMinFy(xi, omegaY,lambda);   /** (b_x,b_y)= argMax f(x,y) s.t. supp(x) \in omegaX, supp(y) \in omegaY*/
//            psiY = identifyDirection(by, s);
//            yi = projectionOnVector(by, psiY);
//            funcValue = function.getFuncValue(xi, yi,lambda);
//        }        
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
    
    public int[] Set2IntArray(Set<Integer> x){
    	int[] data = new int[x.size()];
    	int i = 0;
    	for(int j:x){
    		data[i] = j;
    		i += 1;
    	}
    	return data;
    }

    public Set<Integer> IntArray2Set(int[] x){
    	Set<Integer> mySet = new HashSet<Integer>();
    	for(int i:x){
    		mySet.add(i);
    	}
    	return mySet;
    }
    
    
//    private boolean run() {
//        long startTime = System.nanoTime();
//        double[] x0 = new double[n];
//        double[] y0 = new double[p];
//        Arrays.fill(x0, 0.0D);
//        Arrays.fill(y0, 0.0D);
////        int cnt = 0;
////        for(int j:graph.trueSubGraphNodes){
////        	x0[j] = 1;
////        	cnt += 1;
////        	if(cnt > 5)
////        		break; 
////        }
////        cnt = 0;
////        for(int j:graph.trueFeas){
////        	y0[j] = 1;
////        	cnt += 1;
////        	if(cnt > 1)
////        		break; 
////        }
//        
//        List<double[]> xy0 = CalcInitilvalues(k, s);
//        System.arraycopy(xy0.get(0), 0, x0, 0, x0.length);
//        System.arraycopy(xy0.get(1), 0, y0, 0, y0.length);
//        
////        xi = Stat.getRandomVector(n,k);
////        yi = Stat.getRandomVector(p,s);
//        double[] xi = new double[n];
//        double[] yi = new double[p];
//        System.arraycopy(x0, 0, xi, 0, x0.length);
//        System.arraycopy(y0, 0, yi, 0, y0.length);
////        System.out.println(yi.toString());
//        List<Double> fValues = new ArrayList<>();
//        for (int i = 0; i < this.t; i++) { // t iterations
//            if (verboseLevel > 0) {
//                System.out.println("\n\nS2GraphMPDense k="+k+" lambda="+lambda+" ------------iteration: " + i + "------------");
//            }
////            /** fixed x */
////            Arrays.fill(xi, 0.0D);
////            for(int k:graph.trueSubGraphNodes){
////            	xi[k]=1.0D;        	
////            }
////            /** fixed y */
////          Arrays.fill(yi, 0.0D);
////          for(int m:graph.trueFeas){
////        	yi[m]=1.0D;        	
//////        	break;
////          }
//            fValues.add(function.getFuncValue(xi, yi,lambda));
//            double[] gradientFx = function.getGradientX(xi, yi,lambda);
//            double[] gradientFy = function.getGradientY(xi, yi,lambda);
//            if (verboseLevel > 0) {
//                System.out.println("supp(X): " + Stat.getSupp(xi).length + " " + Stat.supp(xi).toString());
//                System.out.println("supp(Y): " + Stat.getSupp(yi).length + " " + Stat.supp(yi).toString());
//                System.out.println("functionValue: " + function.getFuncValue(xi, yi));
//                System.out.println("X: " + Doubles.join(" ", xi));
//                System.out.println("Y: " + Doubles.join(" ", yi));
//                System.out.println("gradientFx: " + Doubles.join(" ", gradientFx));
//                System.out.println("gradientFy: " + Doubles.join(" ", gradientFy));
//                System.out.println("gradientFx: " +Stat.printWithIndex(gradientFx));
//            }
//            gradientFx = Stat.normalizeGradient(xi, gradientFx, n);
//            gradientFy = Stat.normalizeGradient(yi, gradientFy, p);
//            if (verboseLevel > 0) {
//                System.out.println("Normalized gradientFx: " + Doubles.join(" ", gradientFx));
//                System.out.println("Normalized gradientFy: " + Doubles.join(" ", gradientFy));
//                System.out.println("gradientFx: " +Stat.printWithIndex(gradientFx));
//            }
//            /**Algorithm S2-GraphMPDense No. is the line number in algorithm*/
//            /**5. Gamma_x head approximation**/
//            HeadApprox head = new HeadApprox(graph.edges, graph.edgeCosts, gradientFx, k, g, B, graph.trueSubGraphNodes, false);
//            if(verboseLevel > 0){
//            	System.out.println("S2GraphMP head projection resutl: "+ArrayUtils.toString(head.bestForest.nodesInF));
//            }
//            Set<Integer> gammaX = new HashSet<>(head.bestForest.nodesInF);
//            /**6. Gamma_y = argmax ||\delta_y f(x,y)|| : |R|<2s, identify the direction*/
//            Set<Integer> gammaY = identifyDirection(gradientFy, 2 * s);
//            /**7. Omiga_x*/
//            Set<Integer> omegaX = Sets.union(gammaX, Stat.supp(xi));
//            /**8. Omiga_y*/
//            Set<Integer> omegaY = Sets.union(gammaY, Stat.supp(yi));
//
//            if (verboseLevel > 0) {
//                //Collections.sort(Stat.getArray(gammaX));
//                System.out.println("gammaX: " + gammaX.toString());
//                System.out.println("OmigaX: " + omegaX.toString());
//                System.out.println("head  : " + ArrayUtils.toString(head.bestForest.nodesInF));
//                System.out.println("suppX : " + Stat.supp(xi).toString());
//                System.out.println("OmigaY: " + omegaY.toString());
//            }
//            /** 9.(b_x,b_y)= argMax f(x,y)*/
//            List<double[]> b = function.getArgMinFxy(xi, yi, omegaX, omegaY,lambda);
////            double[] bby = function.getArgMinFy(xi, omegaY,lambda);
////            double[] bbx = function.getArgMinFx(graph.trueNodesSet, graph.trueFeas, yi,omegaX,lambda);
////            List<double[]> b = new ArrayList<>();
////            b.add(bbx);
////            b.add(yi);
////            b.add(xi);
////            b.add(bby);
//            double[] bx = b.get(0);
//            double[] by = b.get(1);
//            if (verboseLevel > 0) {
//                System.out.println("supp(bx): " + Stat.supp(bx));
//                System.out.println("Intersection: " + Sets.intersection(graph.trueNodesSet, Stat.supp(bx)));
//                System.out.println("supp(by): " + Stat.supp(by).toString());
//                System.out.println("bx: " + Doubles.join(" ", bx));
//                System.out.println("by: " + Doubles.join(" ", by));
//            }
//            /**10.\Psi_x^{i+1} tail approximation of added density constrain, where dense(S)>=gamma*/
//            //DenseProjection tail =new DenseProjection(graph, k, g, B,t,lambdas,bx,gamma);
//            TailApprox tail = new TailApprox(graph.edges, graph.edgeCosts, bx, k, g, B, graph.trueSubGraphNodes,false);
//            if (verboseLevel > 0) {
//                System.out.println("number of tail nodes: " + tail.bestForest.nodesInF.size());
//                System.out.println("tail nodes: " + ArrayUtils.toString(tail.bestForest.nodesInF));
//            }
//            psiX = new HashSet<>(tail.bestForest.nodesInF);
//            /**11. \Psi y = argmax ||\b_y|| : |R|<s*/
//            psiY = identifyDirection(b.get(1), s);
//            //System.out.println("psiY: " + psiY.toString());
//            /** 12. calculate x^{i+1} */
//            xi = projectionOnVector(b.get(0), psiX);
//            /** 13. calculate y^{i+1} */
//            yi = projectionOnVector(b.get(1), psiY);
//            if (verboseLevel > 0) {
//                System.out.println("supp(X_i+1): " + Stat.supp(xi).toString());
//                System.out.println("X_i+1: " + Doubles.join(" ", xi));
//                System.out.println("supp(Y_i+1): " + Stat.supp(yi).toString());
//                System.out.println("number of head nodes : " + head.bestForest.nodesInF.size());
//                System.out.println("head nodes: " + head.bestForest.nodesInF.toString());
//                System.out.println("number of tail nodes : " + tail.bestForest.nodesInF.size());
//                System.out.println("tail nodes: " + tail.bestForest.nodesInF.toString());
//            }
////            Arrays.fill(yi, 0.0D);
////            for(int k:graph.trueFeas){
////            	yi[k]=1.0D;        	
////            }
////            Arrays.fill(xi, 0.0D);
////            for(int k:graph.trueSubGraphNodes){
////            	xi[k]=1.0D;        	
////            }
//            funcValue = function.getFuncValue(xi, yi,lambda);
//        }
//        runTime = (System.nanoTime() - startTime) / 1e9;
//        return true;
//    }
    
    private Set<Integer> identifyDirection(double[] gradientFy, int s) {
        Set<Integer> GammaY = new HashSet<>();
        double[] absGradientFy = new double[gradientFy.length];
        for(int i = 0 ; i < gradientFy.length ; i++){
            absGradientFy[i] = Math.abs(gradientFy[i]);
        }
        ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(absGradientFy);
        Integer[] indexes = arrayIndexComparator.getIndices();
        Arrays.sort(indexes, arrayIndexComparator);
        for (int i = 0; i < s; i++) {
            if(i < indexes.length && gradientFy[indexes[i]] != 0){
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
   
    /** ****************************************************************************************************************************************** */
    public static void TestingS2GraphMPDense(String fileName,String resultFileName, int k, int s, double lambda){
		long startTimeAll = System.nanoTime();
		APDMInputFormat apdm = null;
		MultiDataGraph mgdGraph = null;
		

			apdm=new APDMInputFormat(fileName);
		System.out.println("att:"
				+ ArrayUtils.toString(apdm.data.attributes[0]));
			mgdGraph=new MultiDataGraph(apdm);


		PCAScore func=new PCAScore(mgdGraph.W,mgdGraph.AdjMatrix);				

		int g=1;
		int iterations=10;
		double B = k - g + 0.0D;

		TestSGPursuitDenseDetection s2GraphMPDense = new TestSGPursuitDenseDetection(
				fileName, mgdGraph, k, s, g, B, iterations, func, lambda);

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
		double[] a1 = indicateVector(s2GraphMPDense.psiX, N); // top
																				// 1
																				// results
																				// indicatr
		// vector
		double[] a1C = Matrix.VecSubstract(I, a1);
		double aW1_a = Matrix.dot(Matrix.VecMultiplyMat(a1, adjMatrix), a1C);
		double aw1 = Matrix.dot(Matrix.VecMultiplyMat(a1, adjMatrix), I);
		nCut = aW1_a / aw1;
		System.out.println("X:-" + s2GraphMPDense.psiX.toString());
		System.out.println("Normalized Cut=:" + nCut);
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
    
	public static double[] indicateVector(Set<Integer> a, int n) {
		double[] x = new double[n];
		for (int i : a) {
			x[i] = 1.0;
		}

		return x;
	}


	public static void log(String line){
		FileWriter fstream;
		System.out.println(line);		
		try {
			fstream = new FileWriter("data/DenseGraph/Log/log.txt", true);
			fstream.write(line + "\r\n");
			fstream.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

    public static  void UnitTest1() {
    	long startTimeAll = System.nanoTime();
		double lambda     = 5.00D;
		double in_p       = 0.35;
		double out_p      = 0.15;
		int numFea        = 100;
		int trueSubSize   = 200;
        ExecutorService pool = Executors.newFixedThreadPool(8);
        
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date = new Date();
		log("\n\n========================" + dateFormat.format(date) + "=============================\n\n");
		
		String root="data/DenseGraph/DenseSubgraph_APDM/";
    	for (double r : new double[]{0.05}){// , 0.10, 0.15, 0.20
    		for(double sigmas1:new double[]{0.0316D}){ // , 0.05D, 0.05D, 0.1D, 0.5D, 0.6D, 0.7D, 0.8D, 0.9D
	            for (int cases=0;cases<50;cases++) {
					Double dd = numFea * r;
					int numTueFeat = dd.intValue();
//	            	String fileName="APDM_Dense_subgraph_in_"+in_p+"_out_"+out_p+"_FeasNum_100_trueFeasNum_"+numTueFeat+"_sigmas1_" + sigmas1 + "_case_"+cases+".txt";
	            	String fileName="VaryingNumOfAttributes_APDM_Dense_subgraph_in_"+in_p+"_out_"+out_p+"_TrueSGSize_"+trueSubSize+"_FeasNum_"+numFea+"_trueFeasNum_"+numTueFeat+ "_sigmas1_" + sigmas1 + "_case_"+cases+".txt";
	            	String outFile="outputs/S2GraphMPDense/" + fileName;
	            	int trueSubSize1 = 50;
	                pool.execute(new Thread() {
	                    public void run() {
	                    	TestingS2GraphMPDense(root+fileName,outFile,trueSubSize1, numTueFeat, lambda);                        	
	                    }
	                });
	            }
    		}
        }
        pool.shutdown();
		System.out.println("Total running time: " + (System.nanoTime() - startTimeAll) / 1e9);
	}
    
    public static  void UnitTest2() {
    	long startTimeAll = System.nanoTime();
		double lambda     = 0.00D;
		double p_in       = 0.35;
		double p_out      = 0.15;
		int numFea        = 100;
//		int trueSubSize   = 200;
		int numClusters = 20;
		int clusterSize = 100;
        ExecutorService pool = Executors.newFixedThreadPool(8);
        
		DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
		Date date = new Date();
		log("\n\n========================" + dateFormat.format(date) + "=============================\n\n");
		
		String root="data/DenseGraph/DenseSubgraph_APDM/";
    	for (double r : new double[]{0.05}){// , 0.10, 0.15, 0.20
    		for(double sigmas1:new double[]{0.0316D}){ // , 0.05D, 0.05D, 0.1D, 0.5D, 0.6D, 0.7D, 0.8D, 0.9D
	            for (int instance=0;instance<100;instance++) {
					Double dd = numFea * r;
					int numTueFeat = dd.intValue();
//	            	String fileName="APDM_Dense_subgraph_in_"+in_p+"_out_"+out_p+"_FeasNum_100_trueFeasNum_"+numTueFeat+"_sigmas1_" + sigmas1 + "_case_"+cases+".txt";
//	            	String fileName="APDM_Dense_subgraph_in_"+p_in+"_out_"+p_out+"_TrueSGSize_"+trueSubSize+"_FeasNum_"+numFea+"_trueFeasNum_"+numTueFeat+ "_sigmas1_" + sigmas1 + "_case_"+cases+".txt";
	            	String fileName="APDM_Dense_subgraph_in_"+p_in+"_out_"+p_out+"_numClusters_"+numClusters+"_TrueSGSize_"+clusterSize+"_FeasNum_"+numFea+"_trueFeasNum_"+numTueFeat+ "_sigmas1_" + sigmas1 + "_case_"+instance+".txt";
	            	String outFile="outputs/S2GraphMPDense/Result_"+p_in+"_out_"+p_out+"_FeasNum_100_trueFeasNum_"+numTueFeat+"_case_"+instance+".txt";
//	        		log("\n\n\n\n" + fileName);
//	            	int clusterSize = 50;
	            	int numTueFeat1 = numTueFeat + 15; 
	                pool.execute(new Thread() {
	                    public void run() {
	                    	TestingS2GraphMPDense(root+fileName,outFile,clusterSize, numTueFeat1, lambda);                        	
	                    }
	                });
	            }
    		}
        }
        pool.shutdown();
		System.out.println("Total running time: " + (System.nanoTime() - startTimeAll) / 1e9);
	}
    
	public static void experiments_Realdata() {
		long startTimeAll = System.nanoTime();
    	double lambda = 5.00D;

		String fileName = "data/DenseGraph/GraphMLData/RealData/dblp_top_APDM.txt";
		String outFile = "data/DenseGraph/GraphMLData/RealData/dblp_top_APDM_result.txt";


		TestingS2GraphMPDense(fileName, outFile, 5, 3, lambda);
                      
        System.out.println("Total running time: " + (System.nanoTime() - startTimeAll) / 1e9);
	}    

	public static void experiments_RealdataSets() {

		double lambda = 5.00D;
		String root = "data/DenseGraph/Realdata/apdm/";
		String outroot = "outputs/GAMerResult/";

		for (String fileName : new String[] { "dblp", "imdb",
				"arxivSmall", "genes", "patents", "arxivLarge", "dfb" }) {
			fileName += "_APDM.txt";
			long startTimeAll = System.nanoTime();
			System.out.println(fileName);
			TestingS2GraphMPDense(root + fileName, "", 5, 3, lambda);

			System.out.println("Total running time: "
					+ (System.nanoTime() - startTimeAll) / 1e9);
		}
	}
    public static double mad(double [] autoCorrelationValues){
        double [] tempTable = new double[autoCorrelationValues.length];
        Median m = new Median();
        double medianValue = m.evaluate(autoCorrelationValues);
        for(int i=0 ; i<autoCorrelationValues.length ;i++){
            tempTable[i] = Math.abs(autoCorrelationValues[i] - medianValue);
        }
        return m.evaluate(tempTable); 
    }
    
        

           
    public static void validateAPDMFiles(){
    	String inRoot = "data/DenseGraph/DenseSubgraph_APDM/";
    	File folder = new File(inRoot);
    	File[] listOfFiles = folder.listFiles();
	    for (int i = 0; i < listOfFiles.length; i++) {
	      if (listOfFiles[i].isFile()) {
//	        System.out.println("File " + listOfFiles[i].getName());
			APDMInputFormat apdm = null;
			MultiDataGraph mgdGraph = null;
			String fileName = listOfFiles[i].getName();
			try{
				apdm=new APDMInputFormat(inRoot + fileName);
				mgdGraph=new MultiDataGraph(apdm);
			} 
			catch(Exception ex) {
				System.out.println("!!!!!!!!!!!!!! apdm file reading error..." + fileName);
				File f = new File(inRoot + fileName);
				f.delete();
				System.out.println("!!!!!!!!!!!!!! apdm file deleted..." + fileName);
//			    ex.printStackTrace();
//			    return;
			}   
	        
	      } else if (listOfFiles[i].isDirectory()) {
//	        System.out.println("Directory " + listOfFiles[i].getName());
	      }
	    }    	
    }
    
    public static void main(String[] args){

		experiments_RealdataSets();
   
    }
}
