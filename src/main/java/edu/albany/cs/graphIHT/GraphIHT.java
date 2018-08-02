package edu.albany.cs.graphIHT;

import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.Stat;
import edu.albany.cs.tailApprox.TailApprox;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;

/**
 * Algorithm 1 : Graph-IHT algorithm
 *
 * @author Baojian Zhou (bzhou6@albany.edu)
 */
public class GraphIHT {

    /** graphSize, number of nodes in the input graph */
    private final int graphSize;
    /** edges in our graph, should notice that the graph should be connected. */
    private final ArrayList<Integer[]> edges;
    /** we use identity edge cost 1.0D in our algorithm. */
    private final ArrayList<Double> edgeCosts;
    /** total sparsity of S */
    private final int s;
    /** # of connected components in forest F */
    private final int g;
    /** bound on the total weight w(F) of edges in the */
    private final double B;
    private final double[] yi;
    /** statistic function */
    private final Function function;
    /** parameter delta in our algorithm */
    private final double eta = 0.5D;
    /** connected components of current graph */
    private ConnectedComponents cc;
    /** halting condition */
    private double epsilon = 1e-6;

    private final int[] trueSubGraph;
    private int verboseLevel = 1;

    public int[] supportX;
    /** final vector x that we get in our algorithm. */
    public ArrayRealVector x;
    public HashSet<Integer> resultNodesTail;
    public double funcValueTail = 0.0D;
    public double runTime = 0.0D;
    public ArrayList<Double> fValues;
    /** number of iterations */
    public int iterNum;
    
    public GraphIHT(ArrayList<Integer[]> edges,ArrayList<Double> edgeCosts,int numOfNodes,Function func, double[] yi, int s, int g,
                    double B, int[] trueSubGraph) {
        this.graphSize = numOfNodes;
        this.edgeCosts = edgeCosts;
        this.edges = edges;
        this.yi = yi;
        this.s = s;
        this.g = g;
        this.B = B;
        this.trueSubGraph = trueSubGraph;
        this.function = func;
        this.iterNum=0;
        
      
        if (checkInput()) {
        	x = run();
            
        } else {
            x = null;
            System.out.println("input parameter is invalid. ");
            System.exit(0);
        }
    }

    private ArrayRealVector run() {
        long startTime = System.nanoTime();    	
        ArrayRealVector xi = new ArrayRealVector(Stat.getRandomVector(graphSize, s));
        ArrayRealVector biSave=null;
        fValues = new ArrayList<>();
        fValues.add(function.getFuncValue(xi.toArray(),yi));        
        while (true) {            
             iterNum++;
            /** Gradient for the function. */
           
            double[] gradientFxi = function.getGradientX(xi.toArray(),yi);
            /** TODO check this further */
            gradientFxi = Stat.normalizeGradient(xi.toArray(), gradientFxi,graphSize);
            /** head approximation */
            HeadApprox pcsfHead = new HeadApprox(edges, edgeCosts, gradientFxi, s, g, B, trueSubGraph, false);
            /** get head projection vector. */
            ArrayRealVector projectedHeader = projectVector(gradientFxi, pcsfHead.bestForest.nodesInF);
            /** get yi */
            ArrayRealVector bi = (new ArrayRealVector(xi)).add(projectedHeader.mapMultiply(-eta));
            biSave=bi;
            /** tail approximation */
            TailApprox pcsfTail = new TailApprox(edges, edgeCosts, bi.toArray(), s, g, B, trueSubGraph, false);
            /** get x_{i+1} */
            xi = projectVector(bi.toArray(), pcsfTail.bestForest.nodesInF);
            xi = normalize(xi);
            if(StatUtils.sum(xi.toArray())==0.0D){            	
            	xi = new ArrayRealVector(Stat.getRandomVector(graphSize,1));
            	System.out.println("\n xi all zeros...........\n");
            	System.out.println("---------------------------------\n");
            	System.out.println("s="+s+"\n bi= "+ArrayUtils.toString(biSave.toArray()+"\n"));
            	System.out.println("pcsfTail.bestForest.nodesInF: "+pcsfTail.bestForest.nodesInF.toString()+"\n");
            	System.out.println("---------------------------------\n");
           
            	
            }
            resultNodesTail = new HashSet<>(pcsfTail.bestForest.nodesInF);
            double updatedFuncValue = getFuncTailVal(resultNodesTail);
            if ((updatedFuncValue - fValues.get(fValues.size() - 1)) < epsilon) {
                fValues.add(updatedFuncValue);
                break;
            } else {
                fValues.add(updatedFuncValue);
            }
        }
        
        supportX = Stat.getSupp(xi.toArray());
        runTime = (System.nanoTime() - startTime) / 1e9;
        funcValueTail = getFuncTailVal(resultNodesTail);
        return xi;
    }



    private double getFuncTailVal(HashSet<Integer> nodesInF) {
        double[] tmpX = new double[graphSize];
        Arrays.fill(tmpX, 0.0D);
        for (int i : nodesInF) {
            tmpX[i] = 1.0D;
        }
        return function.getFuncValue(tmpX,yi);
    }

    /**
     * @return true if the input parameters are valid.
     */
    private boolean checkInput() {
        Set<Integer> nodes = new HashSet<>();
        for (Integer[] edge : edges) {
            nodes.add(edge[0]);
            nodes.add(edge[1]);
        }
        /** nodes is inconsistent */
        if (nodes.size() != graphSize) {
            return false;
        }
        ArrayList<ArrayList<Integer>> adj = new ArrayList<>();
        for (int i = 0; i < graphSize; i++) {
            adj.add(new ArrayList<Integer>());
        }
        for (Integer[] edge : this.edges) {
            adj.get(edge[0]).add(edge[1]);
            adj.get(edge[1]).add(edge[0]);
        }
        cc = new ConnectedComponents(adj);
        return cc.connectivity;
    }

    private ArrayRealVector projectVector(double[] x, ArrayList<Integer> projectionSet) {
        ArrayRealVector result = new ArrayRealVector(x);
        if (projectionSet == null) {
            return new ArrayRealVector(x);
        } else {
            for (int i = 0; i < result.getDimension(); i++) {
                if (!projectionSet.contains(i)) {
                    result.setEntry(i, 0.0D);
                }
            }
        }
        return result;
    }

    /**
     * normalize a vector. if x[i] < 0.0, then x[i]:= 0.0 if x[i] > 1.0, then
     * x[i]:= 1.0
     *
     * @param x
     *            input vector
     * @return
     */
    private ArrayRealVector normalize(ArrayRealVector x) {
        ArrayRealVector result = new ArrayRealVector(x);
        for (int i = 0; i < x.getDimension(); i++) {
            if (x.getEntry(i) < 0.0D) {
                result.setEntry(i, 0.0D);
            }
            if (x.getEntry(i) > 1.0D) {
                result.setEntry(i, 1.0D);
            }
        }
        return result;
    }

}
