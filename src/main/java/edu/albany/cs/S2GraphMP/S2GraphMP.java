package edu.albany.cs.S2GraphMP;


import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.Stat;
import edu.albany.cs.tailApprox.TailApprox;

import java.util.*;


/**
 * Algorithm 1 : S2Graph-Mp Aglorithm .
 *
 * @author
 */
public class S2GraphMP {

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
    private final int t;
    /**
     * Function to be optimized.
     */
    private final Function function;

    /**
     * results
     */
    public double[] xi;
    public double[] yi;
    public Set<Integer> psiX;
    public Set<Integer> psiY;
    public double funcValue = -1.0D;
    public double runTime;

    private int verboseLevel = 0;

    public S2GraphMP(Graph graph, int k, int s, int g, Function func) {
        this(graph, k, s, g, k - g + 0.0D, func);
    }

    public S2GraphMP(Graph graph, int k, int s, int g, double B, Function func) {
        this(graph, k, s, g, B, 5, func);
    }

    public S2GraphMP(Graph graph, int k, int s, int g, double B, int t, Function func) {
        this.graph = graph;
        this.n = graph.numOfNodes;
        this.p = graph.numOfFeas;
        this.k = k;
        this.s = s;
        this.g = g;
        this.B = B;
        this.t = t;
        this.function = func;
        // run the algorithm
        if (run()) {
            //System.out.println("S2GraphMP successfully finished.");
        } else {
            //System.out.println("S2GraphMP is failed.");
        }
    }

    private boolean run() {

        long startTime = System.nanoTime();
        xi = Stat.getRandomVector(n);
        yi = Stat.getRandomVector(p);
        List<Double> fValues = new ArrayList<>();
        for (int i = 0; i < this.t; i++) { // t iterations
            if (verboseLevel > 0) {
                System.out.println("\n\n------------iteration: " + i + "------------");
            }
            fValues.add(function.getFuncValue(xi, yi));
            double[] gradientFx = function.getGradientX(xi, yi);
            double[] gradientFy = function.getGradientY(xi, yi);
            if (verboseLevel > 0) {
                System.out.println("supp(X): " + Stat.getSupp(xi).length + " " + Stat.supp(xi).toString());
                System.out.println("supp(Y): " + Stat.getSupp(yi).length + " " + Stat.supp(yi).toString());
                System.out.println("functionValue: " + function.getFuncValue(xi, yi));
                System.out.println("X: " + Doubles.join(" ", xi));
                System.out.println("Y: " + Doubles.join(" ", yi));
                System.out.println("gradientFx: " + Doubles.join(" ", gradientFx));
                System.out.println("gradientFy: " + Doubles.join(" ", gradientFy));
            }
            gradientFx = Stat.normalizeGradient(xi, gradientFx, n);
            gradientFy = Stat.normalizeGradient(yi, gradientFy, p);
            if (verboseLevel > 0) {
                System.out.println("Normalized gradientFx: " + Doubles.join(" ", gradientFx));
                System.out.println("Normalized gradientFy: " + Doubles.join(" ", gradientFy));
            }
            /**Algorithm S2-GraphMP No. is the line number in algorithm*/
            /**5. Gamma_x head approximation H(\delta_x f(x,y))*/
		    HeadApprox head = new HeadApprox(graph.edges, graph.edgeCosts, gradientFx, k, g, B, graph.trueSubGraphNodes, false);
            Set<Integer> gammaX = new HashSet<>(head.bestForest.nodesInF);
            /**6. Gamma_y = argmax ||\delta_y f(x,y)|| : |R|<2s, identify the direction*/
            Set<Integer> gammaY = identifyDirection(gradientFy, 2 * s);
            /**7. Omiga_x*/
            Set<Integer> omegaX = Sets.union(gammaX, Stat.supp(xi));
            /**8. Omiga_y*/
            Set<Integer> omegaY = Sets.union(gammaY, Stat.supp(yi));

            if (verboseLevel > 0) {
                Collections.sort(head.bestForest.nodesInF);
                System.out.println("OmigaX: " + omegaX.toString());
                System.out.println("head  : " + head.bestForest.nodesInF.toString());
                System.out.println("suppX : " + Stat.supp(xi).toString());
                System.out.println("OmigaY: " + omegaY.toString());
            }
            /** 9.(b_x,b_y)= argMax f(x,y)*/
            List<double[]> b = function.getArgMinFxy(omegaX, omegaY);
            double[] bx = b.get(0);
            double[] by = b.get(1);
            if (verboseLevel > 0) {
                System.out.println("supp(bx): " + Stat.supp(bx));
                System.out.println("Intersection: " + Sets.intersection(graph.trueNodesSet, Stat.supp(bx)));
                System.out.println("supp(by): " + Stat.supp(by).toString());
                System.out.println("bx: " + Doubles.join(" ", bx));
                System.out.println("by: " + Doubles.join(" ", by));
            }
            /**10.\Psi_x^{i+1} tail approximation */
            TailApprox tail = new TailApprox(graph.edges, graph.edgeCosts, bx, k, g, B, graph.trueSubGraphNodes, false);
            if (verboseLevel > 0) {
                System.out.println("number of tail nodes: " + tail.bestForest.nodesInF.size());
                System.out.println("tail nodes: " + tail.bestForest.nodesInF.toString());
            }
            psiX = new HashSet<>(tail.bestForest.nodesInF);
            /**11. \Psi y = argmax ||\b_y|| : |R|<s*/
            psiY = identifyDirection(b.get(1), s);
            //System.out.println("psiY: " + psiY.toString());
            /** 12. calculate x^{i+1} */
            xi = projectionOnVector(b.get(0), psiX);
            /** 13. calculate y^{i+1} */
            yi = projectionOnVector(b.get(1), psiY);
            if (verboseLevel > 0) {
                System.out.println("supp(X_i+1): " + Stat.getSupp(xi).toString());
                System.out.println("X_i+1: " + Doubles.join(" ", xi));
                System.out.println("supp(Y_i+1): " + Stat.getSupp(yi).toString());
                System.out.println("number of head nodes : " + head.bestForest.nodesInF.size());
                System.out.println("head nodes: " + head.bestForest.nodesInF.toString());
                System.out.println("number of tail nodes : " + tail.bestForest.nodesInF.size());
                System.out.println("tail nodes: " + tail.bestForest.nodesInF.toString());
            }
            funcValue = function.getFuncValue(xi, yi);
        }
        runTime = (System.nanoTime() - startTime) / 1e9;
        return true;
    }

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
            if(i < indexes.length){
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

}
