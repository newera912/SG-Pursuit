package edu.albany.cs.S2GraphMPDebug;


import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;

import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.graph.BreastCancerGraph;
import edu.albany.cs.graph.CrimesOfChicagoGraph;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.GridGraph;
import edu.albany.cs.graph.YelpGraph;
import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.scoreFuncs.Stat;
import edu.albany.cs.tailApprox.TailApprox;


/**
 * Algorithm 1 : S2Graph-Mp Aglorithm .
 *
 * @author
 */
public class S2GraphMPDebug {

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

    /**
     * For testing only
     */
    private Set<Integer> trueFeatures = null;
    private Set<Integer> trueNodes = null;
    private double[] trueY = null;
    private double[] trueX = null;
    private int verboseLevel = 0;
    private TestS2OnYelp.Data data;

    public S2GraphMPDebug(Graph graph, int k, int s, int g, Function func, TestS2OnYelp.Data data) {
        this(graph, k, s, g, k - g + 0.0D, func);
        this.data = data;
    }

    public S2GraphMPDebug(Graph graph, int k, int s, int g, double B, Function func) {
        this(graph, k, s, g, B, 5, func);
    }

    public S2GraphMPDebug(Graph graph, int k, int s, int g, double B, int t, Function func) {
        this.graph = graph;
        this.n = graph.numOfNodes;
        this.p = graph.numOfFeas;
        this.k = k;
        this.s = s;
        this.g = g;
        this.B = B;
        this.t = t;
        this.function = func;
        if (graph instanceof BreastCancerGraph) {
            trueFeatures = ((BreastCancerGraph) graph).trueFeatures;
            trueY = ((BreastCancerGraph) graph).trueY;
        } else if (graph instanceof CrimesOfChicagoGraph) {
            trueFeatures = ((CrimesOfChicagoGraph) graph).trueFeatures;
            trueNodes = ((CrimesOfChicagoGraph) graph).trueNodes;
            trueX = ((CrimesOfChicagoGraph) graph).trueX;
            trueY = ((CrimesOfChicagoGraph) graph).trueY;
        } else if (graph instanceof GridGraph) {
            trueFeatures = ((GridGraph) graph).trueFeatures;
            trueNodes = ((GridGraph) graph).trueNodes;
            trueX = ((GridGraph) graph).trueX;
            trueY = ((GridGraph) graph).trueY;
        } else if (graph instanceof YelpGraph) {
            trueFeatures = ((YelpGraph) graph).trueFeatures;
            trueNodes = ((YelpGraph) graph).trueNodes;
            trueX = ((YelpGraph) graph).trueX;
            trueY = ((YelpGraph) graph).trueY;
        }
        if (run()) {
            if (verboseLevel == 0) {
                System.out.println("S2GraphMP successfully finished.");
            }
        } else {
            System.out.println("S2GraphMP is failed.");
        }
    }


    private boolean run() {
        long startTime = System.nanoTime();
        double[] xi = new double[n];
        double[] yi = new double[p];
        List<double[]> xy0;
        if (function instanceof ElevatedMeanScan) {
            xy0 = ((ElevatedMeanScan) function).calcInitilvalues(k, s);
            xi = xy0.get(0);
            yi = xy0.get(1);
        }

		// TODO debug initialize x with true nodes
		Arrays.fill(xi, 0.0D);
		// Arrays.fill(yi, 0.0D);
		for (int i : trueNodes) {
			xi[i] = 1.0D;
		}

        System.out.println("finish initial");
        PreRec preRecNodes = new PreRec(Stat.supp(xi), this.trueNodes);
        System.out.println(preRecNodes);
        PreRec preRecFeatures = new PreRec(Stat.supp(yi), this.trueFeatures);
        System.out.println(preRecFeatures);
        //Utils.stop();
        double[] xOld;
        double[] yOld;
        int numOfIter = 0;
        while (true) {
            double[] gradientFx = function.getGradientX(xi, yi);
            double[] gradientFy = function.getGradientY(xi, yi);
            System.out.println("\n\n------------iteration: " + numOfIter + "------------");
            if (verboseLevel > 0) {
                System.out.println("\n\n------------iteration: " + numOfIter + "------------");
                double currentVal = function.getFuncValue(xi, yi);
                System.out.println("function value: " + currentVal);
                System.out.println("xi: " + Doubles.join(" ", xi));
                System.out.println("yi: " + Doubles.join(" ", yi));
                System.out.println("supp(X): " + Stat.supp(xi).size() + " " + Stat.supp(xi));
                Set<Integer> X = new HashSet<>();
                for (int i : ((GridGraph) graph).trueNodes) {
                    if (xi[i] > 0) X.add(i);
                }
                double prec = (X.size() * 1.0) / Stat.supp(xi).size();
                double recall = (X.size() * 1.0) / ((GridGraph) graph).trueNodes.size();
                System.out.println("true nodes " + String.format("(%.2f, %.2f)", prec, recall) + ": " + X);
                System.out.println("supp(Y): " + Stat.supp(yi).size() + " " + Stat.supp(yi));
                Set<Integer> Y = new HashSet<>();
                for (int i : ((GridGraph) graph).trueFeatures) {
                    if (yi[i] > 0) Y.add(i);
                }
                prec = (Y.size() * 1.0) / Stat.supp(yi).size();
                recall = (Y.size() * 1.0) / ((GridGraph) graph).trueFeatures.size();
                System.out.println("true features " + String.format("(%.2f, %.2f)", prec, recall) + ": " + Y);
            }
            gradientFx = Stat.normalizeGradient(xi, gradientFx, n);
            gradientFy = Stat.normalizeGradient(yi, gradientFy, p);
            HeadApprox head = new HeadApprox(graph.edges, graph.edgeCosts, gradientFx, k, g, B, false);
            Set<Integer> gammaX = new HashSet<>(head.bestForest.nodesInF);
            Set<Integer> gammaY = identifyDirection(gradientFy, 2 * s);
            Set<Integer> omegaX = Sets.union(gammaX, Stat.supp(xi));
            Set<Integer> omegaY = Sets.union(gammaY, Stat.supp(yi));
            List<double[]> b = function.getArgMinFxy(xi, yi, omegaX, omegaY);
            double[] bx = b.get(0);
            double[] by = b.get(1);
            TailApprox tail = new TailApprox(graph.edges, graph.edgeCosts, bx, k, g, B, false);
            psiX = new HashSet<>(tail.bestForest.nodesInF); /**10.\Psi_x^{i+1} tail approximation */
            psiY = identifyDirection(by, s); /**11. \Psi y = argmax ||\b_y|| : |R|<s*/
            xOld = xi;
            yOld = yi;
            xi = projectionOnVector(bx, psiX); /** 12. calculate x^{i+1} */
            yi = projectionOnVector(by, psiY); /** 13. calculate y^{i+1} */
            /** to check the halting condition. */
            funcValue = function.getFuncValue(xi, yi);
            double gapX = l2Norm(xOld, xi);
            double gapY = l2Norm(yOld, yi);
            numOfIter++;
            if ((gapX < 1e-3D && gapY < 1e-3D) || numOfIter > t) {
                this.yi = yi;
                this.xi = xi;
                break;
            }
        }
        runTime = (System.nanoTime() - startTime) / 1e9;
        return true;
    }

    /**
     * identify the direction, which means a subset of indices will be returned.
     * The goal is to maximize ||gradient_R||^2 s.t. |R| < s
     *
     * @param gradient the input vector.
     * @param s        the sparsity of that direction. # of nodes in returned subset.
     * @return a subset of indices.
     */
    private Set<Integer> identifyDirection(double[] gradient, int s) {
        if (gradient == null || gradient.length == 0 || s <= 0 || s > gradient.length) {
            System.out.println("gradient is null. ");
            System.exit(0);
        }
        /**if the vector is zero, just return an empty set.*/
        if (StatUtils.sum(gradient) == 0.0D) {
            return new HashSet<>();
        }
        Set<Integer> GammaY = new HashSet<>();
        double[] squareGradient = new double[gradient.length];
        for (int i = 0; i < gradient.length; i++) {
            squareGradient[i] = gradient[i] * gradient[i];
        }
        ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(squareGradient);
        Integer[] indexes = arrayIndexComparator.getIndices();
        Arrays.sort(indexes, arrayIndexComparator);
        for (int i = 0; i < s; i++) {
            if (i < indexes.length && (squareGradient[indexes[i]] > 0.0D)) {
                GammaY.add(indexes[i]);
            }
        }
        return GammaY;
    }

    /**
     * Return the projected vector
     *
     * @param vector the vector to be projected.
     * @param set    the projection set.
     * @return the projected vector.
     */
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
     * get L2 norm of two vectors.
     *
     * @param x1 the first vector.
     * @param x2 the second vector.
     * @return the norm ||x1 - x2||
     */
    private double l2Norm(double[] x1, double[] x2) {
        double norm = 0.0D;
        if (x1 == null || x2 == null || x1.length == 0 || x2.length == 0) {
            return norm;
        }
        if (x1.length != x2.length) {
            System.out.println("error");
            System.exit(0);
        } else {
            ArrayRealVector x1Vector = new ArrayRealVector(x1);
            ArrayRealVector x2Vector = new ArrayRealVector(x2);
            norm = (x1Vector.subtract(x2Vector)).getNorm();
        }
        return norm;
    }

    public static void main(String args[]) {
    }

}
