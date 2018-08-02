package edu.albany.cs.S2GraphMP;


import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.BreastCancerGraph;
import edu.albany.cs.graph.CrimesOfChicagoGraph;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.GridGraph;
import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.scoreFuncs.*;
import edu.albany.cs.tailApprox.TailApprox;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.util.*;


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
    private int verboseLevel = -1;

    public S2GraphMPDebug(Graph graph, int k, int s, int g, Function func) {
        this(graph, k, s, g, k - g + 0.0D, func);
    }

    public S2GraphMPDebug(Graph graph, int k, int s, int g, double B, Function func) {
        this(graph, k, s, g, B, 10, func);
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
        xi = initializeX0(2); //option 4 is the true x.

        if (p == 1) {
            yi = new double[]{1.0D};
        } else {
            yi = initializeY0(2); //option 4 is the true y.
        }

        double[] xOld;
        double[] yOld;
        int numOfIter = 0;
        while (true) {
            //---------------------------------------------------
            /** get gradient of x and y*/
            double[] gradientFx = function.getGradientX(xi, yi);
            double[] gradientFy;
            if (p == 1) {
                gradientFy = new double[]{1.0D};
            } else {
                gradientFy = function.getGradientY(xi, yi);
            }
            //---------------------------------------------------
            if (verboseLevel == 0) {
                System.out.println("initial function value: " + this.function.getFuncValue(xi, yi));
                System.out.println("iteration: " + numOfIter);
                System.out.println("sum(xi): " + StatUtils.sum(xi));
                System.out.println("sum(yi): " + StatUtils.sum(yi));
                System.out.println("gradientFy: " + Doubles.join(" ", gradientFy));
                System.out.println("----abnormal feature gradient ----");
                for (int node = 0; node < p; node++) {
                    if (trueFeatures.contains(node)) {
                        System.out.print(" " + gradientFy[node]);
                    }
                }
                System.out.println();
                System.out.println("----normal feature gradient ----");
                for (int node = 0; node < p; node++) {
                    if (!trueFeatures.contains(node)) {
                        System.out.print(" " + gradientFy[node]);
                    }
                }
                System.out.println();
                Utils.stop();
            }

            if (verboseLevel > 0) {
                System.out.println("\n\n------------iteration: " + numOfIter + "------------");
                double currentVal = function.getFuncValue(xi, yi);
                System.out.println("function value: " + currentVal);
                System.out.println("xi: " + Doubles.join(" ", xi));
                System.out.println("yi: " + Doubles.join(" ", yi));
                System.out.println("supp(X): " + Stat.supp(xi).size() + " " + Stat.supp(xi));
                System.out.println("supp(Y): " + Stat.supp(yi).size() + " " + Stat.supp(yi));
                System.out.println("sum(gradient x): " + StatUtils.sum(gradientFx));
                System.out.println("sum(gradient y): " + StatUtils.sum(gradientFy));
                System.out.println("gradientFy: " + Doubles.join(" ", gradientFy));
                if (trueNodes != null) {
                    System.out.println("----abnormal nodes gradient ----");
                    for (int node : trueNodes) {
                        System.out.print(" " + gradientFx[node]);
                    }
                    System.out.println();
                    System.out.println("----normal nodes gradient ----");
                    for (int k = 0; k < graph.numOfNodes; k++) {
                        if (!trueNodes.contains(k)) {
                            System.out.print(" " + gradientFx[k]);
                        }
                    }
                    System.out.println();
                }
                if (trueFeatures != null) {
                    System.out.println("----abnormal feature gradient ----");
                    for (int node : trueFeatures) {
                        System.out.print(" " + gradientFy[node]);
                    }
                    System.out.println();
                    System.out.println("----normal feature gradient ----");
                    for (int node = 0; node < p; node++) {
                        if (!trueFeatures.contains(node)) {
                            System.out.print(" " + gradientFy[node]);
                        }
                    }
                    System.out.println();
                }
            }

            //--------------------------------------------------------------
            /** do normalization on gradient of x and y*/
            if (function instanceof EMSXYScore || function instanceof ElevatedMeanScan) {
                gradientFx = Stat.normalizeGradient(xi, gradientFx, n);
            }
            gradientFy = Stat.normalizeGradient(yi, gradientFy, p);
            //--------------------------------------------------------------

            if (verboseLevel == 0) {
                System.out.println("after normalization ...");
                System.out.println("sum(gradient x): " + StatUtils.sum(gradientFx));
                System.out.println("sum(gradient y): " + StatUtils.sum(gradientFy));
                System.out.println("the normalized of gradientFy: ");
                System.out.println("----abnormal feature gradient ----");
                for (int node = 0; node < p; node++) {
                    if (trueFeatures.contains(node)) {
                        System.out.print(" " + gradientFy[node]);
                    }
                }
                System.out.println();
                System.out.println("----normal feature gradient ----");
                for (int node = 0; node < p; node++) {
                    if (!trueFeatures.contains(node)) {
                        System.out.print(" " + gradientFy[node]);
                    }
                }
                System.out.println();
                Utils.stop();
            }

            //-------------------------------------------------------------------------------------------------
            /** head projection and maximization .*/
            /**5. Gamma_x head approximation H(\delta_x f(x,y))*/
            if (verboseLevel == 0) {
                System.out.println(" to do head approximation ...");
            }
            HeadApprox head = new HeadApprox(graph.edges, graph.edgeCosts, gradientFx, k, g, B, false);
            if (verboseLevel == 0) {
                System.out.println("head approximation finished ...");
            }
            Set<Integer> gammaX = new HashSet<>(head.bestForest.nodesInF);
            /**6. Gamma_y = argmax ||\delta_y f(x,y)|| : |R|<2s, identify the direction*/
            Set<Integer> gammaY;
            if (p == 1) {
                gammaY = new HashSet<>();
                gammaY.add(0);
            } else {
                gammaY = identifyDirection(gradientFy, 2 * s);
            }
            /**7. , 8. Omega_x and Omega_y*/
            Set<Integer> omegaX = Sets.union(gammaX, Stat.supp(xi));
            Set<Integer> omegaY;
            if (p == 1) {
                omegaY = new HashSet<>();
                omegaY.add(0);
            } else {
                omegaY = Sets.union(gammaY, Stat.supp(yi));
            }
            //-------------------------------------------------------------------------------------------------

            if (verboseLevel == 0) {
                System.out.println("size of omegaX: " + omegaX.size());
                System.out.println("size of omegaY: " + omegaY.size());
            }

            if (verboseLevel > 0) {
                Collections.sort(head.bestForest.nodesInF);
                System.out.println("supp(xi):" + Stat.supp(xi));
                System.out.println("head  : " + head.bestForest.nodesInF.toString());
                System.out.println("gammaX: " + gammaX.toString());
                System.out.println("gammaY: " + gammaY.toString());
                System.out.println("omegaX: " + omegaX.toString());
                System.out.println("omegaY: " + omegaY.toString());
                System.out.println("suppX : " + Stat.supp(xi).toString());
                if (trueNodes != null && trueFeatures != null) {
                    System.out.println("Sx: " + trueNodes);
                    System.out.println("intersection nodes: " + Sets.intersection(omegaX, trueNodes));
                    System.out.println("intersection features: " + Sets.intersection(omegaY, trueFeatures));
                }
                Utils.stop();
            }

            //----------------------------------------------------------------------------------------
            /** to do minimization without constraint */
            /** 9.(b_x,b_y)= argMin f(x,y) s.t. supp(x) \in OmegaX and supp(y) \in OmegaY*/
            List<double[]> b;
            if (function instanceof LogLogisticRegression) {
                List<double[]> initialX0YO = new ArrayList<>();
                initialX0YO.add(xi);
                initialX0YO.add(yi);
                b = ((LogLogisticRegression) function).getArgMinFxy(omegaX, omegaY, initialX0YO);
            } else {
                b = function.getArgMinFxy(omegaX, omegaY);
            }
            double[] bx = b.get(0);
            double[] by;
            if (p == 1) {
                by = new double[]{1.0D};
            } else {
                by = b.get(1);
            }
            //----------------------------------------------------------------------------------------

            if (verboseLevel == 0) {
                System.out.println("-----after argmin ------");
                System.out.println("size of supp(by): " + Stat.supp(by));
                if (trueFeatures != null && trueNodes != null) {
                    System.out.println("**** after arg min f(x,y) ****");
                    System.out.println("supp(bx): \n" + Stat.supp(bx));
                    System.out.println("supp(by): " + Stat.supp(by));
                    System.out.println("bx: " + Doubles.join(" ", bx));
                    System.out.println("by: " + Doubles.join(" ", by));
                    System.out.println("PreRec of features after argmin: " + new PreRec(Stat.supp(by), trueFeatures));
                    System.out.println("PreRec of nodes after agrmin: " + new PreRec(Stat.supp(bx), trueNodes));
                }
            }

            //-----------------------------------------------------------------------------------------
            /** to do tail approximation ...*/
            TailApprox tail = new TailApprox(graph.edges, graph.edgeCosts, bx, k, g, B, false);
            psiX = new HashSet<>(tail.bestForest.nodesInF); /**10.\Psi_x^{i+1} tail approximation */
            psiY = identifyDirection(by, s); /**11. \Psi y = argmax ||\b_y|| : |R|<s*/
            xOld = xi;
            yOld = yi;
            xi = projectionOnVector(bx, psiX); /** 12. calculate x^{i+1} */
            if (p == 1) {
                yi = new double[]{1.0D};
            } else {
                yi = projectionOnVector(by, psiY); /** 13. calculate y^{i+1} */
            }
            //-----------------------------------------------------------------------------------------

            if (verboseLevel == 0) {
                System.out.println("number of nodes: " + psiX.size());
                System.out.println("number of features: " + psiY.size());
                if (trueNodes != null && trueFeatures != null) {
                    PreRec preRec = new PreRec(psiX, trueNodes);
                    System.out.println("preRec of nodes: " + preRec.toString());
                    preRec = new PreRec(psiY, trueFeatures);
                    System.out.println("preRec of features: " + preRec.toString());
                }
                Utils.stop();
            }

            //-------------------------------------------------------------
            /** to check the halting condition. */
            funcValue = function.getFuncValue(xi, yi);
            double gapX = l2Norm(xOld, xi);
            double gapY = l2Norm(yOld, yi);
            numOfIter++;
            if ((gapX < 1e-3D && gapY < 1e-3D) || numOfIter >= t) {
                break;
            }
            //-------------------------------------------------------------
            if (verboseLevel == 0) {
                System.out.println("gapX: " + gapX + ", gapY: " + gapY);
                PreRec preRec = new PreRec(psiX, this.trueNodes);
                System.out.println("nodes: " + preRec + " , features: " + new PreRec(psiY, this.trueFeatures));
                System.out.println("function value: " + this.funcValue);
            }
        }
        runTime = (System.nanoTime() - startTime) / 1e9;
        return true;
    }


    private double[] initializeX0(int option) {
        double[] x0 = new double[n];
        switch (option) {
            case 0:
                if (trueNodes == null) {
                    System.out.println("there is no true nodes ...");
                    System.exit(0);
                }
                /**option 1*/
                for (int i = 0; i < n; i++) {
                    if (trueNodes.contains(i)) {
                        x0[i] = new Random().nextDouble();
                    } else {
                        x0[i] = 0.0D;
                    }
                }
                break;
            case 1:
                /**option 2*/
                for (int i = 0; i < n; i++) {
                    x0[i] = new Random().nextDouble();
                }
                break;
            case 2:
                /**option 3: random select k non-zeros*/
                Set<Integer> indices = new HashSet<>();
                while (indices.size() != k) {
                    indices.add(new Random().nextInt(n));
                }
                Arrays.fill(x0, 0.0D);
                for (int node : indices) {
                    x0[node] = new Random().nextDouble();
                }
                break;
            case 3:
                /**option 4: random 1, 0 vector*/
                x0 = Stat.getRandomVector(n);
                break;
            case 4:
                if (trueX != null) {
                    x0 = trueX;
                } else {
                    System.out.println("trueX is null ...");
                    System.exit(0);
                }
                break;
            case 5:
                Arrays.fill(x0, 0.0D);
                x0[new Random().nextInt(x0.length)] = 1.0D;
                break;
            default:
                System.out.println("error.");
                System.exit(0);
        }
        return x0;
    }

    private double[] initializeY0(int option) {
        double[] y0 = new double[p];
        switch (option) {
            case 0:
                if (trueFeatures == null) {
                    System.out.println("there is no true features. ");
                    System.exit(0);
                }
                y0 = new double[p];
                for (int i = 0; i < p; i++) {
                    if (trueFeatures.contains(i)) {
                        y0[i] = new Random().nextDouble();
                    } else {
                        y0[i] = 0.0D;
                    }
                }
                break;
            case 1:
                for (int i = 0; i < p; i++) {
                    y0[i] = new Random().nextDouble();
                }
                break;
            case 2:
                /**option 3: random select s non-zeros*/
                Set<Integer> indices = new HashSet<>();
                while (indices.size() != s) {
                    indices.add(new Random().nextInt(p));
                }
                Arrays.fill(y0, 0.0D);
                for (int node : indices) {
                    y0[node] = new Random().nextDouble();
                }
                break;
            case 3:
                y0 = Stat.getRandomVector(p);
                break;
            case 4:
                if (trueY != null) {
                    y0 = trueY;
                } else {
                    System.out.println("trueY is null ...");
                    System.exit(0);
                }
                break;
            case 5:
                Arrays.fill(y0, 0.0D);
                y0[new Random().nextInt(y0.length)] = 1.0D;
                break;
            default:
                System.out.println("error.");
                System.exit(0);
        }
        return y0;
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
