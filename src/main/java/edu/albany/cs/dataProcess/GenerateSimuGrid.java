package edu.albany.cs.dataProcess;

import com.google.common.primitives.Doubles;
import edu.albany.cs.base.Edge;
import edu.albany.cs.base.RandomWalk;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.Graphs;
import edu.albany.cs.graph.GridGraph;
import edu.albany.cs.scoreFuncs.GradientDescentOpt;
import edu.albany.cs.scoreFuncs.LogLogisticRegression;
import edu.albany.cs.scoreFuncs.Stat;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.distribution.UniformRealDistribution;
import org.apache.commons.math3.linear.ArrayRealVector;

import java.io.*;
import java.util.*;

/**
 * Created by baojian on 1/14/17.
 * Email: bzhou6@albany.edu
 */
public class GenerateSimuGrid implements Serializable {

    /**
     * number of nodes.
     */
    public int n;
    /**
     * number of features
     */
    public int p;
    /**
     * desired number of true nodes
     */
    public int z;
    /**
     * desired number of true features
     */
    public int r;
    /**
     * settings parameter for vector b
     */
    public double mu;
    /**
     * b vector to control the score function
     */
    public double[] b;
    /**
     * design matrix
     */
    public double[][] W;
    /**
     * The grid graph
     */
    public Graph graph;
    /**
     * The set of abnormal nodes
     */
    public Set<Integer> Sx;
    /**
     * Generate true subgraph using random walk
     */
    public RandomWalk trueSubGraph;
    /**
     * The set of abnormal features
     */
    public Set<Integer> Sy;
    /**
     * The vector of true parameters
     */
    public double[] x;
    /**
     * The vector of true features
     */
    public double[] y;

    transient private final int verboseLevel = 0;


    public GenerateSimuGrid(int n, int p, double z, int r, double mu) {
        this.W = new double[n][p];
        this.n = n;
        this.p = p;
        this.z = (int) z;
        this.r = r;
        this.mu = mu;
        run();
    }

    private void run() {

        graph = generateGridGraph(n);
        trueSubGraph = new RandomWalk(graph.arrayListAdj, z);
        Sx = trueSubGraph.nodesInSubGraph;
        Sy = Graphs.generateRandomSubset(r, p);
        x = generateX(Sx, 1.0D, 10.0D);
        y = getTrueY(Sy);
        algorithm3(Sx, Sy, mu, 1.0D);
        algorithm4(Sy, mu, 1.0D);
        b = generateB();
        LogLogisticRegression logistic = new LogLogisticRegression(W, b);

        double[] y = new double[p];
        Arrays.fill(y, 0.0D);
        for (int i = 0; i < p; i++) {
            if (Sy.contains(i)) {
                y[i] = 1.0D;
            }
        }

        if (verboseLevel > 0) {
            System.out.println("Sx: " + Sx.toString());
            System.out.println("Sy: " + Sy.toString());
            System.out.println("x: " + Doubles.join(" ", x));
            System.out.println("supp(x): " + Stat.supp(x));
            System.out.println("b: " + Doubles.join(" ", b));
            System.out.println("f(x,y): " + logistic.getFuncValue(x, y));
            for (int i = 0; i < p; i++) {
                ArrayRealVector xVector = new ArrayRealVector(x);
                ArrayRealVector wi = new ArrayRealVector(getColumn(i));
                if (Sy.contains(i)) {
                    System.out.println("xTransWi *: " + xVector.dotProduct(wi));
                } else {
                    System.out.println("xTransWi: " + xVector.dotProduct(wi));
                }
            }
            Utils.stop();
        }

        GradientDescentOpt gradientDescent = new GradientDescentOpt(logistic, n, p, 0.01D, 1000000);
        List<double[]> results = gradientDescent.simpleGradientDescent(Sx, Sy, x);

        for (int i = 0; i < p; i++) {
            ArrayRealVector xx = new ArrayRealVector(this.x);
            ArrayRealVector ww = new ArrayRealVector(getColumn(i));
            double val = xx.dotProduct(ww) + b[i];
            if (Sy.contains(i)) {
                System.out.println("*: " + val + "sigmoid funcVal: " + sigmoidFuncVal(val));
            } else {
                System.out.println(": " + val + " sigmoid funcVal: " + sigmoidFuncVal(val));
            }
        }
        System.out.println("supp(x0): " + Stat.supp(this.x));
        System.out.println("f(x0,y): " + logistic.getFuncValue(this.x, y));
        System.out.println("--------------------------------------------------------------------");
        x = results.get(0);
        for (int i = 0; i < p; i++) {
            ArrayRealVector xx = new ArrayRealVector(this.x);
            ArrayRealVector ww = new ArrayRealVector(getColumn(i));
            double val = xx.dotProduct(ww) + b[i];
            if (Sy.contains(i)) {
                System.out.println("*: " + val + "sigmoid funcVal: " + sigmoidFuncVal(val));
            } else {
                System.out.println(": " + val + " sigmoid funcVal: " + sigmoidFuncVal(val));
            }
        }
        System.out.println("supp(x0): " + Stat.supp(this.x));
        System.out.println("f(x*,y): " + logistic.getFuncValue(this.x, y));
    }


    private double[] getTrueY(Set<Integer> Sy) {
        double[] y = new double[p];
        for (int i = 0; i < p; i++) {
            if (Sy.contains(i)) {
                y[i] = 1.0D;
            } else {
                y[i] = 0.0D;
            }
        }
        return y;
    }

    private double[] generateX(Set<Integer> Sx, double lowerBound, double upperBound) {
        UniformRealDistribution uniform = new UniformRealDistribution(lowerBound, upperBound);
        ArrayRealVector x = new ArrayRealVector(n);
        for (int i = 0; i < n; i++) {
            if (!Sx.contains(i)) {
                x.setEntry(i, 0.0D);
            } else {
                x.setEntry(i, uniform.sample());
            }
        }
        return x.toArray();
    }

    private double[] generateB() {
        ArrayRealVector b = new ArrayRealVector(p);
        for (int i = 0; i < p; i++) {
            if (Sy.contains(i)) {
                b.setEntry(i, 0.0D);
            } else {
                b.setEntry(i, 0.0D);
            }
        }
        return b.toArray();
    }

    private void setWColumn(ArrayRealVector col, int index) {
        for (int i = 0; i < n; i++) {
            W[i][index] = col.getEntry(i);
        }
    }

    public double[] getColumn(int index) {
        double[] col = new double[n];
        for (int i = 0; i < n; i++) {
            col[i] = W[i][index];
        }
        return col;
    }

    private void algorithm3(Set<Integer> Sx, Set<Integer> Sy, double mu, double std) {
        NormalDistribution gaussianAbnormal = new NormalDistribution(mu, std);
        NormalDistribution gaussianNormal = new NormalDistribution(-mu, std);
        for (Integer node : Sy) {
            ArrayRealVector r;
            r = new ArrayRealVector(n);
            for (int j = 0; j < n; j++) {
                if (Sx.contains(j)) {
                    r.setEntry(j, gaussianAbnormal.sample());
                } else {
                    r.setEntry(j, gaussianNormal.sample());
                }
            }
            setWColumn(r, node);
        }

    }

    private void algorithm4(Set<Integer> Sy, double mu, double std) {
        NormalDistribution gaussianNormal = new NormalDistribution(-mu, std);
        for (int i = 0; i < p; i++) {
            if (!Sy.contains(i)) {
                ArrayRealVector r = new ArrayRealVector(n);
                for (int k = 0; k < n; k++) {
                    r.setEntry(k, gaussianNormal.sample());
                }
                setWColumn(r, i);
            }
        }
    }

    private double sigmoidFuncVal(double x) {
        return 1.0D / (1.0D + Math.exp(-x));
    }

    private Graph generateGridGraph(int n) {
        return Graphs.generateGridGraph((int) Math.sqrt(n), (int) Math.sqrt(n));
    }

    public static GenerateSimuGrid getGridObj(String filePath) {
        GenerateSimuGrid generateSimuGrid = null;
        try {
            FileInputStream fileIn = new FileInputStream(filePath);
            ObjectInputStream in = new ObjectInputStream(fileIn);
            generateSimuGrid = (GenerateSimuGrid) in.readObject();
            in.close();
            fileIn.close();
        } catch (IOException i) {
            i.printStackTrace();
        } catch (ClassNotFoundException c) {
            System.out.println("GenerateSimuGrid class not found");
            c.printStackTrace();
        }
        return generateSimuGrid;
    }

    public static void testOnSerializedObj() {
        String filePath = "./data/SerializedData/Grid_np_900_10_z_9.0_r_1_mu_1_.ser";
        GenerateSimuGrid generateSimuGrid = GenerateSimuGrid.getGridObj(filePath);
        System.out.println("number of nodes in graph: " + generateSimuGrid.graph.numOfNodes);
        System.out.println("number of true nodes: " + generateSimuGrid.trueSubGraph.numOfNodesInT);
        for (Edge edge : generateSimuGrid.trueSubGraph.subGraph) {
            System.out.println(edge.toString());
        }
    }

    public static void generateGridData(String folder) {
        int n = 900;
        int p = 200;
        int index = 0;
        for (double z : new double[]{0.05D * n}) {
            for (int r : new int[]{(int) (0.05D * p)}) {
                for (double mu : new double[]{1.0D, 2.0D, 3.0D, 4.0D, 5.0D, 10.0D, 15.0D, 20.0D}) {
                    System.out.println("simulation " + (index++));
                    GenerateSimuGrid grid = new GenerateSimuGrid(n, p, z, r, mu);
                    try {
                        String filePath = folder + "Grid_np_" + n + "_" + p + "_z_" + z + "_r_" + r + "_mu_" + mu + "_.ser";
                        FileOutputStream fileOut = new FileOutputStream(filePath);
                        ObjectOutputStream out = new ObjectOutputStream(fileOut);
                        out.writeObject(grid);
                        out.close();
                        fileOut.close();
                        System.out.println("Serialized data is saved in " + filePath);
                        GenerateSimuGrid generateSimuGrid = GenerateSimuGrid.getGridObj(filePath);
                        if (true) {
                            System.out.println("p: " + generateSimuGrid.p);
                            System.out.println("n: " + generateSimuGrid.n);
                            System.out.println("b: " + ArrayUtils.toString(generateSimuGrid.b));
                            System.out.println("x: " + Doubles.join(" ", generateSimuGrid.x));
                        }
                    } catch (IOException i) {
                        i.printStackTrace();
                    }
                }
            }
        }
    }

    public static void main(String args[]) {
        String folder = "./data/SerializedData/";
        if (!new File(folder).exists()) {
            new File(folder).mkdir();
        }
        for (int i = 0; i < 20; i++) {
            folder = "./data/SerializedData/dataset" + i + "/";
            if (!new File(folder).exists()) {
                new File(folder).mkdir();
            }
            generateGridData(folder);
            Utils.stop();
        }
    }
}
