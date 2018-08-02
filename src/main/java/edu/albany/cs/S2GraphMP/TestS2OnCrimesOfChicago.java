package edu.albany.cs.S2GraphMP;

import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.graph.CrimesOfChicagoGraph;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.Graphs;
import edu.albany.cs.scoreFuncs.EMSXYScore;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.Stat;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by baojian bzhou6@albany.edu on 2/8/17.
 */
public class TestS2OnCrimesOfChicago {

    private final String filePath;
    private final String outFilePath;
    private final Data data;
    private final int k;

    private final int s;
    private final double mu;


    public TestS2OnCrimesOfChicago(String filePath, String outFilePath, int k, int s, double mu) {
        this.filePath = filePath;
        this.outFilePath = outFilePath;
        this.mu = mu;
        data = readDataFromFile();
        test(data.graphMatrix);
        this.k = k;
        this.s = s;
        System.out.println(data.toString());
        if (!run()) {
            System.out.println("error");
            System.exit(0);
        }
    }

    private void test(double[][] matrix) {
        double meanAbnormal = 0.0D;
        double meanNormal = 0.0D;
        Set<Integer> random = Graphs.generateRandomSubset(data.trueFeatures.size(), data.p);
        double meanRandom = 0.0D;
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                if (data.trueNodes.contains(i) && data.trueFeatures.contains(j)) {
                    meanAbnormal += matrix[i][j];
                } else {
                    meanNormal += matrix[i][j];
                }
                if (random.contains(j) && data.trueNodes.contains(i)) {
                    meanRandom += matrix[i][j];
                }
            }
        }
        meanAbnormal = meanAbnormal / (data.trueFeatures.size() * 1.0D * data.trueNodes.size());
        meanNormal = meanNormal / (matrix.length * 1.0D * matrix[0].length - data.trueFeatures.size() * 1.0D * data.trueNodes.size());
        meanRandom = meanRandom / (data.trueFeatures.size() * 1.0D * data.trueNodes.size());
        System.out.println("mean of abnormal: " + meanAbnormal);
        System.out.println("mean of normal: " + meanNormal);
        System.out.println("mean of random: " + meanRandom);
        System.out.println("number of abnormal features: " + data.trueFeatures.size());
        System.out.println("number of abnormal nodes: " + data.trueNodes.size());
        //Utils.stop();
    }

    private boolean run() {
        Function func = new EMSXYScore(data.graphMatrix);
        ConnectedComponents cc = new ConnectedComponents(data.getAdj());
        System.out.println("isConnected: " + cc.checkConnectivity());
        Graph graph = new CrimesOfChicagoGraph(data.n, data.p, data.edges, data.edgeCosts, data.trueFeatures, data.trueNodes);
        double trueVal = func.getFuncValue(((CrimesOfChicagoGraph) graph).trueX, ((CrimesOfChicagoGraph) graph).trueY);
        System.out.println("true function value: " + trueVal);
        for (int i = 0; i < 20; i++) {
            Set<Integer> random = Graphs.generateRandomSubset(data.trueFeatures.size(), data.p);
            double[] y = new double[data.p];
            for (int j = 0; j < data.p; j++) {
                if (random.contains(j)) {
                    y[j] = 1.0D;
                } else {
                    y[j] = 0.0D;
                }
            }
            double val = func.getFuncValue(((CrimesOfChicagoGraph) graph).trueX, y);
            System.out.println("random val: " + val);
        }
        for (int i = 0; i < 20; i++) {
            Set<Integer> random = Graphs.generateRandomSubset(100, data.n);
            double[] x = new double[data.n];
            for (int j = 0; j < data.n; j++) {
                if (random.contains(j)) {
                    x[j] = 1.0D;
                } else {
                    x[j] = 0.0D;
                }
            }
            double val = func.getFuncValue(x, ((CrimesOfChicagoGraph) graph).trueY);
            System.out.println("random val: " + val);
        }
        //Utils.stop();
        S2GraphMPDebug s2GraphMP = new S2GraphMPDebug(graph, k, s, 1, func);
        PreRec preRec_nodes = new PreRec(s2GraphMP.psiX, data.trueNodes);
        PreRec preRec_features = new PreRec(s2GraphMP.psiY, data.trueFeatures);
        System.out.println("preRec of nodes: " + preRec_nodes);
        System.out.println("preRec of features: " + preRec_features);
        System.out.println("total running time: " + s2GraphMP.runTime);
        System.out.println("func value: " + s2GraphMP.funcValue);
        System.out.println("true function value: " + func.getFuncValue(((CrimesOfChicagoGraph) graph).trueX, ((CrimesOfChicagoGraph) graph).trueY));
        System.out.println("size of result nodes: " + s2GraphMP.psiX.size());
        System.out.println("size of result features: " + s2GraphMP.psiY.size());
        try {
            FileWriter fileWriter = new FileWriter(outFilePath, true);
            fileWriter.write(mu + " " + preRec_nodes.toFileString() + " " + preRec_features.toFileString() + " " + s2GraphMP.runTime + "\n");
            fileWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return true;
    }

    private Data readDataFromFile() {
        Data data = new Data();
        int index = 1;
        try {
            for (String line : Files.readAllLines(Paths.get(filePath))) {
                if (index == 1) {
                    data.n = Integer.parseInt(line.trim().split(" ")[0]);
                    data.p = Integer.parseInt(line.trim().split(" ")[1]);
                    data.graphMatrix = new double[data.n][data.p];
                    System.out.println(data.n + " " + data.p);
                    index++;
                } else if (index >= 2 && index <= (data.n + 1)) {
                    int i = 0;
                    for (String count : line.trim().split(" ")) {
                        data.graphMatrix[index - 2][i++] = Double.parseDouble(count);
                        if (data.maximalValInMat < Double.parseDouble(count)) {
                            data.maximalValInMat = Double.parseDouble(count);
                        }
                    }
                    index++;
                } else if ((index >= (data.n + 2)) && (index <= (data.n + 2))) {
                    data.numOfEdges = Integer.parseInt(line.trim());
                    data.edges = new ArrayList<>();
                    data.edgeCosts = new ArrayList<>();
                    System.out.println(data.numOfEdges);
                    index++;
                } else if ((index >= (data.n + 3)) && (index <= (data.n + 3 + data.numOfEdges - 1))) {
                    int endpoint0 = Integer.parseInt(line.trim().split(" ")[0]);
                    int endpoint1 = Integer.parseInt(line.trim().split(" ")[1]);
                    Integer[] edge = new Integer[]{endpoint0, endpoint1};
                    data.edges.add(edge);
                    data.edgeCosts.add(1.0D);
                    index++;
                } else if ((index >= (data.n + data.numOfEdges + 3)) && (index <= (data.n + data.numOfEdges + 3))) {
                    data.trueNodes = new HashSet<>();
                    for (String node : line.trim().split(" ")) {
                        data.trueNodes.add(Integer.parseInt(node));
                    }
                    index++;
                } else if ((index >= (data.n + data.numOfEdges + 4)) && (index <= (data.n + data.numOfEdges + 4))) {
                    data.trueFeatures = new HashSet<>();
                    for (String word : line.trim().split(" ")) {
                        data.trueFeatures.add(Integer.parseInt(word));
                    }
                    index++;
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        if (this.mu != 0.0D) {
            data = generateSimulationMatrix(data);
        }
//        double[] col = Matrices.getColumn(data.graphMatrix, 116);
//        data.graphMatrix = new double[data.n][1];
//        for (int i = 0; i < data.n; i++) {
//            data.graphMatrix[i][0] = col[i];
//        }
//        data.trueFeatures = new HashSet<>();
//        data.trueFeatures.add(0);
//        data.p = 1;
        //data = normalizeMatrixMethod1(data);
        //data.graphMatrix = normalizeMatrix(data.graphMatrix);
        return data;
    }

    private Data normalizeMatrixMethod1(Data data) {
        for (int i = 0; i < data.p; i++) {
            double[] columni = Stat.getColumn(i, data.graphMatrix);
            Stat stat = new Stat(columni);
            double mad = stat.mad();
            double mean = stat.mean();
            for (int j = 0; j < columni.length; j++) {
                data.graphMatrix[j][i] = (columni[j] - mean) / mad;
            }
        }
        return data;
    }

    private double[][] normalizeMatrixMethod2(double[][] mat) {
        double minVal = Double.MAX_VALUE;
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                if (mat[i][j] < minVal) {
                    minVal = mat[i][j];
                }
            }
        }
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                if (mat[i][j] < minVal) {
                    mat[i][j] = mat[i][j] - minVal;
                }
            }
        }
        double maxVal = Double.MIN_VALUE;
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                mat[i][j] = mat[i][j] - minVal;
                if (maxVal < mat[i][j]) {
                    maxVal = mat[i][j];
                }
            }
        }
        for (int i = 0; i < mat.length; i++) {
            for (int j = 0; j < mat[0].length; j++) {
                mat[i][j] = mat[i][j] / maxVal;
            }
        }
        return mat;
    }


    private Data generateSimulationMatrix(Data data) {
        int index = 0;
        NormalDistribution normal = new NormalDistribution(0.0D, 0.1D);
        NormalDistribution abnormal = new NormalDistribution(mu, 0.1D);
        System.out.println("to test the simulation data, where mu=" + mu);
        System.out.println("size of true nodes: " + data.trueNodes.size());
        System.out.println("size of true features: " + data.trueFeatures.size());
        for (int i = 0; i < data.n; i++) {
            for (int j = 0; j < data.p; j++) {
                if (data.trueNodes.contains(i) && data.trueFeatures.contains(j)) {
                    data.graphMatrix[i][j] = abnormal.sample() + data.graphMatrix[i][j];
                    index++;
                } else {
                    //data.graphMatrix[i][j] = normal.sample();
                }
            }
        }
        System.out.println("number of abnormal entries: " + index);
        return data;
    }


    class Data {
        int n;
        int p;
        double[][] graphMatrix;
        int numOfEdges;
        ArrayList<Integer[]> edges;
        ArrayList<Double> edgeCosts;
        Set<Integer> trueNodes;
        Set<Integer> trueFeatures;
        double maximalValInMat = 0.0D;

        @Override
        public String toString() {
            return "Data{" +
                    "n=" + n +
                    ", p=" + p +
                    ", graphMatrix=" + graphMatrix.length +
                    ", numOfEdges=" + numOfEdges +
                    ", edges=" + edges.size() +
                    ", edgeCosts=" + edgeCosts.size() +
                    ", trueNodes=" + trueNodes.size() +
                    ", trueFeatures=" + trueFeatures.size() +
                    '}';
        }

        public ArrayList<ArrayList<Integer>> getAdj() {
            ArrayList<ArrayList<Integer>> adj = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                adj.add(new ArrayList<>());
            }
            for (Integer[] edge : edges) {
                adj.get(edge[0]).add(edge[1]);
                adj.get(edge[1]).add(edge[0]);
            }
            return adj;
        }
    }

    public static void test() {
        String rootPath = "data/CrimesOfChicago/graph/";
        String filePath = rootPath + "processed_graph_2010.txt";
        final String outFilePath = "./outputs/CrimesOfChicago/result_2010_real.txt";
        int k = 2500;
        int s = 22;
        //final double[] muArr = new double[]{0.08D, 0.20D, 0.18D, 0.16D, 0.14D, 0.12D, 0.10D, 0.08D, 0.06D, 0.04D, 0.02D, 0.0D};
        final double[] muArr = new double[]{1.0D, 2.0D, 1.0D, 0.5D, 0.3D, 0.2D, 0.1D};
        ExecutorService executor = Executors.newFixedThreadPool(1);
        for (double mu : muArr) {
            executor.execute(new Thread() {
                public void run() {
                    new TestS2OnCrimesOfChicago(filePath, outFilePath, k, s, mu);
                }
            });
            break;
        }
        executor.shutdown();
        try {
            if (!executor.awaitTermination(24, TimeUnit.HOURS)) {
                executor.shutdown();//Cancel currently executing tasks.
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static void main(String args[]) {
        String rootPath = "data/CrimesOfChicago/graph/";
        String filePath = rootPath + "processed_ASSAULT_test_case_1.txt";
        final String outFilePath = "./outputs/CrimesOfChicago/result_2010_single_y.txt";
        int k = 107;//2269
        int s = 16;//
        new TestS2OnCrimesOfChicago(filePath, outFilePath, k, s, 0.0D);
    }

}
