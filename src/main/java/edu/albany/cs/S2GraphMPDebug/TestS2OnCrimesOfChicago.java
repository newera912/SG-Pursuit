package edu.albany.cs.S2GraphMPDebug;

import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.graph.CrimesOfChicagoGraph;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.Graphs;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
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
    private final double lambda = 10.0D;
    private final int s;


    public TestS2OnCrimesOfChicago(String filePath, String outFilePath, String keyword) {
        this.filePath = filePath;
        this.outFilePath = outFilePath;
        data = readDataFromFile();
        this.k = data.trueNodes.size() / 2;
        this.s = data.trueFeatures.size();
        System.out.println(data.toString());
        if (!run(keyword)) {
            System.out.println("error");
            System.exit(0);
        }
    }

    private boolean run(String keyword) {
        Function func = new ElevatedMeanScan(data.graphMatrix, lambda);
        ConnectedComponents cc = new ConnectedComponents(data.getAdj());
        System.out.println("isConnected: " + cc.checkConnectivity());
        Graph graph = new CrimesOfChicagoGraph(data.n, data.p, data.edges, data.edgeCosts, data.trueFeatures, data.trueNodes);
        double trueVal = func.getFuncValue(((CrimesOfChicagoGraph) graph).trueX, ((CrimesOfChicagoGraph) graph).trueY);
        System.out.println("true function value: " + trueVal);
        S2GraphMPDebug s2GraphMP = new S2GraphMPDebug(graph, k, s, 1, func, null);
        PreRec preRec_nodes = new PreRec(s2GraphMP.psiX, data.trueNodes);
        PreRec preRec_features = new PreRec(s2GraphMP.psiY, data.trueFeatures);
        System.out.println("preRec of nodes: " + preRec_nodes);
        System.out.println("preRec of features: " + preRec_features);
        System.out.println("total running time: " + s2GraphMP.runTime);
        System.out.println("func value: " + s2GraphMP.funcValue);
        System.out.println("true function value: " + func.getFuncValue(((CrimesOfChicagoGraph) graph).trueX, ((CrimesOfChicagoGraph) graph).trueY));
        System.out.println("size of result nodes: " + s2GraphMP.psiX.size());
        System.out.println("true nodes: " + ((CrimesOfChicagoGraph) graph).trueNodes);
        System.out.println("result nodes: " + s2GraphMP.psiX);
        System.out.println("size of result features: " + s2GraphMP.psiY.size());
        try {
            FileWriter fileWriter = new FileWriter(outFilePath, true);
            fileWriter.write(keyword + " " + preRec_nodes.toFileString() + " " + preRec_features.toFileString() + " ");
            fileWriter.write(s2GraphMP.psiX.size() + " " + s2GraphMP.psiY.size() + " ");
            fileWriter.write(s2GraphMP.funcValue + " " + s2GraphMP.runTime + "\n");
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


    public static void testTheft() {
        ExecutorService pool = Executors.newFixedThreadPool(5);
        String rootPath = "/home/baojian/Dropbox/expriments/CrimesOfChicago/graph/";
        final String[] keywords = new String[]{"THEFT_test_case_0", "THEFT_test_case_1",
                "THEFT_test_case_2", "THEFT_test_case_3", "THEFT_test_case_4"};
        for (int i : new int[]{0, 1, 2, 3, 4}) {
            pool.execute(new Thread() {
                public void run() {
                    String filePath = rootPath + "processed_" + keywords[i] + ".txt";
                    final String outFilePath = "./outputs/CrimesOfChicago/final_result_S2GraphMP.txt";
                    new TestS2OnCrimesOfChicago(filePath, outFilePath, keywords[i]);
                }
            });
        }
        pool.shutdown();
        try {
            if (!pool.awaitTermination(24, TimeUnit.HOURS)) {
                pool.shutdown();//Cancel currently executing tasks.
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static void testBurglary() {
        ExecutorService pool = Executors.newFixedThreadPool(5);
        String rootPath = "/home/baojian/Dropbox/expriments/CrimesOfChicago/graph/";
        final String[] keywords = new String[]{"BURGLARY_test_case_0", "BURGLARY_test_case_1",
                "BURGLARY_test_case_2", "BURGLARY_test_case_3", "BURGLARY_test_case_4"};
        for (int i : new int[]{0, 1, 2, 3, 4}) {
            pool.execute(new Thread() {
                public void run() {
                    String filePath = rootPath + "processed_" + keywords[i] + ".txt";
                    final String outFilePath = "./outputs/CrimesOfChicago/final_result_S2GraphMP.txt";
                    new TestS2OnCrimesOfChicago(filePath, outFilePath, keywords[i]);
                }
            });
        }
        pool.shutdown();
        try {
            if (!pool.awaitTermination(24, TimeUnit.HOURS)) {
                pool.shutdown();//Cancel currently executing tasks.
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static void testBattery() {
        ExecutorService pool = Executors.newFixedThreadPool(5);
        String rootPath = "/home/baojian/Dropbox/expriments/CrimesOfChicago/graph/";
        final String[] keywords = new String[]{"BATTERY_test_case_0", "BATTERY_test_case_1",
                "BATTERY_test_case_2", "BATTERY_test_case_3", "BATTERY_test_case_4"};
        for (int i : new int[]{0, 1, 2, 3, 4}) {
            pool.execute(new Thread() {
                public void run() {
                    String filePath = rootPath + "processed_" + keywords[i] + ".txt";
                    final String outFilePath = "./outputs/CrimesOfChicago/final_result_S2GraphMP.txt";
                    new TestS2OnCrimesOfChicago(filePath, outFilePath, keywords[i]);
                }
            });
        }
        pool.shutdown();
        try {
            if (!pool.awaitTermination(24, TimeUnit.HOURS)) {
                pool.shutdown();//Cancel currently executing tasks.
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static void wholeTest() {
        ExecutorService pool = Executors.newFixedThreadPool(5);
        String rootPath = "/home/baojian/Dropbox/expriments/CrimesOfChicago/graph/";
        final List<String> filePaths = new ArrayList<>();
        final List<String> keywords = new ArrayList<>();
        final List<Integer> indices = new ArrayList<>();
        int index = 0;
        for (File file : new File(rootPath).listFiles()) {
            filePaths.add(file.getAbsolutePath());
            keywords.add(file.getName().split("_")[1]);
            indices.add(index++);
        }
        for (final Integer i : indices) {
            pool.execute(new Thread() {
                public void run() {
                    String keyword = keywords.get(i);
                    String filePath = filePaths.get(i);
                    final String outFilePath = "./outputs/CrimesOfChicago/final_result_S2GraphMP_fix.txt";
                    new TestS2OnCrimesOfChicago(filePath, outFilePath, keyword);
                }
            });
        }
        pool.shutdown();
        try {
            if (!pool.awaitTermination(24, TimeUnit.HOURS)) {
                pool.shutdown();//Cancel currently executing tasks.
            }
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
    }

    public static void main(String args[]) {
        wholeTest();
    }
}
