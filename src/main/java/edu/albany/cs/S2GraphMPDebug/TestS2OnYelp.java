package edu.albany.cs.S2GraphMPDebug;

import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.YelpGraph;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by baojian on 2/17/17.
 */
public class TestS2OnYelp {


    private final String inputFilePath;
    private final String outFilePath;
    private final Data data;
    private final int k;
    private final int s;
    private final double lambda = 10.0D;
    private final String keyword;

    public TestS2OnYelp(String inputFilePath, String outFilePath, String keyword) {
        this.inputFilePath = inputFilePath;
        this.outFilePath = outFilePath;
        this.data = generateDataFromFile();
        this.k = data.k / 2;
        this.s = 10;
        this.keyword = keyword;
        if (!run()) {
            System.out.println("error...");
            System.exit(0);
        }
    }


    private boolean run() {
        Function func = new ElevatedMeanScan(data.graphMatrix, lambda, data);
        ConnectedComponents cc = new ConnectedComponents(data.getAdj());
        System.out.println("isConnected: " + cc.checkConnectivity());
        Graph graph = new YelpGraph(data.n, data.p, data.edges, data.edgeCosts, data.trueFeatures, data.trueNodes);
        double trueVal = func.getFuncValue(((YelpGraph) graph).trueX, ((YelpGraph) graph).trueY);
        System.out.println("true function value: " + trueVal);
        S2GraphMPDebug s2GraphMP = new S2GraphMPDebug(graph, k, s, 1, func, data);
        PreRec preRec_nodes = new PreRec(s2GraphMP.psiX, data.trueNodes);
        PreRec preRec_features = new PreRec(s2GraphMP.psiY, data.trueFeatures);
        System.out.println("preRec of nodes: " + preRec_nodes);
        System.out.println("preRec of features: " + preRec_features);
        System.out.println("total running time: " + s2GraphMP.runTime);
        System.out.println("func value: " + s2GraphMP.funcValue);
        System.out.println("true function value: " + func.getFuncValue(((YelpGraph) graph).trueX, ((YelpGraph) graph).trueY));
        System.out.println("size of result nodes: " + s2GraphMP.psiX.size());
        System.out.println("true nodes: " + ((YelpGraph) graph).trueNodes);
        System.out.println("result nodes: " + s2GraphMP.psiX);
        System.out.println("size of result features: " + s2GraphMP.psiY.size());
        List<String> words = getSortedWords(s2GraphMP);
        Utils.appendRecord(outFilePath, true, keyword, preRec_nodes, preRec_features, s2GraphMP.psiX.size(),
                s2GraphMP.psiY.size(), s2GraphMP.funcValue, s2GraphMP.runTime, words);
        for (String word : getSortedWords(s2GraphMP)) {
            System.out.print(word + " ");
        }
        return true;
    }

    private List<String> getSortedWords(S2GraphMPDebug s2GraphMP) {
        ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(s2GraphMP.yi);
        Integer[] indexes = arrayIndexComparator.getIndices();
        Arrays.sort(indexes, arrayIndexComparator);
        List<String> words = new ArrayList<>();
        Map<Integer, String> reversedMap = new HashMap<>();
        for (String word : data.wordsDict.keySet()) {
            int index = data.wordsDict.get(word);
            reversedMap.put(index, word);
        }
        for (int i = 0; i < s; i++) {
            words.add(reversedMap.get(indexes[i]));
            System.out.println(s2GraphMP.yi[indexes[i]]);
        }
        return words;
    }


    private Data generateDataFromFile() {
        Data data = new Data();
        try {
            int index = 0;
            int graphMatrixRowIndex = 0;
            for (String eachLine : Files.readAllLines(Paths.get(inputFilePath))) {
                String[] items = eachLine.trim().split(" ");
                if (index == 0) {
                    data.n = Integer.parseInt(items[0]);
                    data.p = Integer.parseInt(items[1]);
                    data.k = Integer.parseInt(items[2]);
                    data.s = Integer.parseInt(items[3]);
                    data.numOfEdges = Integer.parseInt(items[4]);
                    data.graphMatrix = new double[data.n][];
                    data.nodesDict = new HashMap<>();
                    data.wordsDict = new HashMap<>();
                    data.trueNodes = new HashSet<>();
                    data.reverseWordsDict = new HashMap<>();
                    data.trueFeatures = new HashSet<>();
                    data.edgeCosts = new ArrayList<>();
                    data.edges = new ArrayList<>();
                } else if (index >= 1 && index <= data.n) {
                    double[] row = new double[data.p];
                    for (int i = 0; i < data.p; i++) {
                        row[i] = Double.parseDouble(items[i]);
                    }
                    data.graphMatrix[graphMatrixRowIndex++] = row;
                } else if (index >= (data.n + 1) && index <= (2 * data.n)) {
                    data.nodesDict.put(items[0], Integer.parseInt(items[1]));
                } else if (index >= (2 * data.n + 1) && index <= (2 * data.n + data.p)) {
                    data.wordsDict.put(items[0], Integer.parseInt(items[1]));
                    data.reverseWordsDict.put(Integer.parseInt(items[1]), items[0]);
                } else if (index == (2 * data.n + data.p + 1)) {
                    for (String item : items) {
                        data.trueNodes.add(Integer.parseInt(item));
                    }
                } else if (index == (2 * data.n + data.p + 2)) {
                    for (String item : items) {
                        data.trueFeatures.add(Integer.parseInt(item));
                    }
                } else if (index >= (2 * data.n + data.p + 3) &&
                        index <= (2 * data.n + data.p + 3 + data.numOfEdges)) {
                    int endPoint0 = Integer.parseInt(items[0]);
                    int endPoint1 = Integer.parseInt(items[1]);
                    Integer[] edge = new Integer[]{endPoint0, endPoint1};
                    data.edges.add(edge);
                    data.edgeCosts.add(1.0D);
                }
                index++;
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.format("%d %d %d %d %d\n", data.n, data.p, data.k, data.s, data.numOfEdges);
        System.out.println("finish to load data ...");
        return data;
    }

    public static void wholeTest() {
        ExecutorService pool = Executors.newFixedThreadPool(2);
        String username = System.getProperty("user.name");
        String rootPath = "/home/" + username + "/Dropbox/expriments/Yelp/graph/";
        final List<String> filePaths = new ArrayList<>();
        final List<String> keywords = new ArrayList<>();
        final List<Integer> indices = new ArrayList<>();
        int index = 0;
        for (File file : new File(rootPath).listFiles()) {
            filePaths.add(file.getAbsolutePath());
            keywords.add("2014_2015");
            indices.add(index++);
        }
        for (final Integer i : indices) {
            pool.execute(new Thread() {
                public void run() {
                    String keyword = keywords.get(i);
                    String filePath = filePaths.get(i);
                    String outputFilePath = "./outputs/Yelp/result_s2GraphMP_yelp.txt";
                    new TestS2OnYelp(filePath, outputFilePath, keyword);
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

    class Data {
        int n;
        int p;
        int k;
        int s;
        double[][] graphMatrix;
        int numOfEdges;
        Map<String, Integer> nodesDict;
        Map<String, Integer> wordsDict;
        Map<Integer, String> reverseWordsDict;
        ArrayList<Integer[]> edges;
        ArrayList<Double> edgeCosts;
        Set<Integer> trueNodes;
        Set<Integer> trueFeatures;


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
}
