package edu.albany.cs.S2GraphMP;

import edu.albany.cs.base.Constants;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.graph.BreastCancerGraph;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.LogLogisticRegression;
import edu.albany.cs.scoreFuncs.Stat;
import org.apache.commons.math3.stat.StatUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.net.InetAddress;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

/**
 * Created by baojian bzhou6@albany.edu on 2/4/17.
 */
public class TestS2OnBreastCancerData {

    private final int n;
    private int p;
    private final int k;
    private final int s;
    private final double[][] data;
    private Set<Integer> trueFeatureNodes;
    private double[] trueY;
    private Map<Integer, Integer> featuresMap;
    private Map<Integer, Integer> nodesMap;
    private Set<Integer> retainNodes;
    private final Graph graph;
    private final Function func;
    private final int verboseLevel = 0;
    private final int partialK;
    private final String serverName;

    public TestS2OnBreastCancerData(String filePath, int k, int s, int partialK, String serverName) {
        nodesMap = getNodesMap(filePath + "/entrez.csv");
        System.out.println("number of features(genes): " + nodesMap.size());
        ArrayList<Integer[]> edges = getEdges(filePath + "/edge.txt");
        trueFeatureNodes = getTrueFeatures(filePath + "/Y.csv");
        ArrayList<Double> edgeCosts = getEdgeCosts(edges);
        this.data = getDataMatrix(filePath + "/X.csv");
        this.n = data.length;
        this.p = data[0].length;
        this.k = k;
        this.s = s;
        this.partialK = partialK;
        this.serverName = serverName;
        this.graph = new BreastCancerGraph(this.n, this.p, edges, edgeCosts, trueFeatureNodes, trueY, this.partialK);
        this.func = new LogLogisticRegression(this.data);
        if (verboseLevel > 0) {
            System.out.println("number of nodes: " + n + " , number of features: " + p);
            System.out.println("number of edges in gene graph: " + edges.size());
            System.out.println("true features: " + ((BreastCancerGraph) graph).trueFeatures);
            for (int i = 0; i < p; i++) {
                double[] wi = ((LogLogisticRegression) func).getColumn(i);
                if (trueFeatureNodes.contains(i)) {
                    System.out.println("*: " + StatUtils.sum(wi));
                } else {
                    System.out.println(": " + StatUtils.sum(wi));
                }
            }
        }
        if (!Stat.isConnected(edges)) {
            System.exit(0);
        } else {
            System.out.println("the graph is connected.");
        }
        if (!run()) {
            System.out.println("error ...");
            System.exit(0);
        }
    }

    private boolean run() {
        S2GraphMPDebug s2GraphMP = new S2GraphMPDebug(graph, k, s, 10, func);
        PreRec preRec_features = new PreRec(s2GraphMP.psiY, trueFeatureNodes);
        try {
            FileWriter fileWriter = new FileWriter("outputs/Breast-Cancer-Data/result_hostname_" + serverName + ".txt", true);
            double pre = preRec_features.pre;
            double rec = preRec_features.rec;
            double fm = preRec_features.fmeasure;
            fileWriter.write(partialK + " " + pre + " " + rec + " " + fm + "\n");
            fileWriter.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println(preRec_features);
        return true;
    }

    private double[][] getDataMatrix(String filePath) {
        double[][] data = null;
        try {
            List<String> lines = Files.readAllLines(Paths.get(filePath));
            Iterator<String> iterator = lines.iterator();
            int featureIndex = 0;
            while (iterator.hasNext()) {
                String line = iterator.next();
                if (!retainNodes.contains(featureIndex)) {
                    if (verboseLevel > 0) {
                        System.out.println("removed feature: " + line);
                    }
                    iterator.remove();
                }
                featureIndex++;
            }
            data = new double[lines.size()][];
            featureIndex = 0;
            for (String line : lines) {
                int i = 0;
                String[] elements = line.trim().split(",");
                double[] col = new double[elements.length];
                for (String val : elements) {
                    col[i++] = Double.parseDouble(val);
                }
                data[featureIndex] = col;
                featureIndex++;
            }

        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("there are " + data.length + " features in data matrix.");
        System.out.println("there are " + data[0].length + " samples in data matrix.");
        return data;
    }

    private ArrayList<Integer[]> getEdges(String filePath) {
        ArrayList<Integer[]> edges = new ArrayList<>();
        Set<Integer> nodes = new HashSet<>();
        try {
            List<String> lines = Files.readAllLines(Paths.get(filePath));
            for (String line : lines) {
                int endPoint0 = Integer.parseInt(line.trim().split("\t")[0]) - 1;
                int endPoint1 = Integer.parseInt(line.trim().split("\t")[1]) - 1;
                if (endPoint0 == endPoint1) {
                    System.out.println("self loop edge");
                    System.exit(0);
                }
                nodes.add(endPoint0);
                nodes.add(endPoint1);
                edges.add(new Integer[]{endPoint0, endPoint1});
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        System.out.println("number of feature(genes) in edges: " + nodes.size());
        System.out.println("maximum feature id: " + Stat.getMaximum(nodes));
        System.out.println("minimum feature id: " + Stat.getMinimum(nodes));
        System.out.println("number of edges(before deletion): " + edges.size());
        retainNodes = Stat.getLargestCCNodes(edges);
        //To delete the edges which are not in the largest connected components
        Iterator<Integer[]> iterator = edges.iterator();
        int numEdgesDel = 0;
        while (iterator.hasNext()) {
            Integer[] edge = iterator.next();
            if ((!retainNodes.contains(edge[0])) || (!retainNodes.contains(edge[1]))) {
                iterator.remove();
                numEdgesDel++;
            }
        }
        System.out.println("delete " + numEdgesDel + " edges ");
        System.out.println("number of edges after deletion: " + edges.size());
        System.out.println("number of retain features(genes): " + retainNodes.size());
        List<Integer> listNodes = new ArrayList<>();
        for (Integer node : retainNodes) {
            listNodes.add(node);
        }
        assert retainNodes.size() == listNodes.size();
        for (int i = 0; i < edges.size(); i++) {
            Integer[] newEdge = new Integer[2];
            newEdge[0] = listNodes.indexOf(edges.get(i)[0]);
            newEdge[1] = listNodes.indexOf(edges.get(i)[1]);
            if (newEdge[0] == -1 || newEdge[1] == -1) {
                System.out.println(Arrays.toString(newEdge));
                System.exit(0);
            }
            edges.set(i, newEdge);
        }
        return edges;
    }

    private ArrayList<Double> getEdgeCosts(ArrayList<Integer[]> edges) {
        ArrayList<Double> edgeCosts = new ArrayList<>();
        for (Integer[] edge : edges) {
            edgeCosts.add(1.0D);
        }
        return edgeCosts;
    }

    private Set<Integer> getTrueFeatures(String filePath) {
        double[] trueFeatures = new double[p];
        Set<Integer> trueFeatureNodes = new HashSet<>();
        featuresMap = new HashMap<>();
        try {
            List<String> lines = Files.readAllLines(Paths.get(filePath));
            p = lines.size();
            trueFeatures = new double[lines.size()];
            for (int i = 0; i < lines.size(); i++) {
                Integer feature = Integer.parseInt(lines.get(i).trim().split(",")[0]);
                trueFeatures[i] = Double.parseDouble(lines.get(i).trim().split(",")[1]);
                featuresMap.put(feature, i);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        int positive = 0;
        int index = 0;
        for (double feature : trueFeatures) {
            if (feature > 0.0D) {
                positive++;
                trueFeatureNodes.add(index);
            }
            index++;
        }
        double[] y = new double[p];
        Arrays.fill(y, 0.0D);
        for (int node : trueFeatureNodes) {
            y[node] = 1.0D;
        }
        trueY = y;
        System.out.println("number of positive features (metastasis): " + positive);
        return trueFeatureNodes;
    }

    /**
     * This method is tested.
     *
     * @param filePath
     * @return the nodes map
     */
    private Map<Integer, Integer> getNodesMap(String filePath) {
        Map<Integer, Integer> result = new HashMap<>();
        try {
            List<String> lines = Files.readAllLines(Paths.get(filePath));
            for (int i = 0; i < lines.size(); i++) {
                Integer node = Integer.parseInt(lines.get(i).trim());
                result.put(node, i);
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
        return result;
    }


    public static void main(String args[]) {
        try {
            InetAddress ip = InetAddress.getLocalHost();
            System.out.println(ip.getHostAddress());
            System.out.println(ip.getHostName());
            ExecutorService pool = Executors.newFixedThreadPool(6);
            int k = 2000;//total sparsity of x
            int s = 100;//total sparsity of y
            final List<Integer> candidates = new ArrayList<>();
            for (int partialk = 30; partialk < 78; partialk += 5) {
                candidates.add(partialk);
            }
            String serverName = ip.getHostName();
            int index = 0;
            for (final int node : candidates) {
                if (serverName.equals(Constants.serverApdm01) && index % 2 == 0) {
                    pool.execute(new Thread() {
                        public void run() {
                            new TestS2OnBreastCancerData("data/Breast-Cancer-Data/", k, s, node, serverName);
                        }
                    });
                }
                if (serverName.equals(Constants.serverApdm02) && index % 2 == 1) {
                    pool.execute(new Thread() {
                        public void run() {
                            new TestS2OnBreastCancerData("data/Breast-Cancer-Data/", k, s, node, serverName);
                        }
                    });
                }
                if (serverName.equals(Constants.desktop) && index % 2 == 0) {
                    pool.execute(new Thread() {
                        public void run() {
                            new TestS2OnBreastCancerData("data/Breast-Cancer-Data/", k, s, node, serverName);
                        }
                    });
                }
                index++;
            }
            pool.shutdown();
            if (!pool.awaitTermination(24, TimeUnit.HOURS)) {
                pool.shutdown();//Cancel currently executing tasks.
            }

        } catch (IOException e) {
            e.printStackTrace();
        } catch (InterruptedException e) {
            e.printStackTrace();
        }

    }
}
