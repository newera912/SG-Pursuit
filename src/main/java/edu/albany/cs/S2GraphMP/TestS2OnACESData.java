package edu.albany.cs.S2GraphMP;

import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.DisjointSet;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.ACEGeneGraph;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.scoreFuncs.LogLogisticRegression;
import edu.albany.cs.scoreFuncs.Stat;

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.*;

/**
 * Created by baojian bzhou6@albany.edu on 2/4/17.
 */
public class TestS2OnACESData {

    private final int n;
    private final int p;
    private final List<String> patientLabels;
    private final List<Boolean> patientClassLabels;
    private final List<String> geneLabels;
    private final Graph graph;
    private final LogLogisticRegression func;

    private final double[][] expressionData;
    private final int verboseLevel = 0;

    public TestS2OnACESData(String filePath) {
        List<String> lines = getAllLines(filePath);
        int index = 0;
        int n = 0;
        int p = 0;
        double[][] data = null;
        ArrayList<Integer[]> edges = new ArrayList<>();
        ArrayList<Double> edgeCosts = new ArrayList<>();
        patientLabels = new ArrayList<>();
        patientClassLabels = new ArrayList<>();
        geneLabels = new ArrayList<>();
        Set<String> hashSet = new HashSet<>();
        for (String str : lines) {
            if (index == 0) {
                n = Integer.parseInt(str.trim().split(" ")[0]);
                p = Integer.parseInt(str.trim().split(" ")[1]);
                data = new double[n][p];
            } else if (index >= 1 && index <= (p)) {
                patientLabels.add(str.trim().split(" ")[0]);
                if (str.trim().split(" ")[1].equals("True")) {
                    patientClassLabels.add(true);
                } else {
                    patientClassLabels.add(false);
                }
            } else if (index >= (p + 1) && index <= (p + n)) {
                geneLabels.add(str.trim());
            } else if (index >= (p + n + 1) && index <= (p + 2 * n)) {
                int i = 0;
                for (String val : str.trim().split(" ")) {
                    data[index - (p + n + 1)][i++] = Double.parseDouble(val);
                }
            } else if (index >= (p + n + 1)) {
                hashSet.add(str.trim().split(" ")[0]);
                hashSet.add(str.trim().split(" ")[1]);
                int endPoint0 = geneLabels.indexOf(str.trim().split(" ")[0]);
                int endPoint1 = geneLabels.indexOf(str.trim().split(" ")[1]);
                if (endPoint0 == -1 || endPoint1 == -1) {
                    System.out.println("geneLabel does not exist!");
                    System.out.println(str);
                    System.out.println("number of gene labels: " + geneLabels.size());
                    //System.exit(0);
                }
                edges.add(new Integer[]{endPoint0, endPoint1});
                edgeCosts.add(1.0D);
            }
            index++;
        }
        this.n = n;
        this.p = p;
        this.expressionData = data;
        this.graph = new ACEGeneGraph(this.n, this.p, edges, edgeCosts);
        this.func = new LogLogisticRegression(this.expressionData);
        if (verboseLevel == 0) {
            System.out.println("number of nodes: " + n + " , number of features: " + p);
            System.out.println("number of patients: " + patientLabels.size());
            System.out.println("number of genes: " + geneLabels.size());
            System.out.println("number of edges in gene graph: " + edges.size());
            System.out.println("number of nodes in edges: "+ hashSet.size());
        }
        Utils.stop();
        if (!Stat.isConnected(edges)) {
            System.exit(0);
        }
        if (!run()) {
            System.out.println("error ...");
            System.exit(0);
        }
    }

    private boolean run() {
        System.out.println();
        return true;
    }

    private List<String> getAllLines(String filePath) {
        try {
            List<String> allLines = Files.readAllLines(Paths.get(filePath));
            return allLines;
        } catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }


    public static void main(String args[]) {
        String filePath = "/home/baojian/Dropbox/expriments/ACES_Data/data/U133A_combat_DMFS_nwEdgesKEGG_data.txt";
        new TestS2OnACESData(filePath);
    }
}
