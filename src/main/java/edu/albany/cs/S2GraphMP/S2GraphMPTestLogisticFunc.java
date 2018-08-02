package edu.albany.cs.S2GraphMP;

import com.google.common.primitives.Doubles;
import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;
import edu.albany.cs.dataProcess.GenerateSimuGrid;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.LogLogisticRegression;
import edu.albany.cs.scoreFuncs.LogisticRegression;
import edu.albany.cs.scoreFuncs.Stat;
import org.apache.commons.math3.linear.ArrayRealVector;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;


/**
 * Created by baojian on 1/17/17.
 * Email: bzhou6@albany.edu
 */
public class S2GraphMPTestLogisticFunc {

    private final int k;
    private final int s;
    private final String filePath;
    private final int verboseLevel = 0;

    public PreRec preRec_nodes;
    public PreRec preRec_features;


    public S2GraphMPTestLogisticFunc(int k, int s, String filePath) {
        this.k = k;
        this.s = s;
        this.filePath = filePath;
        run();
    }

    public void run() {
        GenerateSimuGrid gridData = GenerateSimuGrid.getGridObj(filePath);
        Graph graph = gridData.graph;
        Function logisticRegression = new LogisticRegression(gridData.W, gridData.b);
        ConnectedComponents cc = new ConnectedComponents(gridData.graph.arrayListAdj);
        int numOfCC = cc.computeCCSubGraph(gridData.trueSubGraph.nodesInSubGraph);
        if (verboseLevel > 0) {
            System.out.println("maximum subgraph size: " + k);
            System.out.println("maximum size of selected features: " + s);
            int trueSize = gridData.trueSubGraph.nodesInSubGraph.size();
            System.out.println("number of true subgraph nodes: " + trueSize);
            System.out.println("true nodes: " + gridData.Sx.toString());
            System.out.println("true features: " + gridData.Sy.toString());
            System.out.println("number of connected components in true subgraph: " + numOfCC);
        }
        S2GraphMP s2GraphMP = new S2GraphMP(graph, k, s, 1, logisticRegression);
        preRec_nodes = new PreRec(s2GraphMP.psiX, gridData.Sx);
        preRec_features = new PreRec(s2GraphMP.psiY, gridData.Sy);
        if (verboseLevel == 0) {
            System.out.println("return nodes: " + s2GraphMP.psiX.toString());
            System.out.println("true nodes: " + gridData.Sx.toString());
            System.out.println(preRec_nodes.toString());
            System.out.println(preRec_features.toString());
        }
    }


    @Override
    public String toString() {
        return "S2GraphMPTestLogisticFunc{" +
                "k=" + k +
                ", s=" + s +
                ", filePath='" + filePath + '\'' +
                ", verboseLevel=" + verboseLevel +
                ", preRec_nodes=" + preRec_nodes +
                ", preRec_features=" + preRec_features +
                '}';
    }

    public static void testSingleFile() {
        String filePath = "./data/SerializedData/Grid_np_100_10_z_20.0_r_5_mu_10_.ser";
        S2GraphMPTestLogisticFunc s2Best = null;
        for (int k = 5; k < 30; k++) {
            for (int s = 2; s < 8; s++) {
                S2GraphMPTestLogisticFunc s2 = new S2GraphMPTestLogisticFunc(k, s, filePath);
                if (s2Best == null) {
                    s2Best = s2;
                } else {
                    double sumFmeasure = s2.preRec_features.fmeasure + s2.preRec_nodes.fmeasure;
                    double sumFmeasure_best = s2Best.preRec_features.fmeasure + s2Best.preRec_nodes.fmeasure;
                    if (sumFmeasure > sumFmeasure_best) {
                        s2Best = s2;
                    }
                }
            }
        }
        System.out.println(s2Best.toString());
    }


    public static void testGroup() {
        try {
            List<PreRec> bestPreRec = null;
            int bestK = 0;
            FileWriter resultFile = new FileWriter("outputs/SerializedData/output.txt");
            for (File file : new File("./data/SerializedData/").listFiles()) {
                GenerateSimuGrid gridData = GenerateSimuGrid.getGridObj("./data/SerializedData/" + file.getName());
                for (int k = 9; k < 15; k++) {
                    System.out.println("-------------------------------------------------------------------");
                    List<PreRec> preRec = testFix(k, 5, "./data/SerializedData/Grid_np_100_10_z_20.0_r_5_mu_20_.ser");
                    if (bestPreRec == null) {
                        bestPreRec = preRec;
                        bestK = k;
                    }
                    if (bestPreRec.get(0).fmeasure < preRec.get(0).fmeasure) {
                        bestPreRec = preRec;
                        bestK = k;
                    }
                }
                resultFile.write(gridData.mu + " " + gridData.Sx.size() + " " + gridData.Sy.size() + " ");
                resultFile.write(bestPreRec.get(0).pre + " " + bestPreRec.get(0).rec + " " + bestPreRec.get(0).fmeasure + "\n");
            }
            resultFile.close();
            System.out.println("======================================");
            System.out.println(" " + bestPreRec.toString());
            System.out.println("bestK: " + bestK);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static List<PreRec> testFix(int k, int s, String filePath) {
        int verboseLevel = 0;
        GenerateSimuGrid gridData = GenerateSimuGrid.getGridObj(filePath);
        Graph graph = gridData.graph;
        Function logisticRegression = new LogLogisticRegression(gridData.W, gridData.b);
        ConnectedComponents cc = new ConnectedComponents(gridData.graph.arrayListAdj);
        int numOfCC = cc.computeCCSubGraph(gridData.Sx);
        if (verboseLevel > 0) {
            System.out.println("maximum subgraph size: " + k);
            System.out.println("maximum size of selected features: " + s);
            int trueSize = gridData.trueSubGraph.nodesInSubGraph.size();
            System.out.println("number of true subgraph nodes: " + trueSize);
            System.out.println("true x: " + Doubles.join(" ", gridData.x));
            System.out.println("size of supp(x): " + Stat.supp(gridData.x).size());
            System.out.println("supp(x): " + Stat.supp(gridData.x));
            System.out.println("true nodes: " + gridData.Sx.toString());
            System.out.println("true features: " + gridData.Sy.toString());
            System.out.println("intercept b: " + Doubles.join(" ", gridData.b));
            System.out.println("number of connected components in true subgraph: " + numOfCC);
            for (int i = 0; i < gridData.p; i++) {
                ArrayRealVector xVector = new ArrayRealVector(gridData.x);
                ArrayRealVector wi = new ArrayRealVector(gridData.getColumn(i));
                if (gridData.Sy.contains(i)) {
                    System.out.println("xTransWi *: " + xVector.dotProduct(wi));
                } else {
                    System.out.println("xTransWi: " + xVector.dotProduct(wi));
                }
            }
            Utils.stop();
        }
        S2GraphMPDebug s2GraphMP = new S2GraphMPDebug(graph, k, s, 1, logisticRegression);
        PreRec preRec_nodes = new PreRec(s2GraphMP.psiX, gridData.Sx);
        PreRec preRec_features = new PreRec(s2GraphMP.psiY, gridData.Sy);
        if (verboseLevel > 0) {
            System.out.println("return nodes: " + s2GraphMP.psiX.toString());
            System.out.println(" true nodes: " + gridData.Sx.toString());
            System.out.println(preRec_nodes.toString());
            System.out.println(preRec_features.toString());
            Utils.stop();
        }
        List<PreRec> list = new ArrayList<>();
        list.add(preRec_nodes);
        list.add(preRec_features);
        return list;
    }

    public static void test_fix_y() {
        double[] muArr = new double[]{0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0};
        List<List<double[]>> results = new ArrayList<>();
        for (int i = 0; i < 20; i++) {
            results.add(new ArrayList<>());
            for (double mu : muArr) {
                List<PreRec> bestPreRec = null;
                int bestK = 0;
                for (int k = 5; k < 18; k++) {
                    System.out.println("--------------------------------------------------------------------------");
                    System.out.println("sparsity: " + k + " feature parameter: " + 5);
                    List<PreRec> preRec = testFix(k, 5, "./data/SerializedData/dataset" + i + "/Grid_np_100_10_z_20.0_r_5_mu_" + mu + "_.ser");
                    if (bestPreRec == null) {
                        bestPreRec = preRec;
                        bestK = k;
                    }
                    if (bestPreRec.get(0).fmeasure < preRec.get(0).fmeasure) {
                        bestPreRec = preRec;
                        bestK = k;
                    }
                }
                String result = "bestK: " + bestK + " " + bestPreRec.toString();
                results.get(i).add(new double[]{bestPreRec.get(0).pre, bestPreRec.get(0).rec, bestPreRec.get(0).fmeasure});
                System.out.println(result);
            }
        }

        for (int i = 0; i < results.size(); i++) {
            System.out.println("---------- result " + i + "-----------");
            for (double[] result : results.get(i)) {
                System.out.println(Doubles.join(" ", result));
            }
        }

        System.out.println("==============================================================");
        List<double[]> averageResult = new ArrayList<>();
        for (int i = 0; i < muArr.length; i++) {
            double pre = 0.0D;
            double rec = 0.0D;
            double fm = 0.0D;
            for (int j = 0; j < results.size(); j++) {
                pre += results.get(j).get(i)[0];
            }
            pre = pre / (results.size() * 1.0D);
            for (int j = 0; j < results.size(); j++) {
                rec += results.get(j).get(i)[1];
            }
            rec = rec / (results.size() * 1.0D);
            for (int j = 0; j < results.size(); j++) {
                fm += results.get(j).get(i)[2];
            }
            fm = fm / (results.size() * 1.0D);
            averageResult.add(new double[]{pre, rec, fm});
            System.out.println(pre + " " + rec + " " + fm + "\n");
        }
    }

    public static void test_fix_x() {
        double[] muArr = new double[]{10.0D};
        List<double[]> results = new ArrayList<>();
        for (double mu : muArr) {
            List<PreRec> bestPreRec = null;
            int bestS = 0;
            for (int k = 9; k < 15; k++) {
                for (int s = 5; s < 6; s++) {
                    System.out.println("number of features: " + s + " true feature parameter: " + 5);
                    List<PreRec> preRec = testFix(k, s, "./data/SerializedData/dataset0/Grid_np_100_10_z_20.0_r_5_mu_" + mu + "_.ser");
                    if (bestPreRec == null) {
                        bestPreRec = preRec;
                        bestS = s;
                    }
                    if (bestPreRec.get(1).fmeasure < preRec.get(1).fmeasure) {
                        bestPreRec = preRec;
                        bestS = s;
                    }
                }
            }
            String result = "bestK: " + bestS + " " + bestPreRec.toString();
            results.add(new double[]{bestPreRec.get(1).pre, bestPreRec.get(1).rec, bestPreRec.get(1).fmeasure});
            System.out.println(result);
        }
    }

    public static void test_random_x_random_y() {
        List<List<PreRec>> results = new ArrayList<>();
        for (int k = 9; k < 15; k++) {
            for (int s = 5; s < 6; s++) {
                System.out.println("===============================k: " + k + "==============================" + s + "=========================");
                List<PreRec> preRec = testFix(k, s, "./data/SerializedData/dataset0/Grid_np_100_10_z_20.0_r_5_mu_20.0_.ser");
                results.add(preRec);
            }
        }
        for (List<PreRec> result : results) {
            System.out.println(result);
        }
    }

    public static void main(String args[]) {
        int debugLevel = 2;
        if (debugLevel == 0) {
            test_fix_y();
        } else if (debugLevel == 1) {
            test_fix_x();
        } else {
            test_random_x_random_y();
        }
    }
}
