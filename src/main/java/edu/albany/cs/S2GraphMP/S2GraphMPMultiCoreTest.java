package edu.albany.cs.S2GraphMP;

import com.google.common.collect.Sets;
import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.graph.MultiDataGraph;
import edu.albany.cs.scoreFuncs.EMSXYScore;
import edu.albany.cs.scoreFuncs.Function;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class S2GraphMPMultiCoreTest {

    public int count = 0;

    public S2GraphMPMultiCoreTest(String fileName, int[] trueFea, String mu, String resultFileName, int k, int s, int g) {
        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(resultFileName, true);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        APDMInputFormat apdm = new APDMInputFormat(fileName);
        MultiDataGraph mgdGraph = new MultiDataGraph(apdm);
        double B = k - g + 0.0D;
        System.out.println("True: " + mgdGraph.trueSubGraphNodes);
        Function func = new EMSXYScore(mgdGraph.W);
        S2GraphMP s2GraphMP = new S2GraphMP(mgdGraph, k, s, g, B, func);
        Set<Integer> intersect = Sets.intersection(s2GraphMP.psiX, mgdGraph.trueNodesSet);
        double precision = (intersect.size() * 1.0D / s2GraphMP.psiX.size() * 1.0D);
        double recall = (intersect.size() * 1.0D / mgdGraph.trueNodesSet.size() * 1.0D);
        Set<Integer> yTrueArrayList= Sets.newHashSet(new Integer[]{0,1,2,3,4});
        Set<Integer> intersect2 = Sets.intersection(s2GraphMP.psiY, yTrueArrayList);
        double precision2 = (intersect2.size() * 1.0D / s2GraphMP.psiY.size() * 1.0D);
        double recall2 = (intersect2.size() * 1.0D / yTrueArrayList.size() * 1.0D);
        System.out.println(fileName.split("-")[6] + " k=" + k + " " + precision + " " + recall + " " + precision2 + " " + recall2);
        try {
            fileWriter.write(precision + " " + recall + " " + precision2 + " " + recall2 + "\n");
            fileWriter.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public S2GraphMPMultiCoreTest(String fileName, int[] trueFea, String mu, String resultFileName, int g) {
        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(resultFileName, true);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        ;
        APDMInputFormat apdm = new APDMInputFormat(fileName);
        MultiDataGraph mgdGraph = new MultiDataGraph(apdm);
        Function func = new EMSXYScore(mgdGraph.W);
        //System.out.println(mgdGraph.W.length+" "+mgdGraph.W[0].length);
        S2GraphMP bestS2GraphMP = null;
        double bestFunValue = Double.NEGATIVE_INFINITY;
        int bestS = 0;
        int bestK = 0;
        for (double ratio = 0.05; ratio < 0.4; ratio += 0.02) {
            int k = (int) Math.round(ratio * mgdGraph.numOfNodes);
            for (double ratioS = 0.1; ratioS < 0.6; ratioS += 0.1) {
                int s = (int) Math.round(ratioS * mgdGraph.numOfFeas);
                double B = k - g + 0.0D;
                S2GraphMP s2GraphMP = new S2GraphMP(mgdGraph, k, s, g, B, func);
                if (bestFunValue < s2GraphMP.funcValue) {
                    bestFunValue = s2GraphMP.funcValue;
                    bestS2GraphMP = s2GraphMP;
                    bestK = k;
                    bestS = s;
                }
            }

        }
        Set<Integer> intersect = Sets.intersection(bestS2GraphMP.psiX, mgdGraph.trueNodesSet);
        double precision = (intersect.size() * 1.0D / bestS2GraphMP.psiX.size() * 1.0D);
        double recall = (intersect.size() * 1.0D / mgdGraph.trueNodesSet.size() * 1.0D);
        //System.out.println();
        Set<Integer> yTrueArrayList= Sets.newHashSet(new Integer[]{0,1,2,3,4});
        Set<Integer> intersect2 = Sets.intersection(bestS2GraphMP.psiY, yTrueArrayList);
        double precision2 = (intersect2.size() * 1.0D / bestS2GraphMP.psiY.size() * 1.0D);
        double recall2 = (intersect2.size() * 1.0D / yTrueArrayList.size() * 1.0D);
        //System.out.println("File"+count+" X:pre=" + precision+"rec=" + recall+" Y: pre=" + precision2+"rec=" + recall2);
        System.out.println(fileName.split("-")[6] + " k=" + bestK + " s=" + bestS + " " + precision + " " + recall + " " + precision2 + " " + recall2 + " FunVal=" + bestFunValue);
        try {
            fileWriter.write(bestK + " " + bestS + " " + precision + " " + recall + " " + precision2 + " " + recall2 + " " + bestS2GraphMP.funcValue + "\n");
            fileWriter.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        //System.out.println(ArrayUtils.toString(mgdGraph.Y));
        //System.out.println(ArrayUtils.toString(W));


    }

    public static void TestSignleFile() {
        int[] trueFeatures = new int[]{0, 1, 2, 3, 4};
        String mu = "10";
        String resultPath = "outputs/S2GraphMPResult/debug_mu_" + mu + ".txt";
        int s = 5; //maximum feature number constrain
        //int k=10;   // sparsity constrain
        int g = 1;
        for (int k : new int[]{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}) {
            System.out.println("s=" + k);
            S2GraphMPMultiCoreTest s2GraphMPTest = new S2GraphMPMultiCoreTest("data/SimulationData/multiGridData_mu10/APDM-GridData-100-precen-0.1-noise_0.0-1.txt", trueFeatures, mu, resultPath, k, s, g);
        }
    }

    public static void TestMultiFile() {
        int[] trueFeatures = new int[]{0, 1, 2, 3, 4};
        long startTimeAll = System.nanoTime();
        ExecutorService pool = Executors.newFixedThreadPool(2);
        //double mu=10.0D;
        int s = 5; //maximum feature number constrain
        int k = 6;   // sparsity constrain
        int g = 1;

        for (int truSize : new int[]{10}) {//,15,20}){
            for (String mu : new String[]{"10.0"}) {//"0.5","1.0","2.0","3.0","5.0","10.0"}){
                System.out.println("--------------mu=" + mu + "-------------");
                //data/SimulationData/truesub_5/multiGridData_mu0.5
                String rootFolder = "data/SimulationData/multiCC/debug_truesub_" + truSize + "/multiGridData_mu" + mu + "/";
                for (File rawFile : new File(rootFolder).listFiles()) {
                    //genSignleAPDMFile(rawFile,rootFolder,outFolder,type);
                    String resultPath = "outputs/S2GraphMPResult/truesub_" + truSize + "_mu_" + mu + ".txt";
                    String fileName = rootFolder + rawFile.getName();
                    System.out.println(fileName);
                    pool.execute(new Thread() {
                        public void run() {

                            S2GraphMPMultiCoreTest s2GraphMPTest = new S2GraphMPMultiCoreTest(fileName, trueFeatures, mu, resultPath, g);

                        }
                    });

                }
            }
        }
        System.out.println("Total running time: " + (System.nanoTime() - startTimeAll) / 1e9);
    }

    public static void main(String[] args) {
        // TODO Auto-generated method stub
        TestMultiFile();
        //TestMultiFile();
    }

}
