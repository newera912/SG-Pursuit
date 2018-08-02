package edu.albany.cs.LassoGraphMP;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.MultiDataGraph;
import edu.albany.cs.scoreFuncs.EMSXYScore;
import edu.albany.cs.scoreFuncs.Stat;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class TestingLassoGraphMP {


    public static int METHOD_VERSION = 2;    

    public static void TestLGraphMPMulti(String fileName, int[] trueFea, String mu, String resultFileName, int g) {
        FileWriter fileWriter = null;
        try {
            fileWriter = new FileWriter(resultFileName, true);
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
        
        APDMInputFormat apdm = new APDMInputFormat(fileName);
        MultiDataGraph mgdGraph = new MultiDataGraph(apdm);
        EMSXYScore func = new EMSXYScore(mgdGraph.W);
        LassoGraphMP bestLassoGraphMP = null;
        double bestFunValue = Double.POSITIVE_INFINITY;
        int bestS = 0;
        int bestK = 0;
        for (double ratio = 0.05; ratio < 0.6; ratio +=0.05) {
            int k = (int) Math.round(ratio * mgdGraph.numOfNodes);
            for (double ratioS = 0.2; ratioS < 0.6; ratioS += 0.1) {
                int s = (int) Math.round(ratioS * mgdGraph.numOfFeas);

                double B = k - g + 0.0D;
                
                LassoGraphMP lassoGraphMP = new LassoGraphMP(mgdGraph, func, k, s, g, B, METHOD_VERSION);
                if (bestFunValue > lassoGraphMP.getFunValue()) {
                    bestFunValue = lassoGraphMP.getFunValue();
                    bestLassoGraphMP = lassoGraphMP;
                    bestK = k;
                    bestS = s;
                }
            }
        }
        int[] intersect = Utils.intersect(Stat.getArray(bestLassoGraphMP.getResultNodes()), mgdGraph.trueSubGraphNodes);
        double precision = (intersect.length * 1.0D / bestLassoGraphMP.getResultNodes().size() * 1.0D);
        double recall = (intersect.length * 1.0D / mgdGraph.trueSubGraphNodes.length * 1.0D);
        
        int[] yTrueArrayList = new int[]{0, 1, 2, 3, 4};
        int[] intersect2 = Utils.intersect(Stat.getArray(bestLassoGraphMP.getResultFeatures()), yTrueArrayList);
        double precision2 = (intersect2.length * 1.0D / bestLassoGraphMP.getResultFeatures().size() * 1.0D);
        double recall2 = (intersect2.length * 1.0D / yTrueArrayList.length * 1.0D);
        
        System.out.println("mu=" + mu + " " + fileName.split("-")[fileName.split("-").length-1] + " k=" + bestK + " s=" + bestS + " " + precision + " " + recall + " " + precision2 + " " + recall2 + " FunVal=" + bestFunValue);
        try {
            fileWriter.write(bestK + " " + bestS+" "+precision + " " + recall + " " + precision2 + " " + recall2 + " " +bestLassoGraphMP.getFunValue() + "\n");
            fileWriter.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
           }

    public static void TestMultiFile() {
        int[] trueFeatures = new int[]{0, 1, 2, 3, 4};
        double[] pre=new double[3];
        double[] rec=new double[3];
        
        long startTimeAll = System.nanoTime();
        ExecutorService pool = Executors.newFixedThreadPool(30);
        //double mu=10.0D;

        int g = 1;

        int iter = 0;
        for (int truSize : new int[]{5,10,15,20}){//5,10,15,20}){//5,10,15,20}) {//5,10,15,20}){
            for (String mu : new String[]{"0.5","1.0","2.0","3.0","5.0","10.0"}){//"0.5","1.0","2.0","3.0","5.0","10.0"}){//,"2.0","3.0","5.0","10.0"}) {//"0.5","1.0","2.0","3.0","5.0","10.0"}){
                System.out.println("--------------mu=" + mu + "-------------");                
                String rootFolder = "data/SimulationData/truesub_" + truSize + "/multiGridData_mu" + mu + "/";
                for (File rawFile : new File(rootFolder).listFiles()) {                	
                	
                    //genSignleAPDMFile(rawFile,rootFolder,outFolder,type);
                    String resultPath = "outputs/lassoGraphMPResult/lasso_GraSP_truesub_" + truSize + "_mu_" + mu + ".txt";
                    String fileName = rootFolder + rawFile.getName();
                    System.out.println(fileName);                    
                    pool.execute(new Thread() {
                        public void run() {
                        	TestLGraphMPMulti(fileName, trueFeatures, mu, resultPath, g);
                        	
                        }
                        
                    });
                    iter++;
                }
            }
        }
        pool.shutdown();
        System.out.println("Total running time: " + (System.nanoTime() - startTimeAll) / 1e9);
       
		
    }

    public static void main(String[] args) {
        // TODO Auto-generated method stub
        TestMultiFile();
        //TestMultiFile();
    }

}
