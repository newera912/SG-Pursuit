package edu.albany.cs.S2GraphMPDebug;

import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.RandomWalk;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.Graphs;
import edu.albany.cs.graph.GridGraph;

import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;
import java.util.Random;

/**
 * Created by baojian bzhou6@albany.edu on 2/8/17.
 */
public class TestS2OnGridData {
    public static void log(String line) {
        FileWriter fstream;
        System.out.println(line);
        try {
            fstream = new FileWriter("data/DenseGraph/Log/log.txt", true);
            fstream.write(line + "\r\n");
            fstream.close();
        } catch (IOException e) {
            // TODO Auto-generated catch block
            e.printStackTrace();
        }
    }

    public TestS2OnGridData(int trueSubgraphSize, int numOfNodes, int numOfFeatures, int numOfAbnormalFeatures) {
        int n = numOfNodes;
        int p = numOfFeatures;
        int k = trueSubgraphSize / 2;
        int s = numOfAbnormalFeatures;
        int g = 1;
        double lambda = 10.0D;
        Graph graph = Graphs.generateGridGraph(numOfNodes, numOfFeatures);
        RandomWalk randomWalk = new RandomWalk(graph.arrayListAdj, trueSubgraphSize);
        Set<Integer> subsetOfFeatures = Graphs.generateRandomSubset(numOfAbnormalFeatures, numOfFeatures);
        GridGraph grid = new GridGraph(numOfNodes, numOfFeatures, graph.edges, graph.edgeCosts, subsetOfFeatures, randomWalk.nodesInSubGraph);
        Function func = new ElevatedMeanScan(getDataMatrix(numOfNodes, numOfFeatures, grid.trueFeatures, grid.trueNodes), lambda);
        S2GraphMPDebug s2GraphMP = new S2GraphMPDebug(grid, k, s, g, func,null);
        PreRec preRec_nodes = new PreRec(s2GraphMP.psiX, grid.trueNodes);
        PreRec preRec_features = new PreRec(s2GraphMP.psiY, grid.trueFeatures);
        log("------------------------------------------------------------");
        log(String.format("(n, p, k, s): (%d, %d, %d, %d); (%.2f, %.2f, %.2f, %.2f). %.2f", n, p, k, s, preRec_nodes.pre, preRec_nodes.rec, preRec_features.pre, preRec_features.rec, s2GraphMP.funcValue));
        log("------------------------------------------------------------");
    }

//	if(true){
//		String line = "****************************\nfeatures: ";
//		for(int i:grid.trueFeatures){
//			line += i + ",";
//		}
//		System.out.println(line);
//		line = "nodes: ";
//		for(int i:grid.trueNodes){
//			line += i + ",";
//		}
//		System.out.println(line);
//	}


//    public TestS2OnGridData(int trueSubgraphSize, int numOfNodes, int numOfFeatures, int numOfAbnormalFeatures, int lowerK, int upperK) {
//    	int n = numOfNodes;
//    	int p = numOfFeatures;
//    	int k = trueSubgraphSize;
//    	int s = numOfAbnormalFeatures;
//        Graph graph = Graphs.generateGridGraph(numOfNodes, numOfFeatures);
//        RandomWalk randomWalk = new RandomWalk(graph.arrayListAdj, trueSubgraphSize);
//        Set<Integer> subsetOfFeatures = Graphs.generateRandomSubset(numOfAbnormalFeatures, numOfFeatures);
//        GridGraph grid = new GridGraph(numOfNodes, numOfFeatures, graph.edges, graph.edgeCosts, subsetOfFeatures, randomWalk.nodesInSubGraph);
//        Function func = new ElevatedMeanScan(getDataMatrix(numOfNodes, numOfFeatures, grid.trueFeatures, grid.trueNodes));
//        S2GraphMPDebug bestS2 = null;
//        int bestK = 0;
//    	int g = 1;
//        for (int k = lowerK; k < upperK; k++) {
//    		if(true){
//    			String line = "****************************\nfeatures: ";
//    			for(int i:grid.trueFeatures){
//    				line += i + ",";
//    			}
//    			System.out.println(line);
//    			line = "nodes: ";
//    			for(int i:grid.trueNodes){
//    				line += i + ",";
//    			}
//    			System.out.println(line);
//    		}
//            S2GraphMPDebug s2GraphMP = new S2GraphMPDebug(grid, k, s, g, func);
//            if (bestS2 == null) {
//                bestS2 = s2GraphMP;
//                bestK = k;
//            }
//            if (bestS2.funcValue > s2GraphMP.funcValue) {
//                bestS2 = s2GraphMP;
//                bestK = k;
//            }
//        }
//        PreRec preRec_nodes = new PreRec(bestS2.psiX, grid.trueNodes);
//        PreRec preRec_features = new PreRec(bestS2.psiY, grid.trueFeatures);
//        log("------------------------------------------------------------");
//        log(String.format("(n, p, k, s): (%d, %d, %d, %d); (%.2f, %.2f, %.2f, %.2f). %.2f", n, p, k, s, preRec_nodes.pre, preRec_nodes.rec, preRec_features.pre, preRec_features.rec);
////        log("mu: " + mu + " , bestK: " + bestK);
////        log("preRec of nodes: " + preRec_nodes);
////        log("preRec of features: " + preRec_features);
////        System.out.println("true function value: " + func.getFuncValue(grid.trueX, grid.trueY));
////        System.out.println("estimated function value: " + bestS2.funcValue);
//        System.out.println("------------------------------------------------------------");
//    }


    private double[][] getDataMatrix(int numOfNodes, int numOfFeatures, Set<Integer> trueFeatures, Set<Integer> trueNodes) {
        double[][] mat = new double[numOfNodes][numOfFeatures];
        NormalDistribution normal = new NormalDistribution(0.0D, 1.0D);
        for (int j = 0; j < numOfFeatures; j++) {
//            double mu = new Random().nextDouble() * 2.0D + 0.6D; 
            double mu = 1.2;
            NormalDistribution abnormal = new NormalDistribution(mu, 1.0D);
            for (int i = 0; i < numOfNodes; i++) {
                if (trueFeatures.contains(j) && trueNodes.contains(i)) {
                    mat[i][j] = abnormal.sample();
                } else {
                    mat[i][j] = normal.sample();
                }
            }
        }
        return mat;
    }

    public static void main(String args[]) {
        int trueSubGraphSize = 15;
        int numOfNodes = 100;
        int numOfFeatures = 50;
        int numOfAbnormalFeatures = 10;
        new TestS2OnGridData(trueSubGraphSize, numOfNodes, numOfFeatures, numOfAbnormalFeatures);
    }
}
