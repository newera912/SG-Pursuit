package edu.albany.cs.S2GraphMP;

import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.RandomWalk;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.Graphs;
import edu.albany.cs.graph.GridGraph;
import edu.albany.cs.scoreFuncs.EMSXYScore;
import edu.albany.cs.scoreFuncs.ElevatedMeanScan;
import edu.albany.cs.scoreFuncs.Function;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.util.Set;
import java.util.Random;

/**
 * Created by baojian bzhou6@albany.edu on 2/8/17.
 */
public class TestS2OnGridData {


    public TestS2OnGridData(int trueSubgraphSize, int numOfNodes, int numOfFeatures, int numOfAbnormalFeatures, double mu, int lowerK, int upperK) {
        Graph graph = Graphs.generateGridGraph(numOfNodes, numOfFeatures);
        RandomWalk randomWalk = new RandomWalk(graph.arrayListAdj, trueSubgraphSize);
        Set<Integer> subsetOfFeatures = Graphs.generateRandomSubset(numOfAbnormalFeatures, numOfFeatures);
        GridGraph grid = new GridGraph(numOfNodes, numOfFeatures, graph.edges, graph.edgeCosts, subsetOfFeatures, randomWalk.nodesInSubGraph);
        Function func = new EMSXYScore(getDataMatrix(numOfNodes, numOfFeatures, grid.trueFeatures, grid.trueNodes, mu));
        S2GraphMPDebug bestS2 = null;
        int bestK = 0;
        for (int k = lowerK; k < upperK; k++) {
            S2GraphMPDebug s2GraphMP = new S2GraphMPDebug(grid, k, 10, 1, func);
            if (bestS2 == null) {
                bestS2 = s2GraphMP;
                bestK = k;
            }
            if (bestS2.funcValue > s2GraphMP.funcValue) {
                bestS2 = s2GraphMP;
                bestK = k;
            }
        }
        PreRec preRec_nodes = new PreRec(bestS2.psiX, grid.trueNodes);
        PreRec preRec_features = new PreRec(bestS2.psiY, grid.trueFeatures);
        System.out.println("------------------------------------------------------------");
        System.out.println("mu: " + mu + " , bestK: " + bestK);
        System.out.println("preRec of nodes: " + preRec_nodes);
        System.out.println("preRec of features: " + preRec_features);
        System.out.println("true function value: " + func.getFuncValue(grid.trueX, grid.trueY));
        System.out.println("estimated function value: " + bestS2.funcValue);
        System.out.println("------------------------------------------------------------");
    }

    private double[][] getDataMatrix(int numOfNodes, int numOfFeatures, Set<Integer> trueFeatures, Set<Integer> trueNodes, double mu) {
        double[][] mat = new double[numOfNodes][numOfFeatures];
        NormalDistribution normal = new NormalDistribution(0.0D, 1.0D);
        double[] muArr = new double[]{2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 18.0, 20.0};
        int index = 0;
        for (int j = 0; j < numOfFeatures; j++) {
            if(trueFeatures.contains(j)){
                mu = muArr[index++];
            }
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
        for (double mu = 0.0D; mu < 10.0; mu += 1.0D) {
            int trueSubGraphSize = 15;
            int numOfNodes = 100;
            int numOfFeatures = 50;
            int numOfAbnormalFeatures = 10;
            int lowerK = 5;
            int upperK = 15;
            new TestS2OnGridData(trueSubGraphSize, numOfNodes, numOfFeatures, numOfAbnormalFeatures, mu, lowerK, upperK);
        }
    }
}
