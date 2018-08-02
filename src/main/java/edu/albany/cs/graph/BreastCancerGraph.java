package edu.albany.cs.graph;

import org.apache.commons.math3.stat.StatUtils;

import java.util.*;


/**
 * Created by baojian bzhou6@albany.edu on 2/4/17.
 */
public class BreastCancerGraph extends Graph {

    public Set<Integer> trueFeatures;
    public double[] trueY;
    public int partialK;
    private final int p;

    public BreastCancerGraph(int n, int p, ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, Set<Integer> trueFeatures, double[] trueY, int partialK) {
        super(n, p, edges, edgeCosts);
        this.trueFeatures = trueFeatures;
        this.trueY = trueY;
        this.p = trueY.length;
        this.partialK = partialK;
        if (StatUtils.sum(this.trueY) != 78.0D) {
            System.out.println("sume of true y is " + StatUtils.sum(this.trueY));
            System.out.println("sum of true y is not 78.0");
            System.exit(0);
        }
    }

    public double[] getPartialTrueY(int k) {
        double[] trueY = new double[p];
        Arrays.fill(trueY, 0.0D);
        Set<Integer> randNodes = new HashSet<>();
        List<Integer> listedTrueNodes = new ArrayList<>(trueFeatures);
        Random rand = new Random();
        while (randNodes.size() < k) {
            int randIndex = rand.nextInt(listedTrueNodes.size());
            randNodes.add(listedTrueNodes.get(randIndex));
        }
        for (int node : randNodes) {
            trueY[node] = 1.0D;
        }
        return trueY;
    }

}
