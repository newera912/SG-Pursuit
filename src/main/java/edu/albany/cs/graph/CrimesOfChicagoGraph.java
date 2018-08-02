package edu.albany.cs.graph;


import java.util.*;

/**
 * Created by baojian bzhou6@albany.edu on 2/8/17.
 */
public class CrimesOfChicagoGraph extends Graph {

    public Set<Integer> trueFeatures;
    public Set<Integer> trueNodes;
    public double[] trueY;
    public double[] trueX;
    public int p;
    public int n;

    public CrimesOfChicagoGraph(int n, int p, ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, Set<Integer> trueFeatures, Set<Integer> trueNodes) {
        super(n, p, edges, edgeCosts);
        this.trueFeatures = trueFeatures;
        this.trueNodes = trueNodes;
        this.p = p;
        this.n = n;
        trueX = new double[this.n];
        for (int i = 0; i < this.n; i++) {
            if (this.trueNodes.contains(i)) {
                trueX[i] = 1.0D;
            } else {
                trueX[i] = 0.0D;
            }
        }
        trueY = new double[this.p];
        for (int i = 0; i < this.p; i++) {
            if (this.trueFeatures.contains(i)) {
                trueY[i] = 1.0D;
            } else {
                trueY[i] = 0.0D;
            }
        }
    }
}
