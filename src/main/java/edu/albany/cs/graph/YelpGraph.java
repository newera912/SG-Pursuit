package edu.albany.cs.graph;

import java.util.ArrayList;
import java.util.Set;

/**
 * Created by baojian on 2/17/17.
 */
public class YelpGraph extends Graph {


    public Set<Integer> trueFeatures;
    public Set<Integer> trueNodes;
    public double[] trueY;
    public double[] trueX;
    public int p;
    public int n;

    public YelpGraph(int n, int p, ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts, Set<Integer> trueFeatures, Set<Integer> trueNodes) {
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
