package edu.albany.cs.graph;

import edu.albany.cs.base.APDMInputFormat;

import java.util.ArrayList;

/**
 * Created by baojian bzhou6@albany.edu on 2/4/17.
 */
public class ACEGeneGraph extends Graph {

    public ACEGeneGraph(int n, int p, ArrayList<Integer[]> edges, ArrayList<Double> edgeCosts) {
        super(n, p, edges, edgeCosts);
    }
}
