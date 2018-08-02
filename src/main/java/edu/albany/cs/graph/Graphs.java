package edu.albany.cs.graph;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Created by baojian bzhou6@albany.edu on 2/8/17.
 */
public final class Graphs {

    private Graphs() {
    }

    public static Graph generateGridGraph(int numOfNodes, int numOfFeatures) {
        int length = (int) Math.sqrt(numOfNodes);
        int width = length;
        return generateGridGraph(length, width, numOfFeatures);
    }

    public static Graph generateGridGraph(int length, int width, int p) {
        if (length < width) {
            System.out.println("length should be larger than width. ");
            System.exit(0);
        }
        int n = length * width;
        ArrayList<Integer[]> edges = new ArrayList<>();
        ArrayList<Double> edgeCosts = new ArrayList<>();
        int index = 0;
        for (int i = 0; i < length; i++) {
            for (int j = 0; j < width; j++) {
                if ((index % length) != (length - 1)) {
                    edges.add(new Integer[]{index, index + 1});
                    edgeCosts.add(1.0D);
                    if ((index + length) < n) {
                        edges.add(new Integer[]{index, index + length});
                        edgeCosts.add(1.0D);
                    }
                    //broader
                } else {
                    if ((index + length) < n) {
                        edges.add(new Integer[]{index, index + length});
                        edgeCosts.add(1.0D);
                    }
                }
                index++;
            }
        }
        Graph graph = new GridGraph(n, p, edges, edgeCosts);

        return graph;
    }

    public static Set<Integer> generateRandomSubset(int sizeOfSubset, int dimension) {
        Set<Integer> abnormalFeatures = new HashSet<>();
        Random random = new Random();
        while (true) {
            int feature = random.nextInt(dimension);
            if (!abnormalFeatures.contains(feature) && abnormalFeatures.size() < sizeOfSubset) {
                abnormalFeatures.add(feature);
            }
            if (abnormalFeatures.size() == sizeOfSubset) {
                break;
            }
        }
        return abnormalFeatures;
    }


}
