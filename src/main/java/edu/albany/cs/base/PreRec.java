package edu.albany.cs.base;

import org.apache.commons.lang3.ArrayUtils;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * Given the result subgraph nodes and true subgraph nodes, we calculate the
 * precion, recall and fmeasure.
 * <p>
 * Please see more details of precision and recall:
 * {@link https://en.wikipedia.org/wiki/Precision_and_recall}
 * <p>
 * Please see more details of fmeasure:
 * {@link https://en.wikipedia.org/wiki/F1_score}
 *
 * @author baojian
 */
public class PreRec {

    public double pre;
    public double rec;
    public double fmeasure;
    public String preStr;
    public String recStr;
    public String fmeasureStr;

    public PreRec() {
        run(null, null);
    }

    public PreRec(int[] result, int[] groundTru) {
        run(result, groundTru);
    }

    public PreRec(ArrayList<Integer> nodes, int[] groundTru) {
        int[] result = null;
        for (int i : nodes) {
            result = ArrayUtils.add(result, i);
        }
        run(result, groundTru);
    }

    public PreRec(HashSet<Integer> nodes, int[] groundTru) {
        int[] result = null;
        for (int i : nodes) {
            result = ArrayUtils.add(result, i);
        }
        run(result, groundTru);
    }

    public PreRec(Set<Integer> nodes, Set<Integer> groundTruth) {
        int[] groundTruArr = new int[groundTruth.size()];
        int index = 0;
        for (int node : groundTruth) {
            groundTruArr[index++] = node;
        }
        int[] nodesArr = new int[nodes.size()];
        index = 0;
        for (int node : nodes) {
            nodesArr[index++] = node;
        }
        run(nodesArr, groundTruArr);
    }

    public PreRec(int[] nodes, Set<Integer> groundTruth) {
        int[] groundTruArr = new int[groundTruth.size()];
        int index = 0;
        for (int node : groundTruth) {
            groundTruArr[index++] = node;
        }
        int[] nodesArr = new int[nodes.length];
        index = 0;
        for (int node : nodes) {
            nodesArr[index++] = node;
        }
        run(nodesArr, groundTruArr);
    }

    public PreRec(List<Integer> nodes, Set<Integer> groundTruth) {
        int[] groundTruArr = new int[groundTruth.size()];
        int index = 0;
        for (int node : groundTruth) {
            groundTruArr[index++] = node;
        }
        int[] nodesArr = new int[nodes.size()];
        index = 0;
        for (int node : nodes) {
            nodesArr[index++] = node;
        }
        run(nodesArr, groundTruArr);
    }

    public void run(int[] result, int[] groundTru) {
        if (result == null || groundTru == null || result.length == 0 || groundTru.length == 0) {
            pre = 0.0D;
            rec = 0.0D;
            fmeasure = 0.0D;
        }
        int[] intersect = Utils.intersect(result, groundTru);
        if (intersect == null || intersect.length == 0) {
            pre = 0.0D;
            rec = 0.0D;
            fmeasure = 0.0D;
        } else {
            pre = (intersect.length * 1.0D / result.length * 1.0D);
            rec = (intersect.length * 1.0D / groundTru.length * 1.0D);
            fmeasure = (2.0D * pre * rec) / (pre + rec);
        }
        DecimalFormat df = new DecimalFormat("#.######");
        preStr = df.format(pre);
        recStr = df.format(rec);
        fmeasureStr = df.format(fmeasure);
    }

    @Override
    public String toString() {
        return "PreRec [pre=" + pre + ", rec=" + rec + ", fmeasure=" + fmeasure + "]";
    }

    public String toFileString() {
        return pre + " " + rec + " " + fmeasure;
    }

}
