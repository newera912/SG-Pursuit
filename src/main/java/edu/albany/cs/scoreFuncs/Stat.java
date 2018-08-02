package edu.albany.cs.scoreFuncs;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;

import edu.albany.cs.base.DisjointSet;

public final class Stat {

    private final double[] data;
    private final int size;

    public Stat(double[] data) {
        if (data == null || data.length == 0) {
            this.data = null;
            this.size = 0;
        } else {
            this.data = new double[data.length];
            for (int i = 0; i < data.length; i++) {
                this.data[i] = data[i];
            }
            size = data.length;
        }
    }

    public double mean() {
        double sum = 0.0;
        if (data == null || data.length == 0) {
            return sum;
        }
        for (double a : data) {
            sum += a;
        }
        return sum / size;
    }

    public static double[] getColumn(int index, double[][] matrix) {
        double[] col = new double[matrix.length];
        for (int i = 0; i < matrix.length; i++) {
            col[i] = matrix[i][index];
        }
        return col;
    }

    private double mean(double[] x) {
        double sum = 0.0;
        for (double a : x) {
            sum += a;
        }

        return sum / (x.length + 0.0D);
    }

    double getVariance() {
        double mean = mean();
        double temp = 0;
        for (double a : data)
            temp += (mean - a) * (mean - a);
        return temp / size;
    }

    double getStdDev() {
        return Math.sqrt(getVariance());
    }

    public double median() {
        Arrays.sort(data);
        if (data.length % 2 == 0) {
            return (data[(data.length / 2) - 1] + data[data.length / 2]) / 2.0;
        } else {
            return data[data.length / 2];
        }
    }

    public double mad() {
        double[] abs = new double[data.length];
        for (int i = 0; i < data.length; i++) {
            abs[i] = Math.abs(data[i] - mean(data));
        }
        return mean(abs);
    }


    /**
     * To generate a random vector with size n.
     * Each of the entry in the return vector is 1.0 or 0.0.
     *
     * @param n the returned vector size.
     */
    public static double[] getRandomVector(int n) {
        double[] x0 = new double[n];
        Random rand = new Random();
        Arrays.fill(x0, 0.0D);
        /** to avoid zero vector.*/
        while (StatUtils.sum(x0) == 0) {
            for (int i = 0; i < x0.length; i++) {
                if (rand.nextDouble() < 0.5D) {
                    x0[i] = 1.0D;
                } else {
                    x0[i] = 0.0D;
                }
            }
        }

        return x0;
    }

    public static double[] getRandomVectorWithConstrain(int size, int nonzeroNum) {
        double[] x0 = new double[size];
        Random rand = new Random();
        int[] index = new Random().ints(0, size).distinct().limit(nonzeroNum).toArray();

        Arrays.fill(x0, 0.0D);
        for (int i : index) {
            x0[i] = 1.0D;
        }
        return x0;
    }

    public static double[] getRandomVector(int size, int s) {
        double[] x0 = new double[size];
        int[] index = new Random().ints(0, size).distinct().limit(s).toArray();
        Arrays.fill(x0, 0.0D);
        for (int i : index) {
            x0[i] = new Random().nextDouble();
        }
        //System.out.println("Sum Of Init " + StatUtils.sum(x0));
        return x0;
    }

    /**
     * Minimization Gradient Normalization
     * in order to fit the gradient, we need to normalize the gradient if the
     * domain of function is within [0,1]^n
     * <p>
     * grad_i { = 0       1) if grad_i<0 and x_i=1
     * 2)or if grad_i>0 and x_i =0
     * = grad_i     ,otherwise
     *
     * @param x         input vector x
     * @param gradient  gradient vector
     * @param graphSize the number of nodes
     * @return normGradient the normalized gradient
     */
    public static double[] normalizeGradient(double[] x, double[] gradient, int graphSize) {
        double[] normalizedGradient = new double[graphSize];
        for (int i = 0; i < graphSize; i++) {
            if ((gradient[i] < 0.0D) && (x[i] == 1.0D)) {
                normalizedGradient[i] = 0.0D;
            } else if ((gradient[i] > 0.0D) && (x[i] == 0.0D)) {
                normalizedGradient[i] = 0.0D;
            } else {
                normalizedGradient[i] = gradient[i];
            }
        }
        return normalizedGradient;
    }

    /**
     * get the nodes returned by algorithm
     *
     * @param x array x
     * @return the result nodes
     */
    public static int[] getSupp(double[] x) {
        int[] nodes = null;
        for (int i = 0; i < x.length; i++) {
            if (x[i] != 0.0D) {
                /** get nonzero nodes */
                nodes = ArrayUtils.add(nodes, i);
            }
        }
        if (nodes == null)
            nodes = new int[0];
        Arrays.sort(nodes);
        return nodes;
    }


    public static Set<Integer> getSuppSet(double[] x) {
        Set<Integer> nodes = new HashSet<Integer>();
        for (int i = 0; i < x.length; i++) {
            if (x[i] != 0.0D) {
                /** get nonzero nodes indecies */
                nodes.add(i);
            }
        }

        return nodes;
    }


    public static double[] getIndicateVec(Set<Integer> S, int n) {
        double[] Xs = new double[n];
        Arrays.fill(Xs, 0.0D);
        for (Integer i : S) {
            Xs[i] = 1.0D;
        }

        return Xs;
    }


    public static int[] getArray(Set<Integer> x) {
        int[] result = null;
        for (Integer i : x) {
            result = ArrayUtils.add(result, i);
        }
        if (result == null) {
            result = new int[0];
        }
        return result;
    }

    public static Set<Integer> Array2Set(int[] x) {


        Set<Integer> setX = new HashSet<>();
        if (x.length == 0) {
            return null;
        }
        for (int i = 0; i < x.length; i++) {
            setX.add(x[i]);
        }
        return setX;
    }

    public static String printWithIndex(double[] x) {
        String str = "{ ";
        for (int i = 0; i < x.length; i++) {
            str += i + ":" + Math.round(x[i] * 1000.0D) / 1000.0D + " ";
        }
        str += "}";
        return str;
    }

    public static String printWithIndex(int[] x) {
        String str = "{ ";
        for (int i = 0; i < x.length; i++) {
            str += i + ":" + x[i] + " ";
        }
        str += "}";
        return str;
    }

    public static Set<Integer> supp(double[] x) {

        Set<Integer> nodes = new HashSet<>();
        for (int i = 0; i < x.length; i++) {
            if (x[i] != 0.0D) {
                nodes.add(i);
            }
        }
        return nodes;
    }

    public static double[] sign(double[] x) {
        if (x == null || x.length == 0) {
            return null;
        }
        double[] signVector = new double[x.length];
        for (int i = 0; i < x.length; i++) {
            if (x[i] < 0.0D) {
                signVector[i] = -1.0D;
            } else {
                signVector[i] = 1.0D;
            }
        }
        return signVector;
    }

    public static double getL2Norm(double[] x1, double[] x2) {
        double l2norm = 0.0D;
        if (x1 == null || x1.length == 0.0D) {
            return l2norm;
        }
        for (int j = 0; j < x1.length; j++) {
            l2norm += (x1[j] - x2[j]) * (x1[j] - x2[j]);
        }
        return Math.sqrt(l2norm);
    }


    public static boolean isConnected(ArrayList<Integer[]> edges) {
        DisjointSet<Integer> dis = new DisjointSet<Integer>();
        Set<Integer> nodes = new HashSet<Integer>();

        for (Integer[] edge : edges) {
            nodes.add(edge[0]);
            nodes.add(edge[1]);
        }

        for (Integer node : nodes) {
            dis.makeSet(node);
        }

        for (Integer[] edge : edges) {
            dis.union(edge[0], edge[1]);
        }

        if (dis.numConnectedComponents == 1) {
            return true;
        } else {
            System.out.println("number of connected components is : " + dis.numConnectedComponents);
            System.out.println("number of nodes in graph: " + nodes.size());
            HashMap<Integer, Set<Integer>> cc = dis.getConnectedComponents();
            for (Integer index : cc.keySet()) {
                System.out.println("component[" + index + "]: " + cc.get(index).size());
            }
            return false;
        }
    }

    public static Set<Integer> getLargestCCNodes(ArrayList<Integer[]> edges) {
        Set<Integer> resultNodes = new HashSet<>();
        DisjointSet<Integer> dis = new DisjointSet<>();
        Set<Integer> nodes = new HashSet<>();
        for (Integer[] edge : edges) {
            nodes.add(edge[0]);
            nodes.add(edge[1]);
        }
        for (Integer node : nodes) {
            dis.makeSet(node);
        }
        for (Integer[] edge : edges) {
            dis.union(edge[0], edge[1]);
        }
        HashMap<Integer, Set<Integer>> cc = dis.getConnectedComponents();
        for (Integer index : cc.keySet()) {
            if (resultNodes.size() < cc.get(index).size()) {
                resultNodes = cc.get(index);
            }
        }
        return resultNodes;
    }

    public static int getMaximum(Set<Integer> set) {
        int maximum = Integer.MIN_VALUE;
        for (Integer node : set) {
            if (maximum < node) {
                maximum = node;
            }
        }
        return maximum;
    }

    public static int getMinimum(Set<Integer> set) {
        int minimum = Integer.MAX_VALUE;
        for (Integer node : set) {
            if (minimum > node) {
                minimum = node;
            }
        }
        return minimum;
    }

	public static double[] PRF(int inter, int res, int tru) {
		double[] prf = new double[3];
		prf[0] = 1.0D * inter / res;
		prf[1] = 1.0D * inter / tru;
		if (prf[0] == 0.0D || prf[1] == 0.0D) {
			prf[2] = 0.001D;
		} else {
			prf[2] = 2.0D * prf[0] * prf[1] / (prf[0] + prf[1]);
		}

		return prf;
	}
}