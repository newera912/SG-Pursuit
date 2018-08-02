package edu.albany.cs.S2GraphMPDebug;

import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.base.Matrix;

import edu.albany.cs.scoreFuncs.FuncType;
import org.apache.commons.math3.stat.StatUtils;

import java.util.*;

/**
 * Created by baojian on 2/13/17.
 */
public class ElevatedMeanScan implements Function {

    private final double[][] W;
    private final double[][] WTrans;
    private final int n;
    private final int p;
    private final FuncType funcID;
    private final double lambda;

    private TestS2OnYelp.Data data;

    public ElevatedMeanScan(double[][] W, double lambda, TestS2OnYelp.Data data) {
        this(W, lambda);
        this.data = data;
    }

    public ElevatedMeanScan(double[][] W) {
        this(W, 10.0D);
    }

    public ElevatedMeanScan(double[][] W, double lambda) {
        funcID = FuncType.EMSXYScore;
        if (!checkInput(W)) {
            System.out.println(funcID + " input parameter is invalid.");
        }
        this.W = W;
        this.n = W.length;
        this.p = W[0].length;
        this.lambda = lambda;
        this.WTrans = new double[p][n];
        for (int i = 0; i < p; i++) {
            WTrans[i] = getColumn(i);
        }
    }

    private double getFun(int[] X, int[] Y) {
        double[] x = new double[n];
        double[] y = new double[p];
        Arrays.fill(x, 0.0D);
        for (int i : X) {
            x[i] = 1.0D;
        }
        Arrays.fill(y, 0.0D);
        for (int i : Y) {
            y[i] = 1.0D;
        }
        return getFuncValue(x, y);
    }

    private double getFun(int[] X, ArrayList<Integer> Y) {
        int[] YY = new int[Y.size()];
        int index = 0;
        for (int node : Y) {
            YY[index++] = node;
        }
        return getFun(X, YY);
    }


    List<double[]> calcInitilvalues(int k, int s) {
        double[] x0 = new double[n];
        double[] y0 = new double[p];
        ArrayList<ArrayList<int[]>> res = new ArrayList<>();
        ArrayList<double[]> scores = new ArrayList<>();
        int[] trials = new int[]{2};
        for (int i = 0; i < trials.length; i++) {
            res.add(new ArrayList<>());
            scores.add(new double[p]);
        }
        for (int j = 0; j < p; j++) {
            double[] vals = WTrans[j];
            ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(vals);
            Integer[] indexes = arrayIndexComparator.getIndices();
            Arrays.sort(indexes, arrayIndexComparator);
            int[] rank = new int[k];
            for (int i = 0; i < k; i++) {
                rank[i] = indexes[i];
            }
            int ii = 0;
            double fval = getFun(rank, new int[]{j});
            scores.get(ii)[j] = fval * -1;
            res.get(ii).add(rank);
        }
        double fValue = -1;
        ArrayList<Integer> Y = new ArrayList<>();
        ArrayList<Integer> X = new ArrayList<>();
        for (int ii = 0; ii < trials.length; ii++) {
            ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(scores.get(ii));
            Integer[] indexes = arrayIndexComparator.getIndices();
            Arrays.sort(indexes, arrayIndexComparator);
            for (int jj = 0; jj < s; jj++) {
                ArrayList<Integer> Y1 = new ArrayList<>();
                for (int i = 0; i < jj + 1; i++) {
                    Y1.add(indexes[i]);
                    //System.out.print(" " + data.reverseWordsDict.get(indexes[i].intValue()));
                }
                //System.out.println();
                int[] X1 = res.get(ii).get(indexes[0]);
                double fValue1 = getFun(X1, Y1);
                if (fValue == -1 || fValue > fValue1) {
                    Y.clear();
                    for (int i : Y1) Y.add(i);
                    X.clear();
                    for (int i : X1) X.add(i);
                    fValue = fValue1;
                }
            }
            for (int j = 0; j < s; j++) {
                int[] X1 = res.get(ii).get(indexes[j]);
                ArrayList<Integer> Y1 = new ArrayList<Integer>();
                Y1.add(indexes[j]);
                double fValue1 = getFun(X1, Y1);
                if (fValue == -1 || fValue > fValue1) {
                    Y.clear();
                    for (int i : Y1) Y.add(i);
                    X.clear();
                    for (int i : X1) X.add(i);
                    fValue = fValue1;
                }
            }
        }
        Arrays.fill(y0, 0);
        for (int i : Y) y0[i] = 1;
        Arrays.fill(x0, 0);
        for (int i : X) x0[i] = 1;
        List<double[]> list = new ArrayList<>();
        list.add(x0);
        list.add(y0);
        return list;
    }

    List<double[]> calcInitilvalues_2(int k, int s) {
        double[] x0 = new double[n];
        double[] y0 = new double[p];
        ArrayList<ArrayList<int[]>> res = new ArrayList<>();
        ArrayList<double[]> scores = new ArrayList<>();
        res.add(new ArrayList<>());
        scores.add(new double[p]);
        for (int j = 0; j < p; j++) {
            double[] vals = WTrans[j];
            ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(vals);
            Integer[] indexes = arrayIndexComparator.getIndices();
            Arrays.sort(indexes, arrayIndexComparator);
            int[] rank = new int[k];
            for (int i = 0; i < k; i++) {
                rank[i] = indexes[i];
            }
            int ii = 0;
            double sum = 0.0D;
            for (int i : rank) {
                sum += W[i][j];
            }
            scores.get(ii)[j] = -sum / (rank.length * 1.0D);
            res.get(ii).add(rank);
        }
        double fValue = -1;
        ArrayList<Integer> Y = new ArrayList<>();
        ArrayList<Integer> X = new ArrayList<>();

        ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(scores.get(0));
        Integer[] indexes = arrayIndexComparator.getIndices();
        Arrays.sort(indexes, arrayIndexComparator);
        for (int jj = 0; jj < s; jj++) {
            ArrayList<Integer> Y1 = new ArrayList<>();
            for (int i = 0; i < jj + 1; i++) {
                Y1.add(indexes[i]);
                System.out.print(" " + data.reverseWordsDict.get(indexes[i].intValue()));
            }
            System.out.println();
            int[] X1 = res.get(0).get(indexes[0]);
            double fValue1 = getFun(X1, Y1);
            if (fValue == -1 || fValue > fValue1) {
                Y.clear();
                for (int i : Y1) Y.add(i);
                X.clear();
                for (int i : X1) X.add(i);
                fValue = fValue1;
            }
        }
        for (int j = 0; j < s; j++) {
            int[] X1 = res.get(0).get(indexes[j]);
            ArrayList<Integer> Y1 = new ArrayList<>();
            Y1.add(indexes[j]);
            double fValue1 = getFun(X1, Y1);
            if (fValue == -1 || fValue > fValue1) {
                Y.clear();
                for (int i : Y1) Y.add(i);
                X.clear();
                for (int i : X1) X.add(i);
                fValue = fValue1;
            }
        }
        Arrays.fill(y0, 0);
        for (int i : Y) y0[i] = 1;
        Arrays.fill(x0, 0);
        for (int i : X) x0[i] = 1;
        List<double[]> list = new ArrayList<>();
        list.add(x0);
        list.add(y0);
        return list;
    }

    private double[] getColumn(int j) {
        double[] vals = new double[n];
        for (int i = 0; i < n; i++) {
            vals[i] = W[i][j];
        }
        return vals;
    }

    private boolean checkInput(double[][] W) {
        if (W == null) {
            return false;
        }
        return true;
    }

    @Override
    public FuncType getFuncID() {
        return null;
    }

    @Override
    public double getFuncValue(double[] x, double[] y) {
        double funcValue;
        double xT1 = StatUtils.sum(x);
        double yT1 = StatUtils.sum(y);
        double xTWy = Matrix.dot(Matrix.VecMultiplyMat(x, W), y);
        funcValue = -xTWy / Math.sqrt(xT1) + 0.5D * (lambda * yT1 * yT1);
        return funcValue;
    }

    @Override
    public double[] getGradientX(double[] x, double[] y) {
        if (y.length != p && x.length != n) {
            System.out.println(" GradientX: Error : Invalid parameters ...");
            System.exit(0);
        }
        if (StatUtils.sum(x) == 0) {
            System.out.println("GradientX:Input x vector values are all Zeros !!!");
            System.exit(0);
        }
        if (StatUtils.sum(y) == 0) {
            System.out.println("GradientX(: Input y vector values are all Zeros !!!");
            System.exit(0);
        }
        double[] gradient = new double[n];
        double xT1 = StatUtils.sum(x);
        double[] yW = Matrix.MatMultiplyVec(W, y);
        double[] term1 = Matrix.VarMultiplyVec(1.0D / Math.sqrt(xT1), yW);
        double xTWy = Matrix.dot(x, yW);
        double term2 = 0.5 * xTWy / (Math.sqrt(xT1) * xT1);
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = -(term1[i] - term2);
            if (!Double.isFinite(gradient[i])) {
                System.out.println("gradientX error. ");
                System.exit(0);
            }
        }
        return gradient;
    }

    @Override
    public double[] getGradientY(double[] x, double[] y) {
        if (y.length != p && x.length != n) {
            System.out.println(" GradientY: Error : Invalid parameters ...");
            System.exit(0);
        }
        if (StatUtils.sum(x) == 0) {
            System.out.println("GradientY: Input x vector values are all Zeros !!!");
            System.exit(0);
        }
        if (StatUtils.sum(y) == 0) {
            System.out.println("GradientY: Input y vector values are all Zeros !!!");
            System.exit(0);
        }
        double[] gradient = new double[p];
        double xT1 = StatUtils.sum(x);
        double sqrtXT1 = Math.sqrt(xT1);
        double[] xW = Matrix.MatMultiplyVec(WTrans, x);
        double[] term2 = Matrix.VarMultiplyVec(1.0D / sqrtXT1, xW);
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = -term2[i] + lambda * y[i];
            if (!Double.isFinite(gradient[i])) {
                System.out.println("gradientY error. ");
                System.exit(0);
            }
        }
        return gradient;
    }

    @Override
    public List<double[]> getArgMinFxy(double[] x0, double[] y0, Set<Integer> OmegaX, Set<Integer> OmegaY) {
        return new GradientDescentOpt(this, n, p).projectedGradientDescent(x0, y0, OmegaX, OmegaY);
    }

}
