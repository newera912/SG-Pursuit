package edu.albany.cs.S2GraphMPDebug;


import com.google.common.primitives.Doubles;

import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.base.Utils;

import edu.albany.cs.scoreFuncs.*;

import java.util.*;

/**
 * Created by baojian on 1/17/17.
 * Email: bzhou6@albany.edu
 */
public class GradientDescentOpt {
    /**
     * By default, this function is two variables.
     */
    private final Function func;
    private final int n;
    private final int p;
    private final double stepSize;
    private final int maximumIter;

    private int verboseLevel = 0;

    GradientDescentOpt(Function func, int n, int p) {
        this(func, n, p, 0.01D);
    }

    public GradientDescentOpt(Function func, int n, int p, double stepSize) {
        this(func, n, p, stepSize, 1000);
    }

    public GradientDescentOpt(Function func, int n, int p, double stepSize, int maximumIter) {
        this.func = func;
        this.n = n;
        this.p = p;
        this.stepSize = stepSize;
        this.maximumIter = maximumIter;
    }

    public List<double[]> projectedGradientDescent(double[] x0, double[] y0, Set<Integer> OmegaX, Set<Integer> OmegaY) {
        if (func instanceof ElevatedMeanScan) {
            return multiGradientDescent4EMSScore(x0, y0, OmegaX, OmegaY);
        } else if (func instanceof LogLogisticRegression) {
            return multiGradientDescent4LogLogistic(OmegaX, OmegaY);
        } else {
            return null;
        }
    }


    private List<double[]> multiGradientDescent4EMSScore(double[] x0, double[] y0, Set<Integer> OmegaX, Set<Integer> OmegaY) {
        System.out.println("start argmin ....");
        double[] x = new double[n];
        double[] y = new double[p];
        double[] indicatorX = getIndicateVector(OmegaX, n);
        double[] indicatorY = getIndicateVector(OmegaY, p);
        double[] gradientX;
        double[] gradientY;
        System.arraycopy(x0, 0, x, 0, x.length);
        System.arraycopy(y0, 0, y, 0, y.length);
        double[] xOld = new double[x.length];
        double[] yOld = new double[y.length];
        for (int i = 0; i < maximumIter; i++) {
            gradientX = func.getGradientX(x, y);
            gradientY = func.getGradientY(x, y);
            System.arraycopy(x, 0, xOld, 0, x.length);
            System.arraycopy(y, 0, yOld, 0, y.length);
            x = updatedMinimizerX(gradientX, indicatorX, x, 5);
            y = updatedMinimizerY(gradientY, indicatorY, y, 5);
            double diffNormX = getL2Norm(x, xOld);
            double diffNormY = getL2Norm(y, yOld);
            if (diffNormX <= 1e-6 && diffNormY <= 1e-6) {
                break;
            }
            if (i % 100 == 0) {
                System.out.println("processed: " + i);
            }
        }
        List<double[]> minXY = new ArrayList<>();
        minXY.add(x);
        minXY.add(y);
        return minXY;
    }


    private List<double[]> multiGradientDescent4LogLogistic(Set<Integer> OmegaX, Set<Integer> OmegaY) {

        double[] indicatorX = getIndicateVector(OmegaX, n);
        double[] indicatorY = getIndicateVector(OmegaY, p);
        double[] gradientX;
        double[] gradientY;
        List<double[]> list = getInitialXY(indicatorX, indicatorY);
        double[] x = list.get(0);
        double[] y = list.get(1);
        if (verboseLevel > 0) {
            System.out.println("indicator x: " + Doubles.join(" ", indicatorX));
            System.out.println("indicator y: " + Doubles.join(" ", indicatorY));
            System.out.println("x0: " + Doubles.join(" ", x));
            System.out.println("y0:  " + Doubles.join(" ", y));
            System.out.println("supp(x0): " + Stat.supp(x).toString());
            System.out.println("supp(y0): " + Stat.supp(y).toString());
            Utils.stop();
        }
        /**Gradient descent fixed step size*/
        int numOfIter = 0;
        double[] xOld = new double[x.length];
        double[] yOld = new double[y.length];
        while (true) {
            // TODO : minimize simultaneous
            gradientX = func.getGradientX(x, y);
            gradientY = func.getGradientY(x, y);
            for (int j = 0; j < n; j++) {
                xOld[j] = x[j];
                x[j] = (x[j] - stepSize * gradientX[j]) * indicatorX[j];// update x
            }
            for (int j = 0; j < p; j++) {
                yOld[j] = y[j];
                y[j] = (y[j] - stepSize * gradientY[j]) * indicatorY[j];
                if (y[j] < 0.0) {
                    y[j] = 0.0D;
                } else if (y[j] > 1.0D) {
                    y[j] = 1.0D;
                }
            }
            if (verboseLevel > 0) {
                System.out.println("x: " + Doubles.join(" ", x));
                System.out.println("y: " + Doubles.join(" ", y));
                System.out.println("oldx: " + Doubles.join(" ", xOld));
                System.out.println("oldy: " + Doubles.join(" ", yOld));
                System.out.println("gradX: " + Doubles.join(" ", gradientX));
                System.out.println("gradY: " + Doubles.join(" ", gradientY));
                Utils.stop();
            }
            double gapX = getL2Norm(x, xOld);
            double gapY = getL2Norm(y, yOld);
            numOfIter++;
            if ((gapX <= 1e-3 && gapY <= 1e-3) || (numOfIter > maximumIter)) {
                if (verboseLevel == 0) {
                    System.out.println("Error bound satisfied!!!!");
                    System.out.println("number of steps: " + numOfIter);
                    System.out.println("gapX: " + gapX + " gapY: " + gapY);
                    System.out.println("minimizer x: " + Doubles.join(" ", x));
                    System.out.println("minimizer y: " + Doubles.join(" ", y));
                    //Utils.stop();
                }
                break;
            }
            if (numOfIter % 100 == 0) {
                System.out.println(numOfIter + " gapX: " + gapX + " gapY: " + gapY);
            }
        }
        if (numOfIter >= maximumIter) {
            return multiGradientDescent4LogLogistic(OmegaX, OmegaY);
        }
        System.out.println("size of supp(y): " + Stat.supp(y));

        List<double[]> minXY = new ArrayList<>();
        minXY.add(x);
        minXY.add(y);
        return minXY;
    }


    private double getL2Norm(double[] x1, double[] x2) {
        double l2norm = 0.0D;
        if (x1 == null || x1.length == 0.0D) {
            return l2norm;
        }
        for (int j = 0; j < x1.length; j++) {
            l2norm += (x1[j] - x2[j]) * (x1[j] - x2[j]);
        }
        return Math.sqrt(l2norm);
    }

    /**
     * generate a indicator vector.
     *
     * @param S    vector[i] = 1.0 if i \in S; 0 otherwise.
     * @param size the size of the vector.
     * @return the indicator vector.
     */
    private double[] getIndicateVector(Set<Integer> S, int size) {
        if (size <= 0) {
            return null;
        }
        double[] x = new double[size];
        Arrays.fill(x, 0.0D);
        for (int i : S) {
            x[i] = 1.0D;
        }
        return x;
    }

    /**
     * To generate a initial feasible solution.
     *
     * @return a initial <x0,y0> pair.
     */
    private List<double[]> getInitialXY(double[] indicatorX, double[] indicatorY) {
        double[] x = new double[n];
        double[] y = new double[p];
        Random random = new Random();
        for (int i = 0; i < x.length; i++) {
            x[i] = indicatorX[i] * random.nextDouble();
        }
        for (int i = 0; i < y.length; i++) {
            y[i] = indicatorY[i] * 1.0D;
        }
        List<double[]> list = new ArrayList<>();
        list.add(x);
        list.add(y);
        return list;
    }

    private double[] updatedMinimizerY(double[] gradientY, double[] indicatorY, double[] y, int bound) {
        double[] normalizedY = new double[y.length];
        for (int j = 0; j < y.length; j++) {
            normalizedY[j] = (y[j] - stepSize * gradientY[j]) * indicatorY[j];
        }
        /**Projecting y to [0,1]^p*/
        ArrayIndexSort arrayIndexSort = new ArrayIndexSort(normalizedY);
        Integer[] indexes = arrayIndexSort.getIndices();
        Arrays.sort(indexes, arrayIndexSort);
        int cnt = 0;
        for (int j = 0; j < y.length; j++) {
            if (normalizedY[j] <= 0.0) {
                cnt += 1;
                normalizedY[j] = 0.0D;
            } else if (normalizedY[j] > 1.0D) {
                normalizedY[j] = 1.0D;
            }
        }
        if (cnt == y.length) {
            System.out.println("!!!!!Warning: Sigmas1 is too small and all values in the gradient vector are nonpositive!!!");
            for (int i = 0; i < bound; i++) {
                normalizedY[indexes[i]] = 1;
            }
        }
        return normalizedY;
    }

    private double[] updatedMinimizerX(double[] gradientX, double[] indicatorX, double[] x, int bound) {
        double[] normalizedX = new double[x.length];
        for (int j = 0; j < x.length; j++) {
            normalizedX[j] = (x[j] - stepSize * gradientX[j]) * indicatorX[j];
        }
        ArrayIndexSort arrayIndexSort = new ArrayIndexSort(normalizedX);
        Integer[] indexes = arrayIndexSort.getIndices();
        Arrays.sort(indexes, arrayIndexSort);
        int cnt = 0;
        for (int j = 0; j < x.length; j++) {
            if (normalizedX[j] <= 0.0) {
                cnt += 1;
                normalizedX[j] = 0.0D;
            } else if (normalizedX[j] > 1.0D) {
                normalizedX[j] = 1.0D;
            }
        }
        if (cnt == x.length) {
            System.out.println("!!!!!Warning: Sigmas1 is too large and all values in the gradient vector are nonpositive!!!");
            for (int i = 0; i < bound; i++) {
                normalizedX[indexes[i]] = 1;
            }
        }
        return normalizedX;
    }

    public static void main(String args[]) {
        System.out.println("test: " + Doubles.join(" ", new double[]{1, 2, 3, 4, 5}));
        System.out.println((-0.0D <= 0.0D));
    }
}
