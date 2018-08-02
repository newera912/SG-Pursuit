package edu.albany.cs.scoreFuncs;

//import edu.albany.cs.base.;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import com.google.common.primitives.Doubles;

import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;

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

    //debug use only
    private final int verboseLevel = 0;
    private final int debugLevel = 0;

    public GradientDescentOpt(Function func, int n, int p) {
        this(func, n, p, 0.001D);
    }

    public GradientDescentOpt(Function func, int n, int p, double stepSize) {
        this(func, n, p, stepSize, 1000000);
    }

    public GradientDescentOpt(Function func, int n, int p, double stepSize, int maximumIter) {
        this.func = func;
        this.n = n;
        this.p = p;
        this.stepSize = stepSize;
        this.maximumIter = maximumIter;
    }

    public List<double[]> simpleGradientDescent(Set<Integer> OmegaX, Set<Integer> OmegaY, double[] x0) {
        double[] x = new double[n];
        double[] y = new double[p];
        for (int i = 0; i < p; i++) {
            if (OmegaY.contains(i)) {
                y[i] = 1.0D;
            } else {
                y[i] = 0.0D;
            }
        }
        double[] indicatorX = getIndicateVector(OmegaX, x.length);
        double[] gradientX;
        Random random = new Random();
        for (int i = 0; i < x.length; i++) {
            x[i] = random.nextDouble() * 5.0D;
        }
        if (x0 != null) {
            x = x0;
        }
        if (verboseLevel > 0) {
            System.out.println("x0:       " + Doubles.join(" ", x));
            System.out.println("Indicator x: " + Doubles.join(" ", indicatorX));
            System.out.println("supp(x0): " + Stat.supp(x).toString());
        }
        /**Gradient descent fixed step size*/
        double[] xOld = new double[x.length];
        int numOfIter = 0;
        while (true) {
            gradientX = func.getGradientX(x, y);
            for (int j = 0; j < x.length; j++) {
                xOld[j] = x[j];
                x[j] = (x[j] - stepSize * gradientX[j]) * indicatorX[j];// update x
            }
            if (verboseLevel > 0) {
                System.out.println("x: " + ArrayUtils.toString(x));
                System.out.println("y: " + ArrayUtils.toString(y));
                System.out.println("oldx: " + ArrayUtils.toString(xOld));
                System.out.println("gradX: " + ArrayUtils.toString(gradientX));
            }
            double gap = getL2Norm(x, xOld);
            if (func.getFuncValue(x, y) > func.getFuncValue(xOld, y) && verboseLevel > 0) {
                System.out.println("current function value: " + func.getFuncValue(x, y));
                System.out.println("old function value: " + func.getFuncValue(xOld, y));
                System.exit(0);
            }
            numOfIter++;
            if (gap <= 1e-3 || numOfIter > maximumIter) {
                System.out.println("gap" + gap + " is satisfied. the number of iteration: " + numOfIter);
                break;
            }
            if (numOfIter % 100 == 0) {
                double normGradient = new ArrayRealVector(func.getGradientX(x, y)).getNorm();
                System.out.println("number of iteration udpated: " + numOfIter + " updated gap: " + gap + " ; current gradient norm: " + normGradient);
            }
        }
        if (verboseLevel > 0) {
            System.out.println("In argMin x: " + Doubles.join(" ", x));
            System.out.println("In argMin y: " + Doubles.join(" ", y));
        }
        List<double[]> minXY = new ArrayList<>();
        minXY.add(x);
        minXY.add(y);
        return minXY;
    }

    public List<double[]> alternativeGradientDescent4LogLogistic(Set<Integer> OmegaX, Set<Integer> OmegaY, Set<Integer> trueFeatures) {
        double[] indicatorX = getIndicateVector(OmegaX, n);
        double[] indicatorY = getIndicateVector(OmegaY, p);
        List<double[]> list = getInitialXY(indicatorX, indicatorY);
        double[] x = list.get(0);
        double[] y = list.get(1);
        int numOfIter = 0;
        //update y first
        while (true) {
            double[] xOld = Arrays.copyOf(x, x.length);
            double[] yOld = Arrays.copyOf(y, y.length);

            while (true) {
                double[] yPre = Arrays.copyOf(y, y.length);
                double[] gradientY = func.getGradientY(x, y);
                for (int j = 0; j < p; j++) {
                    y[j] = (y[j] - stepSize * gradientY[j]) * indicatorY[j];
                }
                for (int j = 0; j < p; j++) {
                    if (y[j] < 0.0) {
                        y[j] = 0.0D;
                    } else if (y[j] > 1.0D) {
                        y[j] = 1.0D;
                    }
                }
                double gapY = getL2Norm(y, yPre);
                if ((gapY <= 1e-3) || (numOfIter > maximumIter)) {
                    System.out.println("gap" + gapY + " is satisfied. the number of iteration: " + numOfIter);
                    break;
                }
                if (numOfIter % 100 == 0) {
                    double normGradient = new ArrayRealVector(func.getGradientY(x, y)).getNorm();
                    System.out.println(new PreRec(Stat.supp(y), trueFeatures));
                    System.out.println("number of iteration udpated: " + numOfIter + " updated gapY: " + gapY + " ; current gradient norm: " + normGradient);
                }
                numOfIter++;
            }
            numOfIter = 0;
            //then update x
            while (true) {
                double[] xPre = Arrays.copyOf(x, x.length);
                double[] gradientX = func.getGradientX(xPre, y);
                for (int j = 0; j < n; j++) {
                    x[j] = (xPre[j] - stepSize * gradientX[j]) * indicatorX[j];// update x
                }
                double gapX = getL2Norm(x, xPre);
                if ((gapX <= 1e-6D) || (numOfIter > maximumIter)) {
                    break;
                }
                numOfIter++;
            }

            if ((getL2Norm(x, xOld) <= 1e-6 && getL2Norm(y, yOld) <= 1e-6)) {
                System.out.println("norm of x: " + getL2Norm(x, xOld) + " norm of y: " + getL2Norm(y, yOld));
                break;
            } else {
                System.out.println("norm of x: " + getL2Norm(x, xOld) + " norm of y: " + getL2Norm(y, yOld));
            }
        }
        List<double[]> minXY = new ArrayList<>();
        minXY.add(x);
        minXY.add(y);
        return minXY;
    }

    public List<double[]> multiGradientDescent4LogLogistic_with_initial(Set<Integer> OmegaX, Set<Integer> OmegaY, List<double[]> initial) {
        System.out.println("------------ begin to run multi gradient descent with initial -------------");
        double[] indicatorX = getIndicateVector(OmegaX, n);
        double[] indicatorY = getIndicateVector(OmegaY, p);
        double[] gradientX;
        double[] gradientY;
        double[] x = initial.get(0);
        double[] y = initial.get(1);
        /**Gradient descent fixed step size*/
        int numOfIter = 0;
        double[] xOld = new double[x.length];
        double[] yOld = new double[y.length];
        while (true) {
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
            double gapX = getL2Norm(x, xOld);
            double gapY = getL2Norm(y, yOld);
            numOfIter++;
            if ((gapX <= 1e-3 && gapY <= 1e-3) || (numOfIter > maximumIter)) {
                System.out.println("number of iteration: " + numOfIter);
                break;
            }
            if (numOfIter % 100 == 0) {
                double normGradient = new ArrayRealVector(func.getGradientX(x, y)).getNorm();
                System.out.println("number of iteration udpated: " + numOfIter + " updated gapX: " + gapX + " ,gapY: " + gapY + " ; current gradient norm: " + normGradient);
            }

        }
        if (numOfIter >= maximumIter) {
            return multiGradientDescent4LogLogistic(OmegaX, OmegaY);
        }
        List<double[]> minXY = new ArrayList<>();
        minXY.add(x);
        minXY.add(y);
        return minXY;
    }


    public List<double[]> multiGradientDescent4EMSScore(Set<Integer> OmegaX, Set<Integer> OmegaY, int verboseLevel) {
        if (verboseLevel == 0) {
            System.out.println("start to do argmin for EMS score...");
        }
        double[] indicatorX = getIndicateVector(OmegaX, n);
        double[] indicatorY = getIndicateVector(OmegaY, p);
        double[] gradientX;
        double[] gradientY;
        List<double[]> list = getInitialXY(indicatorX, indicatorY);
        double[] x = list.get(0);
        double[] y = list.get(1);

        if (verboseLevel == 0) {
            System.out.println("indicator x: " + Doubles.join(" ", indicatorX));
            System.out.println("indicator y: " + Doubles.join(" ", indicatorY));
            System.out.println("x0: " + Doubles.join(" ", x));
            System.out.println("y0: " + Doubles.join(" ", y));
            System.out.println("supp(x0):  " + Stat.supp(x).toString());
            System.out.println("supp(y0):  " + Stat.supp(y).toString());
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
                if (x[j] < 0.0D) {
                    x[j] = 0.0D;
                } else if (x[j] > 1.0D) {
                    x[j] = 1.0D;
                }
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
            if (verboseLevel == 0) {
                System.out.println("x:  " + Doubles.join(" ", x));
                System.out.println("y:  " + Doubles.join(" ", y));
                System.out.println("oldx: " + Doubles.join(" ", xOld));
                System.out.println("oldy: " + Doubles.join(" ", yOld));
                System.out.println("gradX: " + Doubles.join(" ", gradientX));
                System.out.println("gradY: " + Doubles.join(" ", gradientY));
                Utils.stop();
            }
            double gapX = getL2Norm(x, xOld);
            double gapY = getL2Norm(y, yOld);
            numOfIter++;
            if ((gapX <= 1e-2 && gapY <= 1e-2) || (numOfIter > maximumIter) || (func.getFuncValue(xOld, yOld) - func.getFuncValue(x, y)) < 1e-8D) {
                if (verboseLevel == 0) {
                    System.out.println("Error bound satisfied!!!! ");
                    System.out.println("number of steps: " + numOfIter);
                    System.out.println("gapX: " + gapX + " gapY: " + gapY);
                    System.out.println("minimizer x: " + Doubles.join(" ", x));
                    System.out.println("minimizer y: " + Doubles.join(" ", y));
                    Utils.stop();
                }
                break;
            }
            if (verboseLevel < 0) {
                if (numOfIter % 100 == 0) {
                    System.out.println(numOfIter + " gapX: " + gapX + " gapY: " + gapY);
                }
            }
        }
        if (numOfIter >= maximumIter) {
            System.out.println("recursive...");
            return multiGradientDescent4LogLogistic(OmegaX, OmegaY);
        }
        List<double[]> minXY = new ArrayList<>();
        minXY.add(x);
        minXY.add(y);
        return minXY;
    }


    public List<double[]> multiGradientDescent4LogLogistic(Set<Integer> OmegaX, Set<Integer> OmegaY) {

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

//    public List<double[]> multiGradientDescent(Set<Integer> OmegaX, Set<Integer> OmegaY) {
//        double[] x = new double[n];
//        double[] y = new double[p];
//        double[] indicatorX = getIndicateVector(OmegaX, n);
//        double[] indicatorY = getIndicateVector(OmegaY, p);
//        double[] gradientX;
//        double[] gradientY;
//        Random random = new Random();
//        /** initialize x,y with a small number [0.0,1.0] */
//        for (int i = 0; i < x.length; i++) {
//            x[i] = indicatorX[i] * random.nextDouble();
//        }
//        for (int i = 0; i < y.length; i++) {
//            y[i] = indicatorY[i] * random.nextDouble();
//        }
//        if (verboseLevel > 0) {
//            System.out.println("x0:       " + Doubles.join(" ", x));
//            System.out.println("supp(x0): " + Stat.supp(x).toString());
//            System.out.println("supp(y0): " + Stat.supp(y).toString());
//            System.out.println("Ix: " + Doubles.join(" ", indicatorX));
//            System.out.println("Iy: " + Doubles.join(" ", indicatorY));
//            Utils.stop();
//        }
//        /**Gradient descent fixed step size*/
//        double[] xOld = new double[n];
//        double[] yOld = new double[p];
//        for (int i = 0; i < maximumIter; i++) {
//            gradientX = func.getGradientX(x, y);
//            gradientY = func.getGradientY(x, y);
//            System.arraycopy(x, 0, xOld, 0, x.length);
//            System.arraycopy(y, 0, yOld, 0, y.length);
//            for (int j = 0; j < x.length; j++) {
//                x[j] = (x[j] - stepSize * gradientX[j]) * indicatorX[j];// update x
//            }
//            y = updatedMinimizerY(gradientY, indicatorY, y);
//            if (verboseLevel > 0) {
//                System.out.println("x: " + ArrayUtils.toString(x));
//                System.out.println("y: " + ArrayUtils.toString(y));
//                System.out.println("oldx: " + ArrayUtils.toString(xOld));
//                System.out.println("oldy: " + ArrayUtils.toString(yOld));
//                System.out.println("gradX: " + ArrayUtils.toString(gradientX));
//                System.out.println("gradY: " + ArrayUtils.toString(gradientY));
//            }
//            double diffNormX = getL2Norm(x, xOld);
//            double diffNormY = getL2Norm(y, yOld);
//            if (diffNormX <= 1e-6 && diffNormY <= 1e-6) {
//                break;
//            }
//        }
//        if (verboseLevel > 0) {
//            System.out.println("In argMin x: " + Doubles.join(" ", x));
//            System.out.println("In argMin y: " + Doubles.join(" ", y));
//        }
//
//        List<double[]> minXY = new ArrayList<>();
//        minXY.add(x);
//        minXY.add(y);
//        return minXY;
//    }

    public void showstat(double[] x, String title, double lowerbound, double upperbound) {
        String line = "";
        int n = 0;
        int n1 = 0;
        for (int j = 0; j < x.length; j++) {
            if (x[j] != 0 && x[j] > lowerbound && x[j] < upperbound) {
//        		line += ", (" + j + ", " + x[j] + ")";
                line += "," + j;
                n += 1;
                if (j >= 300 && j < 400) {
                    n1 += 1;
                }
            }
        }
        System.out.println(title + line);
        System.out.println("In total: " + n + ", " + n1);
    }

    public void showstat(double[] x, String title) {
        String line = "";
        for (int j = 0; j < x.length; j++) {
            if (x[j] != 0) {
                line += ", (" + j + ", " + x[j] + ")";
            }
        }
        System.out.println(title + line);
    }


    public List<double[]> multiGradientDescent(double[] x0, double[] y0, Set<Integer> OmegaX, Set<Integer> OmegaY, double lambda) {
        double[] x = new double[n];
        double[] y = new double[p];
        double[] indicatorX = getIndicateVector(OmegaX, n);
        double[] indicatorY = getIndicateVector(OmegaY, p);
        double[] gradientX;
        double[] gradientY;
        System.arraycopy(x0, 0, x, 0, x.length);
        System.arraycopy(y0, 0, y, 0, y.length);
        double[] xOld = new double[n];
        double[] yOld = new double[p];
        for (int i = 0; i < maximumIter; i++) {
            int t = 0;
            do {
                gradientX = func.getGradientX(x, y, lambda, t * 0.01);
                t += 1;
            } while (StatUtils.min(gradientX) > 0);

            t = 0;
            do {
                gradientY = func.getGradientY(x, y, lambda, t * 0.01);
                t += 1;
            } while (StatUtils.min(gradientY) > 0);
            System.arraycopy(x, 0, xOld, 0, x.length);
            System.arraycopy(y, 0, yOld, 0, y.length);
            x = updatedMinimizerX(gradientX, indicatorX, x, 5);
            y = updatedMinimizerY(gradientY, indicatorY, y, 5);
            double diffNormX = getL2Norm(x, xOld);
            double diffNormY = getL2Norm(y, yOld);
            if (diffNormX <= 1e-6 && diffNormY <= 1e-6) {
                break;
            }
        }
        List<double[]> minXY = new ArrayList<>();
        minXY.add(x);
        minXY.add(y);
        return minXY;
    }

    public double[] l2normalization(double[] x) {
        double[] normalizedx = new double[x.length];
        double l2norm = 0;
        for (int j = 0; j < x.length; j++) {
            l2norm += x[j] * x[j];
        }
        l2norm = Math.sqrt(l2norm);
        for (int j = 0; j < x.length; j++) {
            normalizedx[j] = x[j] / l2norm;
        }
        return normalizedx;
    }

//    public double[] sigGradientDescent(Set<Integer> Sx, int[] Sy, Set<Integer> OmegaX, Set<Integer> OmegaY, double lambda) {
//    	
//    }

//    public List<double[]> multiGradientDescent(Set<Integer> OmegaX, Set<Integer> OmegaY, double lambda) {
//        double[] x = new double[n];
//        double[] y = new double[p];
//        double[] indicatorX = getIndicateVector(OmegaX, n);
//        double[] gradientX;
//        double[] indicatorY = getIndicateVector(OmegaY, p);
//        double[] gradientY;
////        Random random = new Random();
//        /** initialize x,y with a small number [0.0,1.0] */
////        for (int i = 0; i < x.length; i++) {
////            x[i] = indicatorX[i] * random.nextDouble();
////        }
////        Arrays.fill(x, 0.0D);
////        x[388] = 1;
////        x[Sx[1]] = 1;
////        for(int i:Sx){
////        	x[i] = 1;
////        	break; 
////        }
////        for(int i:Sy){
////        	y[i] = 1;
////        	break; 
////        }
//        /**Gradient descent fixed step size*/
//        double[] xOld = new double[n];
//        double[] yOld = new double[p];
//        for (int i = 0; i < maximumIter; i++) {
//            gradientX = func.getGradientX(x, y,lambda);
//            gradientY = func.getGradientY(x, y,lambda);
//            System.arraycopy(x, 0, xOld, 0, x.length); 
//            System.arraycopy(y, 0, yOld, 0, y.length); 
//            x = updatedMinimizerX(gradientX, indicatorX, x, OmegaX.size());
//            y = updatedMinimizerY(gradientY, indicatorY, y, OmegaY.size());
//            gradientX = func.getGradientX(x, y, lambda);
//            gradientY = func.getGradientY(x, y, lambda);
//            System.arraycopy(x, 0, xOld, 0, x.length);
//            System.arraycopy(y, 0, yOld, 0, y.length);
//            for (int j = 0; j < x.length; j++) {
//                x[j] = (x[j] - stepSize * gradientX[j]) * indicatorX[j];// update x
//            }
//            y = updatedMinimizerY(gradientY, indicatorY, y);
//            if (verboseLevel > 0) {
//                System.out.println("y: " + ArrayUtils.toString(x));
//                System.out.println("oldy: " + ArrayUtils.toString(xOld));
//                System.out.println("gradY: " + ArrayUtils.toString(gradientX));
//            }
//            
//            double diffNormX = getL2Norm(x, xOld);
//            if (diffNormX <= 1e-6) {
//                break;
//            }
//        }
//        if (verboseLevel > 0) {
//            System.out.println("In argMin x: " + Doubles.join(" ", x));
//            System.out.println("In argMin y: " + Doubles.join(" ", y));
//        }
//        return x;
//    }    

    public double[] sigGradientDescent_y(double[] x, Set<Integer> OmegaY, double lambda, int bound) {
        double[] y = new double[p];
        double[] indicatorY = getIndicateVector(OmegaY, p);
        double[] gradientY;
        Random random = new Random();
        /** initialize x,y with a small number [0.0,1.0] */
        for (int i = 0; i < y.length; i++) {
            y[i] = indicatorY[i] * random.nextDouble();
        }
        double[] yOld = new double[p];
        for (int i = 0; i < maximumIter; i++) {
            gradientY = func.getGradientY(x, y, lambda);
            System.arraycopy(y, 0, yOld, 0, y.length);
            y = updatedMinimizerY(gradientY, indicatorY, y, bound);
            double diffNormY = getL2Norm(y, yOld);
            if (diffNormY <= 1e-6) {
                break;
            }
        }
        return y;
    }


    public double[] sigGradientDescent_x(Set<Integer> Sx, int[] Sy, double[] y, Set<Integer> OmegaX, double lambda) {
        double[] x = new double[n];
        double[] indicatorX = getIndicateVector(OmegaX, n);
        double[] gradientX;
        Random random = new Random();
        /** initialize x,y with a small number [0.0,1.0] */
        for (int i = 0; i < x.length; i++) {
            x[i] = indicatorX[i] * random.nextDouble();
        }
        Arrays.fill(x, 0.0D);
//        x[388] = 1;
//        x[Sx[1]] = 1;
        for (int i : Sx) {
            x[i] = 1;
            break;
        }
//        y = l2normalization(y);
        if (verboseLevel > 0) {
            System.out.println("supp(x0): " + Stat.supp(x).toString());
            System.out.println("Ix: " + Doubles.join(" ", indicatorX));
            Utils.stop();
        }
        /**Gradient descent fixed step size*/
        double[] xOld = new double[n];
        double[] yOld = new double[p];
        for (int i = 0; i < maximumIter; i++) {
            gradientX = func.getGradientX(x, y, lambda);
            System.arraycopy(x, 0, xOld, 0, x.length);
            x = updatedMinimizerX(gradientX, indicatorX, x, OmegaX.size());
            double diffNormX = getL2Norm(x, xOld);
            if (diffNormX <= 1e-6) {
                break;
            }
        }
        if (verboseLevel > 0) {
            System.out.println("In argMin x: " + Doubles.join(" ", x));
            System.out.println("In argMin y: " + Doubles.join(" ", y));
        }
        return x;
    }


    private double[] getTrueY0(Set<Integer> OmegaY, int debugLevel, double[] y0) {
        double[] y = new double[p];
        if (debugLevel > 0) {
            for (int i = 0; i < p; i++) {
                if (OmegaY.contains(i)) {
                    y[i] = 1.0D;
                } else {
                    y[i] = 0.0D;
                }
            }
            return y;
        } else {
            return y0;
        }
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
            y[i] = indicatorY[i] * random.nextDouble();
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
			// System.out.println("!!!!!Warning: Sigmas1 is too small and all values in the gradient vector are nonpositive!!!");
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
			// System.out.println("!!!!!Warning: Sigmas1 is too large and all values in the gradient vector are nonpositive!!!");
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
