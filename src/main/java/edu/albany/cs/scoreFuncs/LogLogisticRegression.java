package edu.albany.cs.scoreFuncs;

import com.google.common.primitives.Doubles;

import edu.albany.cs.base.ArrayIndexComparator;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.base.Utils;

import org.apache.commons.math3.linear.ArrayRealVector;
import org.apache.commons.math3.stat.StatUtils;

import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.*;

/**
 * Created by baojian on 1/21/17.
 * Email: bzhou6@albany.edu
 */
public class LogLogisticRegression implements Function {

    private final double[][] W;
    private final double[] b;
    private final int n;
    private final int p;
    private int sparsity_p;
    private final double lambdaX = 0.0D;
    private final double lambdaY = 0.0D;
    private final int verboseLevel = 0;

    private int[] trueFea = new int[]{
            256, 129, 257, 258, 259, 5, 136, 14, 143, 16, 19, 147, 23, 24, 280, 281, 282, 27, 283,
            284, 29, 157, 285, 286, 31, 287, 288, 289, 34, 290, 35, 291, 164, 292, 293, 294, 168,
            42, 171, 173, 46, 174, 177, 53, 183, 57, 187, 65, 196, 197, 75, 83, 86, 87, 91, 219,
            94, 95, 223, 96, 97, 225, 98, 227, 229, 102, 103, 105, 106, 237, 240, 116, 118, 248,
            250, 123, 254, 255};
    private Set<Integer> trueFeatures = new HashSet<>();

    public LogLogisticRegression(double[][] W) {
        this(W, null);
    } 

    public LogLogisticRegression(double[][] W, double[] b) {
        this(W, b, 0);
    }

    public LogLogisticRegression(double[][] W, double[] b, int sparsity_p) {
        if (b == null || b.length == 0) {
            this.b = new double[W[0].length];
            Arrays.fill(this.b, 0.0D);
        } else {
            this.b = b;
        }
        for (int f : trueFea) {
            trueFeatures.add(f);
        }
        this.W = W;
        this.sparsity_p = sparsity_p;
        this.n = this.W.length;
        this.p = this.W[0].length;
    }

    @Override
    public FuncType getFuncID() {
        return FuncType.LogLogisticRegression;
    }

    @Override
    public double getFuncValue(double[] x, double[] y) {
        double funcVal = 0.0D;
        if (x == null || y == null || x.length == 0 || y.length == 0) {
            return funcVal;
        }
        double sumObjective = 0.0D;
        ArrayRealVector xVector = new ArrayRealVector(x);
        ArrayRealVector yVector = new ArrayRealVector(y);
        for (int i = 0; i < p; i++) {
            ArrayRealVector wiVector = new ArrayRealVector(getColumn(i));
            double xTransWi = xVector.dotProduct(wiVector);
            double logSigmoidFuncVal = logSigmoidFuncVal(xTransWi);
            double term1 = y[i] * logSigmoidFuncVal;
            double term2 = (1.0D - y[i]) * (-xTransWi + logSigmoidFuncVal);
            sumObjective += (term1 + term2);
        }
        double regularizer = lambdaX * xVector.dotProduct(xVector) + lambdaY * yVector.dotProduct(yVector);
        funcVal = -sumObjective + regularizer;
        if (!Double.isFinite(funcVal)) {
            System.out.println("It is " + funcVal);
            System.exit(0);
        }
        return funcVal;
    }


    @Override
    public double[] getGradientX(double[] x, double[] y) {
        double[] gradient = new double[n];
        ArrayRealVector xVector = new ArrayRealVector(x);
        ArrayRealVector funcValVector = new ArrayRealVector(y);
        double[] inputValues = new double[p];
        for (int i = 0; i < p; i++) {
            ArrayRealVector wiVector = new ArrayRealVector(getColumn(i));
            double xTransWi = xVector.dotProduct(wiVector);
            funcValVector.setEntry(i, sigmoidFuncVal(xTransWi + b[i]));
            inputValues[i] = wiVector.dotProduct(xVector) + b[i];
        }
        if (verboseLevel > 0) {
            System.out.println("funcValVector: " + funcValVector.toString());
            System.out.println("input values: " + Doubles.join(" ", inputValues));
            Utils.stop();
        }
        for (int i = 0; i < n; i++) {
            double term1 = 0.0D;
            double term2 = 0.0D;
            double term3 = 2.0D * lambdaX * x[i];
            for (int j = 0; j < p; j++) {
                term1 += (1.0D - funcValVector.getEntry(j)) * W[i][j] * y[j];
                term2 += (funcValVector.getEntry(j)) * W[i][j] * (1.0D - y[j]);
            }
            gradient[i] = -term1 + term2 + term3;
            if (!Double.isFinite(gradient[i])) {
                System.out.println("it is " + gradient[i]);
                System.exit(0);
            }
        }
        return gradient;
    }

    @Override
    public double[] getGradientY(double[] x, double[] y) {
        double[] gradient = new double[p];
        ArrayRealVector xVector = new ArrayRealVector(x);
        for (int i = 0; i < p; i++) {
            ArrayRealVector wiVector = new ArrayRealVector(getColumn(i));
            double xTransWi = xVector.dotProduct(wiVector);
            double logSigmoidFuncVal = logSigmoidFuncVal(xTransWi);
            double term1 = -logSigmoidFuncVal;
            double term2 = +(-xTransWi + logSigmoidFuncVal);
            double term3 = 2.0D * lambdaY * y[i];
            gradient[i] = term1 + term2 + term3;
            if (!Doubles.isFinite(gradient[i])) {
                System.out.println("error: gradient is infinite...");
                System.exit(0);
            }
            if (!Double.isFinite(gradient[i])) {
                System.out.println("it is " + gradient[i]);
                System.exit(0);
            }
        }
        if (verboseLevel > 0) {
            System.out.println("x: " + Doubles.join(" ", x));
            System.out.println("y: " + Doubles.join(" ", y));
            Utils.stop();
        }
        return gradient;
    }

    @Override
    public List<double[]> getArgMinFxy(Set<Integer> OmegaX, Set<Integer> OmegaY) {
        System.out.println("pre rec fmeasure[OmegaY]: " + new PreRec(OmegaY, trueFeatures));
        System.out.println("size of OmegaX: " + OmegaX.size());
        System.out.println("size of OmegaY: " + OmegaY.size());
        //Utils.stop();
        //List<double[]> initial = getInitial(OmegaX, OmegaY, sparsity_p);
        //return new GradientDescentOpt(this, n, p, 0.001D, 10000).multiGradientDescent4LogLogistic_with_initial(OmegaX, OmegaY, initial);
        //return new GradientDescentOpt(this, n, p, 0.001D, 10000).multiGradientDescent4LogLogistic(OmegaX, OmegaY);
        return new GradientDescentOpt(this, n, p, 0.001D, 10000).alternativeGradientDescent4LogLogistic(OmegaX, OmegaY, trueFeatures);
    }

    public List<double[]> getArgMinFxy(Set<Integer> OmegaX, Set<Integer> OmegaY, List<double[]> initial) {
        System.out.println("sum of x: " + StatUtils.sum(initial.get(0)));
        System.out.println("sum of y: " + StatUtils.sum(initial.get(1)));
        return new GradientDescentOpt(this, n, p, 0.001D, 10000).multiGradientDescent4LogLogistic_with_initial(OmegaX, OmegaY, initial);
    }

    private List<double[]> getInitial(Set<Integer> OmegaX, Set<Integer> OmegaY, int topK) {
        double[] y = getInitial_Y0(OmegaY, topK);
        System.out.println("supp(y): " + Stat.supp(y));
        System.out.println("pre rec fmeasure: " + new PreRec(Stat.supp(y), trueFeatures));
        System.out.println("size of OmegaY: " + OmegaY.size());
        Utils.stop();

        GradientDescentOpt gradientDescentOpt = new GradientDescentOpt(this, n, p, 0.001D, 1000);
        double[] x = gradientDescentOpt.simpleGradientDescent(OmegaX, Stat.supp(y), null).get(0);
        List<double[]> results = new ArrayList<>();
        results.add(x);
        results.add(y);
        return results;
    }

    private double[] getInitial_Y0(Set<Integer> OmegaY, int topK) {
        double[] y0 = new double[p];
        double[] center = new double[n];
        for (int i = 0; i < n; i++) {
            center[i] = new Stat(W[i]).median();
        }
        ArrayRealVector centerVector = new ArrayRealVector(center);
        double[] distances = new double[p];
        for (int i = 0; i < p; i++) {
            distances[i] = centerVector.subtract(new ArrayRealVector(getColumn(i))).getNorm();
        }
        ArrayIndexComparator arrayIndexComparator = new ArrayIndexComparator(distances);
        Integer[] indexes = arrayIndexComparator.indexes;
        Arrays.sort(indexes, arrayIndexComparator);
        Arrays.fill(y0, 0.0D);
        for (int i = 0; i < topK; i++) {
            y0[indexes[i]] = 1.0D;
        }
        return y0;
    }

    @Override
    public List<double[]> getArgMinFx(double[] yi) {
        return null;
    }

    @Override
    public List<double[]> getArgMinFy(double[] xi) {
        return null;
    }

    public double[] getColumn(int index) {
        double[] col = new double[n];
        for (int i = 0; i < n; i++) {
            col[i] = W[i][index];
        }
        return col;
    }

    public static double sigmoidFuncVal(double x) {
        return 1.0D / (1.0D + Math.exp(-x));
    }

    public static double logSigmoidFuncVal(double x) {
        if (-x > 100.0D) {
            return x;
        }
        if (-x < -100.0D) {
            return 0.0D;
        }
        return -Math.log(1.0D + Math.exp(-x));
    }

    public static void main(String args[]) {

        for (double x = -10000000.0D; x < 10000000; x += 100000) {
            System.out.println(x + " " + logSigmoidFuncVal(x));
        }
        System.out.println("--------------------------");
        for (double x = -10000000.0D; x < 10000000; x += 100000) {
            System.out.println(x + " " + sigmoidFuncVal(x));
        }


        BigDecimal bigDecimal = new BigDecimal("-122222222222222222222222222222222222222222222222222222222222222222.1111");
        BigDecimal bb = new BigDecimal("0.123145235");
        System.out.println(Math.exp(-100000));
        System.out.println(Math.log(0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001));
        System.out.println(bigDecimal.abs().toPlainString());
        System.out.println(bigDecimal.abs().toEngineeringString());
        System.out.println(bigDecimal.abs().toString());
        System.out.println(bigDecimal.divide(bb, 100, RoundingMode.CEILING).toString());
        System.out.println(bigDecimal.divide(bb, 100, RoundingMode.CEILING).toPlainString());
        BigDecimal one = new BigDecimal(1.0D);
        BigDecimal bigX = new BigDecimal(10.0D);
        BigDecimal sum = one.add(bigX);
        BigDecimal result = one.divide(sum, 1000, RoundingMode.CEILING);
        System.out.println(result.toString());

    }

    @Override
    public double getFuncValue(double[] x, double[] y, double lambda) {
        // TODO Auto-generated method stub
        return 0;
    }

    @Override
    public double[] getGradientX(double[] x, double[] y, double lambda) {
        // TODO Auto-generated method stub
        return null;
    }

    @Override
    public double[] getGradientY(double[] x, double[] y, double lambda) {
        // TODO Auto-generated method stub
        return null;
    }

	@Override
	public List<double[]> getArgMinFxy(double[] x0, double[] y0,
			Set<Integer> OmegaX, Set<Integer> OmegaY, double lambda) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] getArgMinFy(double[] xi, Set<Integer> OmegaY, double lambda) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] getArgMinFx(Set<Integer> trueNodesSet, int[] trueFeas,
			double[] yi, Set<Integer> omegaX, double lambda) {
		// TODO Auto-generated method stub
		return null;
	}
//    @Override
//    public List<double[]> getArgMinFxy(Set<Integer> OmegaX,
//                                       Set<Integer> OmegaY, double lambda) {
//        // TODO Auto-generated method stub
//        return null;
//    }

	@Override
	public double[] getGradientX(double[] x, double[] y, double lambda, double d) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] getGradientY(double[] x, double[] y, double lambda, double d) {
		// TODO Auto-generated method stub
		return null;
	}
}
