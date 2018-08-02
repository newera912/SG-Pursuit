package edu.albany.cs.scoreFuncs;

import org.apache.commons.math3.linear.ArrayRealVector;

import java.util.List;
import java.util.Set;

/**
 * Created by baojian on 1/16/17.
 * Email: bzhou6@albany.edu
 */
public class LogisticRegression implements Function {

    /**
     * The number of nodes
     */
    private final int n;
    /**
     * The number of features
     */
    private final int p;
    private final double[] b;
    private final double[][] W;

    public LogisticRegression(double[][] W, double[] b) {
        this.b = b;
        this.n = W.length;
        this.p = W[0].length;
        this.W = W;
    }

    @Override
    public FuncType getFuncID() {
        return FuncType.LogisticRegression;
    }

    @Override
    public double getFuncValue(double[] x, double[] y) {
        ArrayRealVector funcValVector = getFuncValVector(x);
        ArrayRealVector vectorX = new ArrayRealVector(x);
        ArrayRealVector vectorY = new ArrayRealVector(y);
        double norm = funcValVector.subtract(vectorY).getNorm();
        double regularizer = 0.5D * (vectorX.dotProduct(vectorX));
        double funcVal = norm * norm + regularizer;
        return funcVal;
    }

    @Override
    public double[] getGradientX(double[] x, double[] y) {
        ArrayRealVector funcValVector = getFuncValVector(x);
        ArrayRealVector term = new ArrayRealVector(new double[p]);
        for (int i = 0; i < p; i++) {
            double val1 = funcValVector.getEntry(i) * (1.0D - funcValVector.getEntry(i));
            double val2 = funcValVector.getEntry(i) - y[i];
            term.setEntry(i, val1 * val2);
        }

        double[] gradient = new double[n];
        for (int i = 0; i < n; i++) {
            /** gradient + regularizer term*/
            gradient[i] = 2.0D * (new ArrayRealVector(W[i]).dotProduct(term)) + x[i];
        }

        return gradient;
    }

    @Override
    public double[] getGradientY(double[] x, double[] y) {
        ArrayRealVector funcValVector = getFuncValVector(x);
        ArrayRealVector gradient = funcValVector.subtract(new ArrayRealVector(y));
        return gradient.mapMultiplyToSelf(-2.0D).toArray();
    }

    @Override
    public List<double[]> getArgMinFxy( Set<Integer> OmegaX, Set<Integer> OmegaY) {
        return new GradientDescentOpt(this, n, p, 0.01D, 1000000).simpleGradientDescent(OmegaX, OmegaY, null);
    }

    @Override
    public List<double[]> getArgMinFx(double[] yi) {
        return null;
    }

    @Override
    public List<double[]> getArgMinFy(double[] xi) {
        return null;
    }

    private double getSigmoidFuncVal(double x) {
        return 1.0D / (1.0D + Math.exp(-x));
    }

    private double[] getColumn(int index) {
        double[] col = new double[n];
        for (int i = 0; i < n; i++) {
            col[i] = W[i][index];
        }
        return col;
    }

    private ArrayRealVector getFuncValVector(double[] x) {
        ArrayRealVector funcValVector = new ArrayRealVector(new double[p]);
        for (int i = 0; i < p; i++) {
            ArrayRealVector wi = new ArrayRealVector(getColumn(i));
            double wTransX = wi.dotProduct(new ArrayRealVector(x));
            //System.out.println("b[i]: "+b[i]);
            //System.out.println("wTransX:"+wTransX);
            funcValVector.setEntry(i, getSigmoidFuncVal(wTransX + b[i]));
        }
        //System.out.println("funcValVector: "+ Doubles.join(" ",funcValVector.toArray()));
        //Utils.stop();
        return funcValVector;
    }

    public static void test() {

    }

    public static void main(String args[]) {
        LogisticRegression.test();
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
