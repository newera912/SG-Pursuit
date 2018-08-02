package edu.albany.cs.scoreFuncs;

import org.apache.commons.math3.linear.ArrayRealVector;

import java.util.List;
import java.util.Set;

/**
 * This extended logistic loss is convex function.
 * Created by baojian bzhou6@albany.edu on 1/31/17.
 */
public class ExtendedLogisticLoss implements Function {

    private final double[][] W;
    private final double C;
    private final int n;
    private final int p;

    public ExtendedLogisticLoss(double[][] W, int n, int p, double C) {
        this.n = n;
        this.p = p;
        this.W = W;
        this.C = C;
        assert W.length == n;
        assert W[0].length == p;
    }

    @Override
    public FuncType getFuncID() {
        return FuncType.LogisticLoss;
    }

    @Override
    public double getFuncValue(double[] x, double[] y) {
        double funcVal = 0.0D;
        for (int i = 0; i < n; i++) {
            funcVal += 0.0D;
        }
        return 0;
    }

    @Override
    public double[] getGradientX(double[] x, double[] y) {
        return new double[0];
    }

    @Override
    public double[] getGradientY(double[] x, double[] y) {
        return new double[0];
    }

    @Override
    public List<double[]> getArgMinFxy(Set<Integer> OmegaX, Set<Integer> OmegaY) {
        return null;
    }

    @Override
    public List<double[]> getArgMinFx(double[] yi) {
        return null;
    }

    @Override
    public List<double[]> getArgMinFy(double[] xi) {
        return null;
    }

    private double extension_func_d(double[] x, double[] wi, double yi) {
        if (yi == 0.0D || yi == 1.0D) {
            return func_d(x, wi, yi);
        } else {
            return func_psi(x, wi, yi);
        }
    }

    private double func_d(double[] x, double[] wi, double yi) {
        double normX = new ArrayRealVector(x).getNorm();
        double term1 = 0.5D * normX * normX;
        double term2 = C * logistic_loss(wi, x, yi);
        return term1 + term2;
    }

    private double func_psi(double[] x, double[] wi, double yi) {
        return 0.0D;
    }

    private double logistic_loss(double[] wi, double[] x, double yi) {
        ArrayRealVector wiVector = new ArrayRealVector(wi);
        ArrayRealVector xVector = new ArrayRealVector(x);
        double xTransWi = wiVector.dotProduct(xVector);
        return 0.0D;
    }

    private double[] getColumn(int index) {
        double[] col = new double[n];
        for (int i = 0; i < n; i++) {
            col[i] = W[i][index];
        }
        return col;
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

	public List<double[]> getArgMinFxy(Set<Integer> OmegaX,
			Set<Integer> OmegaY, double lambda) {
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
