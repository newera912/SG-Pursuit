package edu.albany.cs.scoreFuncs;

import edu.albany.cs.base.Matrix;
import org.apache.commons.math3.stat.StatUtils;

import java.util.List;
import java.util.Set;

/**
 * Created by baojian on 2/13/17.
 */
public class ElevatedMeanScan implements Function {
    private final double[][] W;
    private final int numNodes;
    private final int numFea;
    private final FuncType funcID;

    public ElevatedMeanScan(double[][] W) {
        funcID = FuncType.EMSXYScore;
        if (!checkInput(W)) {
            System.out.println(funcID + " input parameter is invalid.");
        }
        this.W = W;
        this.numNodes = W.length;
        this.numFea = W[0].length;
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
        funcValue = -(xTWy * xTWy) / (xT1) + (0.5D * yT1 * yT1);
        return funcValue;
    }

    @Override
    public double[] getGradientX(double[] x, double[] y) {
        if (y.length != numFea && x.length != numNodes) {
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
        double[] gradient = new double[numNodes];
        double xT1 = StatUtils.sum(x);
        double[] yW = Matrix.MatMultiplyVec(W, y);
        double[] term1 = Matrix.VarMultiplyVec(2.0D / xT1, yW);
        double xTWy = Matrix.dot(x, yW);
        double term2 = xTWy / (xT1 * xT1);
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = -xTWy * (term1[i] - term2);
        }
        return gradient;
    }

    @Override
    public double[] getGradientY(double[] x, double[] y) {
        if (y.length != numFea && x.length != numNodes) {
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
        double[] gradient = new double[numFea];
        double xT1 = StatUtils.sum(x);
        double[] xW = Matrix.MatMultiplyVec(Matrix.transpose(W), x);
        double xTWy = Matrix.dot(xW, y);
        double[] term2 = Matrix.VarMultiplyVec((2.0D * xTWy) / xT1, xW);
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = -term2[i] + y[i];
        }
        return gradient;
    }

    @Override
    public List<double[]> getArgMinFxy(Set<Integer> OmegaX, Set<Integer> OmegaY) {
        return new GradientDescentOpt(this, numNodes, numFea, 0.01D, 1000000).multiGradientDescent4EMSScore(OmegaX, OmegaY, -1);
    }

    @Override
    public List<double[]> getArgMinFx(double[] yi) {
        return null;
    }

    @Override
    public List<double[]> getArgMinFy(double[] xi) {
        return null;
    }

    @Override
    public double getFuncValue(double[] x, double[] y, double lambda) {
        return 0;
    }

    @Override
    public double[] getGradientX(double[] x, double[] y, double lambda) {
        return new double[0];
    }

    @Override
    public double[] getGradientY(double[] x, double[] y, double lambda) {
        return new double[0];
    }

    @Override
    public List<double[]> getArgMinFxy(double[] x0, double[] y0, Set<Integer> OmegaX, Set<Integer> OmegaY, double lambda) {
        return null;
    }

    @Override
    public double[] getArgMinFy(double[] xi, Set<Integer> OmegaY, double lambda) {
        return new double[0];
    }

    @Override
    public double[] getArgMinFx(Set<Integer> trueNodesSet, int[] trueFeas, double[] yi, Set<Integer> omegaX, double lambda) {
        return new double[0];
    }

    @Override
    public double[] getGradientX(double[] x, double[] y, double lambda, double d) {
        return new double[0];
    }

    @Override
    public double[] getGradientY(double[] x, double[] y, double lambda, double d) {
        return new double[0];
    }
}
