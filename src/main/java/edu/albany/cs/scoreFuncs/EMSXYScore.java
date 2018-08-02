package edu.albany.cs.scoreFuncs;


import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.Matrix;
import edu.albany.cs.base.Utils;
import edu.albany.cs.graph.MultiDataGraph;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;

import java.io.*;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/*************************************************************\
 *max f(x,y)=-(x^T W y)^2 / (x^T 1 y^T 1)
 *
 \*************************************************************/
public class EMSXYScore implements Function {


    /**
     * weight matrix
     */
    private final double[][] W;
    /**
     * dimension of a vector
     */
    private final int numNodes;
    private final int numFea;
    private final FuncType funcID;

    public EMSXYScore(double[][] W) {
        funcID = FuncType.EMSXYScore;
        if (!checkInput(W)) {
            System.out.println(funcID + " input parameter is invalid.");
        }

        this.W = W;
        this.numNodes = W.length;
        this.numFea = W[0].length;
    }

    private boolean checkInput(double[][] W) {
        if (W == null)
            return false;
        return true;
    }

    /**
     * delta_x f(x,y)=-x^TWy.{2*W.Y/(x^T1*y^T1)-x^TWy/((x^T1)^2*y^T1)}
     */
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
        double yT1 = StatUtils.sum(y);
        double[] yW = Matrix.MatMultiplyVec(W, y);
        double[] term1 = Matrix.VarMultiplyVec(2.0D / (xT1 * yT1), yW); //(2.0/xT1*yT1)*W.Y
        double xTWy = Matrix.dot(x, yW); //xTWy
        double term2 = xTWy / (xT1 * xT1 * yT1);
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = -xTWy * (term1[i] - term2);
        }
        return gradient;
    }

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
        double yT1 = StatUtils.sum(y);
        double[] xW = Matrix.MatMultiplyVec(Matrix.transpose(W), x);
        double[] term1 = Matrix.VarMultiplyVec(2.0D / (xT1 * yT1), xW);
        double xTWy = Matrix.dot(xW, y);
        double term2 = xTWy / (xT1 * yT1 * yT1);
        for (int i = 0; i < gradient.length; i++) {
            gradient[i] = -xTWy * (term1[i] - term2);
        }
        return gradient;
    }


    @Override
    public List<double[]> getArgMinFx(double[] yi) {
        return null;
    }

    @Override
    public List<double[]> getArgMinFy(double[] xi) {
        return null;
    }

    public double getFuncValue(double[] x, double[] y) {
        double funcValue;
        if (x.length != W.length || y.length != W[0].length) {
            System.out.println("Error : Invalid parameters ...");
            System.exit(0);
        }
        double xT1 = StatUtils.sum(x);
        double yT1 = StatUtils.sum(y);
        double xTWy = Matrix.dot(Matrix.VecMultiplyMat(x, W), y);
        funcValue = -xTWy * xTWy / (xT1 * yT1);
        return funcValue;
    }

    public static void test() {
        int[] x_true = new int[]{54, 53, 43, 64, 63, 44, 45, 55, 56, 46};
        int[] y_true = new int[]{0, 1, 2, 3, 4};
        double[] x = new double[100];
        double[] y = new double[10];
        Arrays.fill(x, 0.0D);
        Arrays.fill(y, 0.0D);

        for (int i : x_true) {
            x[i] = 1.0D;
        }

        for (int i : y_true) {
            y[i] = 1.0D;
        }
        APDMInputFormat apdm = new APDMInputFormat("data/SimulationData/truesub_10/multiGridData_mu10.0/TestScoreFunction_APDM-GridData-100-precen-0.1-noise_0.0-0.txt");
        MultiDataGraph mgdGraph = new MultiDataGraph(apdm);

        System.out.println(ArrayUtils.toString(apdm.data.trueSubGraphNodes));
        EMSXYScore func2 = new EMSXYScore(mgdGraph.W);

        System.out.println("fun:" + func2.getFuncValue(x, y));
        System.out.println("Graient X:" + ArrayUtils.toString(func2.getGradientX(x, y)));
        System.out.println("Graient Y:" + ArrayUtils.toString(func2.getGradientY(x, y)));

        Function func3 = new EMSXYScore(mgdGraph.W);
        System.out.println("fun:" + func3.getFuncValue(x, y));
        System.out.println("Graient X:" + ArrayUtils.toString(func3.getGradientX(x, y)));
        System.out.println("Graient Y:" + ArrayUtils.toString(func3.getGradientY(x, y)));
    }

    @Override
    public FuncType getFuncID() {
        return funcID;
    }


    public static void unittest1() {
        double[] xi = null;
        double[] yi = null;
        Set<Integer> omegaX = new HashSet<>();
        Set<Integer> omegaY = new HashSet<>();
        String fileName = "";
        double lambda = 0;
        try {
            FileInputStream istream = new FileInputStream("testcases/getArgMinFxy/test-case1.dat");
            ObjectInputStream ois = new ObjectInputStream(istream);
            xi = (double[]) ois.readObject();
            yi = (double[]) ois.readObject();
            omegaX = Utils.intArray2Set((int[]) ois.readObject());
            omegaY = Utils.intArray2Set((int[]) ois.readObject());
            lambda = (double) ois.readObject();
            fileName = (String) ois.readObject();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
        APDMInputFormat apdm = new APDMInputFormat(fileName);
        MultiDataGraph mgdGraph = new MultiDataGraph(apdm);
        PCAScore func = new PCAScore(mgdGraph.W, mgdGraph.AdjMatrix);
        List<double[]> minXY = func.getArgMinFxy(xi, yi, omegaX, omegaY, lambda);
    }

    public static void main(String[] args) {
        unittest1();
    }


    private void storeTestCaseData(double[] xi, double[] yi, Set<Integer> omegaX, Set<Integer> omegaY, double lambda, String fileName, String outputfilename) {
        try {
            FileOutputStream fos = new FileOutputStream(outputfilename);
            ObjectOutputStream oos = new ObjectOutputStream(fos);
            oos.writeObject(xi); // write MenuArray to ObjectOutputStream
            oos.writeObject(yi);
            oos.writeObject(Utils.set2IntArray(omegaX));
            oos.writeObject(Utils.set2IntArray(omegaY));
            oos.writeObject(lambda);
            oos.writeObject(fileName);
            oos.close();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
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
        return new GradientDescentOpt(this, numNodes, numFea).multiGradientDescent(x0, y0, OmegaX, OmegaY, lambda);
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
    public List<double[]> getArgMinFxy( Set<Integer> OmegaX, Set<Integer> OmegaY) {
        return new GradientDescentOpt(this, numNodes, numFea, 0.01D, 1000000).multiGradientDescent4EMSScore(OmegaX, OmegaY, -1);
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

