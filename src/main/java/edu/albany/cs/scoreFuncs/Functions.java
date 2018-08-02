package edu.albany.cs.scoreFuncs;

import org.apache.commons.math3.linear.ArrayRealVector;

/**
 * Created by baojian bzhou6@albany.edu on 2/5/17.
 */
public final class Functions {

    private Functions() {
    }

    public static double sigmoid(double x) {
        return 1.0D / (1.0D + Math.exp(-x));
    }

    public static double l2norm(double[] x) {
        return new ArrayRealVector(x).getNorm();
    }

    public static double l2norm(double[] x1, double[] x2) {
        ArrayRealVector x1Vector = new ArrayRealVector(x1);
        ArrayRealVector x2Vector = new ArrayRealVector(x2);
        ArrayRealVector diff = x1Vector.subtract(x2Vector);
        return diff.getNorm();
    }

    public static void main(String args[]) {
        System.out.println(Functions.sigmoid(0.0D));
        System.out.println(Functions.sigmoid(10000.0D));
        System.out.println(Functions.sigmoid(-10000.0D));
    }

}
