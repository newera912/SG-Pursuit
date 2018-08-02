package edu.albany.cs.base;

import org.apache.commons.lang3.ArrayUtils;


public final class Matrix {
    private final int n;
    private final int m;
    private final double[][] mat;

    public Matrix(double[][] mat) {
        this.mat = mat;
        this.n = mat.length;
        this.m = mat[0].length;
    }

    public double[] getColumn(int colIndex) {
        double[] col = new double[n];
        for (int i = 0; i < n; i++) {
            col[i] = mat[i][colIndex];
        }
        return col;
    }

    public void setColumn(int colIndex, double[] col) {
        for (int i = 0; i < n; i++) {
            mat[i][colIndex] = col[i];
        }
    }

    // return a random m-by-n matrix with values between 0 and 1
    public static double[][] random(int m, int n) {
        double[][] a = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                a[i][j] = Math.random();
        return a;
    }

    // return n-by-n identity matrix I
    public static double[][] identity(int n) {
        double[][] a = new double[n][n];
        for (int i = 0; i < n; i++)
            a[i][i] = 1.0D;
        return a;
    }

    // return n-by-1 identity matrix I
    public static double[] identityVector(int n) {
        double[] a = new double[n];
        for (int i = 0; i < n; i++)
            a[i] = 1.0D;
        return a;
    }


    // return x^T y
    public static double dot(double[] x, double[] y) {
        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double sum = 0.0D;
        for (int i = 0; i < x.length; i++)
            sum += x[i] * y[i];
        return sum;
    }

    public static double[] elemWisePro(double[] x, double[] y) {

        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double[] c = new double[x.length];
        for (int i = 0; i < x.length; i++)
            c[i] = x[i] * y[i];
        return c;
    }

    public static double[][] VecOutterPro(double[] x, double[] y) {

        if (x == null || y == null) throw new RuntimeException("Illegal vector dimensions.");
        double[][] c = new double[x.length][y.length];
        for (int i = 0; i < x.length; i++) {
            for (int j = 0; j < y.length; j++) {
                c[i][j] = x[i] * y[j];
            }
        }
        return c;
    }

    public static double[] VecSubstract(double[] x, double[] y) {

        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double[] c = new double[x.length];
        for (int i = 0; i < x.length; i++)
            c[i] = x[i] - y[i];
        return c;
    }

    public static double[] VecAdd(double[] x, double[] y) {

        if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double[] c = new double[x.length];
        for (int i = 0; i < x.length; i++)
            c[i] = x[i] + y[i];
        return c;
    }

    // return x^T y
    public static double[] VarMultiplyVec(double a, double[] x) {
        //if (x.length != y.length) throw new RuntimeException("Illegal vector dimensions.");
        double[] ax = new double[x.length];
        for (int i = 0; i < x.length; i++)
            ax[i] = a * x[i];
        return ax;
    }

    // return B = A^T
    public static double[][] transpose(double[][] a) {
        int m = a.length;
        int n = a[0].length;
        double[][] b = new double[n][m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                b[j][i] = a[i][j];
        return b;
    }

    // return c = a + b
    public static double[][] add(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] + b[i][j];
        return c;
    }

    // return c = a - b
    public static double[][] subtract(double[][] a, double[][] b) {
        int m = a.length;
        int n = a[0].length;
        double[][] c = new double[m][n];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                c[i][j] = a[i][j] - b[i][j];
        return c;
    }

    // return c = a * b
    public static double[][] multiply(double[][] a, double[][] b) {
        int m1 = a.length;
        int n1 = a[0].length;
        int m2 = b.length;
        int n2 = b[0].length;
        if (n1 != m2) throw new RuntimeException("Illegal matrix dimensions.");
        double[][] c = new double[m1][n2];
        for (int i = 0; i < m1; i++)
            for (int j = 0; j < n2; j++)
                for (int k = 0; k < n1; k++)
                    c[i][j] += a[i][k] * b[k][j];
        return c;
    }

    // matrix-vector multiplication (y = A * x)
    public static double[] MatMultiplyVec(double[][] a, double[] x) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != n) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[m];
        for (int i = 0; i < m; i++)
            for (int j = 0; j < n; j++)
                y[i] += a[i][j] * x[j];
        return y;
    }


    // vector-matrix multiplication (y = x^T A)
    public static double[] VecMultiplyMat(double[] x, double[][] a) {
        int m = a.length;
        int n = a[0].length;
        if (x.length != m) throw new RuntimeException("Illegal matrix dimensions.");
        double[] y = new double[n];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                y[j] += a[i][j] * x[i];
        return y;
    }

    // vector-matrix multiplication (y = x^T A)
    public static double[][] VarMultiplyMat(double x, double[][] a) {
        int m = a.length;
        int n = a[0].length;

        double[][] y = new double[n][m];
        for (int j = 0; j < n; j++)
            for (int i = 0; i < m; i++)
                y[j][i] = a[i][j] * x;
        return y;
    }

    public static void main(String args[]) {
        double[][] a = {{1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}};
        double[] b = new double[]{2.0, 2.0, 2.0, 2.0};
        System.out.println(ArrayUtils.toString(Matrix.VecMultiplyMat(b, a)));
        System.out.println(ArrayUtils.toString(Matrix.VecSubstract(b, a[0])));
    }
}    