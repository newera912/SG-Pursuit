package edu.albany.cs.base;

/**
 * Created by baojian bzhou6@albany.edu on 2/13/17.
 */
public final class Matrices {
    private Matrices() {

    }

    public static double[] getColumn(double[][] mat, int colIndex) {
        double[] col = new double[mat.length];
        for (int i = 0; i < mat.length; i++) {
            col[i] = mat[i][colIndex];
        }
        return col;
    }
}
