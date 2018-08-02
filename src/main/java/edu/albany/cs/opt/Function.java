package edu.albany.cs.opt;

/**
 * Created by baojian bzhou6@albany.edu on 1/30/17.
 */
public class Function {

    private final double[][] w;
    private final double[] x;
    private final double[] y;

    public Function(double[][] w, double[] x, double[] y) {
        this.w = w;
        this.x = x;
        this.y = y;
    }

    public double getFuncValue(double z) {
        return (-1 + z * (-4 + z * (-1.0D + z)));
    }

}
