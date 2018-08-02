package edu.albany.cs.opt;


/**
 * % provide the equation you want to solve with R.H.S = 0 form.
 * % Write the L.H.S by using inline function
 * % Give initial guesses.
 * % Solves it by method of bisection.
 * % A very simple code. But may come handy
 * Created by baojian bzhou6@albany.edu on 1/28/17.
 */
public class BisectionMethod {

    private final Function func;
    private double lowerBound;
    private double upperBound;
    private final int maximumIter = 1000000;
    private final double epsilon = 1e-8D;
    private final int verboseLevel = 0;

    public double root;

    public BisectionMethod(Function func, double lowerBound, double upperBound) {
        this.func = func;
        this.lowerBound = lowerBound;
        this.upperBound = upperBound;
        if (!run()) {
            System.out.println("something is wrong.");
            System.exit(0);
        }
    }

    private boolean run() {
        double middle = Double.NaN;
        double diff = upperBound - lowerBound;
        if (func.getFuncValue(lowerBound) * func.getFuncValue(upperBound) >= 0.0D) {
            System.out.println("f(lowerBound): " + func.getFuncValue(lowerBound));
            System.out.println("f(upperBound): " + func.getFuncValue(upperBound));
            System.out.println("improper lower bound or upper bound.");
            System.exit(0);
        }
        int numOfIter = 0;
        while ((Math.abs(upperBound - lowerBound) > epsilon) && (numOfIter < maximumIter)) {
            if (verboseLevel == 0) {
                System.out.println("updated diff: " + diff);
            }
            middle = (lowerBound + upperBound) / 2.0D;
            if (func.getFuncValue(lowerBound) * func.getFuncValue(middle) < 0.0D) {
                upperBound = middle;
                diff = upperBound - lowerBound;
            } else {
                lowerBound = middle;
                diff = upperBound - lowerBound;
            }
            numOfIter++;
        }
        if (numOfIter >= maximumIter || Double.isNaN(middle)) {
            return false;
        }
        root = middle;
        return true;
    }

    public static void main(String args[]) {
        Function func = new Function(null, null, null);
        BisectionMethod bi = new BisectionMethod(func, 0.0D, 3.0D);
        System.out.println(bi.root);
        System.out.println("f(x): " + func.getFuncValue(2.6510934066027403));
    }
}
