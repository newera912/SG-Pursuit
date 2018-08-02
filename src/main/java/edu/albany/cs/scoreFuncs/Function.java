package edu.albany.cs.scoreFuncs;

import java.util.List;
import java.util.Set;


/**
 * Inferface of a general function, this kind of functions is differential.
 *
 * @author Baojian bzhou@albany.edu
 */
public interface Function {

    /**
     * @return the function ID
     */
    FuncType getFuncID();

    double getFuncValue(double[] x, double[] y);

    double[] getGradientX(double[] x, double[] y);

    double[] getGradientY(double[] x, double[] y);

    List<double[]> getArgMinFxy(Set<Integer> OmegaX, Set<Integer> OmegaY);

    List<double[]> getArgMinFx(double[] yi);

    List<double[]> getArgMinFy(double[] xi);


    /**
     * Demnse function
     */
    double getFuncValue(double[] x, double[] y, double lambda);

    double[] getGradientX(double[] x, double[] y, double lambda);

    double[] getGradientY(double[] x, double[] y, double lambda);

    List<double[]> getArgMinFxy(double[] x0, double[] y0, Set<Integer> OmegaX, Set<Integer> OmegaY, double lambda);

    double[] getArgMinFy(double[] xi, Set<Integer> OmegaY, double lambda);

    double[] getArgMinFx(Set<Integer> trueNodesSet, int[] trueFeas,
                         double[] yi, Set<Integer> omegaX, double lambda);

    double[] getGradientX(double[] x, double[] y, double lambda, double d);

    double[] getGradientY(double[] x, double[] y, double lambda, double d);
}
