package edu.albany.cs.S2GraphMPDebug;

import edu.albany.cs.scoreFuncs.FuncType;

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

    List<double[]> getArgMinFxy(double[] xi, double[] yi, Set<Integer> OmegaX, Set<Integer> OmegaY);

}
