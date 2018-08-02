package edu.albany.cs.LassoGraphMP;

import com.google.common.collect.Sets;
import edu.albany.cs.base.ArrayIndexSort;
import edu.albany.cs.graph.MultiDataGraph;
import edu.albany.cs.graphIHT.GraphIHT;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.Stat;
import org.apache.commons.math3.linear.ArrayRealVector;

import java.util.Arrays;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

/**
 * Created by baojian on 1/14/17.
 * Email: bzhou6@albany.edu
 */
public class LassoGraphMP {

    //Graph information
    private final MultiDataGraph graph;
    private final Function func;
    // maximum subgraph size
    private final int k;
    // maximum size of selected features
    private final int s;
    // number of connected components
    private final int g;
    // the total budget
    private final double B;
    // default number of iterations
    private final int t = 5;

    //The estimated subset of nodes and subset of features.
    private Set<Integer> resultNodes;
    private Set<Integer> resultFeatures;

    private double funValue;

    
    public LassoGraphMP(MultiDataGraph graph, Function func, int k, int s, int g, double B, int version) {
        this.graph = graph;
        this.func = func;
        this.s = s;
        this.k = k;
        this.g = g;
        this.B = B;       
		
        if (!run(version)) {
            System.out.println("Lasso GraphMP fails.");
        }

    }


    private boolean run(int version) {
        double[] x;
        double[] y = initializeX_Random(graph.numOfFeas, s);
        int numIteration = 0;
        while (true) {
          
            x = argMinFx(y); 
//            Arrays.fill(x, 0.0D);
//            for(int i:graph.trueSubGraphNodes){
//            	x[i]=1.0D;
//            }
            /** use the first version : projected gradient descent. */
            if (version == 1) {
                y = projectedGradientDescent(x);
            /** use the second version : gradient support pursuit. */
            } else {
                y = gradientSupportPursuit(x, s);
            }
           
            numIteration++;
            
            //halting condition can be added here.
            if (numIteration == t) {
                break;
            }
        }

        resultNodes = this.support(x);
        resultFeatures = this.support(y);
        funValue = func.getFuncValue(x, y);
        return true;
    }

    /**
     * Use Graph-IHT to solve the problem of Line 5
     */
    private double[] argMinFx(double[] yi) {
        GraphIHT graphIHT = new GraphIHT(graph.edges, graph.edgeCosts, graph.numOfNodes, func, yi, k, g, B, graph.trueSubGraphNodes);
        return graphIHT.x.toArray();
    }

    /**
     * Use projected gradient descent to solve the problem of Line 6.
     */
    private double[] projectedGradientDescent(double[] xi) {
        //TODO: change yt to be a feasible solution. Done!
        ArrayRealVector yt = new ArrayRealVector(initializeX_Random(graph.numOfFeas, s));        
        double stepSize = 0.01D;
        int iter = 0;
        while (true) {        	
            ArrayRealVector gradientXt = new ArrayRealVector(func.getGradientY(xi, yt.toArray()));
            ArrayRealVector ytPlus1 = yt.add(gradientXt.mapMultiply(-stepSize));
            //This is the projection step.
            double[] absYtPlus1 = new double[ytPlus1.getDimension()];
            for (int i = 0; i < ytPlus1.getDimension(); i++) {
                if (ytPlus1.getEntry(i) >= 0.0D) {
                    absYtPlus1[i] = ytPlus1.getEntry(i);
                } else {
                    absYtPlus1[i] = 0.0D;
                }
            }

            ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(absYtPlus1);
            Integer[] indexes = arrayIndexComparator.getIndices();
            Arrays.sort(indexes, arrayIndexComparator);
            Set<Integer> selectedIndices = new HashSet<>();
            for (int i = 0; i < s; i++) {
                selectedIndices.add(indexes[i]);
            }
            for (int i = 0; i < ytPlus1.getDimension(); i++) {
                if (!selectedIndices.contains(i)) {
                    ytPlus1.setEntry(i, 0.0D);
                }
                if (ytPlus1.getEntry(i) > 1.0D) {
                    ytPlus1.setEntry(i, 1.0D);
                }
                if (ytPlus1.getEntry(i) < 0.0D) {
                    ytPlus1.setEntry(i, 0.0D);
                }
            }
            double gap = Math.abs(func.getFuncValue(xi, yt.toArray()) - func.getFuncValue(xi, ytPlus1.toArray()));
            yt = ytPlus1;

            if (gap < 1e-3D) {
                System.out.println("Gap < 0.001");
                break;
            }
            if (iter > 300) {
                break;
            }

            iter++;
        }
        return yt.toArray();
    }


    public double[] gradientSupportPursuit(double[] xi, int s) {
        //ArrayRealVector yHat = new ArrayRealVector(new double[this.graph.numOfFeas]);
        ArrayRealVector yHat = new ArrayRealVector(Stat.getRandomVector(this.graph.numOfFeas, s));        
        while (true) {
            // compute local gradient
            double[] localGradient = func.getGradientY(xi, yHat.toArray());            
            // identify directions
            Set<Integer> Z = identifyDirection(localGradient, 2 * s);            
            // merge supports
            Set<Integer> I = Sets.union(Z, support(yHat.toArray()));             
            //\haty_t+1
            ArrayRealVector yHattPlus1=new ArrayRealVector(graph.numOfFeas);
            // minimize over supports
            double[] b = argMin(xi, I);            
            Set<Integer> direction = identifyDirection(b, s);            
            for (int i = 0; i < b.length; i++) {
                if (!direction.contains(i)) {
                	yHattPlus1.setEntry(i, 0.0D);
                } else if(b[i]<0.0D) {  /** a.alim added   Project to [0,1] */
                	yHattPlus1.setEntry(i, 0.0D);
                }else if(b[i]>1.0D) {
                	yHattPlus1.setEntry(i, 1.0D);
                }else{
                	yHattPlus1.setEntry(i, b[i]);
                }
            }           
            /** a.alim Added terminate condition*/
            double gap = Math.abs(func.getFuncValue(xi, yHat.toArray()) - func.getFuncValue(xi, yHattPlus1.toArray()));
            yHat=yHattPlus1;
            //Utils.stop();
            if (gap < 1e-3D) {
                break;
            }
            
        }
        return yHat.toArray();
    }

    private Set<Integer> identifyDirection(double[] gradient, int s) {
        double[] absGradient = new double[gradient.length];
        for(int i = 0 ; i < gradient.length ;i++){
            absGradient[i] = Math.abs(gradient[i]);
        }
        ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(absGradient);
        Integer[] indexes = arrayIndexComparator.getIndices();
        Arrays.sort(indexes, arrayIndexComparator);
        Set<Integer> Z = new HashSet<>();        
        for (int i = 0; i < s; i++) {
            Z.add(indexes[i]);
        }
        return Z;
    }

    /**
     * Find the minimizer using gradient descent.
     *
     * @param I
     * @return
     */
    private double[] argMin(double[] xi, Set<Integer> I) {
        ArrayRealVector minimizer;

        double stepSize = 0.001D;

        int iter=0;
        int MAX_iter=1000;
        //initialize a feasible solution.
        ArrayRealVector yt = new ArrayRealVector(new double[graph.numOfFeas]);
        Random random = new Random();
        for (int i = 0; i < yt.getDimension(); i++) {
            if (!I.contains(i)) {  //a.alim added "!"
                yt.setEntry(i, 0.0D);
            } else {
                yt.setEntry(i, random.nextDouble());
            }
        }
        ArrayRealVector ytPlus1;
        while (true) {            
            double[] gradientYt = func.getGradientY(xi, yt.toArray());           
            ytPlus1 = yt.subtract(new ArrayRealVector(gradientYt).mapMultiplyToSelf(stepSize));
            for (int i = 0; i < graph.numOfFeas; i++) {
                if (!I.contains(i)) {
                    ytPlus1.setEntry(i, 0.0D);
                }else{
                    if(ytPlus1.getEntry(i) > 1.0D){
                        ytPlus1.setEntry(i,1.0D);
                    }
                    if(ytPlus1.getEntry(i) < 0.0D){
                        ytPlus1.setEntry(i,0.0D);
                    }
                }
            }
            double gap = Math.abs(func.getFuncValue(xi, ytPlus1.toArray()) - func.getFuncValue(xi, yt.toArray()));           
            yt=ytPlus1;  //a.alim added
            if (gap < 1e-3D) {
                break;
            }
            if(iter>MAX_iter){
            	break;
            }
            iter++;
        }
        //Utils.stop();
        minimizer = ytPlus1;
        return minimizer.toArray();
    }


    private Set<Integer> support(double[] x) {
        Set<Integer> supp = new HashSet<>();
        for (int i = 0; i < x.length; i++) {
            if (x[i] != 0.0D) {
                supp.add(i);
            }
        }
        return supp;
    }


    public Set<Integer> getResultNodes() {
        return resultNodes;
    }

    public Set<Integer> getResultFeatures() {
        return resultFeatures;
    }

    public double getFunValue() {
        return funValue;
    }

    private double[] initializeX_Random(int size, int s) {
        return Stat.getRandomVector(size, s);
    }

    public static void main(String args[]) {

    }
}