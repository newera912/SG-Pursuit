package edu.albany.cs.S2GraphMP;

import com.google.common.collect.Sets;
import com.google.common.primitives.Doubles;
import edu.albany.cs.base.ArrayIndexComparator;
import edu.albany.cs.base.PreRec;
import edu.albany.cs.dataProcess.GenerateSimuGrid;
import edu.albany.cs.scoreFuncs.Function;
import edu.albany.cs.scoreFuncs.LogLogisticRegression;
import edu.albany.cs.scoreFuncs.Stat;
import org.apache.commons.math3.linear.ArrayRealVector;

import java.io.FileWriter;
import java.io.IOException;
import java.util.*;

/**
 * Created by baojian bzhou6@albany.edu on 1/28/17.
 */
public class TestArgMin {

    public static void testOptimalCase() {
        GenerateSimuGrid gridData = GenerateSimuGrid.getGridObj("./data/SerializedData/dataset0/Grid_np_100_10_z_20.0_r_5_mu_20.0_.ser");
        Function logisticRegression = new LogLogisticRegression(gridData.W, gridData.b);
        double[] xStar = gridData.x;
        double[] yStar = gridData.y;
        System.out.println("optimal value: f(x,y): " + logisticRegression.getFuncValue(xStar, yStar));
        Set<Integer> omegaSx = gridData.Sx;
        Set<Integer> omegaSy = gridData.Sy;
        List<double[]> results = logisticRegression.getArgMinFxy(omegaSx, omegaSy);
        System.out.println("found minimum value: " + logisticRegression.getFuncValue(results.get(0), results.get(1)));
    }

    public static List<PreRec> testCase1(int k, String filePath) {
        System.out.println("----------------------------------------------------------------------------------------");
        GenerateSimuGrid gridData = GenerateSimuGrid.getGridObj(filePath);
        Function func = new LogLogisticRegression(gridData.W, gridData.b, gridData.Sy.size() + k);
        double[] xStar = gridData.x;
        double[] yStar = gridData.y;
        System.out.println("optimal value: f(x,y): " + func.getFuncValue(xStar, yStar));
        System.out.println("f(0,y): " + func.getFuncValue(new double[gridData.n], gridData.y));
        System.out.println("f(0,0): " + func.getFuncValue(new double[gridData.n], new double[gridData.p]));
        System.out.println("gradientX(0,y): " + Doubles.join(" ", func.getGradientX(new double[gridData.n], gridData.y)));
        Set<Integer> omegaSx = Sets.union(gridData.Sx, generateRandomNormalNodes(gridData.Sx, 10, gridData.n));
        Set<Integer> omegaSy = Sets.union(gridData.Sy, generateRandomNormalNodes(gridData.Sy, 10, gridData.p));
        omegaSx = new HashSet<>();
        omegaSy = new HashSet<>();
        for (int i = 0; i < gridData.n; i++) {
            omegaSx.add(i);
        }
        for (int i = 0; i < gridData.p; i++) {
            omegaSy.add(i);
        }
        System.out.println("true Sy: " + gridData.Sy);
        List<double[]> results = func.getArgMinFxy(omegaSx, omegaSy);
        System.out.println("-----------------------");
        double minVal = func.getFuncValue(results.get(0), results.get(1));
        double optimalVal = func.getFuncValue(gridData.x, gridData.y);
        double[] x = results.get(0);
        double[] y = results.get(1);
        System.out.println("found x: " + Doubles.join(" ", x));
        System.out.println("found y: " + Doubles.join(" ", y));
        System.out.println("found minimum value: " + minVal + " vs optimal value: " + optimalVal);
        System.out.println("true Sx: " + gridData.Sx);
        System.out.println("supp(x): " + Stat.supp(results.get(0)));
        System.out.println(new PreRec(Stat.supp(results.get(0)), gridData.Sx));
        System.out.println("supp(y): " + Stat.supp(results.get(1)));
        System.out.println("y: " + Doubles.join(" ", results.get(1)));
        System.out.println("true Sy: " + gridData.Sy);
        System.out.println(new PreRec(Stat.supp(results.get(1)), gridData.Sy));
        System.out.println("||gradientX||: " + new ArrayRealVector(func.getGradientX(x, y)).getNorm());
        System.out.println("||gradientY||: " + new ArrayRealVector(func.getGradientY(x, y)).getNorm());
        System.out.println("gradientY: " + Doubles.join(" ", func.getGradientY(x, y)));
        System.out.println("norm(x-x*): " + Stat.getL2Norm(x, gridData.x));
        for (int i = 0; i < gridData.n; i++) {
            if (gridData.Sx.contains(i)) {
                System.out.print("[" + i + "]" + gridData.x[i] + " ");
            }
        }
        System.out.println();
        Set<Integer> suppX = Stat.supp(results.get(0));
        for (int i = 0; i < gridData.n; i++) {
            if (suppX.contains(i)) {
                System.out.print("[" + i + "]" + x[i] + " ");
            }
        }
        System.out.println();
        if (Stat.getL2Norm(x, gridData.x) > 0.005D) {
            System.out.println(Stat.getL2Norm(x, gridData.x) + "bound is too large.");
            //Utils.stop();
            //return true;
        }

        PreRec preRec_nodes = new PreRec(getTopKNodes(gridData.Sx.size(), x), gridData.Sx);
        PreRec preRec_features = new PreRec(Stat.supp(results.get(1)), gridData.Sy);
        List<PreRec> res = new ArrayList<>();
        res.add(preRec_nodes);
        res.add(preRec_features);
        System.out.println(res);
        return res;
    }

    public static Set<Integer> generateRandomNormalNodes(Set<Integer> abnormalNodes, int size, int dimension) {
        Set<Integer> randomNormalNodes = new HashSet<>();
        Random random = new Random();
        while (true) {
            if (randomNormalNodes.size() == size) {
                return randomNormalNodes;
            }
            int randomNode = random.nextInt(dimension);
            if (!abnormalNodes.contains(randomNode)) {
                randomNormalNodes.add(randomNode);
            }
        }
    }


    public static Set<Integer> getTopKNodes(int k, double[] x) {
        for (int i = 0; i < x.length; i++) {
            x[i] = Math.abs(x[i]);
        }
        Set<Integer> results = new HashSet<>();
        ArrayIndexComparator arrayIndexComparator = new ArrayIndexComparator(x);
        Integer[] indexes = arrayIndexComparator.indexes;
        Arrays.sort(indexes, arrayIndexComparator);
        for (int i = 0; i < k; i++) {
            results.add(indexes[i]);
        }
        return results;
    }


    public static double[] getAverage(List<List<PreRec>> results) {
        double pre_1 = 0.0D;
        double rec_1 = 0.0D;
        double fmeasure_1 = 0.0D;
        double pre_2 = 0.0D;
        double rec_2 = 0.0D;
        double fmeasure_2 = 0.0D;
        for (List<PreRec> result : results) {
            pre_1 += result.get(0).pre;
            rec_1 += result.get(0).rec;
            fmeasure_1 += result.get(0).fmeasure;
            pre_2 += result.get(1).pre;
            rec_2 += result.get(1).rec;
            fmeasure_2 += result.get(1).fmeasure;
        }
        double[] result = new double[]{pre_1, rec_1, fmeasure_1, pre_2, rec_2, fmeasure_2};
        for (int i = 0; i < result.length; i++) {
            result[i] = result[i] / (results.size() * 1.0D);
        }
        return result;
    }

    public static void main(String args[]) throws IOException {

        for (int k = 15; k < 20; k++) {

            for (double mu : new double[]{1.0D, 2.0D, 3.0D, 4.0D, 5.0D, 10.0D, 15.0D, 20.0D}) {
                FileWriter fileWriter = new FileWriter("./outputs/SerializedData/output_900_200_0.05_apdm02.txt", true);
                List<List<PreRec>> results = new ArrayList<>();
                for (int i = 0; i < 10; i++) {
                    String filePath = "./data/SerializedData/dataset0/Grid_np_900_200_z_45.0_r_10_mu_" + mu + "_.ser";
                    results.add(testCase1(k, filePath));
                }
                for (List<PreRec> result : results) {
                    System.out.println(result);
                }
                fileWriter.write(mu + " " + k + " " + Doubles.join(" ", getAverage(results)) + "\n");
                fileWriter.close();
                System.out.println(Doubles.join(" ", getAverage(results)));
                //Utils.stop();
            }
        }

    }

}
