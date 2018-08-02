package edu.albany.cs.base;

import edu.albany.cs.headApprox.HeadApprox;
import edu.albany.cs.tailApprox.TailApprox;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/**
 * Created by baojian on 1/14/17.
 * Email: bzhou6@albany.edu
 */
public class HeadTailApprox {

    private final int sparsity;
    private final int height;
    private final int width;
    private final int length;
    private final int numCC;
    private final double[] dataVector;

    private ArrayList<Integer[]> edges;
    private ArrayList<Double> costs;

    private double[] result_head_vector;
    private double[] result_tail_vector;
    private int[] result_head_indices;
    private int[] result_tail_indices;

    public HeadTailApprox(int sparsity, int height, int width, int length, int numCC, double[] dataVector) {
        this.sparsity = sparsity;
        this.height = height;
        this.width = width;
        this.length = length;
        this.numCC = numCC;
        this.dataVector = dataVector;
        getEdgesAndCosts();
    }

    public void runHeadApprox() {
        HeadApprox headApprox = new HeadApprox(edges, costs, dataVector,sparsity,numCC, sparsity*1.0D - 1.0D,null, false);
        double[] resultHeadVector = new double[dataVector.length];
        Arrays.fill(resultHeadVector,0.0D);
        result_head_indices = new int[headApprox.bestForest.nodesInF.size()];
        int index = 0;
        for(Integer node : headApprox.bestForest.nodesInF){
            resultHeadVector[node] = dataVector[node];
            result_head_indices[index++] = node;
        }
        result_head_vector = resultHeadVector;
    }

    private void getEdgesAndCosts(){
        edges = new ArrayList<Integer[]>();
        costs = new ArrayList<Double>();
        int index = 0;
        for(int i = 0 ; i < length ; i ++){
            for(int j = 0 ; j < width ; j ++){
                if( index % length != (length - 1)){
                    edges.add(new Integer[]{index,index+1});
                    costs.add(1.0D);
                }else{
                    edges.add(new Integer[]{index,index+1});
                    costs.add(1.0D);
                    if((index + length) >= length*width ){
                    }else{
                        edges.add(new Integer[]{index,index+length});
                        costs.add(1.0D);
                    }
                }
                index++;
            }
        }

        index = 0;
        for(int i = 0 ; i < height ; i++){
            for(int j = 0 ; j < length*width ; j++){
                if((index + length*width) >= height*length*width){
                }else{
                    edges.add(new Integer[]{index, index+length*width});
                    costs.add(1.0D);
                    index ++;
                }
            }
        }
    }

    public void runTailApprox() {
        TailApprox tailApprox = new TailApprox(edges, costs, dataVector,sparsity,numCC, sparsity*1.0D - 1.0D,null, false);
        double[] resultTailVector = new double[dataVector.length];
        Arrays.fill(resultTailVector,0.0D);
        result_tail_indices = new int[tailApprox.bestForest.nodesInF.size()];
        int index = 0;
        for(Integer node : tailApprox.bestForest.nodesInF){
            resultTailVector[node] = dataVector[node];
            result_tail_indices[index++] = node;
        }
        result_tail_vector = resultTailVector;
    }

    public double[] getResultHeadVector() {
        return result_head_vector;
    }

    public int[] getResultHeadIndices() {
        return result_head_indices;
    }

    public double[] getResultTailVector() {
        return result_tail_vector;
    }

    public int[] getResultTailIndices() {
        return result_tail_indices;
    }

    public static void main(String args[]) {
        int height = 3;
        int width = 60;
        int length = 80;
        int n = height*width*length;
        int sparsity = (int)(0.3D*n);
        int numCC = 3;
        double[] dataVector = new double[n];
        Random random = new Random();
        for(int i = 0 ; i < n;i++){
            dataVector[i] = random.nextDouble();
        }
        HeadTailApprox headTailApprox = new HeadTailApprox(sparsity, height, width, length, numCC, dataVector);
        headTailApprox.runHeadApprox();
        headTailApprox.runTailApprox();
        System.out.println(headTailApprox.getResultHeadIndices().length);
        System.out.println(headTailApprox.getResultTailIndices().length);


    }
}
