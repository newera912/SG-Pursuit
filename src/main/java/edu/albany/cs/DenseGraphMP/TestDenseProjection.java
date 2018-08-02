package edu.albany.cs.DenseGraphMP;

import com.google.common.collect.Sets;
import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.graph.Graph;
import edu.albany.cs.graph.MultiDataGraph;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Arrays;
import java.util.HashSet;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class TestDenseProjection {
	
	public static void SingleFileTest(double mu, int k){
		FileWriter fileWriter=null;
		try {
			fileWriter=new FileWriter("outputs/S2GraphMPDense/debug_DenseProjection_results_new.txt", true);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		};
		APDMInputFormat apdm=new APDMInputFormat("data/DenseGraph/Dense_APDM/cluster-10/APDM_r_5_Cluster_10_in_0.35_out_0.05_case_0.txt");
		MultiDataGraph mgdGraph=new MultiDataGraph(apdm);
		double[] lambdas={1.0};//3.0,5.0,7.0,10.0,15,20.0,25.0};
		double gamma=3; //density threshold
		double in_p=0.35;
		double out_p=0.05;
		int t=3;
		int g=1;
		double B=k-g+0.0D;
		double[] gradientFx=new double[apdm.data.numNodes];
		Arrays.fill(gradientFx, 0.0D);
        Set<Integer> trueSet=new HashSet<Integer>();
        
		 for(int q=0;q<gradientFx.length;q++){
         	if(q>=400 && q<500){  
         		trueSet.add(q);
         		gradientFx[q]=new NormalDistribution(mu, 1.0D).sample();					
         	}else{
         		gradientFx[q]=new NormalDistribution(0, 1.0D).sample();
         	}
         }
		 Graph graph=mgdGraph;
		 
		DenseProjection head = new DenseProjection(graph, k, g, B,t,lambdas,gradientFx,gamma);
		System.out.println("\n\nResults>>> "+Sets.intersection(trueSet, head.bestS).size()+" "+head.bestS.size()+" "+trueSet.size());
        System.out.println("mu="+mu+" k="+k+" prec="+Sets.intersection(trueSet, head.bestS).size()*1.0D/head.bestS.size()+" rec="+Sets.intersection(trueSet, head.bestS).size()*1.0D/trueSet.size()+"\n\n");
        try {
			fileWriter.write("mu="+mu+" k="+k+" prec="+Sets.intersection(trueSet, head.bestS).size()*1.0D/head.bestS.size()+" rec="+Sets.intersection(trueSet, head.bestS).size()*1.0D/trueSet.size()+"\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			fileWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		ExecutorService pool = Executors.newFixedThreadPool(1);
		
		
		for(double mu:new double[]{100.0D}){//,3.0,2.0,1.0,0.5}){
			
		 for(int k:new int[]{50}){//,50,70,90,100,120,150}){
			 
         //Graph graph, int k, int g, double B, int t,double[] lambdas,double[] b,double gamma
			 pool.execute(new Thread() {
                 public void run() {
                	 SingleFileTest(mu,k); 
                 }
                 
             });
                 }
		 
		}
		pool.shutdown();
	}

}
