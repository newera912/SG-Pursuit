package edu.albany.cs.DenseGraphMP;


import com.google.common.collect.Sets;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.graph.MultiDataGraph;
import edu.albany.cs.scoreFuncs.PCAScore;
import edu.albany.cs.scoreFuncs.Stat;

import org.apache.commons.lang3.ArrayUtils;

import java.io.FileWriter;
import java.io.IOException;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class TestS2GraphMPDense {
	
	public int count=0;	
	
	public static void TestingS2GraphMPDense(String fileName,String resultFileName,double lambda){
		FileWriter fileWriter=null;
		try {
			fileWriter=new FileWriter(resultFileName, true);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		};
		System.out.println(fileName);
		APDMInputFormat apdm=new APDMInputFormat(fileName);
		
		MultiDataGraph mgdGraph=new MultiDataGraph(apdm);
//		for(int e=0;e<10;e++){
//			System.out.println(ArrayUtils.toString(mgdGraph.edges.get(e)));
//		}
		System.out.println("True Nodes: "+ArrayUtils.toString(mgdGraph.trueFeatures));
		System.out.println("True Nodes: "+ArrayUtils.toString(mgdGraph.trueSubGraphNodes));
		
		PCAScore func=new PCAScore(mgdGraph.W,mgdGraph.AdjMatrix);				
		S2GraphMPDense bestS2GraphMPDense=null;
		double bestFunValue=Double.POSITIVE_INFINITY;
		int bestS=0;
		int bestK=0;
		int g=1;
		int t=3;
		for(double ratio=0.05;ratio<0.06;ratio+=0.02){
			int k=(int)Math.round(ratio*mgdGraph.numOfNodes);
			for(double ratioS=0.5;ratioS<0.6;ratioS+=0.1){
				int s=(int)Math.round(ratioS*mgdGraph.numOfFeas);
				double B = k - g + 0.0D;
				//Graph graph, int k, int s, int g, double B, int t, Function func,double gamma,double[] lambdas
				S2GraphMPDense s2GraphMPDense=new S2GraphMPDense(mgdGraph,k,s, g, B,t, func,lambda);
				if(bestFunValue>s2GraphMPDense.funcValue){
					bestFunValue=s2GraphMPDense.funcValue;
					bestS2GraphMPDense=s2GraphMPDense;
					bestK=k;
					bestS=s;
				}
			}
		}

		Set<Integer> intersect = Sets.intersection(bestS2GraphMPDense.psiX, mgdGraph.trueNodesSet);
		double precision = (intersect.size() * 1.0D / bestS2GraphMPDense.psiX.size() * 1.0D);
		double recall = (intersect.size() * 1.0D / mgdGraph.trueNodesSet.size() * 1.0D);		
		Set<Integer> yTrueArrayList= Stat.Array2Set(mgdGraph.trueFeatures);
		Set<Integer> intersect2 = Sets.intersection(bestS2GraphMPDense.psiY, yTrueArrayList);
		double precision2 = (intersect2.size() * 1.0D / bestS2GraphMPDense.psiY.size() * 1.0D);
		double recall2 = (intersect2.size() * 1.0D / yTrueArrayList.size() * 1.0D);
		//System.out.println("File"+count+" X:pre=" + precision+"rec=" + recall+" Y: pre=" + precision2+"rec=" + recall2);
		System.out.println(fileName.split("_")[11]+" k="+bestK+" s="+bestS+" "+precision+" " + recall+" " + precision2+" " + recall2 +" FunVal="+bestFunValue);
		try {
			fileWriter.write(precision+" " + recall+" " + precision2+" " + recall2+" "+bestK+" "+bestS +"\n");
			fileWriter.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
	}
	
	
	public static void TestMultiFile(){
		//int[] trueFeatures= new int[]{0,1,2,3,4};
		long startTimeAll = System.nanoTime();
		double lambda=0.0D;
		
		double in_p=0.35;
		double out_p=0.05;
		
	
        ExecutorService pool = Executors.newFixedThreadPool(1);
        
		///int g = 1;
        String root="data/DenseGraph/DenseSubgraph_APDM/";
    	for (int r : new int[]{5}){           
            for (int cases=0;cases<1;cases++) {
            	String fileName="APDM_Dense_subgraph_in_"+in_p+"_out_"+out_p+"_trueFeasNum_"+r+"_case_"+cases+".txt";				
            	String outFile="outputs/S2GraphMPDense/Result_"+in_p+"_out_"+out_p+"_trueFeasNum_"+r+"_case_"+cases+".txt";				
                
            	System.out.println(fileName);                    
                pool.execute(new Thread() {
                    public void run() {
                    	TestingS2GraphMPDense(root+fileName,outFile,lambda);                        	
                    }
                    
                });
                
            }
        }
        
        pool.shutdown();
	
		System.out.println("Total running time: " + (System.nanoTime() - startTimeAll) / 1e9);
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		TestMultiFile();
		//TestMultiFile();
	}

}
