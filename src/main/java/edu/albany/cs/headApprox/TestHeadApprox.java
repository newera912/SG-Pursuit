package edu.albany.cs.headApprox;

import edu.albany.cs.base.APDMInputFormat;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.random.RandomDataGenerator;

import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class TestHeadApprox {

	public static void testSingleFile(String fileName, String output,
			int trueSubSize, double edgeWeight, double noiseLevel) {
		FileWriter outFile = null;

		try {
			outFile = new FileWriter(output, true);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		APDMInputFormat apdm = new APDMInputFormat(fileName);
		int g = 1;
		ArrayList<Double> edges0_9 = new ArrayList<Double>();
		for (int q = 0; q < apdm.data.edgeCosts.size(); q++) {
			edges0_9.add(edgeWeight);
		}
		
		
		
		/**Abnormal noise ---------------------------------------------*/		
		int nosieNodeNum = (int) Math.ceil(noiseLevel * trueSubSize);
		int[] alreadyDone = null;
		for (int k = 0; k < nosieNodeNum; k++) {

			while (true) {
				int valueToFind = new Random().nextInt(trueSubSize);
				if (!ArrayUtils.contains(alreadyDone, valueToFind)) {
					apdm.data.PValue[apdm.data.trueSubGraphNodes[valueToFind]] = 0.0;
					alreadyDone = ArrayUtils.add(alreadyDone, valueToFind);
					break;
				}
			}
		}
		
		int[] trueNodes = null;
		for (int node : apdm.data.trueSubGraphNodes) {
			trueNodes = ArrayUtils.add(trueNodes, node);
		}
		
		/** Normal noise */
		int normalNoise= ((int) Math.ceil(noiseLevel * (apdm.data.numNodes-trueSubSize)));
		int[] normalNodes = ArrayUtils.removeElements(new RandomDataGenerator().nextPermutation(apdm.data.PValue.length, apdm.data.PValue.length), trueNodes);
		alreadyDone = null;
        for (int j = 0; j < normalNoise; j++) {
            while (true) {
                int valueToFind = new Random().nextInt(normalNodes.length);
                if (!ArrayUtils.contains(alreadyDone, valueToFind)) {                	
                	apdm.data.PValue[normalNodes[valueToFind]] = 1.0;                	
                    alreadyDone = ArrayUtils.add(alreadyDone, valueToFind);
                    break;
                }
            }
        }
         System.out.println(" Normal noise number :" + alreadyDone.length);
        //System.out.println("\nNum of Nosie abnodes     :" + nosieNodeNum);
		//System.out.println("Num of Nosie normal-nodes:" + normalNoise);
		int KK=(int) Math.ceil(0.5 * trueSubSize);
         
		double bestPre=0.0D;
		double bestRec=0.0D;
		double bestFscore=0.0D;
		double bestK=0;
		for (int k = KK; k < KK+1; k+=50) {
			k=KK;
			double B = k - g + 0.0D;
			HeadApprox tail = new HeadApprox(apdm.data.intEdges, edges0_9,
					apdm.data.PValue, k, g, B, false);
			double intersect = intersection(tail.bestForest.nodesInF, trueNodes).length * 1.0D;
			double pre = intersect * 1.0D / tail.bestForest.nodesInF.size();
			double recall = intersect * 1.0D / trueNodes.length;
			
			String result = "k=" + k + " #nodes: "
					+ tail.bestForest.nodesInF.size() + " pre="
					+ Math.round(pre * 100D) / 100.0D + " rec="
					+ Math.round(recall * 100D) / 100.0D;

			System.out.println(result);
		if(f1Score(pre, recall)>bestFscore){
			bestFscore=f1Score(pre, recall);
			bestPre=pre;
			bestRec=recall;
			bestK=k;
		}

		}
		try {
			outFile.write(bestK  + " " + Math.round(bestPre * 100D) / 100.0D + " "+ Math.round(bestRec * 100D) / 100.0D + " "+ Math.round(bestFscore * 100D) / 100.0D + "\n");
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		try {
			outFile.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String outRoot = "outputs/TailGridTest/";
		ExecutorService pool = Executors.newFixedThreadPool(1);
		for(double edgeWeight:new double[]{1.0} ){
			for(double noise:new double[]{0.05,0.1,0.15,0.2,0.25}){
				for (int graphSize : new int[] { 2500 }) {// 100,400,900,1600,2500}){
					for (int truSize : new int[] { 125,250 }) {
						String root = "data/SimulationData/SimuGridTailTest/Grid-2500/truesub_"	+ truSize;
						System.out.println("\nNum of Nosie abnodes     :" + truSize*noise);
						System.out.println("Num of Nosie normal-nodes:" + (graphSize-truSize)*noise);
						for (int i = 0; i < 100; i++) {
		
							String fileName = root + "/APDM-GridData-" + graphSize
									+ "trueSub_" + truSize + "_" + i + ".txt";
							String outFile = "H50%Nosie_EdgeWeight_"+edgeWeight+"_Noise_"+String.valueOf(noise)+"_trueSize_" + truSize + "_"+graphSize+".txt";
							String outputFileName = outRoot + outFile;
							// int numTrueNodes, int numOfNodes,int numFestures, double
							// noiseLevel,String outPutFileName, boolean flag
							 pool.execute(new Thread() {
				                 public void run() {
				                	 testSingleFile(fileName, outputFileName, truSize,edgeWeight,noise);
				                 }
				                 
				             });
						}
					}
		
				}
			}
		}
		pool.shutdown();
	}
	
	public static double f1Score(double pre,double rec){
		if(rec==0.0D || pre==0.0D){
			return 0.001;
		}else{
			return 2.0D*pre*rec/(pre+rec);
		}
	}
	public static int[] intersection(ArrayList<Integer> nums1, int[] nums2) {
		HashSet<Integer> set1 = new HashSet<Integer>();
		for (int i : nums1) {
			set1.add(i);
		}

		HashSet<Integer> set2 = new HashSet<Integer>();
		for (int i : nums2) {
			if (set1.contains(i)) {
				set2.add(i);
			}
		}

		int[] result = new int[set2.size()];
		int i = 0;
		for (int n : set2) {
			result[i++] = n;
		}

		return result;
	}
}
