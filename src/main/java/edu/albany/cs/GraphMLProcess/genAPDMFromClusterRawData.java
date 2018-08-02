package edu.albany.cs.GraphMLProcess;


import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.Edge;



public class genAPDMFromClusterRawData {
	
	public static void genAPDMFromClusterData(String fileName,String root,String subDir,int r,int numNodes,int numFea,int clusterNumber){
		//bio   arxiv99-03 dblp_top  patent_91-95
		
		String apdmFile = root + "Dense_APDM/cluster" + numNodes + "/APDM_r_"
				+ r + "_" + fileName;
		apdmFile = root + "Dense_APDM/GAMer/APDM_r_" + r + "_" + fileName;
		// String dataFile=root+fileName;
		
	    //i9.graph.gamer.base.Parameter.numberOfAtts = myGraph.getNumberOfAtts();
	    
	    int numOfNodes=numNodes;
	    int numFestures=numFea;
	    double[][] data=new double[numFestures][numOfNodes];
	    HashMap<Integer, Integer> originId2NodeID= new HashMap<Integer, Integer>();
	    ArrayList<Edge> edges= new ArrayList<Edge>();
	    ArrayList<Edge> truedges= new ArrayList<Edge>();
	    
	    Set<String> edgeSet = new HashSet<>();
		String usedAlgorithm = "NULL";
		
		String dataSource = "MultiVarDataset";
		String[] attributeNames=new String[numFestures];
		for(int i=0;i<numFestures;i++){
			attributeNames[i]="att"+i;
		}
		
		
		/** random select the ground truth features and cluster */
		Random random = new Random();
		String truthFea="";
		int[] trueFeat=new int[r];
		ArrayList<Integer> list = new ArrayList<Integer>();
        for (int i=1; i<numFestures; i++) {
            list.add(new Integer(i));
        }
        Collections.shuffle(list);
        for (int i=0; i<r; i++) {
        	truthFea+=list.get(i)+" ";
        	trueFeat[i]=list.get(i);
            //System.out.println(list.get(i));
        }
		
		int edgeID=0;
		int clusterID=0;
		int randomCluster=random.nextInt(clusterNumber);
		ArrayList<Integer> abnormalNodes=new ArrayList<Integer>();
		
		
		/** Read raw data file, get edge and truth nodes **/
		try {
			for (String eachLine : Files.readAllLines(Paths.get(root+subDir,fileName))) {
				
				/**Read Edge info */
				if (eachLine.split(" ").length == 2) {
					int end0 = Integer.parseInt(eachLine.split(" ")[0]);
					int end1 = Integer.parseInt(eachLine.split(" ")[1]);
					edges.add(new Edge(end0, end1, edgeID, 1.0D));
					edgeID++;					
					if(end0<end1){
						edgeSet.add(end0+"_"+end1);	
					}else{
						edgeSet.add(end1+"_"+end0);
					}
					}
				/**Read Cluster nodes */
				if (eachLine.split(" ").length > 2) {
					if(clusterID==randomCluster){
						System.out.println("Random Cluster No. "+clusterNumber+"/"+randomCluster);
						for(String node:eachLine.split(" ")){
							abnormalNodes.add(Integer.parseInt(node));
						}
					}
					clusterID++;
				}
				if (eachLine.length()<1) {
					continue;
				}
				
			}
		} catch (NumberFormatException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.out.println(abnormalNodes);
		int count=0;
		for(String edg:edgeSet){
			String[] e=edg.split("_");	
			
			if(abnormalNodes.contains(Integer.parseInt(e[0])) && abnormalNodes.contains(Integer.parseInt(e[1]))){
				//System.out.println(e[0]+" "+e[1]);
				truedges.add(new Edge(Integer.parseInt(e[0]), Integer.parseInt(e[1]), count++, 1.0D));
				
			}
		}
		
		int[] trueNodes = null;		
		for (Edge e : truedges) {			
			if (!ArrayUtils.contains(trueNodes, e.i)) {
				trueNodes = ArrayUtils.add(trueNodes, e.i);
			}
			if (!ArrayUtils.contains(trueNodes, e.j)) {
				trueNodes = ArrayUtils.add(trueNodes, e.j);
			}
		}
		//System.out.println(ArrayUtils.toString(trueNodes));
		
		/**DEfine \mu_j ,where j \in S_y, draw from the uniform distribution [0,1]*/
		double[] mu=new double[trueFeat.length];
		for(int j=0;j<mu.length;j++){
			mu[j]=Math.random();
		}
		for(int i=0;i<data.length;i++){
			for(int j=0;j<data[0].length;j++){
				if(Arrays.asList(ArrayUtils.toObject(trueNodes)).contains(j) && Arrays.asList(ArrayUtils.toObject(trueFeat)).contains(i)){
					data[i][j]=new NormalDistribution(mu[ArrayUtils.indexOf(trueFeat, i)], 0.001D).sample();					
				}else{
					data[i][j]=new NormalDistribution(0, 1.0D).sample();
				}
			}
		}
		//GenerateSingleGrid g = new GenerateSingleGrid(numOfNodes);
		
	    System.out.println("Generatign APDM file..... edges size="+edges.size()+" True Edges "+truedges.size());
	    APDMInputFormat.generateAPDMFile(usedAlgorithm, dataSource, edges, truedges,apdmFile,attributeNames,numFestures,truthFea,data);
	    System.out.println("Done.....");
	}
	
	public static void genAPDMFromGAMerClusterData(String fileName,
			String root, String outroot, int r, int numNodes, int numFea,
			int clusterNumber) {
		// bio arxiv99-03 dblp_top patent_91-95

		String apdmFile = outroot + "/APDM_" + fileName;
		// apdmFile = root + "Dense_APDM/GAMer/APDM_r_" + r + "_" + fileName;
		// String dataFile=root+fileName;

		// i9.graph.gamer.base.Parameter.numberOfAtts =
		// myGraph.getNumberOfAtts();

		int numOfNodes = numNodes;
		int numFestures = numFea;
		double[][] data = new double[numFestures][numOfNodes];
		HashMap<Integer, Integer> originId2NodeID = new HashMap<Integer, Integer>();
		ArrayList<Edge> edges = new ArrayList<Edge>();
		ArrayList<Edge> truedges = new ArrayList<Edge>();

		Set<String> edgeSet = new HashSet<>();
		String usedAlgorithm = "NULL";

		String dataSource = "MultiVarDataset";
		String[] attributeNames = new String[numFestures];
		for (int i = 0; i < numFestures; i++) {
			attributeNames[i] = "att" + i;
		}

		/** random select the ground truth features and cluster */
		Random random = new Random();
		String truthFea = "";
		int[] trueFeat = new int[r];
		ArrayList<Integer> list = new ArrayList<Integer>();
		for (int i = 1; i < numFestures; i++) {
			list.add(new Integer(i));
		}
		Collections.shuffle(list);
		for (int i = 0; i < r; i++) {
			truthFea += list.get(i) + " ";
			trueFeat[i] = list.get(i);
			// System.out.println(list.get(i));
		}

		int edgeID = 0;
		int clusterID = 0;
		int randomCluster = random.nextInt(clusterNumber);
		ArrayList<Integer> abnormalNodes = new ArrayList<Integer>();

		/** Read raw data file, get edge and truth nodes **/
		try {
			for (String eachLine : Files
					.readAllLines(Paths.get(root, fileName))) {

				/** Read Edge info */
				if (eachLine.split(" ").length == 2) {
					int end0 = Integer.parseInt(eachLine.split(" ")[0]);
					int end1 = Integer.parseInt(eachLine.split(" ")[1]);
					edges.add(new Edge(end0, end1, edgeID, 1.0D));
					edgeID++;
					if (end0 < end1) {
						edgeSet.add(end0 + "_" + end1);
					} else {
						edgeSet.add(end1 + "_" + end0);
					}
				}
				/** Read Cluster nodes */
				if (eachLine.split(" ").length > 2) {
					if (clusterID == randomCluster) {
						System.out.println("Random Cluster No. "
								+ clusterNumber + "/" + randomCluster);
						for (String node : eachLine.split(" ")) {
							abnormalNodes.add(Integer.parseInt(node));
						}
					}
					clusterID++;
				}
				if (eachLine.length() < 1) {
					continue;
				}

			}
		} catch (NumberFormatException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
		System.out.println(abnormalNodes);
		int count = 0;
		for (String edg : edgeSet) {
			String[] e = edg.split("_");

			if (abnormalNodes.contains(Integer.parseInt(e[0]))
					&& abnormalNodes.contains(Integer.parseInt(e[1]))) {
				// System.out.println(e[0]+" "+e[1]);
				truedges.add(new Edge(Integer.parseInt(e[0]), Integer
						.parseInt(e[1]), count++, 1.0D));

			}
		}

		int[] trueNodes = null;
		for (Edge e : truedges) {
			if (!ArrayUtils.contains(trueNodes, e.i)) {
				trueNodes = ArrayUtils.add(trueNodes, e.i);
			}
			if (!ArrayUtils.contains(trueNodes, e.j)) {
				trueNodes = ArrayUtils.add(trueNodes, e.j);
			}
		}
		// System.out.println(ArrayUtils.toString(trueNodes));

		/**
		 * DEfine \mu_j ,where j \in S_y, draw from the uniform distribution
		 * [0,1]
		 */
		double[] mu = new double[trueFeat.length];
		for (int j = 0; j < mu.length; j++) {
			mu[j] = Math.random();
		}
		for (int i = 0; i < data.length; i++) {
			for (int j = 0; j < data[0].length; j++) {
				if (Arrays.asList(ArrayUtils.toObject(trueNodes)).contains(j)
						&& Arrays.asList(ArrayUtils.toObject(trueFeat))
								.contains(i)) {
					data[i][j] = new NormalDistribution(mu[ArrayUtils.indexOf(
							trueFeat, i)], 0.001D).sample();
				} else {
					data[i][j] = new NormalDistribution(0, 1.0D).sample();
				}
			}
		}
		// GenerateSingleGrid g = new GenerateSingleGrid(numOfNodes);

		System.out.println("Generatign APDM file..... edges size="
				+ edges.size() + " True Edges " + truedges.size());
		APDMInputFormat
				.generateAPDMFile(usedAlgorithm, dataSource, edges, truedges,
						apdmFile, attributeNames, numFestures, truthFea, data);
		System.out.println("Done.....");
	}
	public static void genDenseData(){
		// int N = 2000;
		int[] m = new int[] { 20 };
		int numFea = 20;
		int[] r = new int[] { 5 };
		double p_in = 0.6;
		double p_out=0.01;
		// , 2500,3000, 3500, 4000, 4500, 5000, 5500, 6000
		for (int N : new int[] { 300, 500, 1000, 1200, 1500, 1700, 2000, 2500,
				3000, 3500, 4000, 4500, 5000, 5500, 6000 }) {
		String root="data/DenseGraph/";
		String subDir = "raw/rawData" + N + "/";
		
			for(int mm:m){
				for(int rr:r){
					for (int i = 0; i < 1; i++) {
						String fileName = "Cluster_" + mm + "_in_" + p_in
								+ "_out_" + p_out + "_15_case_" + i + ".txt";
						genAPDMFromClusterData(fileName,root,subDir,rr,N,numFea,mm);
					}
				}
			}
		}
		
	}

	public static void genGAMerDenseData() {
		// int N = 2000;
		String outroot = "data/DenseGraph/Dense_APDM/GAMer/";
		String root = "data/DenseGraph/raw/rawGAMer/";
		int numFea = 20;
		// int cluster=10;
		for (File rawFile : new File(root).listFiles()) {
			System.out.println(rawFile.getName());
			String fileName = rawFile.getName();
			int N = Integer.valueOf(rawFile.getName().split("_")[0]);
			int mm = Integer.valueOf(rawFile.getName().split("_")[2]);
			genAPDMFromGAMerClusterData(fileName, root, outroot, 5, N, numFea,
					mm);

		}

	}
	public static void main(String args[]) throws FileNotFoundException{
		genGAMerDenseData();
		
	}
}
