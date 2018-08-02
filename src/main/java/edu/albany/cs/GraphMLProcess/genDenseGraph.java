package edu.albany.cs.GraphMLProcess;


import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Random;
import java.util.Set;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.distribution.NormalDistribution;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.Edge;

public class genDenseGraph {
	
	public static void genAPDMFromClusterData(String outfileName,String root,int r,int numNodes,int trueSubSize,int numFea,double p_in,double p_out, double sigmas1){
		
		
		String apdmFile=root+outfileName;	
	    
	    int numOfNodes=numNodes;
	    int numFestures=numFea;
	    double[][] data=new double[numFestures][numNodes];
	    
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
		int trueEdgeID=0;
		int clusterID=0;
		int randomIdx=random.nextInt((int)Math.floor(numNodes/trueSubSize)-1);
		int startIdx=randomIdx*trueSubSize;
		int endIdx=(randomIdx+1)*trueSubSize-1;
		System.out.println("True subnode ID "+startIdx+" ~ "+ endIdx);
		for(int i=0;i<numNodes;i++){ 
			for(int j=i+1;j<numNodes;j++){
				if((i>=startIdx && i<=endIdx)){							
					if((j>=startIdx && j<=endIdx)){  /** edges inside of the dense sub_graph*/
						double chance=new BinomialDistribution(1,p_in).sample();
						if(chance==1.0D){
							edges.add(new Edge(i, j, edgeID, 1.0D));
							truedges.add(new Edge(i, j, trueEdgeID, 1.0D));
							edgeID++;
							trueEdgeID++;
						}
						
					}else{   /** edges from dense graph to outside**/
						double chance=new BinomialDistribution(1,p_out).sample();
						if(chance==1.0D){
							edges.add(new Edge(i, j, edgeID, 1.0D));
							edgeID++;
						}
					}
					
				}else{ //edges outside of the dense subgraph
					double chance=new BinomialDistribution(1,p_out).sample();
					if(chance==1.0D){
						edges.add(new Edge(i, j, edgeID, 1.0D));
						edgeID++;
					}
				}
			}
		}
		ArrayList<ArrayList<Integer>> graphAdjList = new ArrayList<>();
		
		int[][] graphAdj = new int[numNodes][];
		
		for (int i = 0; i < numNodes; i++) {
			graphAdjList.add(new ArrayList<>());			
		}

		for (Edge edge : edges) {
			if (!graphAdjList.get(edge.i).contains(edge.j)) {
				graphAdjList.get(edge.i).add(edge.j);
				
			}
			if (!graphAdjList.get(edge.j).contains(edge.i)) {
				graphAdjList.get(edge.j).add(edge.i);
				
			}
		}
		ConnectedComponents cc = new ConnectedComponents(graphAdjList);
	    System.out.println("is connected: " + cc.checkConnectivity());
		
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
		System.out.println("true nodes: "+trueNodes.length+" " +ArrayUtils.toString(trueNodes));
		//Utils.stop();
		
		/**DEfine \mu_j ,where j \in S_y, draw from the uniform distribution [0,1]*/
		double[] mu=new double[trueFeat.length];
		for(int j=0;j<mu.length;j++){
			mu[j]=Math.random();
		}
		for(int i=0;i<data.length;i++){
			for(int j=0;j<data[0].length;j++){
				if(Arrays.asList(ArrayUtils.toObject(trueNodes)).contains(j) && Arrays.asList(ArrayUtils.toObject(trueFeat)).contains(i)){
					data[i][j]=new NormalDistribution(mu[ArrayUtils.indexOf(trueFeat, i)], sigmas1).sample();					
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
	
	public static void genDenseGAMerData() {
		int N = 500;
		int trueSubSize = 35;
		int numFea = 20;
		double[] ratios = new double[] { 0.25 };// , 0.10, 0.15, 0.20
		double p_in = 0.35;
		double p_out = 0.1;

		String root = "data/DenseGraph/Dense_APDM/GAMerSingle/";

		for (double sigmas1 : new double[] { 0.0316D }) {// , 0.1D, 0.5D, 0.6D,
															// 0.7D, 0.8D, 0.9D
			for (double r : ratios) {
				Double dd = numFea * r;
				int numTueFeat = dd.intValue();
				for (int i = 0; i < 5; i++) {
					String fileName = "APDM_Dense_subgraph_in_" + p_in
							+ "_out_" + p_out + "_FeasNum_" + numFea
							+ "_trueFeasNum_" + numTueFeat + "_sigmas1_"
							+ sigmas1 + "_case_" + i + ".txt";
					genAPDMFromClusterData(fileName, root, numTueFeat, N,
							trueSubSize, numFea, p_in, p_out, sigmas1);
				}
			}
		}
	}
	
	public static void genAPDMFromClusterData_VaryingNumOfAttributes(String outfileName,String root,int r,int numClusters,int clusterSize,int numFea,double p_in,double p_out, double sigmas1){
		String apdmFile=root+outfileName;	
	    int numOfNodes=numClusters * clusterSize;
	    int numFestures=numFea;
	    ArrayList<Edge> edges= new ArrayList<Edge>();
		String usedAlgorithm = "NULL";
		String dataSource = "MultiVarDataset";
		String[] attributeNames=new String[numFestures];
		for(int i=0;i<numFestures;i++){
			attributeNames[i]="att"+i;
		}
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
        }
		ArrayList<ArrayList<Integer>> graphAdjList = new ArrayList<>();
		for (int i = 0; i < numOfNodes; i++) {
			graphAdjList.add(new ArrayList<>());			
		}

		int edgeID=0;
        for(int i=0; i<numOfNodes; i++){
        	for(int j=i+1; j<numOfNodes; j++){
        		double chance = 0;
        		if(Math.floor(i/(clusterSize*1.0)) == Math.floor(j/(clusterSize*1.0))){
					chance=new BinomialDistribution(1,p_in).sample();
        		}else{
					chance=new BinomialDistribution(1,p_out).sample();
        		}
				if(chance==1.0D){
					edges.add(new Edge(i, j, edgeID, 1.0D));
					graphAdjList.get(i).add(j);
					graphAdjList.get(j).add(i);
					edgeID++;
				}
        	}
        }
		ConnectedComponents cc = new ConnectedComponents(graphAdjList);
	    System.out.println("is connected: " + cc.checkConnectivity());
	    
        // randomly pick one cluster as an anomalous cluster. cluster id starts from 0. 
	    ArrayList<Edge> truedges= new ArrayList<Edge>();
        int trueClusterId = (int)(Math.random() * numClusters); 
        System.out.println("ID: " + trueClusterId + "," + numClusters);
        int trueEdgeID = 0;
        int[] trueNodes = new int[clusterSize];
        for(int i=0; i<clusterSize; i++){
        	int nid = trueClusterId * clusterSize + i;
        	trueNodes[i] = nid;
        }

        for(int i=0; i<clusterSize; i++){
        	for(int j=i+1; j<clusterSize; j++){
        		if(graphAdjList.get(trueNodes[i]).indexOf(trueNodes[j]) != -1){
            		trueEdgeID += 1;
                	truedges.add(new Edge(trueNodes[i], trueNodes[j], trueEdgeID, 1.0D));
            		trueEdgeID += 1;
                	truedges.add(new Edge(trueNodes[j], trueNodes[i], trueEdgeID, 1.0D));
        		}
        	}
        }
    	
		/**DEfine \mu_j ,where j \in S_y, draw from the uniform distribution [0,1]*/
		double[] mu=new double[trueFeat.length];
		for(int j=0;j<mu.length;j++){
			mu[j]=Math.random();
		}
	    double[][] data=new double[numFestures][numOfNodes];
		for(int i=0;i<numFestures;i++){
			for(int j=0;j<numOfNodes;j++){
				if(Arrays.asList(ArrayUtils.toObject(trueNodes)).contains(j) && Arrays.asList(ArrayUtils.toObject(trueFeat)).contains(i)){
					data[i][j]=new NormalDistribution(mu[ArrayUtils.indexOf(trueFeat, i)], sigmas1).sample();					
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
	
	public static void genAPDMFromClusterData1(String outfileName,String root,int r,int numClusters,int clusterSize,int numFea,double p_in,double p_out, double sigmas1){
		String apdmFile=root+outfileName;	
	    int numOfNodes=numClusters * clusterSize;
	    int numFestures=numFea;
	    ArrayList<Edge> edges= new ArrayList<Edge>();
		String usedAlgorithm = "NULL";
		String dataSource = "MultiVarDataset";
		String[] attributeNames=new String[numFestures];
		for(int i=0;i<numFestures;i++){
			attributeNames[i]="att"+i;
		}
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
        }
		ArrayList<ArrayList<Integer>> graphAdjList = new ArrayList<>();
		for (int i = 0; i < numOfNodes; i++) {
			graphAdjList.add(new ArrayList<>());			
		}

		int edgeID=0;
        for(int i=0; i<numOfNodes; i++){
        	for(int j=i+1; j<numOfNodes; j++){
        		double chance = 0;
        		if(Math.floor(i/(clusterSize*1.0)) == Math.floor(j/(clusterSize*1.0))){
					chance=new BinomialDistribution(1,p_in).sample();
        		}else{
					chance=new BinomialDistribution(1,p_out).sample();
        		}
				if(chance==1.0D){
					edges.add(new Edge(i, j, edgeID, 1.0D));
					graphAdjList.get(i).add(j);
					graphAdjList.get(j).add(i);
					edgeID++;
				}
        	}
        }
		ConnectedComponents cc = new ConnectedComponents(graphAdjList);
	    System.out.println("is connected: " + cc.checkConnectivity());
	    
        // randomly pick one cluster as an anomalous cluster. cluster id starts from 0. 
	    ArrayList<Edge> truedges= new ArrayList<Edge>();
        int trueClusterId = (int)(Math.random() * numClusters); 
        System.out.println("ID: " + trueClusterId + "," + numClusters);
        int trueEdgeID = 0;
        int[] trueNodes = new int[clusterSize];
        for(int i=0; i<clusterSize; i++){
        	int nid = trueClusterId * clusterSize + i;
        	trueNodes[i] = nid;
        }

        for(int i=0; i<clusterSize; i++){
        	for(int j=i+1; j<clusterSize; j++){
        		if(graphAdjList.get(trueNodes[i]).indexOf(trueNodes[j]) != -1){
            		trueEdgeID += 1;
                	truedges.add(new Edge(trueNodes[i], trueNodes[j], trueEdgeID, 1.0D));
            		trueEdgeID += 1;
                	truedges.add(new Edge(trueNodes[j], trueNodes[i], trueEdgeID, 1.0D));
        		}
        	}
        }
    	
		/**DEfine \mu_j ,where j \in S_y, draw from the uniform distribution [0,1]*/
		double[] mu=new double[trueFeat.length];
		for(int j=0;j<mu.length;j++){
			mu[j]=Math.random();
		}
	    double[][] data=new double[numFestures][numOfNodes];
		for(int i=0;i<numFestures;i++){
			for(int j=0;j<numOfNodes;j++){
				if(Arrays.asList(ArrayUtils.toObject(trueNodes)).contains(j) && Arrays.asList(ArrayUtils.toObject(trueFeat)).contains(i)){
					data[i][j]=new NormalDistribution(mu[ArrayUtils.indexOf(trueFeat, i)], sigmas1).sample();					
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
	
	
	public static void genAPDMFromClusterData_VaryingClusterSizes(String outfileName,String root,int r,int numClusters,int clusterSize_lower, int clusterSize_upper,int numFea,double p_in,double p_out, double sigmas1){
		String apdmFile=root+outfileName;	
		int[] clusterSizes = new int[numClusters];
		int numOfNodes = 0;
		int[][] range = new int[numClusters][2];
		for(int i=0; i<numClusters; i++){
			int cs = (int)(Math.random() * (clusterSize_upper - clusterSize_lower)) + clusterSize_lower;
			clusterSizes[i] = cs;
			range[i][0] = numOfNodes;
			numOfNodes += cs;
			range[i][1] = numOfNodes-1;
		}
	    int numFestures=numFea;
	    ArrayList<Edge> edges= new ArrayList<Edge>();
		String usedAlgorithm = "NULL";
		String dataSource = "MultiVarDataset";
		String[] attributeNames=new String[numFestures];
		for(int i=0;i<numFestures;i++){
			attributeNames[i]="att"+i;
		}
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
        }
		ArrayList<ArrayList<Integer>> graphAdjList = new ArrayList<>();
		for (int i = 0; i < numOfNodes; i++) {
			graphAdjList.add(new ArrayList<>());			
		}

		int edgeID=0;
        for(int i=0; i<numOfNodes; i++){
        	for(int j=i+1; j<numOfNodes; j++){
        		double chance=-1;
        		for(int c=0; c<numClusters; c++){
	        		if(i > range[c][0] && i < range[c][1] && j > range[c][0] && j < range[c][1]){
	        			chance=new BinomialDistribution(1,p_in).sample();
	        			break;
	        		}
        		}
        		if(chance==-1){
        			chance=new BinomialDistribution(1,p_out).sample();
        		}
				if(chance==1.0D){
					edges.add(new Edge(i, j, edgeID, 1.0D));
					graphAdjList.get(i).add(j);
					graphAdjList.get(j).add(i);
					edgeID++;
				}
        	}
       	}

		ConnectedComponents cc = new ConnectedComponents(graphAdjList);
	    System.out.println("is connected: " + cc.checkConnectivity());
	    
        // randomly pick one cluster as an anomalous cluster. cluster id starts from 0. 
	    ArrayList<Edge> truedges= new ArrayList<Edge>();
        int trueClusterId = (int)(Math.random() * numClusters); 
        System.out.println("ID: " + trueClusterId + "," + numClusters);
        int trueEdgeID = 0;
        int trueClusterSize = range[trueClusterId][1] - range[trueClusterId][0] + 1;
        int[] trueNodes = new int[trueClusterSize];
        for(int i=0; i<trueClusterSize; i++){
        	trueNodes[i] = range[trueClusterId][0] + i;
        }
        for(int i=0; i<trueClusterSize; i++){
        	for(int j=i+1; j<trueClusterSize; j++){
        		if(graphAdjList.get(trueNodes[i]).indexOf(trueNodes[j]) != -1){
            		trueEdgeID += 1;
                	truedges.add(new Edge(trueNodes[i], trueNodes[j], trueEdgeID, 1.0D));
            		trueEdgeID += 1;
                	truedges.add(new Edge(trueNodes[j], trueNodes[i], trueEdgeID, 1.0D));
        		}
        	}
        }
    	
		/**DEfine \mu_j ,where j \in S_y, draw from the uniform distribution [0,1]*/
		double[] mu=new double[trueFeat.length];
		for(int j=0;j<mu.length;j++){
			mu[j]=Math.random();
		}
	    double[][] data=new double[numFestures][numOfNodes];
		for(int i=0;i<numFestures;i++){
			for(int j=0;j<numOfNodes;j++){
				if(Arrays.asList(ArrayUtils.toObject(trueNodes)).contains(j) && Arrays.asList(ArrayUtils.toObject(trueFeat)).contains(i)){
					data[i][j]=new NormalDistribution(mu[ArrayUtils.indexOf(trueFeat, i)], sigmas1).sample();					
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
	
	

	
	/*
	 * 1 anomalous cluster, number of normal clusters: {10, 15, 20, 25, 30, 40}, same size 100 nodes. 
	 * 5 anomalous attributes, number of normal attributes: {100}
	 */
	public static void genDenseData_VaryingNumOfClusters(int[] S, int thresholds, boolean validateFile){
		int numFea=100;
		double p_in =0.35;
		double p_out=0.10; 
//		int clusterSize = 50; 
		double sigmas1 = 0.0316D;
//		int numTrueFeat = 5;
        ExecutorService pool = Executors.newFixedThreadPool(thresholds);
		String root="data/DenseGraph/DenseSubgraph_APDM/";
		for(int clusterSize:new int[]{15, 20, 25, 50, 100}){
			for(int numTrueFeat:new int[]{5, 10}){
				for(int numClusters:S){
					for(int instance=0;instance<50;instance++){
						String fileName="VaryingNumOfClusters_APDM_Dense_subgraph_in_"+p_in+"_out_"+p_out+"_numClusters_"+numClusters+"_TrueSGSize_"+clusterSize+"_FeasNum_"+numFea+"_trueFeasNum_"+numTrueFeat+"_sigmas1_" + sigmas1 + "_case_"+instance+".txt";
						File f = new File(root + fileName);
						if(f.exists()){
							try{
								if(validateFile) {
									APDMInputFormat apdm=new APDMInputFormat(root + fileName);
									continue; 
								}
							} catch(Exception ex) {
								System.out.println("delting: " + root + fileName);
								f.delete();
					            ex.printStackTrace();
					        }   
						}
		                pool.execute(new Thread() {
		                    public void run() {
		                    	genAPDMFromClusterData_VaryingNumOfAttributes(fileName,root,numTrueFeat,numClusters,clusterSize,numFea,p_in,p_out,sigmas1);  
		                    }
		                });
					}
				}
			}
		}
        pool.shutdown();
	}
	
	
	/*
	 * 1 anomalous cluster, 9 normal clusters, variable cluster sizes: {(50, 100), (50, 150), (50, 200), (50, 300), (50, 400)}. 
	 * 5 anomalous attributes, number of normal attributes: 100.
	 */
	public static void genDenseData_VaryingClusterSizes(int[] S, int thresholds, boolean validateFile){
		int numFea=100;
		double p_in =0.35;
//		double p_out=0.10;
		int clusterSize = 100; 
		double sigmas1 = 0.0316D;
//		int numTrueFeat = 5;
		int numClusters = 10;
        ExecutorService pool = Executors.newFixedThreadPool(thresholds);
        int clusterSize_lower = 15;
		String root="data/DenseGraph/DenseSubgraph_APDM/";
		for(double p_out:new double[]{0.10, 0.15}){
			for(int numTrueFeat:new int[]{5, 10}){
				for(int clusterSize_upper:S){
					for(int instance=0;instance<50;instance++){
						String fileName="VaryingClusterSizes_APDM_Dense_subgraph_in_"+p_in+"_out_"+p_out+"_numClusters_"+numClusters+"_clusterSizeLower_"+clusterSize_lower+"_clusterSizeUpper_"+clusterSize_upper+"_FeasNum_"+numFea+"_trueFeasNum_"+numTrueFeat+"_sigmas1_" + sigmas1 + "_case_"+instance+".txt";
						File f = new File(root + fileName);
						if(f.exists()){
							try{
								if(validateFile) {
									APDMInputFormat apdm=new APDMInputFormat(root + fileName);
									continue; 
								}
							} catch(Exception ex) {
								System.out.println("delting: " + root + fileName);
								f.delete();
					            ex.printStackTrace();
					        }   
						}
						pool.execute(new Thread() {
		                    public void run() {
		                    	genAPDMFromClusterData_VaryingClusterSizes(fileName,root,numTrueFeat,numClusters,clusterSize_lower,clusterSize_upper,numFea,p_in,p_out,sigmas1);
		                    }
		                });
					}
				}
			}
		}
        pool.shutdown();
	}

	
	/*
	 * 1 anomalous cluster, 9 normal clusters, same size 100 nodes. 
	 * 5 anomalous attributes, number of normal attributes: {10, 20, 40, 80, 100, 500}
	 */
	public static void genDenseData_VaryingNumOfAttributes(int[]  S, int thresholds, boolean validateFile){
		double p_in =0.35;
		double p_out=0.10;
		int numClusters = 20; 
//		int clusterSize = 15; 
		double sigmas1 = 0.0316D;
//		int numTrueFeat = 5;
        ExecutorService pool = Executors.newFixedThreadPool(thresholds);
		String root="data/DenseGraph/DenseSubgraph_APDM/";
		for(int clusterSize:new int[]{15, 20, 25, 50, 100}){
			for(int numTrueFeat:new int[]{5, 10}){
				for(int numFea:S){
					for(int instance=0;instance<50;instance++){
						String fileName="VaryingNumOfAttributes_APDM_Dense_subgraph_in_"+p_in+"_out_"+p_out+"_numClusters_"+numClusters+"_TrueSGSize_"+clusterSize+"_FeasNum_"+numFea+"_trueFeasNum_"+numTrueFeat+"_sigmas1_" + sigmas1 + "_case_"+instance+".txt";
						File f = new File(root + fileName);
//						System.out.println(f.exists());
						if(f.exists()){
							try{
								if(validateFile){
//									APDMInputFormat apdm=new APDMInputFormat(root + fileName);
									continue;
								}
							} catch(Exception ex) {
								System.out.println("delting: " + root + fileName);
								f.delete();
					            ex.printStackTrace();
					        }   
						}
		                pool.execute(new Thread() {
		                    public void run() {
		                    	genAPDMFromClusterData_VaryingNumOfAttributes(fileName,root,numTrueFeat,numClusters,clusterSize,numFea,p_in,p_out,sigmas1);  
		                    }
		                });
					}
				}
			}
		}
        pool.shutdown();
	}

		
	public static void genDenseData(){
		int N=1000;
		int trueSubSize=100;
		int numFea=100;
		double[] ratios=new double[]{0.05};// , 0.10, 0.15, 0.20
		double p_in = 0.35;
		double p_out = 0.01;
		
		String root="data/DenseGraph/DenseSubgraph_APDM/";
		
		for(double sigmas1:new double[]{0.0316D}){//, 0.1D, 0.5D, 0.6D, 0.7D, 0.8D, 0.9D
			for(double r:ratios){
				Double dd = numFea * r;
				int numTueFeat = dd.intValue();
				for(int i=0;i<10;i++){
					String fileName="APDM_Dense_subgraph_in_"+p_in+"_out_"+p_out+"_FeasNum_"+numFea+"_trueFeasNum_"+numTueFeat+ "_sigmas1_" + sigmas1 + "_case_"+i+".txt";
					genAPDMFromClusterData(fileName,root,numTueFeat,N,trueSubSize,numFea,p_in,p_out,sigmas1);
				}
			}
		}
	}
	public static void main(String args[]) throws FileNotFoundException{
		//genDenseData();	
//		genDenseData_VaryingNumOfAttributes(new int[]{100}, 8, false);		//, 500
		genDenseData_VaryingNumOfAttributes(new int[]{20, 40, 80, 100, 500}, 8, true);		//, 
		
		genDenseData_VaryingClusterSizes(new int[]{100, 150, 200, 300, 400}, 8, true);
//		genDenseData_VaryingClusterSizes(new int[]{400}, 1, false);

		genDenseData_VaryingNumOfClusters(new int[]{10, 15, 20, 23, 24, 25}, 8, true);
//		genDenseData_VaryingNumOfClusters(new int[]{30}, 5, false);

//		genDenseGAMerData();
	}
}
