package edu.albany.cs.GraphMLProcess;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.ConnectedComponents;
import edu.albany.cs.base.Edge;
import edu.albany.cs.base.RandomWalk;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.distribution.NormalDistribution;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class GenerateGridData {
	
	
	
	public GenerateGridData(String[] args) {

	}

	
		
	public void generateMultiGridData(int numTrueNodes, int numOfNodes,int numFestures,int[] trueFeat,double mu, double noiseLevel,
			String outPutFileName, boolean flag) throws IOException {

		// usedAlgorithm, dataSource, edges,tredges, PValue, fileName
		String usedAlgorithm = "NULL";
		String dataSource = "MultiVarDataset";
		String[] attributeNames=new String[]{"W1..WP"};
		GenerateSingleGrid g = new GenerateSingleGrid(numOfNodes);
		ArrayList<Edge> treEdges = new RandomWalk(g.adj, numTrueNodes,1000).subGraph;
		String trueFea="";
		double[][] PValue = new double[numFestures][g.numOfNodes];
		int[] trueNodes = null;
		//Arrays.fill(PValue, 1.0);
		for (Edge e : treEdges) {
			
			if (!ArrayUtils.contains(trueNodes, e.i)) {
				trueNodes = ArrayUtils.add(trueNodes, e.i);
			}
			if (!ArrayUtils.contains(trueNodes, e.j)) {
				trueNodes = ArrayUtils.add(trueNodes, e.j);
			}
		}
		System.out.println(ArrayUtils.toString(trueNodes)+" "+PValue.length+" "+PValue[0].length);
		System.out.println(ArrayUtils.toString(trueFeat));
		trueFea = ArrayUtils.toString(trueFeat).replace("{", "").replace("}", "").replace(",", " ");
		for(int i=0;i<PValue.length;i++){//10
			for(int j=0;j<PValue[0].length;j++){//100
				if(Arrays.asList(ArrayUtils.toObject(trueNodes)).contains(j) && Arrays.asList(ArrayUtils.toObject(trueFeat)).contains(i)){
					PValue[i][j]=new NormalDistribution(mu, 1.0D).sample();					
				}else{
					PValue[i][j]=new NormalDistribution(0, 1.0D).sample();
				}
			}
		}

		APDMInputFormat.generateAPDMFile(usedAlgorithm, dataSource, g.edges, treEdges,outPutFileName,attributeNames,numFestures,trueFea,PValue);
	}

	public void generateGridData(int numTrueNodes, int numOfNodes,String outPutFileName) throws IOException {

		// usedAlgorithm, dataSource, edges,tredges, PValue, fileName
		String usedAlgorithm = "NULL";
		String dataSource = "GridDataset";
		String[] attributeNames=new String[]{"W1..WP"};
		GenerateSingleGrid g = new GenerateSingleGrid(numOfNodes);
		ArrayList<Edge> treEdges = new RandomWalk(g.adj, numTrueNodes,1000).subGraph;
		String trueFea="";
		double[] PValue = new double[g.numOfNodes];
		double[] avgValue = new double[g.numOfNodes];
		int[] trueNodes = null;
		//Arrays.fill(PValue, 1.0);
		for (Edge e : treEdges) {
			
			if (!ArrayUtils.contains(trueNodes, e.i)) {
				trueNodes = ArrayUtils.add(trueNodes, e.i);
			}
			if (!ArrayUtils.contains(trueNodes, e.j)) {
				trueNodes = ArrayUtils.add(trueNodes, e.j);
			}
		}
		System.out.println(ArrayUtils.toString(trueNodes)+" "+PValue.length+" "+PValue.length);
		
		
		for(int i=0;i<PValue.length;i++){//100			
				if(Arrays.asList(ArrayUtils.toObject(trueNodes)).contains(i) ){
					PValue[i]=1.0D;					
				}else{
					PValue[i]=0.0D;
				}
			
		}
		//String usedAlgorithm, String dataSource, ArrayList<Edge> edges, double[] PValue,
		//double[] counts, double[] averValue, ArrayList<Edge> trueSubGraphEdges, String fileName
		APDMInputFormat.generateAPDMFile(usedAlgorithm, dataSource, g.edges,PValue,null,null, treEdges,outPutFileName);
	}
	
	public void genGridData() throws IOException {
		//int graphSize = 100;
		//int numOfTrueNodes = 10;
		ExecutorService pool = Executors.newFixedThreadPool(40);
		double noiseLevel = 0.0D;
		//int[] trueFeat=new Random().ints(0, numOfFeatures).distinct().limit(numOfTrueFeas).toArray();
		//double mu=0.5D;
		for(int graphSize:new int[]{2500}){//100,400,900,1600,2500}){			
				for(int truSize:new int[]{125,250}){	
					for(int i=0;i<100;i++){
						String root = "data/SimulationData/SimuGridTailTest/Grid-2500/truesub_"+truSize;
						if (!new File(root).exists()) {
							new File(root).mkdirs();
						}
						String fileName="/APDM-GridData-" + graphSize
								+ "trueSub_" +truSize + "_"+i+".txt";
						
						String outputFileName = root+fileName;
						//int numTrueNodes, int numOfNodes,int numFestures, double noiseLevel,String outPutFileName, boolean flag
						 pool.execute(new Thread() {
			                 public void run() {
						
			                	 try {
									generateGridData(truSize, graphSize, outputFileName);
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
			                 }
			                 
			             });
						
						//testTrueSubGraph(outputFileName);
					}
				}
			
		}
		pool.shutdown();
	}
	public void genMultiGridData() throws IOException {
		//int graphSize = 100;
		int numOfTrueNodes = 10;
		int numOfFeatures=20;
		int numOfTrueFeas=8;
		double noiseLevel = 0.0D;
		int[] trueFeat=new Random().ints(0, numOfFeatures).distinct().limit(numOfTrueFeas).toArray();
		//double mu=0.5D;
		for(int graphSize:new int[]{400}){//100,400,900,1600,2500}){
			for(double mu:new double[]{0.5,1.0,3.0,5.0,10.0}){
				for(int truSize:new int[]{10,20,40,80}){	
					for(int i=0;i<100;i++){
						String root = "data/SimulationData/MultiGrid-400/truesub_"+truSize+"/multiGridData_mu"+mu;
						if (!new File(root).exists()) {
							new File(root).mkdirs();
						}
						String fileName="/APDM-GridData-" + graphSize
								+ "-precen-0.1-noise_" + noiseLevel + "-"+i+".txt";
						
						String outputFileName = root+fileName;
						//int numTrueNodes, int numOfNodes,int numFestures, double noiseLevel,String outPutFileName, boolean flag
						generateMultiGridData(numOfTrueNodes, graphSize,numOfFeatures,trueFeat,mu, noiseLevel, outputFileName, false);
						
						testTrueSubGraph(outputFileName);
					}
				}
			}
		}
	}
	
	public void testTrueSubGraph(String fileName) {
		APDMInputFormat apdm = new APDMInputFormat(fileName);
		System.out.println(apdm.data.trueSubGraphNodes.length);
		ConnectedComponents cc = apdm.data.cc;
		System.out.println(cc.computeCCSubGraph(apdm.data.trueSubGraphNodes));
	}

	public static void main(String args[]) throws IOException {
		new GenerateGridData(args).genGridData();
	}
}
