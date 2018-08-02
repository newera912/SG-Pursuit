package edu.albany.cs.graph;


import java.util.ArrayList;
import java.util.HashSet;

import edu.albany.cs.base.APDMInputFormat;

public class MultiDataGraph extends Graph {
	
	public  double[][] W = null;
	public  int numOfFeas;
	public  int numOfNodes;	
	public  ArrayList<Double> edgeCosts;
	public  HashSet<Integer> nodes;
	public  ArrayList<Integer[]> edges;
	public int[] trueFeatures;

	public MultiDataGraph(APDMInputFormat apdm) {
		super(apdm);				
		this.W= new double[apdm.data.numNodes][apdm.data.numFeas];		
		this.numOfFeas=apdm.data.numFeas;
		this.numOfNodes=apdm.data.numNodes;
		this.trueFeatures=apdm.data.trueFeas;
		this.edgeCosts = new ArrayList<>();
		this.edges=new ArrayList<>();
		for (int i = 0; i < apdm.data.edgeCosts.size(); i++) {
			edgeCosts.add(1.0D);
		}
		for(int[] edge:apdm.data.edges.keySet()){			
			Integer[] edgeInt = {edge[0],edge[1]};
			this.edges.add(edgeInt);
		}
		this.nodes = new HashSet<>();
		for (Integer[] edge : this.edges) {
			this.nodes.add(edge[0]);
			this.nodes.add(edge[1]);
		}
		for(int i=0;i<W.length;i++){
			for(int j=0;j<W[0].length;j++){
				W[i][j]=apdm.data.attributes[j][i];
			}
		}
		

	}
	public MultiDataGraph(double[][] W, int num_nodes, int num_fea,
			ArrayList<Double> edgeCost, ArrayList<Integer[]> edges) {
		super(num_nodes, num_fea, edges, edgeCost);

		// public Graph(int n, int p, ArrayList<Integer[]> edges,
		// ArrayList<Double> edgeCosts) {
		this.W = W;
		this.numOfFeas = num_fea;
		this.numOfNodes = num_nodes;
		this.trueFeatures = null;
		this.edgeCosts = edgeCost;
		this.edges = new ArrayList<>();

		for (Integer[] edge : edges) {
			Integer[] edgeInt = {edge[0], edge[1]};
			this.edges.add(edgeInt);
		}
		this.nodes = new HashSet<>();
		for (Integer[] edge : this.edges) {
			this.nodes.add(edge[0]);
			this.nodes.add(edge[1]);
		}

	}

	

}
