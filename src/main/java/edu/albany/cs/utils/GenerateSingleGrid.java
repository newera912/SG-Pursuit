package edu.albany.cs.utils;

import edu.albany.cs.base.Edge;
import org.apache.commons.lang3.ArrayUtils;

import java.util.ArrayList;


public class GenerateSingleGrid {
	
	public final int numOfNodes ;
	public final ArrayList<ArrayList<Integer>> adj ;
	public final double[][] adjMatrix ;
	public final int[][] adjInt ;
	public final ArrayList<Edge> edges ;
	
	public GenerateSingleGrid(int N){
		numOfNodes = N ;
		adj = this.generateGraph((int)Math.sqrt(N)) ;
		adjMatrix = this.generateGraph((int)Math.sqrt(N),1.0) ;
		adjInt = this.generateGraph((int)Math.sqrt(N), false) ;
		edges = new ArrayList<Edge>() ;
		int count = 0 ;
		double weight = 1.0 ;
		for(int i = 0 ; i < adj.size() ; i++){
			for( int j = 0 ; j < adj.get(i).size() ; j ++){
				edges.add(new Edge(i,adj.get(i).get(j),count++,weight)) ;
			}
		}
	}
	
	private double[][] generateGraph(int N,double weight){
		double[][] graph = new double[N*N][N*N] ;
		for(int k=0;k<graph.length;k++){
			for(int j=0;j<graph.length;j++){
				graph[k][j] = -1 ;
			}
		}
		int[][] Nodes = new int[N][N] ;
		int count = 0;
		for(int i=0;i<Nodes.length;i++)
			for(int j=0;j<Nodes.length;j++){
				Nodes[i][j] = count++ ;
			}
		count = 0 ;
		int i = 1 ;
		for(int k=1;k<=N*N;k++){
			count++ ;
			if(count>N){
				count = 1 ;
				i++ ;
			}
			int j = count;
			if ((i-1) < 1)
				graph[k-1][ Nodes[N-1][j-1] ]= weight ;
			else graph[k-1] [ Nodes[i-2][j-1] ] = weight ;
			if((i+1)>N)
				graph[k-1][ Nodes[0][j-1] ] = weight ;
			else graph[k-1][ Nodes[i][j-1] ] = weight ;
			if((j-1)<1)
				graph[k-1][ Nodes[i-1][N-1] ] = weight ;
			else graph[k-1][ Nodes[i-1][j-2] ] = weight ;
			if((j+1)>N)
				graph[k-1][ Nodes[i-1][0] ] = weight ;
			else graph[k-1][ Nodes[i-1][j] ] = weight ;
		}
		return graph ;
	}
	
	private ArrayList<ArrayList<Integer>> generateGraph(int N){
		double[][] graph = new double[N*N][N*N] ;
		for(int k=0;k<graph.length;k++){
			for(int j=0;j<graph.length;j++){
				graph[k][j] = -1 ;
			}
		}
		int[][] Nodes = new int[N][N] ;
		int count = 0;
		for(int i=0;i<Nodes.length;i++)
			for(int j=0;j<Nodes.length;j++){
				Nodes[i][j] = count++ ;
			}
		count = 0 ;
		int i = 1 ;
		for(int k=1;k<=N*N;k++){
			count++ ;
			if(count>N){
				count = 1 ;
				i++ ;
			}
			int j = count;
			if ((i-1) < 1)
				graph[k-1][ Nodes[N-1][j-1] ]= 1.0 ;
			else graph[k-1] [ Nodes[i-2][j-1] ] = 1.0 ;
			if((i+1)>N)
				graph[k-1][ Nodes[0][j-1] ] = 1.0 ;
			else graph[k-1][ Nodes[i][j-1] ] = 1.0 ;
			if((j-1)<1)
				graph[k-1][ Nodes[i-1][N-1] ] = 1.0 ;
			else graph[k-1][ Nodes[i-1][j-2] ] = 1.0 ;
			if((j+1)>N)
				graph[k-1][ Nodes[i-1][0] ] = 1.0 ;
			else graph[k-1][ Nodes[i-1][j] ] = 1.0 ;
		}
		
		ArrayList<ArrayList<Integer>> adj = new ArrayList<ArrayList<Integer>>() ;
		
		for(int j = 0 ; j < graph.length ; j ++){
			ArrayList<Integer> arr = new ArrayList<Integer>() ;
			for( int k = 0 ; k < graph[j].length ; k ++){
				if(graph[j][k] > 0){
					arr.add(k) ;	
				}
			}
			adj.add(arr) ;
		}
		return adj ;
	}
	
	private int[][] generateGraph(int N, boolean flag){
		double[][] graph = new double[N*N][N*N] ;
		for(int k=0;k<graph.length;k++){
			for(int j=0;j<graph.length;j++){
				graph[k][j] = -1 ;
			}
		}
		int[][] Nodes = new int[N][N] ;
		int count = 0;
		for(int i=0;i<Nodes.length;i++)
			for(int j=0;j<Nodes.length;j++){
				Nodes[i][j] = count++ ;
			}
		count = 0 ;
		int i = 1 ;
		for(int k=1;k<=N*N;k++){
			count++ ;
			if(count>N){
				count = 1 ;
				i++ ;
			}
			int j = count;
			if ((i-1) < 1)
				graph[k-1][ Nodes[N-1][j-1] ]= 1.0 ;
			else graph[k-1] [ Nodes[i-2][j-1] ] = 1.0 ;
			if((i+1)>N)
				graph[k-1][ Nodes[0][j-1] ] = 1.0 ;
			else graph[k-1][ Nodes[i][j-1] ] = 1.0 ;
			if((j-1)<1)
				graph[k-1][ Nodes[i-1][N-1] ] = 1.0 ;
			else graph[k-1][ Nodes[i-1][j-2] ] = 1.0 ;
			if((j+1)>N)
				graph[k-1][ Nodes[i-1][0] ] = 1.0 ;
			else graph[k-1][ Nodes[i-1][j] ] = 1.0 ;
		}
		int[][] adj = new int[N*N][] ;
		for(int j = 0 ; j < graph.length ; j ++){
			int[] arr = null ;
			for( int k = 0 ; k < graph[j].length ; k ++){
				if(graph[j][k] > 0){
					arr = ArrayUtils.add(arr, k) ;	
				}
			}
			adj[j] = arr ;
		}
		return adj ;
	}
	
	public static void main(String args[]){
	}
	
}
