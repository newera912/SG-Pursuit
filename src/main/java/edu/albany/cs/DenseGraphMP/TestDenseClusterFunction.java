package edu.albany.cs.DenseGraphMP;

import java.util.ArrayList;

import edu.albany.cs.base.APDMInputFormat;

public class TestDenseClusterFunction {
	public static void experiments_RealdataSets() {

		double lambda = 5.00D;
		String root = "data/DenseGraph/Realdata/apdm/";
		String outroot = "outputs/new_testing.txt";
		// daasets :"dblp", "genes", "imdb", "dfb", "arxivSmall"
		for (String fileName : new String[]{"dblp"}) {
			fileName += "_norm_APDM.txt";
			long startTimeAll = System.nanoTime();
			System.out.println(fileName);
			APDMInputFormat apdm = new APDMInputFormat(root + fileName);
			double[][] data = new double[apdm.data.numNodes][apdm.data.numFeas+1];

			for (int i = 0; i < data.length; i++) {
				for (int j = 0; j < data[0].length; j++) {
					if (j == 0) {
						data[i][j] = i;
					}else {
					data[i][j] = apdm.data.attributes[j-1][i];
				}
				}
			}

			ArrayList<Integer[]> edgeList = new ArrayList<>();
			for (int[] edge : apdm.data.edges.keySet()) {
				Integer[] edgeInt = {edge[0], edge[1]};
				edgeList.add(edgeInt);
			}
			// double[][] data, ArrayList<Integer[]> edgeList, int s, int k,
			// double lambda, int numCluster, String resultFile
			DenseClusterDetection.SGPursuitDenseClusterDetection(data,
					edgeList, 5, 10, lambda, 10, outroot);

			System.out.println("Total running time: "
					+ (System.nanoTime() - startTimeAll) / 1e9);
		}
	}
	
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		experiments_RealdataSets();
	}

}
