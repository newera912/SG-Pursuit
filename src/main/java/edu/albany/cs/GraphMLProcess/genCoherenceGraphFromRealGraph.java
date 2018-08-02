package edu.albany.cs.GraphMLProcess;

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

public class genCoherenceGraphFromRealGraph {
	public static void genAPDMFromRealClusterData(String fileName,
			String apdmFile, int r, int numFea, int clusterNumber, double sigma) {
		int numFestures = numFea;
		double[][] data = null;
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
			for (String eachLine : Files.readAllLines(Paths.get(fileName))) {

				/** Read Edge info */
				if (eachLine.split(" ").length < 2) {
					int numOfNodes = Integer.parseInt(eachLine);
					data = new double[numFestures][numOfNodes];

				} else if (eachLine.split(" ").length == 2) {
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
		System.out.println("True nodes:" + abnormalNodes);
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
							trueFeat, i)], sigma).sample();
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

	public static void main(String args[]) {
		String root = "data/CrimesOfChicago/";
		String fileName = "graph/CrimeOfChicagoClusterGraph.txt";
		String outputFile = root
				+ "APDM_simu/CrimeOfChicagoClusterGraph_APDM.txt";
		int numFeas = 10;
		int r = 5;
		int clusterNumber = 20;
		double sigma = 0.001D;

		genAPDMFromRealClusterData(root + fileName, outputFile, r, numFeas,
				clusterNumber, sigma);
	}
}
