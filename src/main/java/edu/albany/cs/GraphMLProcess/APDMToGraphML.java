package edu.albany.cs.GraphMLProcess;

import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

import org.apache.commons.lang3.ArrayUtils;

import edu.albany.cs.base.APDMInputFormat;

public class APDMToGraphML {

    // public APDMToGraphMLP(String fileName,String root,String subDir,int r,int
    // numNodes,int numFea,int clusterNumber){
    // //bio arxiv99-03 dblp_top patent_91-95
    // String name="patent1_91-95";
    // String gmlFile="data/GraphMLData/RealData/"+name+".graphml";
    // String apdmFile="data/GraphMLData/RealData/"+name+"_APDM.txt";
    // APDMInputFormat apdm = new APDMInputFormat(apdmFile);
    //
    // }

	private static void genGraphMLPFromAPDM(String fileName, String root,
			int r, int n, int numFea, int mm) {
		// TODO Auto-generated method stub
	String apdmFile = root + "Dense_APDM/cluster" + n + "/cluster-" + mm
		+ "/APDM_r_" + r + "_" + fileName;
	String gmlFileName = "GraphML_" + n + "_r_" + r + "_"
		+ fileName.substring(0, fileName.length() - 3) + "graphml";
	String graphMLPFile = root + "Dense_GraphMLP/" + gmlFileName;
	String graphMLPropFile = root + "Dense_GraphMLP/GraphML_" + n + "_r_"
		+ r + "_" + fileName.substring(0, fileName.length() - 3)
		+ "properties";

	int n_min = 5; // minimum number of nodes returned clique.
	double gamma_min = 0.5; // return clique gamm(s)>=gamma_min
	int s_min = 1; // minimum feature
	double w_max = 0.01;

	int param_a = 1;
	int param_b = 1;
	int param_c = 1;
	double r_dim = 0.5;
	double r_obj = 0.5;

	System.out.println(graphMLPFile);
		FileWriter output = null;
		FileWriter outProp = null;
		try {
	    output = new FileWriter(graphMLPFile, false);
	    outProp = new FileWriter(graphMLPropFile, false);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	/** Write the property file */
		try {
	    outProp.write("filename = " + graphMLPFile + "\n \n");
	    outProp.write("n_min = " + n_min + "\n");
	    outProp.write("gamma_min = " + gamma_min + "\n");
	    outProp.write("s_min = " + s_min + "\n");
	    outProp.write("w_max = " + w_max + "\n\n");

	    outProp.write("param_a = " + param_a + "\n");
	    outProp.write("param_b = " + param_b + "\n");
	    outProp.write("param_c = " + param_c + "\n");
	    outProp.write("r_obj = " + r_obj + "\n");
	    outProp.write("r_dim = " + r_dim + "\n");

			outProp.close();
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}
	APDMInputFormat apdm = new APDMInputFormat(apdmFile);
	// System.out.println("Data "+ArrayUtils.toString(apdm.data.attributes.length));

	// Graph myGraph = new i9.graph.gamer.graph.Graph(numFea);
	// ArrayList<Node> nodeList=new ArrayList<Node>();

	// UndirectedSparseGraph<Node, ALink> g = new UndirectedSparseGraph();
	// GraphMLWriter<Node, ALink> writer = new GraphMLWriter<Node, ALink>();

	String head = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		+ "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" "
		+ "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
		+ "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns "
		+ "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
		+ "<graph id=\"G\" edgedefault=\"undirected\">\n";

		try {
			output.write(head);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	int nodeID = 0;
	for (int i = 0; i < apdm.data.numNodes; i++) {
	    String nodeStr = "<node id=\"" + nodeID + "\" ";
	    for (int j = 0; j < apdm.data.numFeas; j++) {
		if (j != apdm.data.numFeas - 1) {
		    nodeStr += "a" + j + "=\"" + apdm.data.attributes[j][i]
			    + "\" ";
		} else {
		    nodeStr += "a" + j + "=\"" + apdm.data.attributes[j][i]
			    + "\"";
				}

			}
	    nodeStr += "/>\n";
			try {
				output.write(nodeStr);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			nodeID++;

		}

	// <edge source="14615" target="25271" weight="2.0"/>
	for (int[] edge : apdm.data.edges.keySet()) {
	    // myGraph.addEdge(edge[0], edge[1]);
	    String edgeStr = "<edge source=\"" + edge[0] + "\" target=\""
		    + edge[1] + "\"/>\n";
	    // edgeStr+="<edge source=\""+edge[1]+"\" target=\""+edge[0]+"\" weight=\"1.0\"/>\n";
			try {
				output.write(edgeStr);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		try {
			output.write("</graph>\n </graphml>\n");
			output.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	// i9.graph.gamer.graph.Graph myGraph = null;

	// writer.
	// reader.load(filename, g);
	// OutputStream out=new FileOutputStream(graphMLPFile);
	// /GraphMLWriter.outputGraph(myGraph, out);
	// i9.graph.gamer.base.Parameter.numberOfAtts =
	// myGraph.getNumberOfAtts();

	}

    public static String APDM2GraphML(APDMInputFormat apdm, String outputFile,
	    int minNumNode, double min_density, int minNumFea,
	    double max_cluster_fea,
	    double a, double b, double c, double dimRedundancy, double objRedundancy) {
		System.out.println(ArrayUtils.toString(apdm.data.attributes[0]));
	String graphMLPFile = outputFile + ".graphml";
	String graphMLPropFile = outputFile + ".properties";
	int n_min = minNumNode; // minimum number of nodes returned clique.
	double gamma_min = min_density; // return clique gamm(s)>=gamma_min
	int s_min = minNumFea; // minimum feature
	double w_max = max_cluster_fea; // the maximal extent of a cluster in the attribute space
	// Q(C)=density(O)^a . |O|^b . |S|^c
	double param_a = a; //
	double param_b = b; //
	double param_c = c; //
	double r_dim = dimRedundancy;  //
	double r_obj = objRedundancy;
	FileWriter output = null;
	FileWriter outProp = null;
	try {
	    output = new FileWriter(graphMLPFile, false);
	    outProp = new FileWriter(graphMLPropFile, false);
	} catch (IOException e) {
	    e.printStackTrace();
	}
	/** Write the property file */
	try {
	    outProp.write("filename = " + graphMLPFile + "\n \n");
	    outProp.write("n_min = " + n_min + "\n");
	    outProp.write("gamma_min = " + gamma_min + "\n");
	    outProp.write("s_min = " + s_min + "\n");
	    outProp.write("w_max = " + w_max + "\n\n");

	    outProp.write("param_a = " + param_a + "\n");
	    outProp.write("param_b = " + param_b + "\n");
	    outProp.write("param_c = " + param_c + "\n");
	    outProp.write("r_obj = " + r_obj + "\n");
	    outProp.write("r_dim = " + r_dim + "\n");

	    outProp.close();
	} catch (IOException e1) {
	    // TODO Auto-generated catch block
	    e1.printStackTrace();
	}
	
	/**Head of GraphMLFiles*/
	String head = "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n"
		+ "<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\" "
		+ "xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" "
		+ "xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns "
		+ "http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n"
		+ "<graph id=\"G\" edgedefault=\"undirected\">\n";

	try {
	    output.write(head);
	} catch (IOException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}

	int nodeID = 0;
	for (int i = 0; i < apdm.data.numNodes; i++) {
	    String nodeStr = "<node id=\"" + nodeID + "\" ";
	    for (int j = 0; j < apdm.data.numFeas; j++) {
		if (j != apdm.data.numFeas - 1) {

		    nodeStr += "a" + j + "=\"" + apdm.data.attributes[j][i]
			    + "\" ";
		} else {
		    nodeStr += "a" + j + "=\"" + apdm.data.attributes[j][i]
			    + "\"";
				}

	    }
	    nodeStr += "/>\n";
	    try {
		output.write(nodeStr);
	    } catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
	    }
	    nodeID++;
	}

	// <edge source="14615" target="25271" weight="2.0"/>
	for (int[] edge : apdm.data.edges.keySet()) {

	    // myGraph.addEdge(edge[0], edge[1]);
	    String edgeStr = "<edge source=\"" + edge[0] + "\" target=\""
		    + edge[1] + "\"/>\n";
	    // edgeStr+="<edge source=\""+edge[1]+"\" target=\""+edge[0]+"\" weight=\"1.0\"/>\n";
	    try {
		output.write(edgeStr);
	    } catch (IOException e) {
		// TODO Auto-generated catch block
		e.printStackTrace();
			}
	}
	try {
	    output.write("</graph>\n </graphml>\n");
	    output.close();
	} catch (IOException e) {
	    // TODO Auto-generated catch block
	    e.printStackTrace();
	}
	return graphMLPropFile;
	}

    public static void main(String args[]) throws FileNotFoundException {
	int N = 500;
	int[] m = new int[] { 10 };
	int numFea = 10;
	int[] r = new int[] { 2, 3, 5 };
	double p_in = 0.25;
	double p_out = 0.01;
	String root = "data/DenseGraph/";
	for (int mm : m) {
	    for (int rr : r) {
		for (int i = 0; i < 10; i++) {
		    String fileName = "Cluster_" + mm + "_in_" + p_in + "_out_"
			    + p_out + "_case_" + i + ".txt";
		    System.out.println(fileName);
		    genGraphMLPFromAPDM(fileName, root, rr, N, numFea, mm);
		}
	    }
	}

    }

}
