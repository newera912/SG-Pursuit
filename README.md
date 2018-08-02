-------------------------------------------------------------------------------------------------------
--------                  Welcome to the readme.txt of SGPursuit Algorithm.                    --------
-------------------------------------------------------------------------------------------------------


------------------------------------------ General Information ----------------------------------------

This section contains the general information about our code and dataset. 
This folder contains our proposed algorithm, SGPursuit, which is under review of ICDM-2017. 

In this folder, there are  three sub-folders:
1. code: contains our SGPursuit algorithm and the baseline methods used in our paper.
2. dataset: has 3 datasets that we used in our experiments.
3. dataProcessing: contains the code and data for preprocessing of the real data.

!!!Important: Before you go, make sure that these 3 folders are in the same parent folder.
I assume that your put these three sub-folders into {your/root/path} directory.
-------------------------------------------------------------------------------------------------------


-------------------------------------------------------------------------------------------------------
--------                        How to run our SGPursuit Algorithm?                            --------
-------------------------------------------------------------------------------------------------------

Firstly, go to {your/root/path}/code/SGPursuit/ directory. There are 2 ways to run SGPursuit. 

1. use mvn exec:java command:
mvn exec:java -Dexec.mainClass="edu.albany.cs.SGPursuit.TestOnCrimesOfChicago"  for Crimes of chicago
mvn exec:java -Dexec.mainClass="edu.albany.cs.SGPursuit.TestOnYelp"  		for Yelp data
mvn exec:java -Dexec.mainClass="edu.albany.cs.DenseSubGraph.DenseSubGraphUnitTest" for simulation data
mvn exec:java -Dexec.mainClass="edu.albany.cs.DenseSubGraph.MultiClusterSGPursuit.RealdataSetsExperiment" for Real dataset data
!!!Important: You should install maven and Java-1.8 before you run this program.
>>>[Datasets]: 
[Simualtion Dataset Main Folder] data/DenseGraph/DenseSubgraph_APDM/
[Real world Dataset(GraphML files and apdm files) ] data/DenseGraph/Realdata/apdm/

2. use Intellij IDE:
You can use Intellij IDE to import our code as a maven project. 
After you set up the project, then you can find TestOnCrimesOfChicago.java, TestOnYelp.java, and
DenseSubGraphUnitTest.java. Finally, you can run it.


-------------------------------------------------------------------------------------------------------
--------                         How to run the baseline methods?                              --------
-------------------------------------------------------------------------------------------------------

Baseline methods are FocusCO,GAMer,SODA,AMEN. 
SODA and AMEN are in code main folder.
FocusCO and GAMer are in "{SGPursuit project}/edu/albany/cs/baselines/".

1. AMEN: go to the folder {your/root/path}/code/AMEN. Make sure you have matlab environment.
Then just type the following command to go to maltab environment:
                >>matlab -nodesktop
                >>test_amen_on_crimes_of_chicago
2. SODA: go to {your/root/path}/code/SODA. Make sure you have python and cvxopt (and solvers).
                >>python soda_crimesOfchicago.py
3. FocusCO:
[Simualtion Dataset Main Folder] /SG-Pursuit/data/DenseGraph/DenseSubgraph_APDM/ 
[Sub foldors] 1)VaryingNumOfAttributes 2)VaryingNumOfClusters 3)VaryingClusterSizes  
Requirements: 1. A recent version of Matlab
	          2.Java 6+
There are two programs in this baseline method. The first is a matlab script which learns a distance 
metric and reweighs the input graph. The second is a java program which extracts communities & outliers 
from the reweighted graph. We use a matlab proxy tool to call matlab code from the java project. 

To run FocusCO use maven:
"mvn exec:java -Dexec.mainClass="edu.albany.cs.baselines.FocusCO.FocusCOMatlab.TestFocusCO.focusCOTest1VaryNumOfAtt"
this expeiment function coressponding to the Figure 4 in our paper.(or you can test other experiment 
function in same class: focusCOTest1VaryNumOfAtt(),focusCOTest2VaryNumOfCluster(),
focusCOTest3VaryClusterSize(). Note*: you can download the original code from github: 
https://github.com/phanein/focused-clustering. 

4. GAMer:
This method provided by the authors in page: "http://dme.rwth-aachen.de/en/gamer". 
They provided a executable jar file. You can find our test code in 
"edu.albany.cs.baselines.GAMer". 
4.1 Simualtion dataset experiment
[Simualtion Dataset Main Folder] /SG-Pursuit/data/DenseGraph/DenseSubgraph_APDM/ 
[Sub foldors] 1)VaryingNumOfAttributes 2)VaryingNumOfClusters 3)VaryingClusterSizes 
To run GAMer test code use maven:
"mvn exec:java -Dexec.mainClass="edu.albany.cs.baselines.GAMer.TestGAMer.test1VaryNumOfAtt" 
this expeiment resutls coressponding to the Figure 4 in our paper.(or you can test other experiment 
function in same class: test2VaryNumOfCluster(),test3VaryClusterSize() 
4.2 Real dataset experiemnt
[Real world Dataset(GraphML files) ] /SG-Pursuit/data/DenseGraph/Realdata/. 
this expeiment resutls coressponding to the Table 4 in our paper.
To run GAMer realdata test code use maven:
"mvn exec:java -Dexec.mainClass="edu.albany.cs.baselines.GAMer.TestGAMerRealDataset.testRealdatasets" 
-------------------------------------------------------------------------------------------------------
---------                               Data processing                                       ---------
-------------------------------------------------------------------------------------------------------
In case that you are interested in our data preprocessing. The folder dataProcessing contains the 
raw data and the pre-processing code. The pre-processing code is not documented very well. However,
you can email me (bzhou6@albany.edu / aalimu@albany.edu) if you have any questions.
