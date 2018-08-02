package edu.albany.cs.base;

public final class Constants {

    //prevent instantiation
    private Constants() {
    }

    public static final String serverApdm01 = "apdm01-OptiPlex-3046";
    public static final String serverApdm02 = "apdm02-OptiPlex-3046";
    public static final String desktop = "baojian-OptiPlex-9020";

    //define default alphaMax value be 0.15D ;
    public static final double AlphaMax = 0.15D;

    //Input data folders
    public static final String Input_Root_Folder = "./data";
    public static final String Input_Simulation_Folder = Input_Root_Folder + "/" + "SimulationData/GridData/SingleCom";
    public static final String Input_WaterNetwork_Folder = Input_Root_Folder + "/" + "WaterData/source_12500";

    //Output data folders
    public static final String Output_Root_Folder = "./output";
    public static final String Output_Simulation_Folder = Output_Root_Folder + "/" + "SimulationData/";
    public static final String Output_WaterNetwork_Folder = "./output";
    public static final String Output_Simulation_IHT_Folder = Output_Root_Folder + "/" + "SimulationData_IHT/";
    public static final String Output_Simulation_AmCoSamp_Folder = Output_Root_Folder + "/" + "SimulationData_AMCoSamp/";

    public enum DataSet {
        WaterNetwork,
        Transportation,
        CrimeOfChicago,
        CitHepPh,
        Weibo
    }

    //different algorithms
    public enum Algorithms {
        FuesdLasso,
        EventTree,
        LTSS,
        AdditiveScan,
        DepthFirstScan,
        GraphLaplacian,
        EdgeLasso,
        AMCoSamp,
        GraphMP,
        GraphIHT,
        GraphIHTP
    }
}