package edu.albany.cs.scoreFuncs;

/**
 * Types of different score functions.
 *
 * @author baojian bzhou6@albany.edu
 */
public enum FuncType {
	EMSXYScore, Unknown,  LogisticRegression,DensityProjectScore,LogLogisticRegression,PCAScore, LogisticLoss;


	public static FuncType defaultFuncType() {
		return FuncType.Unknown;
	}
}
