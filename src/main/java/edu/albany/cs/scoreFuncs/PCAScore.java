package edu.albany.cs.scoreFuncs;


import java.io.FileInputStream;
import java.io.ObjectInputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.apache.commons.math3.stat.StatUtils;

import com.google.common.collect.Sets;

import edu.albany.cs.base.APDMInputFormat;
import edu.albany.cs.base.Matrix;
import edu.albany.cs.graph.MultiDataGraph;

/*************************************************************\
 *min f(x,y)= \sigma_{1..n} (w_iy - xWy/1^Tx)^2 - lambda*xAx/1^Tx
 *
 \*************************************************************/
public class PCAScore implements Function {


    /**
     * weight matrix
     */
    private final double[][] W;
    
    /**
     * adjacent Matrix
     */
    private final double[][] A; 
    
    
    /**
     * dimension of a vector
     */
    private final int numNodes;
    private final int numFea;
    private final FuncType funcID;
	// double sigmas1 = 0.1D;
	double sigmas1 = 0.1D;
	double sigmas2 = 1.0D;
    private int verboseLevel = 0;

    public PCAScore(double[][] W,double[][] A) {
        funcID = FuncType.PCAScore;
        if (!checkInput(W)) {
            System.out.println(funcID + " input parameter is invalid.");
        }
        this.W = W;
        this.A=A;
        this.numNodes = W.length;
        this.numFea = W[0].length;
    }

    private boolean checkInput(double[][] W) {

        if (W == null)
            return false;
        return true;

    }
    
    public double[] median_WTx(double[] x){
    	int p = W[0].length;
    	double sum = StatUtils.sum(x);
    	double[] medians = new double[p];
    	ArrayList<Integer> S = new ArrayList<Integer>();
    	for(int i = 0; i < x.length; i++){
    		if(x[i] != 0){
    			S.add(i);
    		}
    	}
    	for(int j=0; j<p; j++){
    		double[] vals = new double[S.size()];
    		for(int k=0; k<S.size(); k++){
    			vals[k] = x[S.get(k)] * W[S.get(k)][j];
    		}
    		medians[j] = StatUtils.percentile(vals, 50) * sum;
    	}
    	return medians;
    }

//  if (y.length != numFea && x.length != numNodes) {
//  System.out.println(" GradientX: Error : Invalid parameters ...");
//  System.exit(0);
//}
//if(StatUtils.sum(x)==0){        	
//	System.out.println("GradientX:Input x vector values are all Zeros !!!");
//  System.exit(0);
//}
//if(StatUtils.sum(y)==0){        	
//	System.out.println("GradientX(: Input y vector values are all Zeros !!!");
//  System.exit(0);
//}
    public double[] getGradientX(double[] x, double[] y, double lambda) {
    	return getGradientX(x, y, lambda, 0);
    }

    public double[] getGradientX(double[] x, double[] y, double lambda, double adjust) {
        double[] Ix = new double[numNodes];
        Arrays.fill(Ix, 1.0D);
        double xT1 = Matrix.dot(x, Ix); // 1^Tx
        int ncnt = 0;
        for(int i=0; i<x.length; i++){
        	if(x[i] < 1 && x[i] > 0)
        		ncnt += 1;
        }
        double[] WTx;
		// // TODO for debug
		// ncnt = 1;
        if(ncnt > 0)
        	WTx=Matrix.VecMultiplyMat(x,W);           //Wx
        else
        	WTx = median_WTx(x);
        double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
        double[] comp1 = Matrix.VarMultiplyVec(1.0/(sigmas1 + adjust), squaretimevector(Matrix.subtract(W,Matrix.VecOutterPro(Ix, WTx_1Tx)), y));
        double[] comp2 = Matrix.VarMultiplyVec(1.0/sigmas2, squaretimevector(W, y));
        double[] ATx=Matrix.VecMultiplyMat(x,A);           //Wx
        double[] ATx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, ATx); // Wy/(1^Tx)
        double[] dens_comp = Matrix.VarMultiplyVec(lambda , Matrix.VecSubstract(ATx_1Tx, Matrix.VarMultiplyVec(1/(xT1*xT1) * Matrix.dot(ATx, x), Ix)));
        double[] gradient = Matrix.VecSubstract(Matrix.VecSubstract(comp1, comp2), dens_comp);
        return gradient;
    }    

    
//  System.out.println(Matrix.dot(ATx, x));
//  showstat(ATx_1Tx, "1");
//  showstat(Matrix.VarMultiplyVec(0.5/(xT1*xT1) * Matrix.dot(ATx, x), Ix), "2");
//    showstat(dens_comp, "density: ", 0, +1000);
    
    
    
//    /**
//     * delta_x f(x,y)=-2* \sigma_{1..n} (w_iy - xTWy/1^Tx)*(Wy/1^Tx - xTWy.1/(1^Tx)^2) - \lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2               
//     *               =-2 * \sigma rightside* leftside -\lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2
//     */
//    public double[] getGradientX(double[] x, double[] y,double lambda) {
//        if (y.length != numFea && x.length != numNodes) {
//            System.out.println(" GradientX: Error : Invalid parameters ...");
//            System.exit(0);
//        }
//        if(StatUtils.sum(x)==0){        	
//        	System.out.println("GradientX:Input x vector values are all Zeros !!!");
//            System.exit(0);
//        }
//        if(StatUtils.sum(y)==0){        	
//        	System.out.println("GradientX(: Input y vector values are all Zeros !!!");
//            System.exit(0);
//        }
//        
//        double[] gradient = new double[numNodes];
//        double[] Ix = new double[numNodes];
//        Arrays.fill(Ix, 1.0D);
//       
//       
//        
//        
//        /**
//         *  leftside= (w_iy - xTWy/1^Tx)             
//         *  rightside= (Wy/1^Tx - xTWy.1/(1^Tx)^2) 
//         *  densityTerm= -\lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2
//         * */
//        double x1T = Matrix.dot(x, Ix);        
//        double[] Wy=Matrix.MatMultiplyVec(W, y);           //Wy
//        double[] Wy_1Tx=Matrix.VarMultiplyVec(1.0D/x1T, Wy); // Wy/(1^Tx)
//        double xTWy=Matrix.dot(Matrix.VecMultiplyMat(x, W),y); 
//        
//        double[] rightSide = Matrix.VecSubstract(Wy_1Tx, Matrix.VarMultiplyVec(1.0D*xTWy/(x1T*x1T), Ix));
//        
//        /** density term */ 
//        double xTAx=0.0D;
//		double[] Ax1Tx=new double[x.length];
//		double[] xAx_1=new double[x.length];
//        double[] denseTerm=new double[x.length]; //density Term  -\lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2
//		Ax1Tx=Matrix.VarMultiplyVec(-2.0D*x1T,Matrix.MatMultiplyVec(A, x)); //2*Ax.1Tx
//		xTAx=Matrix.dot(Matrix.VecMultiplyMat(x, A), x);
//		xAx_1=Matrix.VarMultiplyVec(xTAx,Ix);								//x^TAx.1^T
//		
//		denseTerm=Matrix.VarMultiplyVec(-lambda/(x1T*x1T), Matrix.VecSubstract(Ax1Tx,xAx_1));
//		double leftSide=0.0D;
//		
//		/**\sigma_{i=0,..n-1}  (w_iy - xTWy/1^Tx)*/
//		for(int i=0;i<numNodes;i++){
//	        	leftSide+=-2.0D*(Matrix.dot(W[i], y)-xTWy/x1T);
//	    }
//           
//        gradient = Matrix.VecAdd(denseTerm,Matrix.VarMultiplyVec(leftSide,rightSide));  
//        	
//        return gradient;
//    }

    
    
//    //\delta_x f(x,y)=2* \sigma_{1..n} (w_i - W^Ty/1^Tx)(w_i - W^Ty/1^Tx)^T.y
//    public double[] getGradientY(double[] x, double[] y,double lambda) {
//
//        if (y.length != numFea && x.length != numNodes) {
//        	System.out.println(" GradientY: Error : Invalid parameters ...");
//            System.exit(0);
//        }
//        if(StatUtils.sum(x)==0){        	
//        	System.out.println("GradientY: Input x vector values are all Zeros !!!");
//            System.exit(0);
//        }
//        if(StatUtils.sum(y)==0){        	
//        	System.out.println("GradientY: Input y vector values are all Zeros !!!"); 
//        	System.exit(0);        	
//        }
//        
//        double[] gradient = new double[numNodes];
//        double[] Ix = new double[numNodes];    
//
//        Arrays.fill(Ix, 1.0D);       
//
//        double xT1 = Matrix.dot(x, Ix);
//        
//        double[] WTx=Matrix.VecMultiplyMat(x,W);           //Wx
//        double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
//        
//        /**
//         *  (w_i - W^Tx/1^Tx).(w_i - W^Tx/1^Tx)^T
//         * 
//         * */
//        //double[] leftSide = Matrix.VecSubstract(Wy_1Tx, Matrix.VarMultiplyVec(xTWy/(xT1*xT1), Ix));
//        for(int i=0;i<numFea;i++){
//        	double[] wi_BarW=Matrix.VecSubstract(W[i], WTx_1Tx);
//        	double[][] Sx=Matrix.VecOutterPro(wi_BarW,wi_BarW);
//        	//term1 Wy/1^Tx
//        	gradient = Matrix.VarMultiplyVec(2.0D * x[i],Matrix.MatMultiplyVec(Sx,y));  
//        	
//        }
//        return gradient;
//    }
    

    public double[] squaretimevector(double[][] M, double[] v){
    	double[] out = new double[M.length];
    	for(int i = 0; i < M.length; i++){
    		out[i] = 0;
    		for(int j = 0; j < v.length; j++){
    			out[i] += M[i][j] * M[i][j] * v[j];
    		}
    	}
    	return out;
    }

    public double[] getGradientY(double[] x, double[] y, double lambda) {
    	return getGradientY(x, y, lambda, 0);
    }
    
    //\delta_x f(x,y)=2* \sigma_{1..n} (w_i - W^Ty/1^Tx)(w_i - W^Ty/1^Tx)^T.y
    public double[] getGradientY(double[] x, double[] y, double lambda, double adjust) {
        if (y.length != numFea && x.length != numNodes) {
        	System.out.println(" GradientY: Error : Invalid parameters ...");
            System.exit(0);
        }
        if(StatUtils.sum(x)==0){        	
        	System.out.println("GradientY: Input x vector values are all Zeros !!!");
            System.exit(0);
        }
        if(StatUtils.sum(y)==0){        	
        	System.out.println("GradientY: Input y vector values are all Zeros !!!"); 
        	System.exit(0);        	
        }

        double[] Ix = new double[numNodes];    
        Arrays.fill(Ix, 1.0D);       
        double[][] W_transpose = Matrix.transpose(W);
        
        double xT1 = Matrix.dot(x, Ix); // 1^Tx

        int ncnt = 0;
        for(int i=0; i<x.length; i++){
        	if(x[i] < 1 && x[i] > 0)
        		ncnt += 1;
        }
        double[] WTx;
		// // TODO for debug
		// ncnt = 1;
        if(ncnt > 0)
        	WTx=Matrix.VecMultiplyMat(x,W);           //Wx
        else
        	WTx = median_WTx(x);
        
//        double[] WTx=Matrix.VecMultiplyMat(x,W);           //Wx
//        double[] WTx = median_WTx(x);
        double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
        double[] comp1 = Matrix.VarMultiplyVec(1.0/(sigmas1 + adjust), squaretimevector(Matrix.subtract(W_transpose,Matrix.VecOutterPro(WTx_1Tx, Ix)), x));
        double[] comp2 = Matrix.VarMultiplyVec(1.0/sigmas2, squaretimevector(W_transpose, x));
        double[] gradient = Matrix.VecSubstract(comp1, comp2);
        
//        ArrayIndexSort arrayIndexComparator = new ArrayIndexSort(gradient);
//        Integer[] indexes = arrayIndexComparator.getIndices();
//        Arrays.sort(indexes, arrayIndexComparator);
        
        return gradient;
    }    
    

//    public double[] getGradientY(double[] x, double[] y,double lambda) {
//
//        if (y.length != numFea && x.length != numNodes) {
//        	System.out.println(" GradientY: Error : Invalid parameters ...");
//            System.exit(0);
//        }
//        if(StatUtils.sum(x)==0){        	
//        	System.out.println("GradientY: Input x vector values are all Zeros !!!");
//            System.exit(0);
//        }
//        if(StatUtils.sum(y)==0){        	
//        	System.out.println("GradientY: Input y vector values are all Zeros !!!"); 
//        	System.exit(0);        	
//        }
//        
//
//        double[] gradient = new double[numFea];
//        double[] Ix = new double[numNodes];    
//
//        Arrays.fill(Ix, 1.0D);       
//
//        double xT1 = Matrix.dot(x, Ix);
//        
//        double[] WTx=Matrix.VecMultiplyMat(x,W);           //Wx
//        double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
//               
//        /**
//         *  2.0*\sum_i x_i (w_i - W^Tx/1^Tx).(w_i - W^Tx/1^Tx)^T.y
//         * 
//         * */
//        double[] Sxy=new double[numFea];
//        Arrays.fill(Sxy, 0.0D);  
//        for(int i=0;i<numNodes;i++){
//        	double[] wi_BarW=Matrix.VecSubstract(W[i], WTx_1Tx);
//        	double[][] Sx_i=Matrix.VecOutterPro(wi_BarW,wi_BarW);
//        	//term1 S(x)y
//        	if(x[i]!=0.0D){
//        		Sxy = Matrix.VecAdd(Sxy,Matrix.VarMultiplyVec(x[i],Matrix.MatMultiplyVec(Sx_i,y)));
//        	}else{
//        		continue;
//        	} 
//        	
//        }
//        gradient=Matrix.VarMultiplyVec(2.0D,Sxy); //2*S(x).y
//        return gradient;
//    }    
    
    @Override
    public List<double[]> getArgMinFxy(double[] x0, double[] y0, Set<Integer> OmegaX, Set<Integer> OmegaY,double lambda) {
        return new GradientDescentOpt(this,numNodes,numFea,0.1,100).multiGradientDescent(x0, y0, OmegaX, OmegaY, lambda);
    }

    @Override
    public double[] getArgMinFx(Set<Integer> Sx, int[] Sy, double[] yi, Set<Integer> OmegaX, double lambda) {
        return new GradientDescentOpt(this,numNodes,numFea,0.01).sigGradientDescent_x(Sx, Sy, yi, OmegaX, lambda);
    }

    @Override
    public double[] getArgMinFy(double[] xi, Set<Integer> OmegaY, double lambda) {
        return new GradientDescentOpt(this,numNodes,numFea,0.01).sigGradientDescent_y(xi, OmegaY, lambda, 5);
    }
    

    /** f(x,y)= \sigma_{1..n} ((w_i - xW/1^Tx)y)^2        */
    public double getFuncValue(double[] x, double[] y,double lambda) {

    	
        double[] Ix = new double[numNodes];
        Arrays.fill(Ix, 1.0D);
        double xT1 = Matrix.dot(x, Ix); // 1^Tx
        double[] WTx=Matrix.VecMultiplyMat(x,W);           //Wx
        double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
        double[] comp1 = Matrix.VarMultiplyVec(1.0/sigmas1, squaretimevector(Matrix.subtract(W,Matrix.VecOutterPro(Ix, WTx_1Tx)), y));
        double[] comp2 = Matrix.VarMultiplyVec(1.0/sigmas2, squaretimevector(W, y));
        double[] diff = Matrix.VecSubstract(comp1, comp2);
        double funvalue = Matrix.dot(diff, x);
        if(lambda > 0){
            double[] ATx=Matrix.VecMultiplyMat(x,A);           //Wx
            double[] ATx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, ATx); // Wy/(1^Tx)
            funvalue = funvalue - lambda * Matrix.dot(ATx_1Tx, x);
        }
        return funvalue;
        
//        double funcValue = 0.0D;
//        
//        if (x.length != W.length || y.length != W[0].length) {
//            new IllegalArgumentException("Error : Invalid parameters ...");
//            System.exit(0);
//        }
//        
//        /** Dense term*/
//        double xTAx=0.0D;		
//		double x1T=StatUtils.sum(x); //1^Tx		
//		double denseTerm=0.0D;           //term2  \lambda (x^TAx) / 1^T x
//		xTAx=Matrix.dot(Matrix.VecMultiplyMat(x, A), x);		
//		denseTerm=lambda*xTAx/x1T;
//		
//        double[] Ix = new double[numNodes];
//        Arrays.fill(Ix, 1.0D);       
//
//        double[] WTx = Matrix.MatMultiplyVec(Matrix.transpose(W), x);
//        double[] BarW=Matrix.VarMultiplyVec(1/x1T, WTx); 
//        for(int i=0;i<numNodes;i++){
//        	double temp=Matrix.dot(Matrix.VecSubstract(W[i], BarW),y); //((w_i - xW/1^Tx)
//        	funcValue+=temp*temp*x[i];
//        }
//        funcValue=funcValue+denseTerm;
//        return funcValue;
    }

//    @Override
//    public FuncType getFuncID() {
//        return funcID;
//    }
//
//	@Override
//	public double getFuncValue(double[] x, double[] y) {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public double[] getGradientX(double[] x, double[] y) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public double[] getGradientY(double[] x, double[] y) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public List<double[]> getArgMinFxy(Set<Integer> OmegaX, Set<Integer> OmegaY) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public List<double[]> getArgMinFy(double[] xi) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public List<double[]> getArgMinFx(double[] yi) {
//=======

//    private int verboseLevel = 1;
//
//    public PCAScore(double[][] W,double[][] A) {
//        funcID = FuncType.PCAScore;
//        if (!checkInput(W)) {
//            System.out.println(funcID + " input parameter is invalid.");
//        }
//
//        this.W = W;
//        this.A=A;
//        this.numNodes = W.length;
//        this.numFea = W[0].length;
//    }
//
//    private boolean checkInput(double[][] W) {
//
//        if (W == null)
//            return false;
//        return true;
//
//    }
//
//    /**
//     * delta_x f(x,y)=-2* \sigma_{1..n} (w_iy - xTWy/1^Tx)*(Wy/1^Tx - xTWy.1/(1^Tx)^2) - \lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2               
//     *               =-2 * \sigma rightside* leftside -\lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2
//     */
//    public double[] getGradientX(double[] x, double[] y,double lambda) {
//        if (y.length != numFea && x.length != numNodes) {
//            System.out.println(" GradientX: Error : Invalid parameters ...");
//            System.exit(0);
//        }
//        if(StatUtils.sum(x)==0){        	
//        	System.out.println("GradientX:Input x vector values are all Zeros !!!");
//            System.exit(0);
//        }
//        if(StatUtils.sum(y)==0){        	
//        	System.out.println("GradientX(: Input y vector values are all Zeros !!!");
//            System.exit(0);
//        }
//        
//        double[] gradient = new double[numNodes];
//        double[] Ix = new double[numNodes];
//        Arrays.fill(Ix, 1.0D);        
//        
//        /**
//         *  leftside= (w_iy - xTWy/1^Tx)             
//         *  rightside= (Wy/1^Tx - xTWy.1/(1^Tx)^2) 
//         *  densityTerm= -\lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2
//         * */
//        
//        /**rightside= (Wy/1^Tx - xTWy.1/(1^Tx)^2)  */
//        double x1T = Matrix.dot(x, Ix);        
//        double[] Wy=Matrix.MatMultiplyVec(W, y);           //Wy
//        double[] Wy_1Tx=Matrix.VarMultiplyVec(1.0D/x1T, Wy); // Wy/(1^Tx)
//        double xTWy=Matrix.dot(Matrix.VecMultiplyMat(x, W),y); 
//        
//        double[] rightSide = Matrix.VecSubstract(Wy_1Tx, Matrix.VarMultiplyVec(1.0D*xTWy/(x1T*x1T), Ix));
//        
//        /** density term */ 
//        double xTAx=0.0D;
//		double[] Ax1Tx=new double[x.length];
//		double[] xAx_1=new double[x.length];
//        double[] denseTerm=new double[x.length]; //density Term  -\lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2
//		Ax1Tx=Matrix.VarMultiplyVec(-2.0D*x1T,Matrix.MatMultiplyVec(A, x)); //2*Ax.1Tx
//		xTAx=Matrix.dot(Matrix.VecMultiplyMat(x, A), x);
//		xAx_1=Matrix.VarMultiplyVec(xTAx,Ix);								//x^TAx.1^T
//		
//		denseTerm=Matrix.VarMultiplyVec(-lambda/(x1T*x1T), Matrix.VecSubstract(Ax1Tx,xAx_1));
//		double leftSide=0.0D;
//		
//		/**leftside: -2.0D * \sigma_{i=0,..n-1}  (w_iy - xTWy/1^Tx)*/
//		for(int i=0;i<numNodes;i++){
//	        	leftSide+=-2.0D*(Matrix.dot(W[i], y)-xTWy/x1T);
//	    }
//           
//        gradient = Matrix.VecAdd(denseTerm,Matrix.VarMultiplyVec(leftSide,rightSide));  
//        	
//        
//        
//        
//        return gradient;
//    }
//
//    //\delta_x f(x,y)=2* \sigma_{1..n} (w_i - W^Ty/1^Tx)(w_i - W^Ty/1^Tx)^T.y
//    public double[] getGradientY(double[] x, double[] y,double lambda) {
//
//        if (y.length != numFea && x.length != numNodes) {
//        	System.out.println(" GradientY: Error : Invalid parameters ...");
//            System.exit(0);
//        }
//        if(StatUtils.sum(x)==0){        	
//        	System.out.println("GradientY: Input x vector values are all Zeros !!!");
//            System.exit(0);
//        }
//        if(StatUtils.sum(y)==0){        	
//        	System.out.println("GradientY: Input y vector values are all Zeros !!!"); 
//        	System.exit(0);        	
//        }
//        
//
//        double[] gradient = new double[numFea];
//        double[] Ix = new double[numNodes];    
//
//        Arrays.fill(Ix, 1.0D);       
//
//        double xT1 = Matrix.dot(x, Ix);
//        
//        double[] WTx=Matrix.VecMultiplyMat(x,W);           //Wx
//        double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
//               
//        
//        /**
//         *  2.0*\sigma_{i=1..n} x_i*(w_i - W^Tx/1^Tx).(w_i - W^Tx/1^Tx)^T.y
//         * 
//         * */
//        double[][] Sx=new double[numFea][numFea];
//        for(int i=0;i<numNodes;i++){
//        	double[] wi_BarW=Matrix.VecSubstract(W[i], WTx_1Tx);
//        	double[][] Sx_i=Matrix.VecOutterPro(wi_BarW,wi_BarW);
//        	//term1 S(x)y
//        	if(x[i]!=0.0D){
//        		//Sxy = Matrix.VarMultiplyVec(x[i],Matrix.VecAdd(Sxy,Matrix.MatMultiplyVec(Sx_i,y)));
//        		//\sigma_{i=1..n} x_i*(w_i - W^Tx/1^Tx).(w_i - W^Tx/1^Tx)^T
//        		Sx = Matrix.add(Sx,Matrix.VarMultiplyMat(x[i], Sx_i)); 
//        	}else{
//        		continue;
//        	} 
//        	
//        }
//        gradient=Matrix.VarMultiplyVec(2.0D,Matrix.MatMultiplyVec(Sx, y) ); //2*S(x).y
//        return gradient;
//    }
//
//    @Override
//    public List<double[]> getArgMinFxy(Set<Integer> OmegaX, Set<Integer> OmegaY,double lambda) {
//        return new GradientDescentOpt(this,numNodes,numFea,0.1).multiGradientDescent(OmegaX,OmegaY, lambda);
//    }

    @Override
    public List<double[]> getArgMinFx(double[] yi) {
        return null;
    }

    @Override
    public List<double[]> getArgMinFy(double[] xi) {
        return null;
    }

	@Override
	public FuncType getFuncID() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double getFuncValue(double[] x, double[] y) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double[] getGradientX(double[] x, double[] y) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double[] getGradientY(double[] x, double[] y) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<double[]> getArgMinFxy( Set<Integer> OmegaX, Set<Integer> OmegaY) {
		// TODO Auto-generated method stub
		return null;
	}
    public static int[] Set2IntArray(Set<Integer> x){
    	int[] data = new int[x.size()];
    	int i = 0;
    	for(int j:x){
    		data[i] = j;
    		i += 1;
    	}
    	return data;
    }

    public static Set<Integer> IntArray2Set(int[] x){
    	Set<Integer> mySet = new HashSet<Integer>();
    	for(int i:x){
    		mySet.add(i);
    	}
    	return mySet;
    }
	
    public static Set<Integer> GetNonzeroEntries(double[] x){
    	Set<Integer> S = new HashSet<Integer>();
    	for(int i=0; i<x.length; i++){
    		if(x[i] != 0){
    			S.add(i);
    		}
    	}
    	return S;
    }

    public static void printarray(double[] x){
    	System.out.println(array2string(x));
    }
    public static String array2string(double[] x){
    	String line = "";
    	for(double i:x){
    		line += "," + i;
    	}
    	return line; 
    }
    
    public static void showstat(double[] x, String title){
        String line = "";
        for(int j = 0; j < x.length; j++){
        	if(x[j] != 0){
        		line += ", (" + j + ", " + x[j] + ")";
        	}
        }
        System.out.println(title + line);
    }

    public static void showstat(double[] x, String title, double lowerbound, double upperbound){
        String line = "";
        int n = 0;
        int n1 = 0;
        for(int j = 0; j < x.length; j++){
        	if(x[j] != 0 && x[j] > lowerbound && x[j] < upperbound){
        		line += ", (" + j + ", " + x[j] + ")";
//        		line += "," + j;
        		n += 1;
        		if(j >= 300 && j < 400){
        			n1 += 1;
        		}
        	}
        }
        System.out.println(title + line);
        System.out.println("In total: " + n + ", " + n1);
    } 

    public static void CalcStat(List<double[]> minXY, MultiDataGraph mgdGraph){
//    	showstat(minXY.get(0), "X:", -1000, 1000);
//    	showstat(minXY.get(1), "Y:", -1000, 1000);
    	Set<Integer> S = GetNonzeroEntries(minXY.get(0));
    	Set<Integer> R = GetNonzeroEntries(minXY.get(1));
		Set<Integer> intersect = Sets.intersection(S, mgdGraph.trueNodesSet);
		double precision = (intersect.size() * 1.0D / S.size() * 1.0D);
		double recall = (intersect.size() * 1.0D / mgdGraph.trueNodesSet.size() * 1.0D);
		Set<Integer> yTrueArrayList= Stat.Array2Set(mgdGraph.trueFeatures);
		Set<Integer> intersect2 = Sets.intersection(R, yTrueArrayList);
		double precision2 = (intersect2.size() * 1.0D / R.size() * 1.0D);
		double recall2 = (intersect2.size() * 1.0D / yTrueArrayList.size() * 1.0D);
		System.out.println("x: precision: " + precision + ", recall:" + recall + "\r\ny: precision:  " + precision2+", recall " + recall2);
    }
    
    /*
     *  OUTPUT: 
     *  x: precision: 1.0, recall:0.94
	 *	y: precision: 1.0, recall 0.4
     */
    public static void unittest1(){
    	double[] xi = null;
    	double[] yi = null;
    	Set<Integer> omegaX = new HashSet<Integer>();    	
    	Set<Integer> omegaY = new HashSet<Integer>();    	
    	String fileName = "";
    	double lambda = 0; 
        try {
            FileInputStream istream = new FileInputStream("testcases/getArgMinFxy/test-case1.dat");
            ObjectInputStream ois = new ObjectInputStream(istream);
            xi = (double []) ois.readObject(); 
            yi = (double []) ois.readObject(); 
            omegaX = IntArray2Set((int[]) ois.readObject());
            omegaY = IntArray2Set((int[]) ois.readObject());
            lambda = (double) ois.readObject(); 
            fileName = (String) ois.readObject();
        } catch(Exception ex) {
            ex.printStackTrace();
        }
		APDMInputFormat apdm=new APDMInputFormat(fileName);
		MultiDataGraph mgdGraph=new MultiDataGraph(apdm);
    	PCAScore func=new PCAScore(mgdGraph.W,mgdGraph.AdjMatrix);	
    	List<double[]> minXY = func.getArgMinFxy(xi, yi, omegaX, omegaY, lambda);
    	CalcStat(minXY, mgdGraph);
    	
    }
    public static void main(String[] args){
    	unittest1();
    }
	
    
//    /** f(x,y)= \sigma_{1..n} ((w_i - xW/1^Tx)y)^2        */
//    public double getFuncValue(double[] x, double[] y,double lambda) {
//
//        double funcValue = 0.0D;
//        
//        if (x.length != W.length || y.length != W[0].length) {
//            new IllegalArgumentException("Error : Invalid parameters ...");
//            System.exit(0);
//        }
//        
//        /** Dense term*/
//        double xTAx=0.0D;		
//		double x1T=StatUtils.sum(x); //1^Tx		
//		double denseTerm=0.0D;           //term2  \lambda (x^TAx) / 1^T x
//		xTAx=Matrix.dot(Matrix.VecMultiplyMat(x, A), x);		
//		denseTerm=lambda*xTAx/x1T;
//		
//        double[] Ix = new double[numNodes];
//        Arrays.fill(Ix, 1.0D);     
//        
//        double xT1 = Matrix.dot(x, Ix);
//        
//        double[] WTx=Matrix.MatMultiplyVec(Matrix.transpose(W), x);          //Wx
//        double[] WTx_1Tx=Matrix.VarMultiplyVec(1.0D/xT1, WTx); // Wy/(1^Tx)
//        
//        double[][] Sx=new double[numFea][numFea];
//        for(int i=0;i<numNodes;i++){
//        	double[] wi_BarW=Matrix.VecSubstract(W[i], WTx_1Tx);
//        	double[][] Sx_i=Matrix.VecOutterPro(wi_BarW,wi_BarW);
//        	//term1 S(x)y
//        	if(x[i]!=0.0D){
//        		//Sxy = Matrix.VarMultiplyVec(x[i],Matrix.VecAdd(Sxy,Matrix.MatMultiplyVec(Sx_i,y)));
//        		//\sigma_{i=1..n} x_i*(w_i - W^Tx/1^Tx).(w_i - W^Tx/1^Tx)^T
//        		Sx = Matrix.add(Sx,Matrix.VarMultiplyMat(x[i], Sx_i)); 
//        	}else{
//        		continue;
//        	} 
//        	
//        }
//        
//        funcValue=Matrix.dot(Matrix.VecMultiplyMat(y, Sx),y); //yS(s)y
//        
//        funcValue=funcValue+denseTerm;
//        return funcValue;
//    }

//    @Override
//    public FuncType getFuncID() {
//        return funcID;
//    }
//
//	@Override
//	public double getFuncValue(double[] x, double[] y) {
//		// TODO Auto-generated method stub
//		return 0;
//	}
//
//	@Override
//	public double[] getGradientX(double[] x, double[] y) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public double[] getGradientY(double[] x, double[] y) {
//		// TODO Auto-generated method stub
//		return null;
//	}
//
//	@Override
//	public List<double[]> getArgMinFxy(Set<Integer> OmegaX, Set<Integer> OmegaY) {
//		// TODO Auto-generated method stub
//		return null;
//	}

    
}
