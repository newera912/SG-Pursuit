package edu.albany.cs.scoreFuncs;


import edu.albany.cs.DenseGraphMP.FunctionX;
import edu.albany.cs.base.Matrix;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.StatUtils;

import java.math.BigDecimal;
import java.util.*;

/*************************************************************\
 *argmin_x = -(x^TAx) / 1^T x
 * 
\*************************************************************/
public class DensityProjectScore implements FunctionX {

	
	/** Adjacent Matrix */
	private final double[][] A;
	
	/** p(b)=b_s b=supp(x) */
	private double[] b;
	
	private double lambda;
	
	/** dimension of a vector */
	private final int numNodes;	
	private final FuncType funcID;
	
	private int verboseLevel = 1;
	
	public DensityProjectScore(double[] b,double lambda,double[][] A){
		funcID = FuncType.DensityProjectScore;
		if (!checkInput(A)) {
			System.out.println(funcID + " input parameter is invalid.");
		}
		this.b=b;
		
		this.lambda=lambda;
		this.A= A;
		this.numNodes = A.length;
		
	}

	private boolean checkInput(double[][] W) {

		if (W ==null) 
			return false;		
		return true;
		
	}

	//\delta_x f(x)=-(b.b) - \lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2
	public double[] getGradient(double[] x) {
		if (x.length!=numNodes) {
			new IllegalArgumentException("Error : Invalid parameters ...");
			System.exit(0);
		}
		double[] gradient = new double[numNodes];
		double[] Ix=Matrix.identityVector(numNodes);

		double xTAx=0.0D;
		double[] Ax1Tx=new double[x.length];
		double[] xAx_1=new double[x.length];
		//double[] b_b=new double[x.length];
		double x1T=StatUtils.sum(x); //1^Tx
//		double[] term1=new double[x.length]; //term1 (b.b)
		double[] term2=new double[x.length]; //term2  \lambda (2A^Tx 1^Tx - x^TAx.1^T)/( 1^Tx)^2

//		term1=Matrix.elemWisePro(b,b);       // Term1 : (b.b)
		Ax1Tx=Matrix.VarMultiplyVec( 2.0D*x1T,Matrix.MatMultiplyVec(A, x)); //2*Ax.1Tx
		xTAx=Matrix.dot(Matrix.VecMultiplyMat(x, A), x);
		xAx_1=Matrix.VarMultiplyVec(xTAx,Ix);
		term2=Matrix.VecSubstract(Ax1Tx,xAx_1);
//		gradient=Matrix.VecAdd(Matrix.VarMultiplyVec(-1.0D,term1), Matrix.VarMultiplyVec(-lambda/(x1T*x1T), term2));
		gradient = Matrix.VarMultiplyVec(-1/(x1T*x1T), term2);
		return gradient;
	}

	

	/** - x^T*(b \dot b) - \lambda (x^TAx) / 1^T x */
	public double getFuncValue(double[] x) {
		
		double funcValue = 0.0D;
		
		if(x.length!=A.length){
			new IllegalArgumentException("Error : Invalid parameters ...");
			System.exit(0);
		}
			
				
		double xTAx=0.0D;		
		double x1T=StatUtils.sum(x); //1^Tx
		double term1=0.0D;           //term1  x^T*(b \dot b)
		double term2=0.0D;           //term2  \lambda (x^TAx) / 1^T x
		
		
		term1=Matrix.dot(x, Matrix.elemWisePro(b, b));
		xTAx=Matrix.dot(Matrix.VecMultiplyMat(x, A), x);
		
//		term2=lambda*xTAx/x1T;
		term2=xTAx/x1T;
		
		//System.out.println("||| x.bb =" +term1 +" lambda*den(S)="+term2);
//		funcValue=-term1-term2;
		funcValue=-term2;
	
		return funcValue;
	}

	@Override
	public FuncType getFuncID() {
		return funcID;
	}


	public double[] getGradient(int[] S) {
		double[] x = new double[numNodes];
		Arrays.fill(x, 0.0D);
		for (int i : S) {
			x[i] = 1.0D;
		}
		return getGradient(x);
	}


	public double getFuncValue(int[] S) {
		double[] x = new double[numNodes];
		Arrays.fill(x, 0.0D);
		for (int i : S) {
			x[i] = 1.0D;
		}
		return getFuncValue(x);
	}

	public double[] getArgMinFx(ArrayList<Integer> S, boolean debug) {
		/**
		 * as our objective function is f(S), we do min_{S} -f(S), this is
		 * equivalent to maximize f(S)
		 */
		double[] x =Stat.getRandomVector(A.length);	
		
		double[] gradX=new double[x.length];
		double[] Ix=Stat.getIndicateVec(new HashSet(S), x.length);
		//System.out.println("Ix"+Stat.supp(Ix));		
		int maxIteration=10000;
//		double stepSize=0.01D;
		double stepSize = 1.0D;
		Random random= new Random();
		
		if (verboseLevel > 0 && debug == true) {
			System.out.println("x0:       "+ArrayUtils.toString(x));
			System.out.println("supp(x0): "+ArrayUtils.toString(getSupportNodes(x)));					
		}
		/*Gradient descent fixed step size**/
		double[] xOld=new double[numNodes];
		double[] gradXOld=new double[x.length];
		
		x=Matrix.elemWisePro(x, Ix);
		
		int endIteration=0;
		double gap=0;
		for(int i=0;i<maxIteration;i++){
				
			gradX=getGradient(x);				
			System.arraycopy( x, 0, xOld, 0, x.length );	
			System.arraycopy( gradX, 0, gradXOld, 0, gradX.length );
			for(int j:S){
				x[j]=(x[j] - stepSize*gradX[j]);	 // update x
				/*Projecting x to [0,1]^n*/
				if(x[j]<0.0D || x[j]==-0.0D){
					x[j]=0.0D;	
				}else if(x[j]>1.0){
					x[j]=1.0D;
				}
			}		
			
			if (debug == true){
				double max = x[0];
				List myList = new ArrayList();
				for(int j = 1; j < x.length; j++){
					max = Math.max(x[j], max);
					if (x[j] > 0 && x[j] < 1){
						myList.add(x[j]);
					}
				}
				System.out.println("Maximum value: " + max + myList);
			}
			
			
			if (verboseLevel > 0 && debug == true) {
				System.out.println("-------------- GD iter= "+i+" -----------------");
				System.out.println("fun(x): "+ArrayUtils.toString(getFuncValue(x)));
				System.out.println("x: "+ArrayUtils.toString(x));
				System.out.println("oldx: "+ArrayUtils.toString(xOld));
				System.out.println("gradX: "+Stat.printWithIndex(gradX));
				
			}		
			
			 gap=getDifferenceNorm(x,xOld);
			if(gap<=1.e-6){
				endIteration=i;
				break;
			}
			
		}		
		if (debug == true){
			System.out.println("Gap="+gap+" "+"-------------- GD iter= "+endIteration+" -----------------");
		}
		x=Matrix.elemWisePro(x, Ix);
		return x;
	}
	
	public double[] getArgMinFx() {
		/**
		 * as our objective function is f(S), we do min_{S} -f(S), this is
		 * equivalent to maximize f(S)
		 */
		return getArgMaxFx(null);
	}

	/**
	 * This is a gradient descent method. It maximizes objective function f(x,y).
	 * 
	 * @return the maximizer of objective function
	 */
	public double[] getArgMaxFx(ArrayList<Integer> S) {
		return getArgMaxFx(S);
	}




	
	public double[] getIndicateVector(int[] S,int indicator) {
		double[] x = new double[numNodes];
		
		
		Arrays.fill(x, 0.0D);
		for (int i : S) {
			x[i] = 1.0D;
		}
		return x;
	}
	
	private int[] getSupportNodes(double[] x) {
		int[] nodes = null;
		for (int i = 0; i < x.length; i++) {
			if (x[i] != 0.0D) {
				/** get nonzero nodes */
				nodes = ArrayUtils.add(nodes, i);
			}
		}
		Arrays.sort(nodes);
		return nodes;
	}


	public List<double[]> getArgMinFx(int[] OmigaX, int[] OmigaY) {
		// TODO Auto-generated method stub
		return null;
	}







	@Override
	public BigDecimal[] getGradientBigDecimal(BigDecimal[] x) {
		// TODO Auto-generated method stub
		return null;
	}

    private double getDifferenceNorm(double[] x1, double[] x2) {
        double l2norm = 0.0D;
        if (x1 == null || x1.length == 0.0D) {
            return l2norm;
        }
        for (int j = 0; j < x1.length; j++) {
            l2norm += (x1[j] - x2[j]) * (x1[j] - x2[j]);
        }
        return Math.sqrt(l2norm);
    }

	@Override
	public double[] getArgMinFx(ArrayList<Integer> S) {
		// TODO Auto-generated method stub
		return null;
	}
}
