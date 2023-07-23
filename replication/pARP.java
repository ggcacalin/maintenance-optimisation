package replication;

import java.util.ArrayList;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class pARP {

	// 12/36 months, 1/3 years. Any difference?
	private static double alpha = 12;
	private static double beta = 2;
	private static int N = 12;
	private static int M = 24;
	private static double[][] x0_opt = new double[N][M + 1];
	private static double[][] x1_opt = new double[N][M + 1];
	private static ArrayList<Integer> ages = new ArrayList<>();

	public static void main(String[] args) {
		// Various parameters and system time
		double phi = -2 * Math.PI / N;
		long start = System.currentTimeMillis();

		// Transition probability definition
		double[][][][] pi0 = new double[N][M + 1][N][M + 1];
		double[][][][] pi1 = new double[N][M + 1][N][M + 1];
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				if (i2 > 0 && i2 < M) {
					pi0[i1][i2][(i1 + 1)%N][i2 + 1] = 1 - pfail(i2+1);
					pi0[i1][i2][(i1 + 1)%N][0] = pfail(i2+1);
				}
				pi1[i1][i2][(i1 + 1)%N][1] = 1 - pfail(1);
				pi1[i1][i2][(i1 + 1)%N][0] = pfail(1);
			}
		}

		// Cost definition
		double delta = 0;
		double cpbar = 10;
		double cfbar = 50;
		double[] cp = new double[N];
		double[] cf = new double[N];
		for (int i1 = 0; i1 < N; i1++) {
			cp[i1] = cpbar * (1 + delta * Math.cos(2 * Math.PI * (i1 + 1) / N + phi));
			cf[i1] = cfbar * (1 + delta * Math.cos(2 * Math.PI * (i1 + 1) / N + phi));
		}

		try {
			solvePARP(pi0, pi1, cp, cf, N, M);
			System.out.println(ages.toString());
		} catch (IloException e) {
			System.out.println("A Cplex exception occured: " + e.getMessage());
			e.printStackTrace();
		}
		System.out.println("Running time: " + (System.currentTimeMillis() - start)/1000.0 + " seconds");
	}

	public static double wcdf(double x) {
		double exponent = -Math.pow(x / alpha, beta);
		return 1 - Math.pow(Math.E, exponent);
	}

	public static double pfail(int i2) {
		double numerator = wcdf(i2) - wcdf(i2 - 1);
		return numerator / (1 - wcdf(i2 - 1));
	}

	public static void solvePARP(double[][][][] pi0, double[][][][] pi1, double[] cp, double[] cf, int N, int M)
			throws IloException {
		// Create the model
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-5);
		// Variables
		IloNumVar[][] x0 = new IloNumVar[N][M + 1];
		IloNumVar[][] x1 = new IloNumVar[N][M + 1];
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				x1[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
				x0[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
			}
		}
		// Objective
		IloNumExpr objExpr = cplex.constant(0);
		for (int i1 = 0; i1 < N; i1++) {
			objExpr = cplex.sum(objExpr, cplex.prod(cf[i1], x1[i1][0]));
			for (int i2 = 1; i2 <= M; i2++) {
				objExpr = cplex.sum(objExpr, cplex.prod(cp[i1], x1[i1][i2]));
			}
		}
		cplex.addMinimize(objExpr);
		// Constraints 7b
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				// Subtraction left side
				IloNumExpr sumActions = cplex.constant(0);
				sumActions = cplex.sum(sumActions, x1[i1][i2]);
				sumActions = cplex.sum(sumActions, x0[i1][i2]);
				// Subtraction right side
				IloNumExpr sumTransitions = cplex.constant(0);
				for (int j1 = 0; j1 < N; j1++) {
					for (int j2 = 0; j2 <= M; j2++) {
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi1[j1][j2][i1][i2], x1[j1][j2]));
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi0[j1][j2][i1][i2], x0[j1][j2]));
					}
				}
				cplex.addEq(cplex.diff(sumActions, sumTransitions), 0);
			}
		}
		// Constraints 16c (Schouten)
		for (int i1 = 0; i1 < N; i1++) {
			cplex.addEq(x0[i1][0], 0);
			cplex.addEq(x0[i1][M], 0);
		}
		// Constraints 7c
		for (int i1 = 0; i1 < N; i1++) {
			IloNumExpr sumAges = cplex.constant(0);
			for (int i2 = 0; i2 <= M; i2++) {
				sumAges = cplex.sum(sumAges, x1[i1][i2]);
				sumAges = cplex.sum(sumAges, x0[i1][i2]);
			}
			cplex.addEq(sumAges, 1.0 / N);
		}
		// Solve the model
		cplex.solve();
		// Query the solution
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + 12*cplex.getObjValue());
			for (int i1 = 0; i1 < N; i1++) {
				int thisAge = M;
				for (int i2 = 0; i2 <= M; i2++) {
					x1_opt[i1][i2] = cplex.getValue(x1[i1][i2]);
					if (i2 > 0 && i2 < thisAge && x1_opt[i1][i2] > 0) {
						ages.add(i2);
						thisAge = i2;
					}
				}
				if (thisAge == M)
					ages.add(M);
			}

		} else {
			System.out.println("No optimal solution found");
		}
		// Close the model
		cplex.close();
	}

}
