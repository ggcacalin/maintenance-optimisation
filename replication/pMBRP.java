package replication;
import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class pMBRP {

	// 12/36 months, 1/3 years. Any difference?
	private static double alpha = 12;
	private static double beta = 2;
	private static int N = 12;
	private static int M = 2 * (int) alpha;
	private static int m = 1;
	private static double[] y_opt = new double[m * N];
	private static double[] t_opt = new double[m * N];

	public static void main(String[] args) {
		// Various parameters and system time
		double phi = -2 * Math.PI / N;
		long start = System.currentTimeMillis();

		// Transition probability definition
		double[][][][] pi0 = new double[m * N][M + 1][m * N][M + 1];
		double[][][][] pi1 = new double[m * N][M + 1][m * N][M + 1];
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				if (i2 > 0 && i2 < M) {
					pi0[i1][i2][(i1 + 1) % (m * N)][i2 + 1] = 1 - pfail(i2+1);
					pi0[i1][i2][(i1 + 1) % (m * N)][0] = pfail(i2+1);
				}
				pi1[i1][i2][(i1 + 1) % (m * N)][1] = 1 - pfail(1);
				pi1[i1][i2][(i1 + 1) % (m * N)][0] = pfail(1);
			}
		}

		// Cost definition
		double delta = 0;
		double cpbar = 10;
		double cfbar = 50;
		double[] cp = new double[m * N];
		double[] cf = new double[m * N];
		for (int multi = 0; multi < m; multi++) {
			for (int i1 = 0; i1 < N; i1++) {
				cp[multi * N + i1] = cpbar * (1 + delta * Math.cos(2 * Math.PI * (i1 + 1) / N + phi));
				cf[multi * N + i1] = cfbar * (1 + delta * Math.cos(2 * Math.PI * (i1 + 1) / N + phi));
			}
		}

		try {
			solvePMBRP(pi0, pi1, cp, cf, N, M, m);
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

	public static void solvePMBRP(double[][][][] pi0, double[][][][] pi1, double[] cp, double[] cf, int N, int M, int m)
			throws IloException {
		// Create the model
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-5);
		// Variables
		IloNumVar[][] x0 = new IloNumVar[m * N][M + 1];
		IloNumVar[][] x1 = new IloNumVar[m * N][M + 1];
		IloNumVar[] y = new IloNumVar[m * N];
		IloNumVar[] t = new IloNumVar[m * N];
		IloNumVar[][] z = new IloNumVar[m * N][M + 1];
		for (int i1 = 0; i1 < m * N; i1++) {
			y[i1] = cplex.boolVar();
			t[i1] = cplex.numVar(1, Double.POSITIVE_INFINITY);
			for (int i2 = 0; i2 <= M; i2++) {
				x1[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
				x0[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
				z[i1][i2] = cplex.boolVar();
			}
		}
		// Objective
		IloNumExpr objExpr = cplex.constant(0);
		for (int i1 = 0; i1 < m * N; i1++) {
			objExpr = cplex.sum(objExpr, cplex.prod(cf[i1], x1[i1][0]));
			for (int i2 = 1; i2 <= M; i2++) {
				objExpr = cplex.sum(objExpr, cplex.prod(cp[i1], x1[i1][i2]));
			}
		}
		cplex.addMinimize(objExpr);
		// Constraints 12b
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				// Subtraction left side
				IloNumExpr sumActions = cplex.constant(0);
				sumActions = cplex.sum(sumActions, x1[i1][i2]);
				sumActions = cplex.sum(sumActions, x0[i1][i2]);
				// Subtraction right side
				IloNumExpr sumTransitions = cplex.constant(0);
				for (int j1 = 0; j1 < m * N; j1++) {
					for (int j2 = 0; j2 <= M; j2++) {
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi1[j1][j2][i1][i2], x1[j1][j2]));
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi0[j1][j2][i1][i2], x0[j1][j2]));
					}
				}
				cplex.addEq(cplex.diff(sumActions, sumTransitions), 0);
			}
		}
		// Constraints 32c (Schouten)
		for (int i1 = 0; i1 < m * N; i1++) {
			cplex.addEq(x0[i1][0], 0);
			cplex.addEq(x0[i1][M], 0);
		}
		
		// Constraints 10c
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.sum(x0[i1][i2], z[i1][i2]), 1);
			}
		}
		// Constraints 10d
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.diff(x1[i1][i2], z[i1][i2]), 0);
			}
		}
		
		// Constraints 12c
		for (int i1 = 0; i1 < m * N; i1++) {
			IloNumExpr sumAges = cplex.constant(0);
			for (int i2 = 0; i2 <= M; i2++) {
				sumAges = cplex.sum(sumAges, x1[i1][i2]);
				sumAges = cplex.sum(sumAges, x0[i1][i2]);
			}
			cplex.addEq(sumAges, 1.0 / (m * N));
		}
		// Constraints 12d
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				cplex.addLe(cplex.diff(z[i1][i2], y[i1]), 0);
			}
		}
		// Constraints 12e
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				for (int j2 = i2 + 1; j2 <= M; j2++) {
					cplex.addLe(cplex.diff(z[i1][i2], z[i1][j2]), 0);
				}
			}
		}
		// Constraints 12f, +1 where j1 and i1 are used to move set to 1..mN instead of 0..mN-1
		for (int j1 = 0; j1 < m * N; j1++) {
			for (int i1 = j1 + 1; i1 < m * N; i1++) {
				IloNumExpr left = cplex.sum(t[i1], cplex.prod(j1 + 1 + m * N, y[j1]));
				cplex.addLe(left, m * N + i1 + 1);
			}
		}
		// Constraints 12g, same observation as 12f
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int j1 = i1 + 1; j1 < m * N; j1++) {
				IloNumExpr left = cplex.sum(t[i1], cplex.prod(j1 + 1, y[j1]));
				cplex.addLe(left, m * N + i1 + 1);
			}
		}
		// Constraints 12h
		for (int i1=0; i1<m*N; i1++) {
			for (int i2=0; i2<=M; i2++) {
				IloNumExpr left = cplex.diff(cplex.prod(M, cplex.diff(y[i1], z[i1][i2])), t[i1]);
				cplex.addLe(left, M-1-i2);
			}
		}
		// Constraints 12i
		for (int i1=0; i1<m*N; i1++) {
			for (int i2=0; i2<=M; i2++) {
				IloNumExpr left = cplex.sum(cplex.prod(M, z[i1][i2]), t[i1]);
				cplex.addLe(left, M+i2);
			}
		}
		// Solve the model
		cplex.solve();
		// Query the solution
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + 12*cplex.getObjValue());
			System.out.print("Months, Ages: ");
			for (int i1 = 0; i1 < m * N; i1++) {
				y_opt[i1] = cplex.getValue(y[i1]);
				t_opt[i1] = cplex.getValue(t[i1]);
				if (y_opt[i1] > 0)
					System.out.print("(" + (i1+1) + ", " + t_opt[i1] + "), ");
			}
		} else {
			System.out.println("No optimal solution found");
		}
		// Close the model
		cplex.close();
	}

}
