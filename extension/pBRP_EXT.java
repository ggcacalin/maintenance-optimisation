package extension;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Scanner;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class pBRP_EXT {

	// 12/36 months in alpha, 1/3 years in m
	private static double alpha = 12;
	private static double beta = 2;
	private static int N = 12;
	private static int M = 2 * (int) alpha;
	private static int m = 1;
	private static double phi = -0.178;
	private static double rho;
	private static double[] y_opt = new double[m * N];
	private static ArrayList<ArrayList<Double>> solutionSet = new ArrayList<>();
	private static ArrayList<Double> solutionEntry = new ArrayList<>();
	private static ArrayList<Integer> months = new ArrayList<>();
	private static ArrayList<Integer> months_mod = new ArrayList<>();

	public static void main(String[] args) throws FileNotFoundException {

		// Transition probability definition
		double[][][][] pi0 = new double[m * N][M + 1][m * N][M + 1];
		double[][][][] pi1 = new double[m * N][M + 1][m * N][M + 1];
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				if (i2 > 0 && i2 < M) {
					pi0[i1][i2][(i1 + 1) % (m * N)][i2 + 1] = 1 - pfail(i2 + 1);
					pi0[i1][i2][(i1 + 1) % (m * N)][0] = pfail(i2 + 1);
				}
				pi1[i1][i2][(i1 + 1) % (m * N)][1] = 1 - pfail(1);
				pi1[i1][i2][(i1 + 1) % (m * N)][0] = pfail(1);
			}
		}

		// Reading energy prices
		double[] allPrices = new double[36];
		double[] price19 = new double[12];
		double[] price20 = new double[12];
		double[] price21 = new double[12];
		Scanner read = new Scanner(new File("energy prices 19-21.txt"));
		int filler = 0;
		while (read.hasNextDouble()) {
			allPrices[filler] = read.nextDouble();
			filler++;
		}

		// Separating prices by year
		for (int i = 0; i < 12; i++) {
			price19[i] = allPrices[i];
			price20[i] = allPrices[i + 12];
			price21[i] = allPrices[i + 24];
		}

		// Cost definition - single years for alpha = 12, all years for alpha = 36
		// Change price array in cp[], cf[] as needed (19/20/21/all)

		double[] cp = new double[(int) alpha];
		double[] cf = new double[(int) alpha];
		double[] P = new double[N];
		double cpbar = 0;
		double cfbar = 0;
		for (int i1 = 0; i1 < N; i1++)
			P[i1] = 114024 + 21480 * Math.cos(2 * Math.PI * (i1 + 1) / N + phi);
		for (int i1 = 0; i1 < alpha; i1++) {
			cp[i1] = 148.2 + 10 * (P[i1 % N] * price21[i1] / 1000);
			cf[i1] = 592.8 + 40 * (P[i1 % N] * price21[i1] / 1000);
			cpbar += cp[i1];
			cfbar += cf[i1];
		}
		cpbar = cpbar / alpha;
		cfbar = cfbar / alpha;
		double[] constantP = new double[(int) alpha];
		double[] constantF = new double[(int) alpha];
		for (int i1 = 0; i1 < alpha; i1++) {
			constantP[i1] = cpbar;
			constantF[i1] = cfbar;
		}
		// Simple solve for different prices

//		try {
//			solvePBRP(pi0, pi1, cp, cf, N, M, m, 1);
//			solvePBRP(pi0, pi1, constantP, constantF, N, M, m, 1);
//		} catch (IloException e) {
//			System.out.println("A Cplex exception occured: " + e.getMessage());
//			e.printStackTrace();
//		}

		// Solving and output printing CM limit
		/*
		 * try { File cmres = new File("CMresults.txt"); cmres.createNewFile();
		 * FileWriter fileWriter = new FileWriter(cmres, true); BufferedWriter buffres =
		 * new BufferedWriter(fileWriter); buffres.write("pBRP CM results: \n"); int
		 * counter = 0; for (double L = 0.025; L >= 0.001; L -= 0.003) {
		 * 
		 * // Model solved - time-varying first, constant after try { solvePBRP(pi0,
		 * pi1, cp, cf, N, M, m, L); buffres.write(solutions.get(counter) + ",");
		 * counter++; solvePBRP(pi0, pi1, constantP, constantF, N, M, m, L);
		 * buffres.write(solutions.get(counter) + " \n"); counter++; } catch
		 * (IloException e) { System.out.println("A Cplex exception occured: " +
		 * e.getMessage()); e.printStackTrace(); } } buffres.close(); } catch
		 * (IOException e) { System.out.println("Printing went awry..."); }
		 */

		// Solving imperfect PM model
		try {
			File ipmres = new File("iPMresults.txt");
			File ipmage = new File("iPMages.txt");
			ipmres.createNewFile();
			ipmage.createNewFile();
			FileWriter resfileWriter = new FileWriter(ipmres, true);
			FileWriter agefileWriter = new FileWriter(ipmage, true);
			BufferedWriter buffres = new BufferedWriter(resfileWriter);
			BufferedWriter buffage = new BufferedWriter(agefileWriter);
			buffres.write("pBRP iPM results: \n");
			buffage.write("pBRP iPM results: \n");
			try {
				double[][][][] pi2 = new double[m * N][M + 1][m * N][M + 1];
				for (rho = 0.1; rho < 1; rho += 0.2) {
					// Transition probability definition
					pi2 = new double[m * N][M + 1][m * N][M + 1];
					for (int i1 = 0; i1 < m * N; i1++) {
						for (int i2 = 1; i2 < M; i2++) {
							// Sequence after iPM: immediately go to i2-i2*rho to lose i2*rho years, then
							// get to next age if successful repair or 0 if not
							pi2[i1][i2][(i1 + 1) % (m * N)][i2 + 1 - (int) Math.round(i2 * rho)] = 1
									- pfail(i2 + 1 - (int) Math.round(i2 * rho));
							pi2[i1][i2][(i1 + 1) % (m * N)][0] = pfail(i2 + 1 - (int) Math.round(i2 * rho));
						}
					}
					solvePBRP_imp(pi0, pi1, pi2, cp, cf, N, M, m, 1);
					buffage.write(months.toString() + "\n");
					buffage.write(months_mod.toString() + "\n");
					months.clear();
					months_mod.clear();
				}
				solutionSet.add(new ArrayList<Double>(solutionEntry));
				solutionEntry.clear();
				for (rho = 0.1; rho < 1; rho += 0.2) {
					// Transition probability definition
					pi2 = new double[m * N][M + 1][m * N][M + 1];
					for (int i1 = 0; i1 < m * N; i1++) {
						for (int i2 = 1; i2 < M; i2++) {
							// Sequence after iPM: immediately go to i2-i2*rho to lose i2*rho years, then
							// get to next age if successful repair or 0 if not
							pi2[i1][i2][(i1 + 1) % (m * N)][i2 + 1 - (int) Math.round(i2 * rho)] = 1
									- pfail(i2 + 1 - (int) Math.round(i2 * rho));
							pi2[i1][i2][(i1 + 1) % (m * N)][0] = pfail(i2 + 1 - (int) Math.round(i2 * rho));
						}
					}
					solvePBRP_imp(pi0, pi1, pi2, constantP, constantF, N, M, m, 1);
					buffage.write(months.toString() + "\n");
					buffage.write(months_mod.toString() + "\n");
					months.clear();
					months_mod.clear();
				}
				solutionSet.add(new ArrayList<Double>(solutionEntry));
				solutionEntry.clear();
				System.out.println("Solutions:");
				for (ArrayList<Double> set : solutionSet) {
					System.out.println(set.toString());
					buffres.write(set.toString() + "\n");
				}

			} catch (IloException e) {
				System.out.println("A Cplex exception occured: " + e.getMessage());
				e.printStackTrace();
			}
			buffres.close();
			buffage.close();
		} catch (IOException e) {
			System.out.println("File writing error: " + e.getMessage());
			e.printStackTrace();
		}

	}

	public static double wcdf(double x) {
		double exponent = -Math.pow(x / alpha, beta);
		return 1 - Math.pow(Math.E, exponent);
	}

	public static double pfail(int i2) {
		double numerator = wcdf(i2) - wcdf(i2 - 1);
		return numerator / (1 - wcdf(i2 - 1));
	}

	public static void solvePBRP(double[][][][] pi0, double[][][][] pi1, double[] cp, double[] cf, int N, int M, int m,
			double L) throws IloException {
		// Create the model
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-5);
		// Variables
		IloNumVar[][] x0 = new IloNumVar[m * N][M + 1];
		IloNumVar[][] x1 = new IloNumVar[m * N][M + 1];
		IloNumVar[] y = new IloNumVar[m * N];
		for (int i1 = 0; i1 < m * N; i1++) {
			y[i1] = cplex.boolVar();
			for (int i2 = 0; i2 <= M; i2++) {
				x1[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
				x0[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
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
		// Constraints 10b
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
		// Constraints 22c (Schouten)
		for (int i1 = 0; i1 < m * N; i1++) {
			cplex.addEq(x0[i1][0], 0);
			cplex.addEq(x0[i1][M], 0);
		}
		// Constraints 10c
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.sum(x0[i1][i2], y[i1]), 1);
			}
		}
		// Constraints 10d
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.diff(x1[i1][i2], y[i1]), 0);
			}
		}
		// Constraints 10e
		for (int i1 = 0; i1 < m * N; i1++) {
			IloNumExpr sumAges = cplex.constant(0);
			for (int i2 = 0; i2 <= M; i2++) {
				sumAges = cplex.sum(sumAges, x1[i1][i2]);
				sumAges = cplex.sum(sumAges, x0[i1][i2]);
			}
			cplex.addEq(sumAges, 1.0 / (m * N));
		}
		// Extension - CM limit
		IloNumExpr sumCM = cplex.constant(0);
		for (int i1 = 0; i1 < m * N; i1++)
			sumCM = cplex.sum(sumCM, x1[i1][0]);
		cplex.addLe(sumCM, L);
		// Solve the model
		cplex.solve();
		// Query the solution
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + 12 * cplex.getObjValue());
			solutionEntry.add(12 * cplex.getObjValue());
			System.out.print("Months: ");
			for (int i1 = 0; i1 < m * N; i1++) {
				y_opt[i1] = cplex.getValue(y[i1]);
				if (y_opt[i1] > 0)
					System.out.print(i1 + 1 + " ");
			}
			System.out.println();

		} else {
			System.out.println("No optimal solution found");
		}
		// Close the model
		cplex.close();
	}

	public static void solvePBRP_imp(double[][][][] pi0, double[][][][] pi1, double[][][][] pi2, double[] cp,
			double[] cf, int N, int M, int m, double L) throws IloException {
		// Create the model
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-5);
		// Variables
		IloNumVar[][] x0 = new IloNumVar[m * N][M + 1];
		IloNumVar[][] x1 = new IloNumVar[m * N][M + 1];
		IloNumVar[][] x2 = new IloNumVar[m * N][M + 1];
		IloNumVar[] y = new IloNumVar[m * N];
		IloNumVar[] v = new IloNumVar[m * N];
		for (int i1 = 0; i1 < m * N; i1++) {
			y[i1] = cplex.boolVar();
			v[i1] = cplex.boolVar();
			for (int i2 = 0; i2 <= M; i2++) {
				x1[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
				x0[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
				x2[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
			}
		}
		// Objective
		IloNumExpr objExpr = cplex.constant(0);
		for (int i1 = 0; i1 < m * N; i1++) {
			objExpr = cplex.sum(objExpr, cplex.prod(cf[i1], x1[i1][0]));
			for (int i2 = 1; i2 <= M; i2++) {
				objExpr = cplex.sum(objExpr, cplex.prod(cp[i1], x1[i1][i2]));
				objExpr = cplex.sum(objExpr, cplex.prod(rho * cp[i1], x2[i1][i2]));
			}
		}
		cplex.addMinimize(objExpr);
		// Constraints 10b
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				// Subtraction left side
				IloNumExpr sumActions = cplex.constant(0);
				sumActions = cplex.sum(sumActions, x1[i1][i2]);
				sumActions = cplex.sum(sumActions, x0[i1][i2]);
				sumActions = cplex.sum(sumActions, x2[i1][i2]);
				// Subtraction right side
				IloNumExpr sumTransitions = cplex.constant(0);
				for (int j1 = 0; j1 < m * N; j1++) {
					for (int j2 = 0; j2 <= M; j2++) {
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi1[j1][j2][i1][i2], x1[j1][j2]));
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi0[j1][j2][i1][i2], x0[j1][j2]));
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi2[j1][j2][i1][i2], x2[j1][j2]));
					}
				}
				cplex.addEq(cplex.diff(sumActions, sumTransitions), 0);
			}
		}
		// Constraints 22c (Schouten)
		for (int i1 = 0; i1 < m * N; i1++) {
			cplex.addEq(x0[i1][0], 0);
			cplex.addEq(x0[i1][M], 0);
			cplex.addEq(x2[i1][0], 0);
			cplex.addEq(x2[i1][M], 0);
		}
		// Constraints 10c
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.sum(x0[i1][i2], y[i1]), 1);
			}
		}
		// Constraints 10d
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.diff(x1[i1][i2], y[i1]), 0);
			}
		}
		// Extension - no iPM allowed in periods when PM performed
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.sum(x2[i1][i2], y[i1]), 1);
			}
		}
		// Constraints 10c for iPM
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.sum(x0[i1][i2], v[i1]), 1);
			}
		}
		// Constraints 10d for iPM
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.sum(x1[i1][i2], v[i1]), 1);
			}
		}
		// Extension - no PM allowed in periods when iPM performed
		for (int i1 = 0; i1 < m * N; i1++) {
			for (int i2 = 1; i2 <= M; i2++) {
				cplex.addLe(cplex.diff(x2[i1][i2], v[i1]), 0);
			}
		}
		// Constraints 10e
		for (int i1 = 0; i1 < m * N; i1++) {
			IloNumExpr sumAges = cplex.constant(0);
			for (int i2 = 0; i2 <= M; i2++) {
				sumAges = cplex.sum(sumAges, x1[i1][i2]);
				sumAges = cplex.sum(sumAges, x0[i1][i2]);
				sumAges = cplex.sum(sumAges, x2[i1][i2]);
			}
			cplex.addEq(sumAges, 1.0 / (m * N));
		}
		// Extension - CM limit
		IloNumExpr sumCM = cplex.constant(0);
		for (int i1 = 0; i1 < m * N; i1++)
			sumCM = cplex.sum(sumCM, x1[i1][0]);
		cplex.addLe(sumCM, L);
		// Solve the model
		cplex.solve();
		// Query the solution
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + 12 * cplex.getObjValue());
			solutionEntry.add(12 * cplex.getObjValue());
			System.out.print("Months PM: ");
			for (int i1 = 0; i1 < m * N; i1++) {
				if (cplex.getValue(y[i1]) > 0) {
					System.out.print(i1 + 1 + " ");
					months.add(i1 + 1);
				}
			}
			System.out.println();
			System.out.print("Months imperfect PM: ");
			for (int i1 = 0; i1 < m * N; i1++) {
				if (cplex.getValue(v[i1]) > 0) {
					System.out.print(i1 + 1 + " ");
					months_mod.add(i1 + 1);
				}
			}
			System.out.println();

		} else {
			System.out.println("No optimal solution found");
		}
		// Close the model
		cplex.close();
	}

}
