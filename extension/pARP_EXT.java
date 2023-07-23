package extension;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.math.RoundingMode;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Scanner;

import ilog.concert.IloException;
import ilog.concert.IloNumExpr;
import ilog.concert.IloNumVar;
import ilog.cplex.IloCplex;

public class pARP_EXT {

	// 12/36 months, 1/3 years. Any difference?
	private static double alpha = 12;
	private static double beta = 2;
	private static int N = 12;
	private static int M = 2 * (int) alpha;
	private static int gamma = 4;
	private static int H = 6;
	private static double rho;
	private static double[][] x0_opt = new double[N][M + 1];
	private static double[][] x1_opt = new double[N][M + 1];
	private static ArrayList<Double> solutions = new ArrayList<>();
	private static ArrayList<Integer> ages = new ArrayList<>();
	private static ArrayList<Double> solutions_mod = new ArrayList<>();
	private static ArrayList<Integer> ages_mod = new ArrayList<>();

	public static void main(String[] args) throws FileNotFoundException, IOException {

		//Handy decimal limiter
		DecimalFormat df = new DecimalFormat("#.###");
		df.setRoundingMode(RoundingMode.CEILING);
		
		// Transition probability definition
		double[][][][] pi0 = new double[N][M + 1][N][M + 1];
		double[][][][] pi1 = new double[N][M + 1][N][M + 1];
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				if (i2 > 0 && i2 < M) {
					pi0[i1][i2][(i1 + 1) % N][i2 + 1] = 1 - pfail(i2 + 1);
					pi0[i1][i2][(i1 + 1) % N][0] = pfail(i2 + 1);
				}
				pi1[i1][i2][(i1 + 1) % N][1] = 1 - pfail(1);
				pi1[i1][i2][(i1 + 1) % N][0] = pfail(1);
			}
		}

		// Modified transition probabilities for limited gearbox availability
		double[][][][][][] pi0_mod = new double[N][M + 1][H + 1][N][M + 1][H + 1];
		double[][][][][][] pi1_mod = new double[N][M + 1][H + 1][N][M + 1][H + 1];
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				for (int h = 0; h <= H; h++) {
					if (i2 > 0 && i2 < M) {
						pi0_mod[i1][i2][h][(i1 + 1) % N][i2 + 1][Math.min(h + indicator_G(i1 + 1, gamma), H)] = 1
								- pfail(i2 + 1);
						pi0_mod[i1][i2][h][(i1 + 1) % N][0][Math.min(h + indicator_G(i1 + 1, gamma), H)] = pfail(
								i2 + 1);
					}
					if (h > 0) {
						pi1_mod[i1][i2][h][(i1 + 1) % N][1][h + indicator_G(i1 + 1, gamma) - 1] = 1 - pfail(1);
						pi1_mod[i1][i2][h][(i1 + 1) % N][0][h + indicator_G(i1 + 1, gamma) - 1] = pfail(1);
					}
				}
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
		double phi = -0.178;
		double[] cp = new double[N];
		double[] cf = new double[N];
		double[] cd = new double[N];
		double[] P = new double[N];
		double cpbar = 0;
		double cfbar = 0;
		double cdbar = 0;
		for (int i1 = 0; i1 < N; i1++) {
			P[i1] = 114024 + 21480 * Math.cos(2 * Math.PI * (i1 + 1) / N + phi);
			cd[i1] = 30 * (P[i1] * (price19[i1] + price20[i1] + price21[i1]) / 3000);
			cp[i1] = 148.2 + 10 * (P[i1] * (price19[i1] + price20[i1] + price21[i1]) / 3000);
			cf[i1] = 592.8 + 40 * (P[i1] * (price19[i1] + price20[i1] + price21[i1]) / 3000);
			cpbar += cp[i1];
			cfbar += cf[i1];
			cdbar += cd[i1];
		}
		cpbar = cpbar / N;
		cfbar = cfbar / N;
		cdbar = cdbar / N;
		double[] constantP = new double[N];
		double[] constantF = new double[N];
		double[] constantD = new double[N];
		for (int i1 = 0; i1 < N; i1++) {
			constantP[i1] = cpbar;
			constantF[i1] = cfbar;
			constantD[i1] = cdbar;
		}
		// Simple solve for different prices

//		try {
//			solvePARP(pi0, pi1, cp, cf, N, M, 1);
//			System.out.println(ages.toString());
//			ages.clear();
//			solvePARP(pi0, pi1, constantP, constantF, N, M, 1);
//			System.out.println(ages.toString());
//			ages.clear();
//			System.out.println(solutions.toString());
//		} catch (IloException e) {
//			System.out.println("A Cplex exception occured: " + e.getMessage());
//			e.printStackTrace();
//		}

		// Solving and output printing CM limit

//		try {
//			File cmres = new File("CMresults.txt");
//			File cmage = new File("CMagesARP.txt");
//			cmres.createNewFile();
//			cmage.createNewFile();
//			FileWriter resfileWriter = new FileWriter(cmres, true);
//			FileWriter agefileWriter = new FileWriter(cmage, false);
//			BufferedWriter buffres = new BufferedWriter(resfileWriter);
//			BufferedWriter buffage = new BufferedWriter(agefileWriter);
//			buffres.write("pARP CM results: \n");
//			int counter = 0;
//			for (double L = 0.025; L >= 0.001; L -= 0.003) {
//
//				// Model solved - time-varying first, constant after
//				try {
//					buffage.write("L=" + L + "\n");
//					solvePARP(pi0, pi1, cp, cf, N, M, L);
//					buffres.write(solutions.get(counter) + ",");
//					buffage.write(ages.toString());
//					ages.clear();
//					buffage.write("\n");
//					counter++;
//					solvePARP(pi0, pi1, constantP, constantF, N, M, L);
//					buffres.write(solutions.get(counter) + " \n");
//					buffage.write(ages.toString());
//					ages.clear();
//					buffage.write("\n");
//					counter++;
//				} catch (IloException e) {
//					System.out.println("A Cplex exception occured: " + e.getMessage());
//					e.printStackTrace();
//				}
//			}
//			buffres.close();
//			buffage.close();
//		} catch (IOException e) {
//			System.out.println("Printing went awry...");
//		}

		// Solving limited gearbox model

//		try {
//			solvePARP_mod(pi0_mod, pi1_mod, cp, cf, cd, N, M, H, 1);
//			System.out.println(ages_mod.toString());
//			ages_mod.clear(); //
//			solvePARP_mod(pi0_mod, pi1_mod, constantP, constantF, constantD, N, M, H, 1);
//			System.out.println(ages_mod.toString());
//		} catch (IloException e) {
//			System.out.println("A Cplex exception occured: " + e.getMessage());
//			e.printStackTrace();
//		}
//		System.out.println(solutions_mod.toString());

		// Solving imperfect PM model

//		try {
//			File ipmres = new File("iPMresults.txt");
//			File ipmage = new File("iPMages.txt");
//			ipmres.createNewFile();
//			ipmage.createNewFile();
//			FileWriter resfileWriter = new FileWriter(ipmres, true);
//			FileWriter agefileWriter = new FileWriter(ipmage, true);
//			BufferedWriter buffres = new BufferedWriter(resfileWriter);
//			BufferedWriter buffage = new BufferedWriter(agefileWriter);
//			buffres.write("pARP iPM results: \n");
//			buffage.write("pARP iPM varying results: \n");
//			try {
//				double[][][][] pi2 = new double[N][M + 1][N][M + 1];
//				for (rho = 0.1; rho < 1; rho += 0.2) {
//					for (int i1 = 0; i1 < N; i1++) {
//						for (int i2 = 1; i2 < M; i2++) {
//							// Sequence after iPM: immediately go to i2-i2*rho to lose i2*rho years, then
//							// get to next age if successful repair or 0 if not
//							pi2[i1][i2][(i1 + 1) % N][i2 + 1 - (int) Math.round(i2 * rho)] = 1
//									- pfail(i2 + 1 - (int) Math.round(i2 * rho));
//							pi2[i1][i2][(i1 + 1) % N][0] = pfail(i2 + 1 - (int) Math.round(i2 * rho));
//						}
//					}
//					solvePARP_imp(pi0, pi1, pi2, cp, cf, N, M, 1);
//					buffage.write("rho=" + rho);
//					buffage.write(ages.toString() + "\n");
//					buffage.write(ages_mod.toString() + "\n");
//					System.out.println(ages.toString());
//					System.out.println(ages_mod.toString());
//					ages.clear();
//					ages_mod.clear();
//				}
//				ArrayList<Double> varyingSol = new ArrayList<>(solutions);
//				solutions.clear();
//				buffage.write("pARP iPM constant results: \n");
//				for (rho = 0.1; rho < 1; rho += 0.2) {
//					for (int i1 = 0; i1 < N; i1++) {
//						for (int i2 = 1; i2 < M; i2++) {
//							// Sequence after iPM: immediately go to i2-i2*rho to lose i2*rho years, then
//							// get to next age if successful repair or 0 if not
//							pi2[i1][i2][(i1 + 1) % N][i2 + 1 - (int) Math.round(i2 * rho)] = 1
//									- pfail(i2 + 1 - (int) Math.round(i2 * rho));
//							pi2[i1][i2][(i1 + 1) % N][0] = pfail(i2 + 1 - (int) Math.round(i2 * rho));
//						}
//					}
//					solvePARP_imp(pi0, pi1, pi2, constantP, constantF, N, M, 0.06);
//					buffage.write("rho=" + rho + "\n");
//					buffage.write(ages.toString() + "\n");
//					buffage.write(ages_mod.toString() + "\n");
//					System.out.println(ages.toString());
//					System.out.println(ages_mod.toString());
//					ages.clear();
//					ages_mod.clear();
//				}
//				buffres.write(varyingSol.toString() + "\n");
//				buffres.write(solutions.toString() + "\n");
//				System.out.println(varyingSol.toString());
//				System.out.println(solutions.toString());
//			} catch (IloException e) {
//				System.out.println("A Cplex exception occured: " + e.getMessage());
//				e.printStackTrace();
//			}
//			buffres.close();
//			buffage.close();
//		} catch (IOException e) {
//			System.out.println("File writing error: " + e.getMessage());
//			e.printStackTrace();
//		}

		// Solving the all-encompasing extension

		try {
//			File bigres = new File("BIGresults.txt");
//			bigres.createNewFile();
//			FileWriter resfileWriter = new FileWriter(bigres, true);
//			BufferedWriter buffres = new BufferedWriter(resfileWriter);
//			File bigage = new File("BIGages.txt");
//			bigage.createNewFile();
//			FileWriter agefileWriter = new FileWriter(bigage, true);
//			BufferedWriter buffage = new BufferedWriter(agefileWriter);
//			buffres.write("L=0.03, varying \n");
//			buffage.write("L=0.03, varying \n");
			for (rho = 0.1; rho < 1; rho += 0.2) {
				for (gamma = 2; gamma < 7; gamma++) {
						solvePARP_ultimate(pi0_mod, pi1_mod, cp, cf, cd, N, M, H, 0.06, gamma, rho);
//						buffres.write(solutions.get(0) + "\n");
//						solutions.remove(0);
//						buffage.write(ages.toString() + "\n");
//						buffage.write(ages_mod.toString() + "\n");
//						buffage.write("-------------------------------------------------------------- \n");
						ages.clear();
						ages_mod.clear();
				}

			}
			ArrayList<Double> varyingSol = new ArrayList<>(solutions);
			solutions.clear();
//			buffres.write("L=0.03, constant \n");
//			buffage.write("L=0.03, constant \n");
			for (rho = 0.1; rho < 1; rho += 0.2) {
				for (gamma = 2; gamma < 7; gamma++) {
						solvePARP_ultimate(pi0_mod, pi1_mod, constantP, constantF, constantD, N, M, H, 0.06, gamma, rho);
//						buffres.write(solutions.get(0) + "\n");
//						solutions.remove(0);
//						buffage.write(ages.toString() + "\n");
//						buffage.write(ages_mod.toString() + "\n");
//						buffage.write("-------------------------------------------------------------- \n");
						ages.clear();
						ages_mod.clear();
				}

			}
			int count = 0;
			for (rho = 0.1; rho < 1; rho += 0.2) {
				for (gamma = 2; gamma < 7; gamma++) {
					System.out.print(df.format(varyingSol.get(count)) + "(" + df.format((solutions.get(count)/varyingSol.get(count)-1)*100) + ") ");
					count++;
				}
				System.out.println();
			}
			//System.out.println(varyingSol.toString());
			//System.out.println(solutions.toString());
//			buffres.close();
//			buffage.close();
		} catch (IloException e) {
			System.out.println("A Cplex exception occured: " + e.getMessage());
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

	public static int indicator_G(int period, int gamma) {
		if (period % gamma == 1 || gamma == 1)
			return 1;
		return 0;
	}

	public static void solvePARP(double[][][][] pi0, double[][][][] pi1, double[] cp, double[] cf, int N, int M,
			double L) throws IloException {
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
		// Extension - CM limit
		IloNumExpr sumCM = cplex.constant(0);
		for (int i1 = 0; i1 < N; i1++)
			sumCM = cplex.sum(sumCM, x1[i1][0]);
		cplex.addLe(sumCM, L);

		// Solve the model
		cplex.solve();
		// Query the solution
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + 12 * cplex.getObjValue());
			solutions.add(12 * cplex.getObjValue());
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
			for (int i1 = 0; i1 < N; i1++) {
				for (int i2 = 0; i2 <= M; i2++) {
					x0_opt[i1][i2] = cplex.getValue(x0[i1][i2]);
				}
			}

		} else {
			System.out.println("No optimal solution found");
		}
		// Close the model
		cplex.close();
	}

	public static void solvePARP_mod(double[][][][][][] pi0, double[][][][][][] pi1, double[] cp, double[] cf,
			double[] cd, int N, int M, int H, double L) throws IloException {
		// Create the model
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-5);
		// Variables
		IloNumVar[][][] x0 = new IloNumVar[N][M + 1][H + 1];
		IloNumVar[][][] x1 = new IloNumVar[N][M + 1][H + 1];
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				for (int h = 0; h <= H; h++) {
					x1[i1][i2][h] = cplex.numVar(0, 1);
					x0[i1][i2][h] = cplex.numVar(0, 1);
				}
			}
		}
		// Objective
		IloNumExpr objExpr = cplex.constant(0);
		for (int h = 1; h <= H; h++) {
			for (int i1 = 0; i1 < N; i1++) {
				objExpr = cplex.sum(objExpr, cplex.prod(cf[i1], x1[i1][0][h]));
				for (int i2 = 1; i2 <= M; i2++) {
					objExpr = cplex.sum(objExpr, cplex.prod(cp[i1], x1[i1][i2][h]));
				}
			}
		}
		// Adding downtime cost when system broken
		for (int i1 = 0; i1 < N; i1++)
			objExpr = cplex.sum(objExpr, cplex.prod(cd[i1], x0[i1][0][0]));
		cplex.addMinimize(objExpr);
		// Constraints 7b
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				for (int h = 0; h <= H; h++) {
					// Subtraction left side
					IloNumExpr sumActions = cplex.constant(0);
					sumActions = cplex.sum(sumActions, x1[i1][i2][h]);
					sumActions = cplex.sum(sumActions, x0[i1][i2][h]);
					// Subtraction right side
					IloNumExpr sumTransitions = cplex.constant(0);
					for (int j1 = 0; j1 < N; j1++) {
						for (int j2 = 0; j2 <= M; j2++) {
							for (int k = 0; k <= H; k++) {
								sumTransitions = cplex.sum(sumTransitions,
										cplex.prod(pi1[j1][j2][k][i1][i2][h], x1[j1][j2][k]));
								sumTransitions = cplex.sum(sumTransitions,
										cplex.prod(pi0[j1][j2][k][i1][i2][h], x0[j1][j2][k]));
							}
						}
					}
					cplex.addEq(cplex.diff(sumActions, sumTransitions), 0);
				}
			}
		}
		// Constraints 16c (Schouten) only apply when at least 1 gearbox available
		for (int i1 = 0; i1 < N; i1++) {
			for (int h = 1; h <= H; h++) {
				cplex.addEq(x0[i1][0][h], 0);
				cplex.addEq(x0[i1][M][h], 0);
			}
		}
		// Constraints 7c
		for (int i1 = 0; i1 < N; i1++) {
			IloNumExpr sumAges = cplex.constant(0);
			for (int i2 = 0; i2 <= M; i2++) {
				for (int h = 0; h <= H; h++) {
					sumAges = cplex.sum(sumAges, x1[i1][i2][h]);
					sumAges = cplex.sum(sumAges, x0[i1][i2][h]);
				}
			}
			cplex.addEq(sumAges, 1.0 / N);
		}
		// Extension - CM limit
		IloNumExpr sumCM = cplex.constant(0);
		for (int i1 = 0; i1 < N; i1++)
			for (int h = 0; h <= H; h++)
				sumCM = cplex.sum(sumCM, x1[i1][0][h]);
		cplex.addLe(sumCM, L);

		// Extension - cannot maintain if no gearboxes
		for (int i1 = 0; i1 < N; i1++)
			for (int i2 = 0; i2 <= M; i2++)
				cplex.addEq(x1[i1][i2][0], 0);

		// Solve the model
		cplex.solve();
		// Query the solution
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + 12 * cplex.getObjValue());
			solutions_mod.add(12 * cplex.getObjValue());
			for (int i1 = 0; i1 < N; i1++) {
				int thisAge = M;
				for (int i2 = 0; i2 <= M; i2++) {
					for (int h = 1; h <= H; h++) {
						if (i2 > 0 && i2 < thisAge && cplex.getValue(x1[i1][i2][h]) > 0) {
							ages_mod.add(i2);
							thisAge = i2;
							break;
						}
					}
				}
				if (thisAge == M)
					ages_mod.add(M);

			}

		} else {
			System.out.println("No optimal solution found");
			System.out.println(cplex.getStatus());
		}
		// Close the model
		cplex.close();
	}

	public static void solvePARP_imp(double[][][][] pi0, double[][][][] pi1, double[][][][] pi2, double[] cp,
			double[] cf, int N, int M, double L) throws IloException {
		// Create the model
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-5);
		// Variables
		IloNumVar[][] x0 = new IloNumVar[N][M + 1];
		IloNumVar[][] x1 = new IloNumVar[N][M + 1];
		IloNumVar[][] x2 = new IloNumVar[N][M + 1];
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				x1[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
				x0[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
				x2[i1][i2] = cplex.numVar(0, Double.POSITIVE_INFINITY);
			}
		}
		// Objective
		IloNumExpr objExpr = cplex.constant(0);
		for (int i1 = 0; i1 < N; i1++) {
			objExpr = cplex.sum(objExpr, cplex.prod(cf[i1], x1[i1][0]));
			for (int i2 = 1; i2 <= M; i2++) {
				objExpr = cplex.sum(objExpr, cplex.prod(cp[i1], x1[i1][i2]));
				objExpr = cplex.sum(objExpr, cplex.prod(rho * cp[i1], x2[i1][i2]));
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
				sumActions = cplex.sum(sumActions, x2[i1][i2]);
				// Subtraction right side
				IloNumExpr sumTransitions = cplex.constant(0);
				for (int j1 = 0; j1 < N; j1++) {
					for (int j2 = 0; j2 <= M; j2++) {
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi1[j1][j2][i1][i2], x1[j1][j2]));
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi0[j1][j2][i1][i2], x0[j1][j2]));
						sumTransitions = cplex.sum(sumTransitions, cplex.prod(pi2[j1][j2][i1][i2], x2[j1][j2]));
					}
				}
				cplex.addEq(cplex.diff(sumActions, sumTransitions), 0);
			}
		}
		// Constraints 16c (Schouten)
		for (int i1 = 0; i1 < N; i1++) {
			cplex.addEq(x0[i1][0], 0);
			cplex.addEq(x0[i1][M], 0);
			cplex.addEq(x2[i1][0], 0);
			cplex.addEq(x2[i1][M], 0);
		}
		// Constraints 7c
		for (int i1 = 0; i1 < N; i1++) {
			IloNumExpr sumAges = cplex.constant(0);
			for (int i2 = 0; i2 <= M; i2++) {
				sumAges = cplex.sum(sumAges, x1[i1][i2]);
				sumAges = cplex.sum(sumAges, x0[i1][i2]);
				sumAges = cplex.sum(sumAges, x2[i1][i2]);
			}
			cplex.addEq(sumAges, 1.0 / N);
		}
		// Extension - CM limit
		IloNumExpr sumCM = cplex.constant(0);
		for (int i1 = 0; i1 < N; i1++)
			sumCM = cplex.sum(sumCM, x1[i1][0]);
		cplex.addLe(sumCM, L);

		// Solve the model
		cplex.solve();
		// Query the solution
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + 12 * cplex.getObjValue());
			solutions.add(12 * cplex.getObjValue());
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
			for (int i1 = 0; i1 < N; i1++) {
				int thisAge = M;
				for (int i2 = 0; i2 <= M; i2++) {
					if (i2 > 0 && i2 < thisAge && cplex.getValue(x2[i1][i2]) > 0) {
						ages_mod.add(i2);
						thisAge = i2;
					}
				}
				if (thisAge == M)
					ages_mod.add(M);
			}
			for (int i1 = 0; i1 < N; i1++) {
				for (int i2 = 0; i2 <= M; i2++) {
					x0_opt[i1][i2] = cplex.getValue(x0[i1][i2]);
				}
			}

		} else {
			System.out.println("No optimal solution found");
		}
		// Close the model
		cplex.close();
	}

	public static void solvePARP_ultimate(double[][][][][][] pi0, double[][][][][][] pi1, double[] cp, double[] cf,
			double[] cd, int N, int M, int H, double L, int gamma, double rho) throws IloException {
		// Inhouse transition probabilities for iPM
		double[][][][][][] pi2_mod = new double[N][M + 1][H + 1][N][M + 1][H + 1];
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 1; i2 < M; i2++) {
				for (int h = 0; h <= H; h++) {
					pi2_mod[i1][i2][h][(i1 + 1) % N][i2 + 1 - (int) Math.round(i2 * rho)][Math
							.min(h + indicator_G(i1 + 1, gamma), H)] = 1 - pfail(i2 + 1 - (int) Math.round(i2 * rho));
					pi2_mod[i1][i2][h][(i1 + 1) % N][0][Math.min(h + indicator_G(i1 + 1, gamma), H)] = pfail(
							i2 + 1 - (int) Math.round(i2 * rho));
				}
			}
		}
		// Create the model
		IloCplex cplex = new IloCplex();
		cplex.setParam(IloCplex.Param.MIP.Tolerances.MIPGap, 1.0e-5);
		// Variables
		IloNumVar[][][] x0 = new IloNumVar[N][M + 1][H + 1];
		IloNumVar[][][] x1 = new IloNumVar[N][M + 1][H + 1];
		IloNumVar[][][] x2 = new IloNumVar[N][M + 1][H + 1];
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				for (int h = 0; h <= H; h++) {
					x1[i1][i2][h] = cplex.numVar(0, 1);
					x0[i1][i2][h] = cplex.numVar(0, 1);
					x2[i1][i2][h] = cplex.numVar(0, 1);
				}
			}
		}
		// Objective
		IloNumExpr objExpr = cplex.constant(0);
		for (int h = 1; h <= H; h++) {
			for (int i1 = 0; i1 < N; i1++) {
				objExpr = cplex.sum(objExpr, cplex.prod(cf[i1], x1[i1][0][h]));
				for (int i2 = 1; i2 <= M; i2++) {
					objExpr = cplex.sum(objExpr, cplex.prod(cp[i1], x1[i1][i2][h]));
					objExpr = cplex.sum(objExpr, cplex.prod(rho * cp[i1], x2[i1][i2][h]));
				}
			}
		}
		// Adding downtime cost when system broken
		for (int i1 = 0; i1 < N; i1++)
			objExpr = cplex.sum(objExpr, cplex.prod(cd[i1], x0[i1][0][0]));
		cplex.addMinimize(objExpr);
		// Constraints 7b
		for (int i1 = 0; i1 < N; i1++) {
			for (int i2 = 0; i2 <= M; i2++) {
				for (int h = 0; h <= H; h++) {
					// Subtraction left side
					IloNumExpr sumActions = cplex.constant(0);
					sumActions = cplex.sum(sumActions, x1[i1][i2][h]);
					sumActions = cplex.sum(sumActions, x0[i1][i2][h]);
					sumActions = cplex.sum(sumActions, x2[i1][i2][h]);
					// Subtraction right side
					IloNumExpr sumTransitions = cplex.constant(0);
					for (int j1 = 0; j1 < N; j1++) {
						for (int j2 = 0; j2 <= M; j2++) {
							for (int k = 0; k <= H; k++) {
								sumTransitions = cplex.sum(sumTransitions,
										cplex.prod(pi1[j1][j2][k][i1][i2][h], x1[j1][j2][k]));
								sumTransitions = cplex.sum(sumTransitions,
										cplex.prod(pi0[j1][j2][k][i1][i2][h], x0[j1][j2][k]));
								sumTransitions = cplex.sum(sumTransitions,
										cplex.prod(pi2_mod[j1][j2][k][i1][i2][h], x2[j1][j2][k]));
							}
						}
					}
					cplex.addEq(cplex.diff(sumActions, sumTransitions), 0);
				}
			}
		}
		// Constraints 16c (Schouten) only apply when at least 1 gearbox available
		for (int i1 = 0; i1 < N; i1++) {
			for (int h = 1; h <= H; h++) {
				cplex.addEq(x0[i1][0][h], 0);
				cplex.addEq(x0[i1][M][h], 0);
				cplex.addEq(x2[i1][0][h], 0);
				cplex.addEq(x2[i1][M][h], 0);
			}
		}
		// Constraints 7c
		for (int i1 = 0; i1 < N; i1++) {
			IloNumExpr sumAges = cplex.constant(0);
			for (int i2 = 0; i2 <= M; i2++) {
				for (int h = 0; h <= H; h++) {
					sumAges = cplex.sum(sumAges, x1[i1][i2][h]);
					sumAges = cplex.sum(sumAges, x0[i1][i2][h]);
					sumAges = cplex.sum(sumAges, x2[i1][i2][h]);
				}
			}
			cplex.addEq(sumAges, 1.0 / N);
		}
		// Extension - CM limit
		IloNumExpr sumCM = cplex.constant(0);
		for (int i1 = 0; i1 < N; i1++)
			for (int h = 0; h <= H; h++)
				sumCM = cplex.sum(sumCM, x1[i1][0][h]);
		cplex.addLe(sumCM, L);

		// Extension - cannot maintain if no gearboxes
		for (int i1 = 0; i1 < N; i1++)
			for (int i2 = 0; i2 <= M; i2++)
				cplex.addEq(x1[i1][i2][0], 0);

		// Solve the model
		cplex.solve();
		// Query the solution
		if (cplex.getStatus() == IloCplex.Status.Optimal) {
			System.out.println("Found optimal solution!");
			System.out.println("Objective = " + 12 * cplex.getObjValue());
			solutions.add(12 * cplex.getObjValue());
			for (int i1 = 0; i1 < N; i1++) {
				int thisAge = M;
				for (int i2 = 0; i2 <= M; i2++) {
					for (int h = 1; h <= H; h++) {
						if (i2 > 0 && i2 < thisAge && cplex.getValue(x1[i1][i2][h]) > 0) {
							ages.add(i2);
							thisAge = i2;
							break;
						}
					}
				}
				if (thisAge == M)
					ages.add(M);
			}
			for (int i1 = 0; i1 < N; i1++) {
				int thisAge = M;
				for (int i2 = 0; i2 <= M; i2++) {
					for (int h = 1; h <= H; h++) {
						if (i2 > 0 && i2 < thisAge && cplex.getValue(x2[i1][i2][h]) > 0) {
							ages_mod.add(i2);
							thisAge = i2;
						}
					}
				}
				if (thisAge == M)
					ages_mod.add(M);
			}

		} else {
			System.out.println("No optimal solution found");
			System.out.println(cplex.getStatus());
		}
		// Close the model
		cplex.close();

	}
}
