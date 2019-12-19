import java.io.FileWriter;
import java.io.PrintWriter;
import java.lang.management.ManagementFactory;
import java.lang.management.ThreadMXBean;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.*;
import java.util.List;
import java.util.Random;
import java.awt.geom.Point2D;

public class TravelingSalesman {
    static int numberOfTrials = 50; //  50 for brute, 1000 for greedy + 1, 50 for greedyX2
    static int MAXINPUTSIZE = 10; // 11 for brute, 50 for greedy + 1, 4100 for greedyX2
    static int MAXVALUE = 100;
    static int VERTICES = 4;

    static String ResultsFolderPath = "/home/ryan/Results/"; // pathname to results folder
    static FileWriter resultsFile;
    static PrintWriter resultsWriter;

    public static void main(String[] args) {
        double totalCost = 0;
        //double[][] testMatrix = new double[VERTICES][VERTICES];
        // very confusing 2d array to store points
        // thinking like a c programmer
        // but later realizing this is what structs are for
        //int[][] testMatrixPoints = new int[VERTICES][2];
        // stores the path from the greedy algorithm i wrote
        List<Integer> path = new ArrayList<Integer>();
        // testing to compare outputs, printing, etc
        for ( VERTICES = 4; VERTICES < 8; VERTICES++) {
            System.out.println("Vertices = " + VERTICES);

            double[][] testMatrix = new double[VERTICES][VERTICES];
            generateRandomCostMatrix(testMatrix, VERTICES, MAXVALUE);
            printCostMatrix(testMatrix, VERTICES);

            totalCost = greedyTSP(path, testMatrix, VERTICES);
            printThePath(path, VERTICES, totalCost);

            // array to store the best tour of the bruteForceTSP
            int[] bestTour = bruteForceTSP(testMatrix);
            double bfTourCost = computeTourCost(bestTour, testMatrix);
            System.out.println("Tour cost for bruteForce is: " + bfTourCost);
            printBFTour(bestTour);
            path.removeAll(path);
        }

/*        double countMatches = 0;
        // testing to compare percentage of costs that are the same
        for (VERTICES = 4; VERTICES < 11; ++VERTICES)
        {
            for (numberOfTrials = 0; numberOfTrials < 100; numberOfTrials++) {

                //System.out.println("Vertices = " + VERTICES);
                double[][] testMatrix = new double[VERTICES][VERTICES];
                generateRandomCostMatrix(testMatrix, VERTICES, MAXVALUE);
                //printCostMatrix(testMatrix, VERTICES);

                totalCost = greedyTSP(path, testMatrix, VERTICES);
                //printThePath(path, VERTICES, totalCost);

                // array to store the best tour of the bruteForceTSP
                int[] bestTour = bruteForceTSP(testMatrix);
                double bfTourCost = computeTourCost(bestTour, testMatrix);
                // System.out.println("Tour cost for bruteForce is: " + bfTourCost);
                //printBFTour(bestTour);
                path.removeAll(path);
                if (totalCost == bfTourCost)
                    countMatches++;

            }
            System.out.println("The percentage of matching tour costs for " + VERTICES + " vertices is " + (countMatches/numberOfTrials * 100) + "%");
            countMatches = 0;
*/
//        }

        // just to see if there is an expected rate
        //System.out.println("The percentage of matching tour costs for " + VERTICES + " vertices is " + (countMatches/numberOfTrials * 100));

        //generateRandomEuclideanCostMatrix(testMatrixPoints, testMatrix, VERTICES, MAXVALUE);
        //printPointMatrix(testMatrixPoints, VERTICES);
        //printCostMatrix(testMatrix, VERTICES);

        // writing to our data files
        /* **********************************************UNCOMMENT ONE******************************************/
 /*       System.out.println("Running first full experiment ...");
        runFullExperiment("BFTSP-Exp1-ThrowAway.txt");
        System.out.println("Running second full experiment ...");
        runFullExperiment("BFTSP-Exp2.txt");
        System.out.println("Running third full experiment ...");
        runFullExperiment("BFTSP-Exp3.txt");
/*
        System.out.println("Running first full experiment ...");
        runFullExperiment("GreedyTSP-Exp1-ThrowAway.txt");
        System.out.println("Running second full experiment ...");
        runFullExperiment("GreedyTSP-Exp2.txt");
        System.out.println("Running third full experiment ...");
        runFullExperiment("GreedyTSP-Exp3.txt");
*/
        // for the doubling but probably wont use due to how slow bruteForce is
//        System.out.println("Running first full experiment ...");
//        runFullExperiment("BFTSPX2-Exp1-ThrowAway.txt");
//        System.out.println("Running second full experiment ...");
//        runFullExperiment("BFTSPX2-Exp2.txt");
//        System.out.println("Running third full experiment ...");
//       runFullExperiment("BFTSPX2-Exp3.txt");
        // for doubling
/*        System.out.println("Running first full experiment ...");
        runFullExperiment("GreedyTSPX2-Exp1-ThrowAway.txt");
        System.out.println("Running second full experiment ...");
        runFullExperiment("GreedyTSP-Exp2.txt");
        System.out.println("Running third full experiment ...");
        runFullExperiment("GreedyTSPX2-Exp3.txt");
*/
        System.out.println("Running first full experiment ...");
        runFullExperiment("GreedyTSPSQR-Exp1-ThrowAway.txt");
        System.out.println("Running second full experiment ...");
        runFullExperiment("GreedyTSPSQR-Exp2.txt");
        System.out.println("Running third full experiment ...");
        runFullExperiment("GreedyTSPSQR-Exp3.txt");
       /* **********************************************UNCOMMENT ONE******************************************/
    }

    static void runFullExperiment(String resultsFileName)
    {
        // making sure the desired files exist or can be created
        try {
            resultsFile = new FileWriter(ResultsFolderPath + resultsFileName);
            resultsWriter = new PrintWriter(resultsFile);
        } catch (Exception e) {
            System.out.println("*****!!!!!  Had a problem opening the results file " + ResultsFolderPath + resultsFileName);
            return;
        }

        ThreadCPUStopWatch BatchStopwatch = new ThreadCPUStopWatch(); // for timing an entire set of trials
        //ThreadCPUStopWatch TrialStopwatch = new ThreadCPUStopWatch(); // for timing an individual trial
        // uncomment this for other tests
//        resultsWriter.println("#Nodes(N)    AverageTime "); // # marks a comment in gnuplot data
//        resultsWriter.flush();

        // for sqr
        resultsWriter.println("#  Nodes(N)       SQR "); // # marks a comment in gnuplot data
        resultsWriter.flush();
        // for each size of input we want to test: in this case incrementing by 1

        for ( int VERTICES = 4; VERTICES <= MAXINPUTSIZE; VERTICES++)
        {
            double tmpsqr = 0;
            for ( int i = 1; i <= numberOfTrials; i++)
            {
                double[][] testMatrix = new double[VERTICES][VERTICES];
                int[][] pointMatrix = new int[VERTICES][2];
                populatePointMatrix(pointMatrix, VERTICES, 100);
                generateRandomEuclideanCostMatrix(pointMatrix, testMatrix, VERTICES, 100);
                System.out.println("tmp sqr is " + tmpsqr);
                tmpsqr += sqrTest(testMatrix, VERTICES);
            }
            double avgSqr = tmpsqr/numberOfTrials;
            System.out.println("avg sqr is " + avgSqr);
            resultsWriter.printf("%12d  %15.2f\n", VERTICES, avgSqr); // might as well make the columns look nice
            resultsWriter.flush();

        }
 /*
        for (int VERTICES = 4; VERTICES <= MAXINPUTSIZE; VERTICES *= 2) {
            /* repeat for desired number of trials (for a specific size of input)... */
 //           System.out.println("Running test for input size " + VERTICES + " ... ");
            // creating new matrix of the appropriate size
 //           double[][] testMatrix = new double[VERTICES][VERTICES];
            // generating random values for the cost matrix
 //           generateRandomCostMatrix(testMatrix, VERTICES, MAXVALUE);
            // will not really be used here but necessary for params
  //          List<Integer> path = new ArrayList<Integer>();
            // will hold total amount of time
            // will reset after each batch of trials
   //         long batchElapsedTime = 0;

            /* force garbage collection before each batch of trials run so it is not included in the time */
/*            System.gc();
            System.out.print("    Running trial batch...");
            BatchStopwatch.start(); // comment this line if timing trials individually

            // run the trials
            for (int trial = 0; trial < numberOfTrials; trial++) {

                //actual beginning of trial
                /* **********************************************UNCOMMENT ONE******************************************/
                //bruteForceTSP(testMatrix);
/*                greedyTSP(path, testMatrix, VERTICES);
                /* **********************************************UNCOMMENT ONE^^^******************************************/
                //System.out.println("....done.");// *** uncomment this line if timing trials individually
//            }

 /*           batchElapsedTime = BatchStopwatch.elapsedTime(); // *** comment this line if timing trials individually

            double averageTimePerTrialInBatch = (double) batchElapsedTime / (double) numberOfTrials; // calculate the average time per trial in this batch
            /* print data for this size of input */

  /*          resultsWriter.printf("%12d  %15.2f\n", VERTICES, (double) averageTimePerTrialInBatch); // might as well make the columns look nice
            resultsWriter.flush();
            System.out.println(" ....done.");
        }
*/
    }

    // creates a cost matrix of size verticesXvertices with an upper bound of maxEdgeCost
    // the diagonal line/no trip = 0
    static void generateRandomCostMatrix(double[][] matrix, int vertices, int maxEdgeCost) {
        //matrix = new int[vertices][vertices];
        for (int i = 0; i < vertices; i++) {
            for (int j = i + 1; j < vertices; j++) {
                if (i != j) {
                    double tmpRand = randomNumGen(maxEdgeCost);
                    matrix[i][j] = tmpRand;
                    matrix[j][i] = tmpRand;
                } else {
                    matrix[i][j] = 0;
                }
            }
        }
    }

    // returns a random number between 0 and maxEdgeCost
    static double randomNumGen(int maxEdgeCost) {
        Random r = new Random();
        return r.nextInt(maxEdgeCost) + 1;
    }

    // prints a matrix. passed vertices for readability purposes
    static void printCostMatrix(double[][] matrix, int vertices) {
        for (int i = 0; i < vertices; i++) {
            for (int j = 0; j < vertices; j++) {
                System.out.printf("%6.2f", matrix[i][j]);
            }
            System.out.print("\n");
        }
    }

    static void generateRandomEuclideanCostMatrix(int[][] matrixPoints, double[][] costMatrix, int vertices, int maxXYValue) {
        populatePointMatrix(matrixPoints, vertices, maxXYValue);

        for (int i = 0; i < vertices; i++) {
            for (int j = i + 1; j < vertices; j++) {
                if (i != j) {
                    double tmpDist = Point2D.distance(matrixPoints[i][0], matrixPoints[i][1], matrixPoints[j][0], matrixPoints[j][1]);
                    costMatrix[i][j] = tmpDist;
                    costMatrix[j][i] = tmpDist;
                } else {
                    costMatrix[i][j] = 0;
                }

            }
        }
    }

    // making random points
    static void populatePointMatrix(int[][] matrixPoints, int vertices, int maxEdgeCost) {
        for (int i = 0; i < vertices; i++) {
            for (int j = 0; j < 2; j++) {
                matrixPoints[i][j] = (int)randomNumGen(maxEdgeCost);
            }
        }
    }

    // prints the point matrix
    static void printPointMatrix(int[][] pointMatrix, int vertices) {
        for (int i = 0; i < vertices; i++) {
            for (int j = 0; j < 2; j++) {
                System.out.printf("%6d", pointMatrix[i][j]);
            }
            System.out.println();
        }
    }

    // takes the shortest path for each point that has not already been used
    static double greedyTSP(List<Integer> path, double[][] costMatrix, int vertices) {
        double totalCost = 0;
        double tmpShortestPath = 0;
        int tmpIndex = 0;
        int curNode = 0;
        // starting location
        //path.add(0);


        while (path.size() < vertices) {

            tmpShortestPath = MAXVALUE + 1;
            tmpIndex = curNode;
            for (int j = 0; j < vertices; j++) {
                // finding the shortest trip for each row and making sure the vertex has not been visited previously
                if (!path.contains(j) && (j != 0)) {
                    if ((costMatrix[tmpIndex][j] > 0) && (costMatrix[tmpIndex][j] <= tmpShortestPath)) {
                        tmpShortestPath = costMatrix[tmpIndex][j];
                        //System.out.println(tmpShortestPath + " j = " + j);
                        curNode = j;
                    }
                }
            }
            //System.out.println("totalCost for greedy " + totalCost);
            path.add(tmpIndex);
        }
        // adding the last trip back to 0
        path.add(0);
        // kept getting a goofy error so had to do this :D
        for (int a = 0; a < vertices; a++) {
            totalCost += costMatrix[path.get(a)][path.get(a + 1)];
        }
        //System.out.println("totalCost is " + totalCost);
        //totalCost += costMatrix[0][curNode];
        return totalCost;
    }

    // pretty prints the path
    static void printThePath(List<Integer> path, int vertices, double cost) {
        System.out.println("The total cost for greedy is: " + cost);
        System.out.print("The path for greedy is: ");
        path.forEach(System.out::print);
        System.out.println();
    }


   /* https://github.com/williamfiset/Algorithms/blob/ccebd26d3b0c9aefec661400e655804d9efda1b9/com/williamfiset/algorithms/graphtheory/TspBruteForce.java#L100*/
    public static int[] bruteForceTSP(double[][] matrix)
    {
        int n = matrix.length;
        int[] permutation = new int[n];
        for ( int i = 0; i < n; i++) permutation[i] = i;

        //copying the current best tour
        int[] bestTour = permutation.clone();
        double bestTourCost = Double.POSITIVE_INFINITY;

        do {

            double tourCost = computeTourCost(permutation, matrix);

            if ( tourCost < bestTourCost) {
                bestTourCost = tourCost;
                bestTour = permutation.clone();
            }
        } while ( nextPermutation(permutation));


        return bestTour;

    }

    // determines tour cost based on the best tour
    public static double computeTourCost(int[] tour, double[][] matrix) {
       double cost = 0;

       for ( int i = 1; i < matrix.length; i++)
       {
           int from = tour[i-1];
           int to = tour[i];
           cost += matrix[from][to];
       }

       int last = tour[matrix.length-1];
       int first = tour[0];
       return cost + matrix[last][first];
    }

    // finds next permutation
    public static boolean nextPermutation(int[] sequence)
    {
        int first = getFirst(sequence);
        if ( first == -1) return false;
        int toSwap = sequence.length -1;
        while ( sequence[first] >= sequence[toSwap]) --toSwap;
        swap(sequence, first++, toSwap);
        toSwap = sequence.length -1;
        while (first < toSwap) swap(sequence, first++, toSwap--);
        return true;
    }

    // returns
    private static int getFirst(int[] sequence) {
        for ( int i = sequence.length - 2; i >= 0; --i) if ( sequence[i] < sequence[i + 1]) return i;

        return -1;
    }

    // swaps values when creating the next permutation
    private static void swap(int[] sequence, int i, int j) {
        int tmp = sequence[i];
        sequence[i] = sequence[j];
        sequence[j] = tmp;
    }

    // prints 0 at the end of the array
    public static void printBFTour(int[] bestTour)
    {
        System.out.print("The best BFTour is: ");
        for ( int i = 0; i < bestTour.length; i++)
        {
            System.out.print(bestTour[i] + " ");
        }
        System.out.println(0);
    }

    public static double sqrTest(double[][] a, int vertices)
    {
        List<Integer> sqrPath = new ArrayList<Integer>();
        int[] bestTour = bruteForceTSP(a);
        printCostMatrix(a, vertices);
        double bestCost = computeTourCost(bestTour, a);
        System.out.println("best is " + bestCost);
        double greedyCost = greedyTSP(sqrPath, a, vertices);
        printThePath(sqrPath, VERTICES,greedyCost);
        System.out.println("greedy is " + greedyCost);
        return bestCost/greedyCost;
    }

}

