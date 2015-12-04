package net.finmath;

public class MatrixUtils {
    public static boolean checkSymmetry(double[][] matrix) {
        for (int row = 0; row < matrix.length; row++) {
            for (int col = 0; col < matrix[row].length; col++) {
                if(Math.abs(matrix[row][col] - matrix[col][row]) > 1e-14) return false;
            }
        }
        return true;
    }

    public static boolean matricesAreEqual(double[][] A, double[][] B, double threshold, int maxViolations) {
        if(A.length != B.length) return false;
        int violations = 0;
        for (int row = 0; row < A.length; row++) {
            if(A[row].length != B[row].length) return false;
            for (int col = 0; col < A[row].length; col++) {
                double a = A[row][col];
                double b = B[row][col];
                if(b == 0 && Math.abs(a) > threshold || b != 0 && Math.abs(a / b - 1.0) > threshold) {
                    violations++;
                }
            }
        }
        return violations <= maxViolations;
    }

    public static String printMatrix(double[][] matrix) {
        StringBuilder stringBuilder = new StringBuilder();
        for (double[] row : matrix) {
            for (double value : row) {
                stringBuilder.append(value).append("\t");
            }
            stringBuilder.append("\n");
        }
        return stringBuilder.toString();
    }
}
