/*
 * Copyright (c) 2018 Javier Orbe
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

package com.javierorbe.math;

import com.javierorbe.math.util.IntPair;
import com.javierorbe.math.util.MathUtils;

import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.BiFunction;
import java.util.function.Consumer;
import java.util.function.Function;

/**
 * Represents a mathematical matrix.
 *
 * @author Javier Orbe
 * @see <a href="https://en.wikipedia.org/wiki/Matrix_(mathematics)" target="_top">Matrix (mathematics) in Wikipedia</a>
 */
public class Matrix {

    private int rows;
    private int cols;
    private double[][] values;

    /**
     * Construct a matrix with the given values.
     *
     * @param values the values for the matrix.
     * @exception InvalidMatrixException if the size of the matrix is not at least 1.
     */
    public Matrix(double[][] values) {
        if (values.length == 0 || values[0].length == 0) {
            throw new InvalidMatrixException("A matrix cannot be of size zero.");
        }

        this.rows = values.length;
        this.cols = values[0].length;
        this.values = values;
    }

    /**
     * Construct a zero matrix of an specified size.
     *
     * @param rows the number of rows.
     * @param cols the number of columns.
     */
    public Matrix(int rows, int cols) {
        this(new double[rows][cols]);
    }

    /**
     * Returns the number of rows in the matrix.
     *
     * @return the number of rows.
     */
    public int getRows() {
        return rows;
    }

    /**
     * Returns the number of columns in the matrix.
     *
     * @return the number of columns.
     */
    public int getColumns() {
        return cols;
    }

    /**
     * Returns a copy of this matrix.
     *
     * @return a copy of this matrix.
     */
    public Matrix copy() {
    	return new Matrix(MathUtils.deepCopy(values));
    }

    /**
     * Returns the element at the specified row and column indices.
     *
     * @param row the row index.
     * @param col the column index.
     * @return the element at the specified indices.
     */
    public double get(int row, int col) {
        return values[row][col];
    }

    /**
     * Set a value at a specified index.
     *
     * @param row the row index.
     * @param col the column index.
     * @param value the value.
     */
    public void set(int row, int col, double value) {
        values[row][col] = value;
    }

    /**
     * Performs the given action for every element in the matrix.
     *
     * @param consumer the action to be performed for each element.
     */
    public void forEach(BiConsumer<Double, IntPair> consumer) {
        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                consumer.accept(values[row][col], new IntPair(row, col));
            }
        }
    }

    /**
     * Performs the given action for every row in the matrix.
     *
     * @param consumer the action to be performed for each row.
     */
    public void forEachRow(Consumer<double[]> consumer) {
        for (int row = 0; row < rows; row++) {
            consumer.accept(values[row]);
        }
    }

    /**
     * Apply a function to every element in this matrix.
     *
     * @param function the function that is applied to every element.
     *                 Gets as a parameter the value of the element being applied.
     */
    public void map(Function<Double, Double> function) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                values[i][j] = function.apply(values[i][j]);
            }
        }
    }
    
    /**
     * Apply a function to every element in this matrix.
     *
     * @param function the function that is applied to every element.
     *                 Gets as parameters the row and column indices of the element being applied.
     * @return this matrix.
     */
    public Matrix map(BiFunction<Integer, Integer, Double> function) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                values[i][j] = function.apply(i, j);
            }
        }

        return this;
    }

    /**
     * Returns a new matrix where the values are the return value of the given function.
     * The given function gets as a parameter every value at the original matrix.
     *
     * @param matrix the original matrix.
     * @param function the function that evaluates every value.
     * @return the new matrix.
     */
    public static Matrix map(Matrix matrix, Function<Double, Double> function) {
        Matrix result = new Matrix(matrix.rows, matrix.cols);

        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.cols; j++) {
                result.values[i][j] = function.apply(matrix.values[i][j]);
            }
        }

        return result;
    }

    /**
     * Randomizes the matrix values with values from 0 to 1.
     */
    public void randomize() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                values[i][j] = Math.random();
            }
        }
    }

    /**
     * Add a scalar to every element in the matrix.
     *
     * @param scalar the scalar
     */
    public void add(double scalar) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                values[i][j] += scalar;
            }
        }
    }

    /**
     * Add a matrix to this matrix.
     *
     * @param matrix the matrix to add to the this matrix.
     * @exception InvalidMatrixException if the matrices are not the same size.
     */
    public void add(Matrix matrix) {
        if (notSameSize(this, matrix)) {
            throw new InvalidMatrixException("The two matrices must be the same size.");
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.values[i][j] += matrix.get(i, j);
            }
        }
    }

    /**
     * Subtract matrix B from matrix A.
     *
     * @param a matrix A
     * @param b matrix B
     * @return the resultant matrix after subtracting B from A
     * @exception InvalidMatrixException if the matrices are not the same size.
     */
    public static Matrix subtract(Matrix a, Matrix b) {
        if (notSameSize(a, b)) {
            throw new InvalidMatrixException("The two matrices must be the same size.");
        }

        Matrix result = new Matrix(a.rows, a.cols);

        for (int i = 0; i < a.rows; i++) {
            for (int j = 0; j < a.cols; j++) {
                result.values[i][j] = a.values[i][j] - b.values[i][j];
            }
        }

        return result;
    }

    /**
     * Multiply the matrix by a scalar.
     *
     * @param scalar the number to multiply by the matrix
     */
    public void multiply(double scalar) {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                this.values[i][j] *= scalar;
            }
        }
    }

    /**
     * Hadamard product.
     * Multiply each element in this matrix with the element in the same position in the given matrix.
     *
     * @param matrix the matrix that will multiply this matrix.
     * @exception InvalidMatrixException if the matrices are not the same size.
     * @see <a href="https://en.wikipedia.org/wiki/Hadamard_product_(matrices)">Hadamard product in Wikipedia</a>
     */
    public void multiply(Matrix matrix) {
        if (notSameSize(this, matrix)) {
            throw new InvalidMatrixException("The two matrices must be the same size.");
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                values[i][j] *= matrix.values[i][j];
            }
        }
    }

    /**
     * Returns the result matrix of multiplying two matrices.
     * The resultant matrix has the same number of rows as the first operand matrix
     * and the same number of columns as the second operand matrix.
     * The order of multiplication matters.
     * The number of columns in the first matrix must equal the number of rows in the second matrix.
     *
     * @param a the first matrix in the multiplication.
     * @param b the second matrix in the multiplication.
     * @return the result of the multiplication.
     * @exception InvalidMatrixException if the number of columns in the first matrix
     *                                   doesn't equal the number of rows in the second matrix.
     */
    public static Matrix multiply(Matrix a, Matrix b) {
        if (a.cols != b.rows) {
            throw new InvalidMatrixException(
                    "The number of columns in the first matrix must equal the number of rows in the second matrix."
            );
        }

        Matrix result = new Matrix(a.rows, b.cols);

        for (int i = 0; i < result.rows; i++) {
            for (int j = 0; j < result.cols; j++) {
                for (int k = 0; k < a.cols; k++) {
                    result.values[i][j] += a.values[i][k] * b.values[k][j];
                }
            }
        }

        return result;
    }

    /**
     * Returns the transpose of a matrix.
     *
     * @param matrix the matrix to transpose.
     * @return the transposed matrix.
     * @see <a href="https://en.wikipedia.org/wiki/Transpose" target="_top">Transpose in Wikipedia</a>
     */
    public static Matrix transpose(Matrix matrix) {
        Matrix result = new Matrix(matrix.cols, matrix.rows);

        for (int i = 0; i < matrix.rows; i++) {
            for (int j = 0; j < matrix.cols; j++) {
                result.values[j][i] = matrix.values[i][j];
            }
        }

        return result;
    }

    /**
     * Get the submatrix of this matrix with a row and a column removed.
     *
     * @param removedRow Index of the removed row.
     * @param removedCol Index of the removed column.
     * @return a submatrix of this matrix.
     */
    private Matrix getSubmatrix(int removedRow, int removedCol) {
        Matrix submatrix = new Matrix(rows - 1, cols - 1);

        int i = 0;
        int j = 0;

        for (int row = 0; row < rows; row++) {
            if (row == removedRow) {
                continue;
            }

            for (int col = 0; col < cols; col++) {
                if (col == removedCol) {
                    continue;
                }

                submatrix.set(i, j, this.get(row, col));
                j++;
            }

            j = 0;
            i++;
        }

        return submatrix;
    }

    /**
     * Returns {@code true} if the matrix is square.
     * If the row count equals the column count, a matrix is square.
     *
     * @return {@code true} if the matrix is square, otherwise {@code false}.
     */
    public boolean isSquare() {
        return rows == cols;
    }

    /**
     * Returns {@code true} if this matrix is an identity matrix.
     *
     * @return {@code true} if this matrix is an identity matrix, otherwise {@code false}.
     * @exception InvalidMatrixException if the matrix is not square.
     */
    public boolean isIdentity() {
        if (!isSquare()) {
            throw new InvalidMatrixException("The matrix must be square.");
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (i == j) {
                    if (values[i][j] != 1) {
                        return false;
                    }
                } else if (values[i][j] != 0) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Get the identity matrix of a specific size.
     *
     * @param size the size of the identity matrix.
     * @return the identity matrix of the specified size.
     */
    public static Matrix getIdentity(int size) {
        Matrix result = new Matrix(size, size);

        for (int i = 0; i < size; i++) {
            result.values[i][i] = 1;
        }

        return result;
    }

    /**
     * Returns {@code true} if the matrix is a zero matrix.
     *
     * @return {@code true} if the matrix is a zero matrix, otherwise {@code false}.
     * @see <a href="https://en.wikipedia.org/wiki/Zero_matrix" target="_top">Zero matrix in Wikipedia</a>
     */
    public boolean isZero() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                if (values[i][j] != 0) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Returns {@code true} if the matrix has inverse.
     *
     * @return {@code true} if the matrix has inverse, otherwise {@code false}.
     */
    public boolean isInvertible() {
        return getDeterminant() != 0;
    }

    /**
     * Calculate the inverse of this matrix.
     *
     * @param inverseFunction the function that will calculate the inverse.
     * @return the inverse of this matrix.
     */
    public Matrix getInverse(InverseFunction inverseFunction) {
        return inverseFunction.function.apply(this);
    }

    /**
     * Get the determinant of this matrix.
     *
     * @return the determinant of this matrix.
     */
    public double getDeterminant() {
        return getDeterminantByLaplace(this);
    }

    /**
     * Calculate the determinant of a matrix using Laplace's formula.
     *
     * @param matrix the matrix.
     * @return the determinant of the matrix.
     * @exception InvalidMatrixException if the matrix is not square.
     * @see <a href="https://en.wikipedia.org/wiki/Laplace_expansion" target="_top">Laplace expansion in Wikipedia</a>
     */
    public static double getDeterminantByLaplace(Matrix matrix) {
        if (!matrix.isSquare()) {
            throw new InvalidMatrixException("The matrix must be square.");
        }

        // subMatrix is square, so rows == cols
        if (matrix.rows == 1) {
            return matrix.get(0, 0);
        } else {
            double sum = 0;

            // Fixed i = 1
            for (int j = 0; j < matrix.cols; j++) {
                double sign = Math.pow(-1, 1 + (j + 1));
                sum += sign * matrix.get(0, j) * matrix.getSubmatrix(0, j).getDeterminant();
            }

            return sum;
        }
    }

    /**
     * Get the (i,j)-minor.
     *
     * @param i Row index.
     * @param j Column index.
     * @return the (i,j)-minor.
     *
     * @see <a href="https://en.wikipedia.org/wiki/Minor_(linear_algebra)" target="_top">Minor (linear algebra) in Wikipedia</a>
     */
    public double getMinor(int i, int j) {
        return getSubmatrix(i, j).getDeterminant();
    }

    /**
     * Get the (i,j) adjugate of a matrix.
     *
     * @param i Row index.
     * @param j Column index.
     * @return the adjugate (i,j) of the matrix.
     */
    public double getAdjugate(int i, int j) {
        int sign = (int) Math.pow(-1, i + j);
        double minor = getMinor(i, j);
        return sign * minor;
    }

    /**
     * Get adjugate of this matrix.
     *
     * @return the adjugate of this matrix.
     */
    public Matrix getAdjugate() {
        return copy().map(this::getAdjugate);
    }

    /**
     * Swap two rows.
     *
     * @param indexRow1 index of the first row.
     * @param indexRow2 index of the second row.
     */
    public void swapRows(int indexRow1, int indexRow2) {
        double[] row1 = values[indexRow1];
        double[] aux = Arrays.copyOf(row1, row1.length);
        values[indexRow1] = values[indexRow2];
        values[indexRow2] = aux;
    }

    /**
     * Swap two columns.
     *
     * @param indexCol1 index of the first column.
     * @param indexCol2 index of the second column.
     */
    public void swapColumns(int indexCol1, int indexCol2) {
        for (int row = 0; row < rows; row++) {
            double aux = values[row][indexCol1];
            values[row][indexCol1] = values[row][indexCol2];
            values[row][indexCol2] = aux;
        }
    }

    /**
     * Subtract to a row another row multiplied by a scalar.
     * row1' = row1 - (scalar * row2)
     *
     * @param row1 The affected row.
     * @param row2 The row to subtract.
     * @param scalar The scalar that multiplies the second row.
     */
    public void subtractRow(int row1, int row2, double scalar) {
        for (int col = 0; col < cols; col++) {
            values[row1][col] -= scalar * values[row2][col];
        }
    }

    /**
     * Subtract to a column another column multiplied by a scalar.
     * col1' = col1 - (scalar * col2)
     *
     * @param col1 The affected column.
     * @param col2 The column to subtract.
     * @param scalar The scalar that multiplies the second column.
     */
    public void subtractColumn(int col1, int col2, double scalar) {
        for (int row = 0; row < rows; row++) {
            values[row][col1] -= scalar * values[row][col2];
        }
    }

    /**
     * Multiply a row of this matrix by a scalar.
     *
     * @param row row to multiply.
     * @param scalar multiplier.
     */
    public void multiplyRow(int row, double scalar) {
        for (int col = 0; col < cols; col++) {
            values[row][col] *= scalar;
        }
    }

    /**
     * Multiply a column of this matrix by a scalar.
     *
     * @param col the column index.
     * @param scalar the multiplier.
     */
    public void multiplyColumn(int col, double scalar) {
        for (int row = 0; row < rows; row++) {
            values[row][col] *= scalar;
        }
    }

    /**
     * Get the row echelon form of this matrix.
     *
     * @param reduced {@code true} if the wanted result is a reduced row echelon form.
     * @return the row echelon form of the matrix.
     * @see <a href="https://en.wikipedia.org/wiki/Row_echelon_form" target="_top">Row echelon form in Wikipedia</a>
     */
    public Matrix getRowEchelon(boolean reduced) {
        Matrix result = copy();

        // Order the rows in descending order by the amount of zeroes in the first elements until the first non zero element.
        List<Map.Entry<Integer, Integer>> sortedList = getSortedZeroCount(result);
        List<Integer> sortedIndices = new ArrayList<>();
        sortedList.forEach(entry -> sortedIndices.add(entry.getKey()));
        for (int i = 0; i < sortedIndices.size(); i++) {
            int rowIndex = sortedIndices.get(i);

            if (rowIndex != i) {
                result.swapRows(rowIndex, i);

                for (int j = 0; j < sortedIndices.size(); j++) {
                    if (sortedIndices.get(j) == i) {
                        sortedIndices.set(j, rowIndex);
                        break;
                    }
                }

                sortedIndices.set(i, i);
            }
        }

        // Go through all the rows. If the result doesn't have to be reduced, the last one is skipped.
        int offset = reduced ? 0 : 1;
        for (int row = 0; row < result.rows - offset; row++) {
            int pivotIndex = getPivotIndex(result.values[row]);

            // If all the elements in the row are zeroes, it means that the matrix is a zero matrix.
            if (pivotIndex == -1) {
                return result;
            }

            // Set the pivot to 1 (multiplying the row by the inverse of the pivot).
            result.multiplyRow(row, 1 / result.values[row][pivotIndex]);

            // Set all the elements below the pivot to zero
            for (int k = row + 1; k < result.rows; k++) {
                result.subtractRow(k, row, result.values[k][pivotIndex]);
            }

            if (reduced) {
                // Set all the elements above the pivot to zero
                for (int k = 0; k < row; k++) {
                    result.subtractRow(k, row, result.values[k][pivotIndex]);
                }
            }
        }

        return result;
    }

    /**
     * Get the index of the pivot in a row.
     * The pivot is the first none zero element.
     *
     * @param row the row.
     * @return the index of the pivot if there is, otherwise {@code -1}.
     */
    private static int getPivotIndex(double[] row) {
        int pivotIndex = -1;

        for (int i = 0; i < row.length; i++) {
            if (row[i] != 0) {
                pivotIndex = i;
                break;
            }
        }

        return pivotIndex;
    }

    // TODO: calcular el rango de una matriz mediante el determinante

    /**
     * Returns the range of this matrix.
     *
     * @return the range of this matrix.
     */
    public int getRange() {
        Matrix echelon = getRowEchelon(false);

        // Get the amount of rows that are not full of zeroes
        int nonZero = rows;

        for (int row = 0; row < rows; row++) {
            int zeroes = 0;
            for (int col = 0; col < cols; col++) {
                if (echelon.values[row][col] == 0) {
                    zeroes++;
                }
            }
            // The row is full of zeroes
            if (zeroes == cols) {
                nonZero--;
            }
        }

        return nonZero;
    }

    /**
     * Convert a column matrix to an array.
     *
     * @return an array with the values of the matrix.
     * @exception InvalidMatrixException if the matrix is not a column matrix.
     */
    public double[] toArray() {
        if (cols != 1) {
            throw new InvalidMatrixException("The matrix must only have a column.");
        }

        double[] result = new double[rows];

        for (int i = 0; i < rows; i++) {
            result[i] = values[i][0];
        }

        return result;
    }

    /**
     * Compares the specified object with this matrix for equality.
     * Returns {@code true} if and only if the specified object is also a matrix,
     * both matrices have the save size, and all corresponding elements in the two
     * matrices are equal.
     *
     * @param o the object to be compared for equality with this matrix.
     * @return {@code true} if the specified object is equal to this matrix, otherwise {@code false}.
     */
    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }

        if (!(o instanceof Matrix)) {
            return false;
        }

        Matrix m = (Matrix) o;

        if (notSameSize(this, m)) {
            return false;
        }

        for (int row = 0; row < rows; row++) {
            for (int col = 0; col < cols; col++) {
                if (Double.compare(values[row][col], m.values[row][col]) != 0) {
                    return false;
                }
            }
        }

        return true;
    }

    /**
     * Returns a representation of this matrix.
     *
     * @return a string representation of this matrix.
     */
    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();

        for (int i = 0; i < rows; i++) {
            builder.append("[");
            for (int j = 0; j < cols; j++) {
                String s = String.valueOf(values[i][j]);
                s = !s.contains(".") ? s : s.replaceAll("0*$", "").replaceAll("\\.$", "");
                builder.append(s).append(j < cols - 1 ? " " : "");
            }

            builder.append("]").append("\n");
        }

        return builder.toString();
    }

    /**
     * Returns {@code true} if the two matrices have the same number of rows and the same number of columns.
     *
     * @param a first matrix.
     * @param b second matrix.
     * @return {@code true} if the the matrices are the same size, otherwise {@code false}.
     */
    private static boolean notSameSize(Matrix a, Matrix b) {
        return a.rows != b.rows || a.cols != b.cols;
    }

    /**
     * Returns the augmented matrix of two matrices.
     *
     * @param a the first matrix.
     * @param b the second matrix.
     * @return the augmented matrix.
     * @exception InvalidMatrixException if the number of rows of both matrices are not the same.
     * @see <a href="https://en.wikipedia.org/wiki/Augmented_matrix" target="_top">Augmented matrix in Wikipedia</a>
     */
    public static Matrix getAugmentedMatrix(Matrix a, Matrix b) {
        if (a.rows != b.rows) {
            throw new InvalidMatrixException("The number of rows of both matrices must be the same.");
        }

        Matrix result = new Matrix(a.rows, a.cols + b.cols);

        for (int row = 0; row < result.rows; row++) {
            for (int col = 0; col < result.cols; col++) {
                if (col < a.cols) {
                    result.values[row][col] = a.values[row][col];
                } else {
                    result.values[row][col] = b.values[row][col - a.cols];
                }
            }
        }

        return result;
    }

    /**
     * Create a column matrix from an array.
     *
     * @param columnValues values for the matrix.
     * @return the column matrix.
     * @exception InvalidMatrixException if there are no values for the matrix.
     */
    public static Matrix getColumnMatrix(double[] columnValues) {
        if (columnValues.length == 0) {
            throw new InvalidMatrixException("A matrix cannot be of size zero.");
        }

        double[][] values = new double[columnValues.length][1];

        for (int i = 0; i < columnValues.length; i++) {
            values[i][0] = columnValues[i];
        }

        return new Matrix(values);
    }

    /**
     * Functions that calculate the inverse of a matrix.
     *
     * @see <a href="https://en.wikipedia.org/wiki/Invertible_matrix" target="_top">Invertible matrix in Wikipedia</a>
     */
    public enum InverseFunction {
        /**
         * Calculate the inverse of a matrix using adjugates.
         */
        ADJUGATES((matrix) -> {
            if (!matrix.isInvertible()) {
                throw new InvalidMatrixException("The matrix must be invertible.");
            }

            Matrix result = Matrix.transpose(matrix.getAdjugate());
            result.multiply((1 / matrix.getDeterminant()));

            return result;
        }),

        /**
         * Calculate the inverse of a matrix using the Gauss method.
         */
        GAUSS((matrix) -> {
            if (!matrix.isInvertible()) {
                throw new InvalidMatrixException("The matrix must be invertible.");
            }

            Matrix initial = matrix.copy();
            Matrix identity = Matrix.getIdentity(matrix.rows);

            // Order the rows in descending order by the amount of zeroes in the first elements until the first non zero element.
            List<Map.Entry<Integer, Integer>> sorted = getSortedZeroCount(initial);
            double[][] sortedValues = new double[initial.rows][initial.cols];
            double[][] identitySwapValues = new double[identity.rows][identity.cols];
            for (int i = 0; i < sorted.size(); i++) {
                int row = sorted.get(i).getKey();
                sortedValues[i] = initial.values[row];
                identitySwapValues[i] = identity.values[row];
            }
            initial.values = sortedValues;
            identity.values = identitySwapValues;

            // Go through all the rows
            for (int row = 0; row < initial.rows; row++) {
                int pivotIndex = getPivotIndex(initial.values[row]);

                // If all the elements in the row are zeroes, it means that the matrix is a zero matrix.
                if (pivotIndex == -1) {
                    return initial;
                }

                // Set the pivot to 1 (multiplying the row by the inverse of the pivot).
                double scalar = 1 / initial.values[row][pivotIndex];
                initial.multiplyRow(row, scalar);
                identity.multiplyRow(row, scalar);

                // Set all the elements below the pivot to zero
                for (int k = row + 1; k < initial.rows; k++) {
                    scalar = initial.values[k][pivotIndex];
                    initial.subtractRow(k, row, scalar);
                    identity.subtractRow(k, row, scalar);
                }

                // Set all the elements above the pivot to zero
                for (int k = 0; k < row; k++) {
                    scalar = initial.values[k][pivotIndex];
                    initial.subtractRow(k, row, scalar);
                    identity.subtractRow(k, row, scalar);
                }
            }

            return identity;
        })
        ;

        /**
         * The default method for getting the inverse of a matrix.
         */
        public static final InverseFunction DEFAULT = InverseFunction.ADJUGATES;

        private Function<Matrix, Matrix> function;

        InverseFunction(Function<Matrix, Matrix> function) {
            this.function = function;
        }
    }

    /**
     * Sort a map entry list in ascending order by its values.
     *
     * @param map the map to sort.
     * @param <K> key type.
     * @param <V> value type.
     * @return a sorted list of map entries.
     */
    private static <K, V extends Comparable<? super V>> List<Map.Entry<K, V>> entriesSortedByValues(Map<K, V> map) {
        List<Map.Entry<K, V>> sortedEntries = new ArrayList<>(map.entrySet());
        sortedEntries.sort(Comparator.comparing(Map.Entry::getValue));
        return sortedEntries;
    }

    /**
     * Get the sorted list of map entries with keys as column indices of the matrix and as values
     * the amount of zeroes in the first elements of that column until the first none zero element.
     *
     * @param matrix the matrix.
     * @return the sorted list.
     */
    private static List<Map.Entry<Integer, Integer>> getSortedZeroCount(Matrix matrix) {
        Map<Integer, Integer> zeroCount = new HashMap<>();
        for (int row = 0; row < matrix.rows; row++) {
            zeroCount.put(row, 0);

            for (int col = 0; col < matrix.cols; col++) {
                double element = matrix.values[row][col];

                if (element == 0) {
                    zeroCount.put(row, zeroCount.get(row) + 1);
                } else {
                    break;
                }
            }
        }

        // Order the rows in descending order by the amount of zeroes in the first elements until the first non zero element.
        return entriesSortedByValues(zeroCount);
    }
}
