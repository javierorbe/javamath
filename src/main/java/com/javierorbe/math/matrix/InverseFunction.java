/*
 * Copyright (c) 2019 Javier Orbe
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

package com.javierorbe.math.matrix;

import java.util.List;
import java.util.Map;
import java.util.function.Function;

import static com.javierorbe.math.util.MathUtils.getPivotIndex;
import static com.javierorbe.math.util.MathUtils.getSortedZeroCount;

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

        Matrix result = matrix.getAdjugate().transpose();
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
        Matrix identity = Matrix.getIdentity(matrix.getRows());

        // Order the rows in descending order by the amount of zeroes in the first elements until the first non zero element.
        List<Map.Entry<Integer, Integer>> sorted = getSortedZeroCount(initial);
        double[][] sortedValues = new double[initial.getRows()][initial.getColumns()];
        double[][] identitySwapValues = new double[identity.getRows()][identity.getColumns()];
        for (int i = 0; i < sorted.size(); i++) {
            int row = sorted.get(i).getKey();
            sortedValues[i] = initial.getRow(row);
            identitySwapValues[i] = identity.getRow(row);
        }
        initial = new Matrix(sortedValues);
        identity = new Matrix(identitySwapValues);

        // Go through all the rows
        for (int row = 0; row < initial.getRows(); row++) {
            int pivotIndex = getPivotIndex(initial.getRow(row));

            // If all the elements in the row are zeroes, it means that the matrix is a zero matrix.
            if (pivotIndex == -1) {
                return initial;
            }

            // Set the pivot to 1 (multiplying the row by the inverse of the pivot).
            double scalar = 1 / initial.get(row, pivotIndex);
            initial.multiplyRow(row, scalar);
            identity.multiplyRow(row, scalar);

            // Set all the elements below the pivot to zero
            for (int k = row + 1; k < initial.getRows(); k++) {
                scalar = initial.get(k, pivotIndex);
                initial.subtractRow(k, row, scalar);
                identity.subtractRow(k, row, scalar);
            }

            // Set all the elements above the pivot to zero
            for (int k = 0; k < row; k++) {
                scalar = initial.get(k, pivotIndex);
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

    public Function<Matrix, Matrix> getFunction() {
        return function;
    }
}