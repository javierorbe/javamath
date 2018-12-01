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

import static com.javierorbe.math.util.MathUtils.trimZeroes;

/**
 * Represents a system of linear equations.
 *
 * @author Javier Orbe
 * @see <a href="https://en.wikipedia.org/wiki/System_of_linear_equations" target="_top">System of linear equations in Wikipedia</a>
 */
public class LinearEquationSystem {

    private Matrix coefficients;
    private Matrix constants;

    /**
     * Construct a system of linear equations.
     *
     * @param coefficients the coefficients.
     * @param constants the constant terms.
     * @exception InvalidLinearEquationSystemException if one or more of the equations is null (all the coefficients are zero).
     */
    public LinearEquationSystem(Matrix coefficients, Matrix constants) {
        coefficients.forEachRow((row) -> {
            for (double element : row) {
                if (element != 0) {
                    return;
                }
            }

            throw new InvalidLinearEquationSystemException("One or more of the equations is null.");
        });

        this.coefficients = coefficients;
        this.constants = constants;
    }

    /**
     * Returns the number of unknowns in the system.
     *
     * @return the number of unknowns.
     */
    public int getUnknownCount() {
        return coefficients.getColumns();
    }

    /**
     * Returns the number of equations in the system.
     *
     * @return the number of equations.
     */
    public int getEquationCount() {
        return coefficients.getRows();
    }

    /**
     * Returns the type of linear equation system.
     *
     * @return the type of linear equation system.
     */
    public LinearEquationSystemType getType() {
        final int unknownCount = coefficients.getColumns();

        Matrix augmented = Matrix.getAugmentedMatrix(coefficients, constants);
        double coefficientsRange = coefficients.getRange();
        double augmentedRange = augmented.getRange();

        if (coefficientsRange == augmentedRange) {
            if (coefficientsRange == unknownCount) {
                return LinearEquationSystemType.UNIQUE_SOLUTION;
            } else {
                return LinearEquationSystemType.MULTIPLE_SOLUTIONS;
            }
        } else {
            return LinearEquationSystemType.INCONSISTENT;
        }
    }

    /**
     * Solve the linear equation system using the Gauss-Jordan elimination.
     *
     * @return a column matrix containing the solution to the system.
     * @see <a href="https://en.wikipedia.org/wiki/Gaussian_elimination" target="_top">Gaussian elimination in Wikipedia</a>
     * @exception InvalidLinearEquationSystemException if the system has none or multiple solutions.
     */
    public Matrix getSolution() {
        LinearEquationSystemType type = getType();

        if (type == LinearEquationSystemType.UNIQUE_SOLUTION) {
            Matrix augmented = Matrix.getAugmentedMatrix(coefficients, constants);
            Matrix reduced = augmented.getRowEchelon(true);
            double[] solutions = new double[getUnknownCount()];

            for (int i = 0; i < solutions.length; i++) {
                double solution = reduced.get(i, getUnknownCount());
                solution /= reduced.get(i, i);
                solutions[i] = solution;
            }

            return Matrix.getColumnMatrix(solutions);
        } else if (type == LinearEquationSystemType.MULTIPLE_SOLUTIONS) {
            throw new InvalidLinearEquationSystemException("The system has multiple solutions.");
        } else {
            throw new InvalidLinearEquationSystemException("The system doesn't have a solution");
        }
    }

    /**
     * Returns {@code true} if the system is a cramer system.
     *
     * @return {@code true} if this is a cramer system, otherwise {@code false}.
     */
    public boolean isCramerSystem() {
        Matrix augmented = Matrix.getAugmentedMatrix(coefficients, constants);
        return getEquationCount() == getUnknownCount() && augmented.isInvertible();
    }

    /**
     * Solve the linear equation system using the cramer matrix resolution.
     *
     * @return a column matrix containing the solution to the system.
     * @exception InvalidLinearEquationSystemException if the system is not a cramer system.
     */
    public Matrix getSolutionUsingCramer() {
        if (!isCramerSystem()) {
            throw new InvalidLinearEquationSystemException("The system is not a cramer system.");
        }

        return Matrix.multiply(coefficients.getInverse(Matrix.InverseFunction.DEFAULT), constants);
    }

    /**
     * Returns a representation of this system.
     *
     * @return a string representation of this system.
     */
    @Override
    public String toString() {
        StringBuilder b = new StringBuilder();

        for (int equation = 0; equation < getEquationCount(); equation++) {
            for (int unknown = 0; unknown < getUnknownCount(); unknown++) {
                double coefficient = coefficients.get(equation, unknown);
                if (coefficient == 0) {
                    b.append("   ").append(getUnknownCount() > 4 ? " " : "");
                } else {
                    if (Math.abs(coefficient) == 1) {
                        b.append(" ");
                        b.append(coefficient > 0 ? "+" : "-");
                    } else {
                        b.append(coefficient > 0 ? "+" : "");
                        b.append(trimZeroes(coefficient));
                    }

                    b.append(getUnknownName(unknown));
                }

                if (unknown - 1 < getUnknownCount()) {
                    b.append(" ");
                }
            }

            b.append(" = ");

            b.append(trimZeroes(constants.get(equation, 0)));

            b.append("\n");
        }

        return b.toString();
    }

    private String getUnknownName(int unknownIndex) {
        if (getUnknownCount() > 4) {
            return "x" + (unknownIndex + 1);
        } else {
            switch (unknownIndex) {
                case 0:
                    return "x";
                case 1:
                    return "y";
                case 2:
                    return "z";
                case 3:
                    return "w";
                default:
                    return "x";
            }
        }
    }

    /**
     * Types of linear equation systems.
     */
    public enum LinearEquationSystemType {
        /**
         * The system has one and only one solution.
         */
        UNIQUE_SOLUTION,

        /**
         * The system has more than one solution.
         */
        MULTIPLE_SOLUTIONS,

        /**
         * The system has no solutions.
         */
        INCONSISTENT
    }
}
