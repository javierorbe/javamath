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

package com.javierorbe.math.util;

import com.javierorbe.math.matrix.Matrix;

import java.util.*;

/**
 * Util math functions.
 *
 * @author Javier Orbe
 */
public final class MathUtils {

	private static final Random random = new Random();

	private MathUtils() {}

	/**
	 * Get a random {@code double} in a range.
	 *
	 * @param min lower range bounds.
	 * @param max upper range bounds.
	 * @return a random number in a range.
	 */
	public static double random(int min, int max) {
		return min + (max - min) * random.nextDouble();
	}

	/**
	 * Get a random {@code int} in a range.
	 *
	 * @param min lower range bounds.
	 * @param max upper range bounds.
	 * @return a random number in a range.
	 */
	public static int randomInt(int min, int max) {
		return random.nextInt((max - min) + 1) + min;
	}

	/**
	 * Generate a random number from a Gaussian distribution.
	 *
	 * @return the random number.
	 * @see <a href="https://en.wikipedia.org/wiki/Normal_distribution" target="_top">Normal distribution in Wikipedia</a>
	 */
	public static double randomGaussian() {
		double x1;
		double x2;
		double rad;

		do {
			x1 = (2 * Math.random()) - 1;
			x2 = (2 * Math.random()) - 1;
			rad = (x1 * x1) + (x2 * x2);
		} while (rad >= 1 || rad == 0);
		
		double c = Math.sqrt((-2 * Math.log(rad)) / rad);
		return x1 * c;
	}

    /**
     * Map a number that is given in some bounds to a value in other bounds.
     *
     * @param number the number.
     * @param inMin the input lower bound.
     * @param inMax the input upper bound.
     * @param outMin the output lower bound.
     * @param outMax the output upper bound.
     * @return the mapped number.
     */
	public static double mapNumber(double number, double inMin, double inMax, double outMin, double outMax) {
		return (((number - inMin) * (outMax - outMin)) / (inMax - inMin)) + outMin;
	}

	/**
	 * Returns a number whose value is limited to the given range.
	 *
	 * @param number the value.
	 * @param min minimum value.
	 * @param max maximum value.
	 * @return the value limited to the range.
	 */
	public static double clamp(double number, double min, double max) {
		if (number > max) {
			return max;
		}

		if (number < min) {
			return min;
		}

		return number;
	}

	/**
	 * Returns a deep copy of a two-dimensional array.
	 *
	 * @param arr the array to copy.
	 * @return the copy of the array.
	 */
	public static double[][] deepCopy(double[][] arr) {
	    if (arr.length == 0) {
	        throw new IllegalArgumentException("The array cannot be empty.");
        }

	    double[][] result = new double[arr.length][arr[0].length];

	    for (int i = 0; i < arr.length; i++) {
	        result[i] = Arrays.copyOf(arr[i], arr[i].length);
        }

        return result;
    }

	/**
	 * Returns the string value of the specified, with the decimal point zeroes trimmed.
	 *
	 * @param number the number to trim.
	 * @return the trimmed number.
	 */
	public static String trimZeroes(double number) {
		String s = String.valueOf(number);
		return !s.contains(".") ? s : s.replaceAll("0*$", "").replaceAll("\\.$", "");
	}

	/**
	 * Sort a map entry list in ascending order by its values.
	 *
	 * @param map the map to sort.
	 * @param <K> key type.
	 * @param <V> value type.
	 * @return a sorted list of map entries.
	 */
	public static <K, V extends Comparable<? super V>> List<Map.Entry<K, V>> entriesSortedByValues(Map<K, V> map) {
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
	public static List<Map.Entry<Integer, Integer>> getSortedZeroCount(Matrix matrix) {
		Map<Integer, Integer> zeroCount = new HashMap<>();
		for (int row = 0; row < matrix.getRows(); row++) {
			zeroCount.put(row, 0);

			for (int col = 0; col < matrix.getColumns(); col++) {
				double element = matrix.get(row, col);

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

	/**
	 * Get the index of the pivot in a row.
	 * The pivot is the first none zero element.
	 *
	 * @param row the row.
	 * @return the index of the pivot if there is, otherwise {@code -1}.
	 */
	public static int getPivotIndex(double[] row) {
		int pivotIndex = -1;

		for (int i = 0; i < row.length; i++) {
			if (row[i] != 0) {
				pivotIndex = i;
				break;
			}
		}

		return pivotIndex;
	}
}
