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

/**
 * Represents a mathematical vector of dimension 2.
 *
 * @author Javier Orbe
 */
public class Vector2 {

    private double x;
    private double y;

    /**
     * Create a vector with the specified coordinates.
     *
     * @param x the x coordinate.
     * @param y the y coordinate.
     */
    public Vector2(double x, double y) {
        this.x = x;
        this.y = y;
    }

    /**
     * Create a vector with its coordinates set to zero.
     */
    public Vector2() {
        this(0, 0);
    }

    /**
     * Returns the value of the x coordinate.
     *
     * @return the value of the x coordinate.
     */
    public double getX() {
        return x;
    }

    /**
     * Set the value of the x coordinate.
     *
     * @param x the new value of the x coordinate.
     */
    public void setX(double x) {
        this.x = x;
    }

    public void addX(double x) {
        this.x += x;
    }

    /**
     * Returns the value of the y coordinate.
     *
     * @return the value of the y coordinate.
     */
    public double getY() {
        return y;
    }

    /**
     * Set the value of the y coordinate.
     *
     * @param y the new value of the y coordinate.
     */
    public void setY(double y) {
        this.y = y;
    }

    public void addY(double y) {
        this.y += y;
    }

    /**
     * Add the values of a vector to this vector.
     *
     * @param vector the vector.
     */
    public void add(Vector2 vector) {
        x += vector.x;
        y += vector.y;
    }

    /**
     * Set the vector values to zero.
     */
    public void zero() {
        x = 0;
        y = 0;
    }

    public void multiply(double scalar) {
        x *= scalar;
        y *= scalar;
    }

    /**
     * Returns the length of the vector squared.
     * @return the length of the vector squared.
     */
    public double squareLength() {
        return x * x + y * y;
    }

    /**
     * Returns the length of the vector.
     * @return the length of the vector.
     */
    public double length() {
        return Math.sqrt(squareLength());
    }

    /**
     * Normalizes the vector.
     */
    public void normalize() {
        double len = length();
        x /= len;
        y /= len;
    }

    /**
     * Randomizes the vector coordinates with values from 0 to 1.
     */
    public void randomize() {
        x = Math.random();
        y = Math.random();
    }

    /**
     * Compares the specified object with this vector for equality.
     * Returns {@code true} if and only if the specified object is also a vector,
     * and both vectors have the same x and y coordinate values.
     *
     * @param o the object to be compared for equality with this vector.
     * @return {@code true} if the specified object is equal to this vector, otherwise {@code false}.
     */
    @Override
    public boolean equals(Object o) {
        if (o == this) {
            return true;
        }

        if (!(o instanceof Vector2)) {
            return false;
        }

        Vector2 v = (Vector2) o;

        return Double.compare(x, v.x) == 0 && Double.compare(y, v.y) == 0;
    }

    /**
     * Returns a representation of this vector.
     *
     * @return a string representation of this vector.
     */
    @Override
    public String toString() {
        return "(" + x + ", " + y + ")";
    }
}
