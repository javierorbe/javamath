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
 * Thrown to indicate that the application has attempted to create or use a matrix
 * that doesn't fulfill the requirements or doesn't have the needed properties.
 *
 * @author Javier Orbe
 */
class InvalidMatrixException extends RuntimeException {

    /**
     * Constructs a {@code InvalidMatrixException} with no detail message.
     */
    InvalidMatrixException() {
        super();
    }

    /**
     * Constructs a {@code InvalidMatrixException} with the specified detail message.
     *
     * @param message the detail message.
     */
    InvalidMatrixException(String message) {
        super(message);
    }
}
