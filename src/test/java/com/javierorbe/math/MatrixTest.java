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

import com.javierorbe.math.util.MathUtils;
import org.jboss.arquillian.container.test.api.Deployment;
import org.jboss.arquillian.junit.Arquillian;
import org.jboss.shrinkwrap.api.ShrinkWrap;
import org.jboss.shrinkwrap.api.asset.EmptyAsset;
import org.jboss.shrinkwrap.api.spec.JavaArchive;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;

import java.util.Arrays;

import static org.junit.Assert.*;

@RunWith(Arquillian.class)
public class MatrixTest {

    @Deployment
    public static JavaArchive createDeployment() {
        return ShrinkWrap.create(JavaArchive.class)
                .addClass(Matrix.class)
                .addAsManifestResource(EmptyAsset.INSTANCE, "beans.xml");
    }

    @Test
    public void copy() {
    }

    @Test
    public void map() {
    }

    @Test
    public void map1() {
    }

    @Test
    public void map2() {
    }

    @Test
    public void randomize() {
    }

    @Test
    public void addScalar() {
        double[][] matrixValues = new double[][]{
                {-3, 2, 1},
                {1,  0, 2},
                {3,  4, 5}
        };
        double scalar = -3;

        Matrix m = new Matrix(MathUtils.deepCopy(matrixValues));
        m.add(scalar);

        m.forEach((value, index) -> Assert.assertEquals(matrixValues[index.i][index.j] + scalar, value.doubleValue(), 0.0001));
    }

    @Test
    public void addMatrix() {
        Matrix a = new Matrix(new double[][]{
                {-3, 2, 1},
                {1,  0, 2},
                {3,  4, 5}
        });
        Matrix b = new Matrix(new double[][]{
                {-3, 6, -3},
                {3,  -4, -1},
                {-2,  4,  6}
        });

        Matrix expected = new Matrix(new double[][] {
                {-6, 8, -2},
                {4, -4, 1},
                {1, 8, 11}
        });

        a.add(b);

        Assert.assertEquals(expected, a);
    }

    @Test
    public void subtract() {
    }

    @Test
    public void multiply() {
    }

    @Test
    public void multiply1() {
    }

    @Test
    public void multiply2() {
    }

    @Test
    public void subtractRow() {
    }

    @Test
    public void subtractRow1() {
    }

    @Test
    public void subtractColumn() {
    }

    @Test
    public void multiplyRow() {
    }

    @Test
    public void multiplyColumn() {
    }

    @Test
    public void swapRows() {
        Matrix m = new Matrix(new double[][]{
                {-3, 2, 1},
                {1,  0, 2},
                {3,  4, 5}
        });

        m.swapRows(0, 1);
        m.swapRows(1, 2);

        Matrix result = new Matrix(new double[][]{
                {1,  0, 2},
                {3,  4, 5},
                {-3, 2, 1}
        });

        Assert.assertEquals(m, result);
    }

    @Test
    public void swapColumns() {
        Matrix m = new Matrix(new double[][]{
                {-3, 2, 1},
                {1,  0, 2},
                {3,  4, 5}
        });

        m.swapColumns(0, 2);
        m.swapColumns(1, 2);

        Matrix result = new Matrix(new double[][]{
                {1, -3, 2},
                {2, 1,  0},
                {5, 3,  4}
        });

        Assert.assertEquals(m, result);
    }

    @Test
    public void isSquare() {
        Matrix m1 = new Matrix(new double[][]{
                {-3, 2, 1},
                {1,  0, 2},
                {3,  4, 5}
        });
        Assert.assertTrue(m1.isSquare());

        Matrix m2 = new Matrix(new double[][]{
                {-3, 2, 1, 5},
                {1,  0, 2, 3},
                {3,  4, 5, 1}
        });
        Assert.assertFalse(m2.isSquare());
    }

    @Test
    public void isIdentity() {
        Matrix m1 = new Matrix(new double[][]{
                {1, 0, 0, 0},
                {0, 1, 0, 0},
                {0, 0, 1, 0},
                {0, 0, 0, 1}
        });
        Assert.assertTrue(m1.isIdentity());

        Matrix m2 = new Matrix(new double[][]{
                {-3, 2, 1},
                {1,  0, 2},
                {3,  4, 5}
        });
        Assert.assertFalse(m2.isIdentity());
    }

    @Test
    public void getIdentity() {
    }

    @Test
    public void isZero() {
        Matrix m1 = new Matrix(new double[][]{
                {0, 0, 0, 0},
                {0, 0, 0, 0},
                {0, 0, 0, 0}
        });
        Assert.assertTrue(m1.isZero());

        Matrix m2 = new Matrix(new double[][]{
                {-3, 2},
                {0,  0},
                {0,  4}
        });
        Assert.assertFalse(m2.isZero());
    }

    @Test
    public void getRowEchelon() {

    }

    @Test
    public void isInvertible() {
    }

    @Test
    public void getInverse() {
    }

    @Test
    public void getMinor() {
    }

    @Test
    public void getAdjugate() {
    }

    @Test
    public void getAdjugateMatrix() {

    }

    @Test
    public void solveSystem() {
    }

    @Test
    public void solveSystemByCramer() {
    }

    @Test
    public void getRange() {
        Matrix m1 = new Matrix(new double[][]{
                {1, 2, -1, 3, -2},
                {2, 1,  0, 1,  1},
                {2, 4, -2, 6, -4},
                {0, 0,  0, 0,  0},
                {5, 4, -1, 5,  0}

        });
        Assert.assertEquals(2, m1.getRange());

        Matrix m2 = new Matrix(new double[][]{
                {2, -1, 0, 7},
                {1,  0, 1, 3},
                {3,  2, 7, 7}
        });
        Assert.assertEquals(2, m2.getRange());

        Matrix m3 = new Matrix(new double[][]{
                {2,  1, 3,  2},
                {3,  2, 5,  1},
                {-1, 1, 0, -7},
                {3, -2, 1, 17},
                {0,  1, 1,  -4}
        });
        Assert.assertEquals(2, m3.getRange());
    }

    @Test
    public void getDeterminant() {
    }

    @Test
    public void getDeterminantByLaplace() {
    }

    @Test
    public void transpose() {

    }

    @Test
    public void toArray() {
    }

    @Test
    public void getColumnMatrix() {
        double[] values = new double[] {1, 2, 3, 4, 5};

        Matrix result = Matrix.getColumnMatrix(values);

        Matrix expected = new Matrix(new double[][]{
                {1},
                {2},
                {3},
                {4},
                {5}
        });

        Assert.assertEquals(expected, result);
    }
}
