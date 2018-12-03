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

import org.jboss.arquillian.container.test.api.Deployment;
import org.jboss.arquillian.junit.Arquillian;
import org.jboss.shrinkwrap.api.ShrinkWrap;
import org.jboss.shrinkwrap.api.asset.EmptyAsset;
import org.jboss.shrinkwrap.api.spec.JavaArchive;
import org.junit.Assert;
import org.junit.Test;
import org.junit.runner.RunWith;

import static org.junit.Assert.*;

@RunWith(Arquillian.class)
public class LinearEquationSystemTest {

    private static final double DOUBLE_DELTA = 0.0001;

    @Deployment
    public static JavaArchive createDeployment() {
        return ShrinkWrap.create(JavaArchive.class)
                .addClass(LinearEquationSystem.class)
                .addAsManifestResource(EmptyAsset.INSTANCE, "beans.xml");
    }

    @Test
    public void getUnknownCount() {
        LinearEquationSystem system = new LinearEquationSystem(
                new Matrix(new double[][]{
                        {0, 1, 1},
                        {-2, -3, -1},
                        {1, 1, 3}
                }),
                Matrix.getColumnMatrix(new double[] {1, 1, 5})
        );

        Assert.assertEquals(3, system.getUnknownCount());
    }

    @Test
    public void getType() {
        LinearEquationSystem system = new LinearEquationSystem(
                new Matrix(new double[][]{
                        {0, 1, 1},
                        {-2, -3, -1},
                        {1, 1, 3}
                }),
                Matrix.getColumnMatrix(new double[] {1, 1, 5})
        );

        Assert.assertEquals(LinearEquationSystem.LinearEquationSystemType.UNIQUE_SOLUTION, system.getType());
    }

    @Test
    public void getSolution() {
        LinearEquationSystem system = new LinearEquationSystem(
                new Matrix(new double[][]{
                        {0, 1, 1},
                        {-2, -3, -1},
                        {1, 1, 3}
                }),
                Matrix.getColumnMatrix(new double[] {1, 1, 5})
        );
        Matrix solution = system.getSolution();

        Assert.assertEquals(0, solution.get(0, 0), DOUBLE_DELTA);
        Assert.assertEquals(-1, solution.get(1, 0), DOUBLE_DELTA);
        Assert.assertEquals(2, solution.get(2, 0), DOUBLE_DELTA);
    }
}
