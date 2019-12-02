/*
 * www.javagl.de - Geom - Geometry utilities
 *
 * Copyright (c) 2013-2016 Marco Hutter - http://www.javagl.de
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */
package de.javagl.geom;

import java.awt.geom.Point2D;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * Simple implementation of an open Catmull-Rom-Spline
 */
public class CatmullRomSpline{

	private static final Function<Point2D, double[]> CONVERTOR = p -> new double[]{
		p.getX(),
		p.getY()
	};


	private static final Function<double[], Point2D> REVERSE_CONVERTOR = v -> {
		if(v.length != 2){
			throw new IllegalArgumentException("Vector of length " + v.length + " cannot be converted to " + Point2D.class.getSimpleName());
		}
		return new Point2D.Double(v[0], v[1]);
	};

	private static int dimensionality(List<double[]> values){
		final var dims = values.stream()
			.mapToInt(v -> v.length)
			.distinct()
			.toArray();

		if(dims.length != 1){
			throw new IllegalArgumentException("Dimensionality must be fixed across all points.");
		}

		return dims[0];
	}

	public static CatmullRomSpline create(
		List<? extends Point2D> points, int stepsPerSegment, double alpha){
		return create(points, stepsPerSegment, alpha, false);
	}

	public static CatmullRomSpline create(
		List<? extends Point2D> points,
		int stepsPerSegment,
		double alpha,
		boolean closed
	){
		return new CatmullRomSpline(
			points.stream().map(CONVERTOR).collect(Collectors.toList()),
			stepsPerSegment,
			alpha,
			closed
		);
	}

	private double alpha;
	private final int dimensionality;
	private final List<double[]> controlPoints;
	private final int stepsPerSegment;
	private final List<double[]> interpolatedPoints;

	private boolean updateRequired = true;
	private final boolean closed;

	private CatmullRomSpline(
		List<double[]> points,
		int stepsPerSegment,
		double alpha,
		boolean closed
	){
		this.dimensionality = dimensionality(points);
		this.stepsPerSegment = stepsPerSegment;
		this.alpha = alpha;
		int numInterpolatedPoints = (points.size() - 1) * stepsPerSegment + 1;
		if(closed){
			numInterpolatedPoints += stepsPerSegment;
			this.controlPoints = createPoints(points.size() + 3);
			this.controlPoints.set(1,
				controlPoints.get(controlPoints.size() - 2));
		}
		else{
			this.controlPoints = createPoints(points.size() + 2);
		}
		this.interpolatedPoints = createPoints(numInterpolatedPoints);
		this.closed = closed;
		updateControlPoints(points);
	}

	private List<double[]> createPoints(int n){
		return IntStream.range(0, n)
			.mapToObj(i -> new double[dimensionality])
			.collect(Collectors.toList());
	}

	public void setInterpolation(double alpha){
		this.alpha = alpha;
		updateRequired = true;
	}

	void updateControlPoint(int index, Point2D point){
		updateControlPoint(index, CONVERTOR.apply(point));
	}

	void updateControlPoint(int index, double[] point){
		int numPoints = controlPoints.size() - (closed ? 3 : 2);
		if(index < 0){
			throw new IndexOutOfBoundsException(
				"Index " + index + " must be positive");
		}
		if(index >= controlPoints.size() - 1){
			throw new IndexOutOfBoundsException(
				"Index was " + index + ", but number of control " +
					"points was " + numPoints);
		}
		double[] cp = controlPoints.get(index + 1);
		System.arraycopy(point, 0, cp, 0, dimensionality);
		updateRequired = true;
	}

	public void updateControlPoints(List<double[]> points){
		int numPoints = controlPoints.size() - (closed ? 3 : 2);
		if(points.size() != numPoints){
			throw new IllegalArgumentException(
				"Expected " + numPoints +
					" points, but got " + points.size());
		}
		for(int j = 0; j < points.size(); j++){
			double[] p = points.get(j);
			double[] cp = controlPoints.get(j + 1);
			System.arraycopy(p, 0, cp, 0, dimensionality);
		}
		updateRequired = true;
	}

	public List<Point2D> getInterpolatedPoints(){
		validatePoints();
		return Collections.unmodifiableList(
			interpolatedPoints.stream()
				.map(REVERSE_CONVERTOR)
				.collect(Collectors.toList())
		);
	}

	private void validatePoints(){
		if(updateRequired){
			updateAdditionalControlPoints();
			updateInterpolatedPoints();
			updateRequired = false;
		}
	}

	private void updateInterpolatedPoints(){
		int numPoints = controlPoints.size() - 2;
		for(int i = 0; i < numPoints - 1; i++){
			int stepsInCurrentSegment = stepsPerSegment;
			int lastStepInSegment = stepsInCurrentSegment;
			if(i == numPoints - 2){
				stepsInCurrentSegment++;
				lastStepInSegment = stepsInCurrentSegment - 1;
			}
			updateInterpolatedPoints(
				i, stepsInCurrentSegment, lastStepInSegment);
		}
	}

	private void copy(double[] src, double[] dst){
		System.arraycopy(src, 0, dst, 0, dimensionality);
	}
	private void updateAdditionalControlPoints(){
		if(closed){
			copy(
				controlPoints.get(controlPoints.size() - 3),
				controlPoints.get(0)
			);
			copy(
				controlPoints.get(2),
				controlPoints.get(controlPoints.size() - 1)
			);
		}
		else{
			double[] p0 = controlPoints.get(1);
			double[] p1 = controlPoints.get(2);
			double[] cp0 = controlPoints.get(0);
			apply(p1, p0, SUBTRACT, cp0);
			apply(p0, cp0, SUBTRACT, cp0);

			double[] py = controlPoints.get(controlPoints.size() - 3);
			double[] pz = controlPoints.get(controlPoints.size() - 2);
			double[] cpz = controlPoints.get(controlPoints.size() - 1);
			apply(pz, py, SUBTRACT, cpz);
			apply(pz, cpz, ADD, cpz);
		}

	}

	private static final DoubleToDoubleBiFunction SUBTRACT = (d0, d1) -> d0 - d1;
	private static final DoubleToDoubleBiFunction ADD = (d0, d1) -> d0 + d1;

	private interface DoubleToDoubleBiFunction{
		double applyAsDouble(double v0, double v1);
	}

	private double[] apply(double[] p0, double[] p1, DoubleToDoubleBiFunction function){
		final var rv = new double[dimensionality];
		apply(p0, p1, function, rv);
		return rv;
	}

	private void apply(double[] p0, double[] p1, DoubleToDoubleBiFunction function, double[] result){
		for(int i = 0; i < dimensionality; i++){
			result[i] = function.applyAsDouble(p0[i], p1[i]);
		}
	}

	private double magSquared(double[] value){
		double sum = 0;
		for(double d: value){
			sum += d * d;
		}
		return sum;
	}

	private void updateInterpolatedPoints(
		int index,
		int stepsInCurrentSegment,
		int lastStepInSegment
	){
		final double[] p0 = controlPoints.get(index + 0);
		final double[] p1 = controlPoints.get(index + 1);
		final double[] p2 = controlPoints.get(index + 2);
		final double[] p3 = controlPoints.get(index + 3);
		double t0 = 0;
		double t1 = 1;
		double t2 = 2;
		double t3 = 3;
		if(alpha != 0.0){
			double exponent = alpha * 0.5;
			double d01 = magSquared(apply(p1, p0, SUBTRACT));
			t1 = t0 + Math.pow(d01, exponent);

			double d12 = magSquared(apply(p2, p1, SUBTRACT));
			t2 = t1 + Math.pow(d12, exponent);

			double d23 = magSquared(apply(p3, p2, SUBTRACT));
			t3 = t2 + Math.pow(d23, exponent);

			// System.out.println("Times "+t0+" "+t1+" "+t2+" "+t3);
		}
		double invStep = 1.0 / lastStepInSegment;
		for(int i = 0; i < stepsInCurrentSegment; i++){
			double t = i * invStep;
			int interpolatedPointIndex = index * stepsPerSegment + i;
			double[] interpolatedPoint = interpolatedPoints.get(interpolatedPointIndex);
			copy(
				interpolate(p0, p1, p2, p3, t0, t1, t2, t3, t1 + t * (t2 - t1)),
				interpolatedPoint
			);
		}
	}

	private static void interpolate(
		Point2D p0,
		Point2D p1,
		Point2D p2,
		Point2D p3,
		double t0,
		double t1,
		double t2,
		double t3,
		double t,
		Point2D out
	){
		final double[] result = interpolate(
			CONVERTOR.apply(p0),
			CONVERTOR.apply(p1),
			CONVERTOR.apply(p2),
			CONVERTOR.apply(p3),
			t0,
			t1,
			t2,
			t3,
			t
		);
		out.setLocation(result[0], result[1]);
	}

	private static double[] interpolate(
		double[] p0,
		double[] p1,
		double[] p2,
		double[] p3,
		double t0,
		double t1,
		double t2,
		double t3,
		double t
	){
		final int dimensionality = dimensionality(Arrays.asList(p0, p1, p2, p3));

		double invDt01 = 1.0 / (t1 - t0);
		double invDt12 = 1.0 / (t2 - t1);
		double invDt23 = 1.0 / (t3 - t2);
		double f01a = (t1 - t) * invDt01;
		double f01b = (t - t0) * invDt01;
		double f12a = (t2 - t) * invDt12;
		double f12b = (t - t1) * invDt12;
		double f23a = (t3 - t) * invDt23;
		double f23b = (t - t2) * invDt23;
		double invDt02 = 1.0 / (t2 - t0);
		double invDt13 = 1.0 / (t3 - t1);
		double f012a = (t2 - t) * invDt02;
		double f012b = (t - t0) * invDt02;
		double f123a = (t3 - t) * invDt13;
		double f123b = (t - t1) * invDt13;

		final double[] rv = new double[dimensionality];

		for(int i = 0; i < rv.length; i++){
			final double v01 = f01a * p0[i] + f01b * p1[i];
			final double v12 = f12a * p1[i] + f12b * p2[i];
			final double v23 = f23a * p2[i] + f23b * p3[i];
			final double v012 = f012a * v01 + f012b * v12;
			final double v123 = f123a * v12 + f123b * v23;
			rv[i] = f12a * v012 + f12b * v123;
		}

		return rv;
	}


}