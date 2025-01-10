/* 
 * This program is part of the following research paper:
 * 
 * [1] M. N. Rahimi and G. Moutsanidis, "Modeling concrete failure with smoothed particle hydrodynamics using the Microplane (M7) constitutive model" (Journal)
 * 
 * It also implements the Total Lagrangian Smoothed Particle Hydrodynamics formulations presented in:
 * 
 * [2] M. N. Rahimi and G. Moutsanidis, “A smoothed particle hydrodynamics approach for phase field modeling of brittle fracture,” 
 *     Computer Methods in Applied Mechanics and Engineering, Aug. 2022. doi: https://doi.org/10.1016/j.cma.2022.115191
 * 
 * [3] M. N. Rahimi and G. Moutsanidis, “Modeling dynamic brittle fracture in functionally graded materials using hyperbolic phase field and smoothed particle hydrodynamics,” 
 *     Computer Methods in Applied Mechanics and Engineering, Nov. 2022. doi: https://doi.org/10.1016/j.cma.2022.115642
 * 
 * [4] M. N. Rahimi and G. Moutsanidis, “An SPH-based FSI framework for phase-field modeling of brittle fracture under extreme hydrodynamic events,” 
 *     Engineering with Computers, Aug. 2023. doi: https://doi.org/10.1007/s00366-023-01857-0
 * 
 * [5] M. N. Rahimi and G. Moutsanidis, “IGA-SPH: Coupling isogeometric analysis with smoothed particle hydrodynamics for air-blast–structure interaction,” 
 *     Engineering with Computers, May 2024. doi: https://doi.org/10.1007/s00366-024-01978-0
 * 
 * If this code is helpfull, whether partially or entirely, in your research, please cite the relevant paper(s).
 * 
 * ============================================================================
 */
#pragma once
#ifndef __Matrix3D_H__
#define __Matrix3D__

#include <math.h>

#ifndef INF
#define INF 1E-20
#endif

class Vec3d; // Forward declaration of Vec3d
typedef unsigned int uint;

class Mat3d
{
public:
	double e[3][3];		//Index[i][j]
public:
	Mat3d() {
		for (int i = 0; i < 3; i++) {
			for (int j = 0; j < 3; j++) {
				e[i][j] = 0.0e0;
			}
		}
	}

	Mat3d(
		double _x11, double _x12, double _x13,
		double _x21, double _x22, double _x23,
		double _x31, double _x32, double _x33)
	{
		e[0][0] = _x11; e[0][1] = _x12; e[0][2] = _x13;
		e[1][0] = _x21; e[1][1] = _x22; e[1][2] = _x23;
		e[2][0] = _x31; e[2][1] = _x32; e[2][2] = _x33;
	};

	Mat3d(Mat3d const& vec)
	{
		e[0][0] = vec.e[0][0]; e[0][1] = vec.e[0][1]; e[0][2] = vec.e[0][2];
		e[1][0] = vec.e[1][0]; e[1][1] = vec.e[1][1]; e[1][2] = vec.e[1][2];
		e[2][0] = vec.e[2][0]; e[2][1] = vec.e[2][1]; e[2][2] = vec.e[2][2];
	};

	Mat3d operator + (const Mat3d& vec) const
	{
		return Mat3d(
			e[0][0] + vec.e[0][0], e[0][1] + vec.e[0][1], e[0][2] + vec.e[0][2],
			e[1][0] + vec.e[1][0], e[1][1] + vec.e[1][1], e[1][2] + vec.e[1][2],
			e[2][0] + vec.e[2][0], e[2][1] + vec.e[2][1], e[2][2] + vec.e[2][2]);
	}

	Mat3d operator - (const Mat3d& vec) const
	{
		return  Mat3d(
			e[0][0] - vec.e[0][0], e[0][1] - vec.e[0][1], e[0][2] - vec.e[0][2],
			e[1][0] - vec.e[1][0], e[1][1] - vec.e[1][1], e[1][2] - vec.e[1][2],
			e[2][0] - vec.e[2][0], e[2][1] - vec.e[2][1], e[2][2] - vec.e[2][2]);
	}

	Mat3d operator * (const Mat3d& vec) const
	{
		return  Mat3d(
			e[0][0] * vec.e[0][0], e[0][1] * vec.e[0][1], e[0][2] * vec.e[0][2],
			e[1][0] * vec.e[1][0], e[1][1] * vec.e[1][1], e[1][2] * vec.e[1][2],
			e[2][0] * vec.e[2][0], e[2][1] * vec.e[2][1], e[2][2] * vec.e[2][2]);
	}

	Mat3d operator / (const Mat3d& vec) const
	{
		return  Mat3d(
			e[0][0] / vec.e[0][0], e[0][1] / vec.e[0][1], e[0][2] / vec.e[0][2],
			e[1][0] / vec.e[1][0], e[1][1] / vec.e[1][1], e[1][2] / vec.e[1][2],
			e[2][0] / vec.e[2][0], e[2][1] / vec.e[2][1], e[2][2] / vec.e[2][2]);
	}

	friend  Mat3d operator + (const Mat3d& vec, double n)
	{
		return Mat3d(
			vec.e[0][0] + n, vec.e[0][1] + n, vec.e[0][2] + n,
			vec.e[1][0] + n, vec.e[1][1] + n, vec.e[1][2] + n,
			vec.e[2][0] + n, vec.e[2][1] + n, vec.e[2][2] + n);
	}

	friend  Mat3d operator - (const Mat3d& vec, double n)
	{
		return Mat3d(
			vec.e[0][0] - n, vec.e[0][1] - n, vec.e[0][2] - n,
			vec.e[1][0] - n, vec.e[1][1] - n, vec.e[1][2] - n,
			vec.e[2][0] - n, vec.e[2][1] - n, vec.e[2][2] - n);
	}

	friend  Mat3d operator * (const Mat3d& vec, double n)
	{
		return Mat3d(
			vec.e[0][0] * n, vec.e[0][1] * n, vec.e[0][2] * n,
			vec.e[1][0] * n, vec.e[1][1] * n, vec.e[1][2] * n,
			vec.e[2][0] * n, vec.e[2][1] * n, vec.e[2][2] * n);
	}

	friend  Mat3d operator / (const Mat3d& vec, double n)
	{
		return Mat3d(
			vec.e[0][0] / n, vec.e[0][1] / n, vec.e[0][2] / n,
			vec.e[1][0] / n, vec.e[1][1] / n, vec.e[1][2] / n,
			vec.e[2][0] / n, vec.e[2][1] / n, vec.e[2][2] / n);
	}

	friend  Mat3d operator + (double n, const Mat3d& vec)
	{
		return Mat3d(
			n + vec.e[0][0], n + vec.e[0][1], n + vec.e[0][2],
			n + vec.e[1][0], n + vec.e[1][1], n + vec.e[1][2],
			n + vec.e[2][0], n + vec.e[2][1], n + vec.e[2][2]);
	}

	friend  Mat3d operator - (double n, const Mat3d& vec)
	{
		return Mat3d(
			n - vec.e[0][0], n - vec.e[0][1], n - vec.e[0][2],
			n - vec.e[1][0], n - vec.e[1][1], n - vec.e[1][2],
			n - vec.e[2][0], n - vec.e[2][1], n - vec.e[2][2]);
	}

	friend  Mat3d operator * (double n, const Mat3d& vec)
	{
		return Mat3d(
			n * vec.e[0][0], n * vec.e[0][1], n * vec.e[0][2],
			n * vec.e[1][0], n * vec.e[1][1], n * vec.e[1][2],
			n * vec.e[2][0], n * vec.e[2][1], n * vec.e[2][2]);
	}

	friend  Mat3d operator / (double n, const Mat3d& vec)
	{
		return Mat3d(
			n / vec.e[0][0], n / vec.e[0][1], n / vec.e[0][2],
			n / vec.e[1][0], n / vec.e[1][1], n / vec.e[1][2],
			n / vec.e[2][0], n / vec.e[2][1], n / vec.e[2][2]);
	}

	friend  Mat3d operator - (const Mat3d& vec) { return 0.0 - vec; }

	//Dot product of matrix with matrix, M3 = M . T
	Mat3d DotMat(const Mat3d& vec) const;

	//Transpose of matrix
	Mat3d transpose() const;

	//Inverse of matrix
	Mat3d inverse() const;

	//Trace of matrix, Ajj
	double trace();
	//Determinant of matrix
	double det() const;
	//Dot product of matrix with vecotr, M3 = M.v
	Vec3d DotVec(const Vec3d& vec) const;

	//Set elements of matrix to zero
	void setZero();
};

#endif
