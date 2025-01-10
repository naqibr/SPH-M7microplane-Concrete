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
#ifndef __VECTOR3D_H__
#define __VECTOR3D__

#include <math.h>

typedef unsigned int uint;

class Mat3d; // Forward declaration of Mat3d

class Vec3d
{
public:
	double e[3];			//component i;

public:
	Vec3d() {
		e[0] = 0.0;
		e[1] = 0.0;
		e[2] = 0.0;
	}

	Vec3d(double _x, double _y, double _z) { e[0] = _x; e[1] = _y; e[2] = _z;};

	Vec3d(Vec3d const& vec) { e[0] = vec.e[0]; e[1] = vec.e[1]; e[2] = vec.e[2];};

	Vec3d operator + (const Vec3d& vec) const { return Vec3d(e[0] + vec.e[0], e[1] + vec.e[1], e[2] + vec.e[2]); }
	Vec3d operator - (const Vec3d& vec) const { return Vec3d(e[0] - vec.e[0], e[1] - vec.e[1], e[2] - vec.e[2]); }
	Vec3d operator * (const Vec3d& vec) const { return Vec3d(e[0] * vec.e[0], e[1] * vec.e[1], e[2] * vec.e[2]); }
	Vec3d operator / (const Vec3d& vec) const { return Vec3d(e[0] / vec.e[0], e[1] / vec.e[1], e[2] / vec.e[2]); }

	friend  Vec3d operator + (const Vec3d& vec, double n) { return Vec3d(vec.e[0] + n, vec.e[1] + n, vec.e[2] + n); }
	friend  Vec3d operator - (const Vec3d& vec, double n) { return Vec3d(vec.e[0] - n, vec.e[1] - n, vec.e[2] - n); }
	friend  Vec3d operator * (const Vec3d& vec, double n) { return Vec3d(vec.e[0] * n, vec.e[1] * n, vec.e[2] * n); }
	friend  Vec3d operator / (const Vec3d& vec, double n) { return Vec3d(vec.e[0] / n, vec.e[1] / n, vec.e[2] / n); }

	friend  Vec3d operator + (double n, const Vec3d& vec) { return Vec3d(n + vec.e[0], n + vec.e[1], n + vec.e[2]); }
	friend  Vec3d operator - (double n, const Vec3d& vec) { return Vec3d(n - vec.e[0], n - vec.e[1], n - vec.e[2]); }
	friend  Vec3d operator * (double n, const Vec3d& vec) { return Vec3d(n * vec.e[0], n * vec.e[1], n * vec.e[2]); }
	friend  Vec3d operator / (double n, const Vec3d& vec) { return Vec3d(n / vec.e[0], n / vec.e[1], n / vec.e[2]); }

	friend  Vec3d operator - (const Vec3d& vec) { return 0.0 - vec  ; }

	double DotVec(const Vec3d& vec) const;
	Mat3d Dyadic(const Vec3d& vec) const;

	Vec3d Normal() const;
	void Normalize();
	double Length() const;

	void setZero();
};

class Vec3i
{
public:
	int e[3];

public:
	Vec3i() {}
	Vec3i(int _x, int _y, int _z) { e[0] = _x; e[1] = _y; e[2] = _z;}
	Vec3i(Vec3i& vec) { e[0] = vec.e[0]; e[1] = vec.e[1]; e[2] = vec.e[2];}

	void setZero() { e[0] = 0; e[1] = 0; e[2] = 0;}
};

#endif
