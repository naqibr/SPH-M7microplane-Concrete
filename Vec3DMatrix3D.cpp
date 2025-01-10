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
#include "Vector3D.h"
#include "Matrix3D.h"
#include <iostream>

double Vec3d::DotVec(const Vec3d& vec) const { return e[0] * vec.e[0] + e[1] * vec.e[1] + e[2] * vec.e[2]; }
Mat3d Vec3d::Dyadic(const Vec3d& vec) const { 
    Mat3d result;
    result.e[0][0] = e[0] * vec.e[0];result.e[0][1] = e[0] * vec.e[1];result.e[0][2] = e[0] * vec.e[2];
    result.e[1][0] = e[1] * vec.e[0];result.e[1][1] = e[1] * vec.e[1];result.e[1][2] = e[1] * vec.e[2];
    result.e[2][0] = e[2] * vec.e[0];result.e[2][1] = e[2] * vec.e[1];result.e[2][2] = e[2] * vec.e[2];
    return result; 
}

Vec3d Vec3d::Normal() const { return  *this / Length(); }
void Vec3d::Normalize() { *this = *this / Length(); }
double Vec3d::Length() const { return sqrtf(e[0] * e[0] + e[1] * e[1] + e[2] * e[2]); }

void Vec3d::setZero() { e[0] = 0.0e0; e[1] = 0.0e0; e[2] = 0.0e0; }



//Dot product of matrix with matrix, M3 = M . T
Mat3d Mat3d::DotMat(const Mat3d& vec) const
{
    return Mat3d(
        vec.e[0][0] * e[0][0] + vec.e[1][0] * e[0][1] + vec.e[2][0] * e[0][2], vec.e[0][1] * e[0][0] + vec.e[1][1] * e[0][1] + vec.e[2][1] * e[0][2], vec.e[0][2] * e[0][0] + vec.e[1][2] * e[0][1] + vec.e[2][2] * e[0][2],
        vec.e[0][0] * e[1][0] + vec.e[1][0] * e[1][1] + vec.e[2][0] * e[1][2], vec.e[0][1] * e[1][0] + vec.e[1][1] * e[1][1] + vec.e[2][1] * e[1][2], vec.e[0][2] * e[1][0] + vec.e[1][2] * e[1][1] + vec.e[2][2] * e[1][2],
        vec.e[0][0] * e[2][0] + vec.e[1][0] * e[2][1] + vec.e[2][0] * e[2][2], vec.e[0][1] * e[2][0] + vec.e[1][1] * e[2][1] + vec.e[2][1] * e[2][2], vec.e[0][2] * e[2][0] + vec.e[1][2] * e[2][1] + vec.e[2][2] * e[2][2]);
}

//Transpose of matrix
Mat3d Mat3d::transpose() const
{
    return Mat3d(
        e[0][0], e[1][0], e[2][0],
        e[0][1], e[1][1], e[2][1],
        e[0][2], e[1][2], e[2][2]);
}

//Inverse of matrix
Mat3d Mat3d::inverse() const
{
    double det = (e[0][2] * e[1][1] * e[2][0] - e[0][1] * e[1][2] * e[2][0] - e[0][2] * e[1][0] * e[2][1] + e[0][0] * e[1][2] * e[2][1] + e[0][1] * e[1][0] * e[2][2] - e[0][0] * e[1][1] * e[2][2]);
    if (fabs(det) < INF) {
        std::cout << "\n\nInverse Doesnt Exist!!\n";
        exit(0);
        return Mat3d();
    } else
    return 1 / det *
        Mat3d(
            e[1][2] * e[2][1] - e[1][1] * e[2][2], 0.0e0 - e[0][2] * e[2][1] + e[0][1] * e[2][2], e[0][2] * e[1][1] - e[0][1] * e[1][2],
            0.0e0 - e[1][2] * e[2][0] + e[1][0] * e[2][2], e[0][2] * e[2][0] - e[0][0] * e[2][2], 0.0e0 - e[0][2] * e[1][0] + e[0][0] * e[1][2],
            e[1][1] * e[2][0] - e[1][0] * e[2][1], 0.0e0 - e[0][1] * e[2][0] + e[0][0] * e[2][1], e[0][1] * e[1][0] - e[0][0] * e[1][1]);
}

//Trace of matrix, Ajj
double Mat3d::trace() {
    return (e[0][0] + e[1][1] + e[2][2]);
}

//Determinant of matrix
double Mat3d::det() const
{
    return 0.0e0 - (e[0][2] * e[1][1] * e[2][0] - e[0][1] * e[1][2] * e[2][0] - e[0][2] * e[1][0] * e[2][1] + e[0][0] * e[1][2] * e[2][1] + e[0][1] * e[1][0] * e[2][2] - e[0][0] * e[1][1] * e[2][2]);
}

//Dot product of matrix with vecotr, M3 = M.v
Vec3d Mat3d::DotVec(const Vec3d& vec) const
{
    return Vec3d(
        e[0][0] * vec.e[0] + e[0][1] * vec.e[1] + e[0][2] * vec.e[2],
        e[1][0] * vec.e[0] + e[1][1] * vec.e[1] + e[1][2] * vec.e[2],
        e[2][0] * vec.e[0] + e[2][1] * vec.e[1] + e[2][2] * vec.e[2]);
}

//Set elements of matrix to zero
void Mat3d::setZero() {
    e[0][0] = 0.0e0; e[0][1] = 0.0e0; e[0][2] = 0.0e0;
    e[1][0] = 0.0e0; e[1][1] = 0.0e0; e[1][2] = 0.0e0;
    e[2][0] = 0.0e0; e[2][1] = 0.0e0; e[2][2] = 0.0e0;
}