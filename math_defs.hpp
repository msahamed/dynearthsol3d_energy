#ifndef DYNEARTHSOL3D_MATH_DEFS_HPP
#define DYNEARTHSOL3D_MATH_DEFS_HPP

#include "constants.hpp"

#ifdef USE_EIGEN
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

// Define Eigen types for 2D and 3D
using Vec2 = Eigen::Vector2d;
using Vec3 = Eigen::Vector3d;
using Mat2 = Eigen::Matrix2d;
using Mat3 = Eigen::Matrix3d;

// Tensor type depends on dimension
#if NDIMS == 2
    using Vector = Vec2;
    using Tensor = Mat2;
#else
    using Vector = Vec3;
    using Tensor = Mat3;
#endif

// Helper to map raw pointers to Eigen types
template<typename T>
using MapVec = Eigen::Map<T>;

template<typename T>
using MapMat = Eigen::Map<T>;

#endif // USE_EIGEN

#endif // DYNEARTHSOL3D_MATH_DEFS_HPP
