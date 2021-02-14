#ifndef ZMATH_QUATERNION_H_
#define ZMATH_QUATERNION_H_

#ifdef _WIN32
#if !defined(_USE_MATH_DEFINES)
#define _USE_MATH_DEFINES  // For M_PI.
#endif                     // !defined(_USE_MATH_DEFINES)
#endif                     // _WIN32

#include "matrix.h"
#include "vector.h"
#include <math.h>

/// ZMATH provides a Quaternion class that utilizes SIMD optimized
/// Matrix and Vector classes.

namespace zmath {

/// @addtogroup ZMATH_quaternion
/// @{
/// @class Quaternion
///
/// @brief Stores a Quaternion of type T and provides a set of utility
/// operations on each Quaternion.
/// @tparam T Type of each element in the Quaternion.
template <class T>
class Quaternion {
 public:
  /// @brief Construct an uninitialized Quaternion.
  inline Quaternion() {}

  /// @brief Construct a Quaternion from a copy.
  /// @param q Quaternion to copy.
  inline Quaternion(const Quaternion<T>& q) {
    s_ = q.s_;
    v_ = q.v_;
  }

  /// @brief Construct a Quaternion using scalar values to initialize each
  /// element.
  ///
  /// @param s1 Scalar component.
  /// @param qs1 First element of the Vector component.
  /// @param qs2 Second element of the Vector component.
  /// @param qs3 Third element of the Vector component.
  inline Quaternion(T s1, T qs1, T qs2, T qs3) {
    s_ = s1;
    v_ = Vector<T, 3>(qs1, qs2, qs3);
  }

  /// @brief Construct a quaternion from a scalar and 3-dimensional Vector.
  ///
  /// @param s1 Scalar component.
  /// @param v1 Vector component.
  inline Quaternion(T s1, const Vector<T, 3>& v1) {
    s_ = s1;
    v_ = v1;
  }

  /// @brief Return the scalar component of the quaternion.
  ///
  /// @return The scalar component
  inline const T& scalar() const { return s_; }

  /// @brief Set the scalar component of the quaternion.
  ///
  /// @param s Scalar component.
  inline void set_scalar(T s) { s_ = s; }

  /// @brief Return the vector component of the quaternion.
  ///
  /// @return The vector component
  inline const Vector<T, 3>& vector() const { return v_; }

  /// @brief Set the vector component of the quaternion.
  ///
  /// @param v Vector component.
  inline void set_vector(const Vector<T, 3>& v) { v_ = v; }

  /// @brief Calculate the inverse Quaternion.
  ///
  /// This calculates the inverse such that <code>(q * q).Inverse()</code>
  /// is the identity.
  ///
  /// @return Quaternion containing the result.
  inline Quaternion<T> Inverse() const { return Quaternion<T>(s_, -v_); }

  /// @brief Add this Quaternion to another Quaternion.
  ///
  /// @param q Quaternion to add.
  /// @return Quaternion containing the result.
  inline Quaternion<T> operator+(const Quaternion<T>& q) const {
    return Quaternion<T>(s_ + q.s_, v_ + q.v_);
  }

  /// @brief Add another Quaternion to this, and store the result.
  ///
  /// @param q Quaternion to add.
  /// @return Quaternion containing the result.
  inline Quaternion<T>& operator+=(const Quaternion<T>& q) {
    s_ += q.s_;
    v_ += q.v_;
    return *this;
  }

  /// @brief Multiply this Quaternion with another Quaternion.
  ///
  /// @note This is equivalent to
  /// <code>FromMatrix(ToMatrix() * q.ToMatrix()).</code>
  /// @param q Quaternion to multiply with.
  /// @return Quaternion containing the result.
  inline Quaternion<T> operator*(const Quaternion<T>& q) const {
    return Quaternion<T>(
        s_ * q.s_ - Vector<T, 3>::DotProduct(v_, q.v_),
        s_ * q.v_ + q.s_ * v_ + Vector<T, 3>::CrossProduct(v_, q.v_));
  }

  /// @brief Multiply this Quaternion by a scalar.
  ///
  /// This conditions the Quaternion to be a rotation <= 180 degrees, then
  /// multiplies the angle of the rotation by a scalar factor.
  ///
  /// If the scalar factor is < 1, the resulting rotation will be on the shorter
  /// of the two paths to the identity orientation, which is often intuitive but
  /// can trip you up if you really did want to take the longer path.
  ///
  /// If the scalar factor is > 1, the resulting rotation will be on the longer
  /// of the two paths to the identity orientation, which can be unintuitive.
  /// For example, you are not guaranteed that (q * 2) * .5 and q * (2 * .5)
  /// are the same orientation, let alone the same quaternion.
  ///
  /// @param s1 Scalar to multiply with.
  /// @return Quaternion containing the result.
  inline Quaternion<T> operator*(T s1) const {
    T angle;
    Vector<T, 3> axis;
    ToAngleAxis(&angle, &axis);
    angle *= s1;
    // The axis coming from ToAngleAxis() is already normalized, but
    // ToAngleAxis may return slightly non-normal axes in unstable cases.
    // It should arguably handle that internally, allowing us to remove
    // the Normalized() here.
    return Quaternion<T>(cos(0.5f * angle),
                         axis.Normalized() * static_cast<T>(sin(0.5f * angle)));
  }

  /// @brief Multiply a Vector by this Quaternion.
  ///
  /// This will rotate the specified vector by the rotation specified by this
  /// Quaternion.
  ///
  /// @param v1 Vector to multiply by this Quaternion.
  /// @return Rotated Vector.
  inline Vector<T, 3> operator*(const Vector<T, 3>& v1) const {
    T ss = s_ + s_;
    return ss * Vector<T, 3>::CrossProduct(v_, v1) + (ss * s_ - 1) * v1 +
           2 * Vector<T, 3>::DotProduct(v_, v1) * v_;
  }

  /// @brief Normalize this quaternion (in-place).
  ///
  /// @return Length of the quaternion.
  inline T Normalize() {
    T length = sqrt(s_ * s_ + Vector<T, 3>::DotProduct(v_, v_));
    T scale = (1 / length);
    s_ *= scale;
    v_ *= scale;
    return length;
  }

  /// @brief Calculate the normalized version of this quaternion.
  ///
  /// @return The normalized quaternion.
  inline Quaternion<T> Normalized() const {
    Quaternion<T> q(*this);
    q.Normalize();
    return q;
  }

  /// @brief Convert this Quaternion to a shortest-path angle and axis.
  ///
  /// The resulting angle-axis is guaranteed to have have angle <= 180 and
  /// represents the same orientation as *this, but it may not convert back to
  /// *this.
  ///
  /// For example, if *this represents "Rotate 350 degrees left", you will
  /// get the angle-axis "Rotate 10 degrees right".
  ///
  /// @param angle Receives the angle, in the range [0, pi].
  /// @param axis Receives the normalized axis.
  inline void ToAngleAxis(T* out_angle, Vector<T, 3>* out_axis) const {
    const Quaternion<T> q = (s_ > 0) ? *this : Quaternion<T>(-s_, -v_);
    q.ToAngleAxisFull(out_angle, out_axis);
  }

  /// @brief Convert this Quaternion to an angle and axis.
  ///
  /// The resulting angle-axis uses the full range of angles supported by
  /// quaternions, and will convert back to the original Quaternion.
  ///
  /// @param angle Receives the angle, in the range [0, 2pi).
  /// @param axis Receives the normalized axis.
  inline void ToAngleAxisFull(T* out_angle, Vector<T, 3>* out_axis) const {
    Vector<T, 3> axis = v_;
    const T axis_length = axis.Normalize();
    if (axis_length == 0) {
      // Normalize has left NaNs in axis.  This happens at angle = 0 and 360.
      // All axes are correct, so any will do.
      *out_axis = Vector<T, 3>(1, 0, 0);
    } else {
      *out_axis = axis;
    }
    *out_angle = 2 * atan2(axis_length, s_);
  }

  /// @brief Convert this Quaternion to 3 Euler Angles.
  ///
  /// @return 3-dimensional Vector where each element is a angle of rotation
  /// (in radians) around the x, y, and z axes.
  inline Vector<T, 3> ToEulerAngles() const {
    Matrix<T, 3> m(ToMatrix());
    T cos2 = m[0] * m[0] + m[1] * m[1];
    if (cos2 < 1e-6f) {
      return Vector<T, 3>(
          0,
          m[2] < 0 ? static_cast<T>(0.5 * M_PI) : static_cast<T>(-0.5 * M_PI),
          -std::atan2(m[3], m[4]));
    } else {
      return Vector<T, 3>(std::atan2(m[5], m[8]),
                          std::atan2(-m[2], std::sqrt(cos2)),
                          std::atan2(m[1], m[0]));
    }
  }

  /// @brief Convert to a 3x3 Matrix.
  ///
  /// @return 3x3 rotation Matrix.
  inline Matrix<T, 3> ToMatrix() const {
    const T x2 = v_[0] * v_[0], y2 = v_[1] * v_[1], z2 = v_[2] * v_[2];
    const T sx = s_ * v_[0], sy = s_ * v_[1], sz = s_ * v_[2];
    const T xz = v_[0] * v_[2], yz = v_[1] * v_[2], xy = v_[0] * v_[1];
    return Matrix<T, 3>(1 - 2 * (y2 + z2), 2 * (xy + sz), 2 * (xz - sy),
                        2 * (xy - sz), 1 - 2 * (x2 + z2), 2 * (sx + yz),
                        2 * (sy + xz), 2 * (yz - sx), 1 - 2 * (x2 + y2));
  }


template <typename T>
Quaternion<T> Quaternion<T>::identity = Quaternion<T>(1, 0, 0, 0);

/// @addtogroup ZMATH_quaternion
/// @{

/// @brief Multiply a Quaternion by a scalar.
///
/// This multiplies the angle of the rotation of the specified Quaternion
/// by a scalar factor.
/// @param s Scalar to multiply with.
/// @param q Quaternion to scale.
/// @return Quaternion containing the result.
///
/// @related Quaternion
template <class T>
inline Quaternion<T> operator*(T s, const Quaternion<T>& q) {
  return q * s;
}
/// @}

}  // namespace zmath
#endif  // ZMATH_QUATERNION_H_


// this is gonna kill me
