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

namespace zmath {

template <class T>
class Quaternion {
 public:
  inline Quaternion(const Quaternion<T>& q) {
    s_ = q.s_;
    v_ = q.v_;
  }

  inline Quaternion(T s1, T qs1, T qs2, T qs3) {
    s_ = s1;
    v_ = Vector<T, 3>(qs1, qs2, qs3);
  }
  inline Quaternion(T s1, const Vector<T, 3>& v1) {
    s_ = s1;
    v_ = v1;
  }

  inline const T& scalar() const { return s_; }

  inline void set_scalar(T s) { s_ = s; }

  inline const Vector<T, 3>& vector() const { return v_; }

  inline void set_vector(const Vector<T, 3>& v) { v_ = v; }

  inline Quaternion<T> Inverse() const { return Quaternion<T>(s_, -v_); }

  inline Quaternion<T> operator+(const Quaternion<T>& q) const {
    return Quaternion<T>(s_ + q.s_, v_ + q.v_);
  }

  inline Quaternion<T>& operator+=(const Quaternion<T>& q) {
    s_ += q.s_;
    v_ += q.v_;
    return *this;
  }

  inline Quaternion<T> operator*(const Quaternion<T>& q) const {
    return Quaternion<T>(
        s_ * q.s_ - Vector<T, 3>::DotProduct(v_, q.v_),
        s_ * q.v_ + q.s_ * v_ + Vector<T, 3>::CrossProduct(v_, q.v_));
  }

  /// @brief Multiply this Quaternion by a scalar.
  ///
  /// This conditions the Quaternion to be a rotation <= 180 degrees, then
  /// multiplies the angle of the rotation by a scalar factor.

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

  inline Vector<T, 3> operator*(const Vector<T, 3>& v1) const {
    T ss = s_ + s_;
    return ss * Vector<T, 3>::CrossProduct(v_, v1) + (ss * s_ - 1) * v1 +
           2 * Vector<T, 3>::DotProduct(v_, v1) * v_;
  }

  inline T Normalize() {
    T length = sqrt(s_ * s_ + Vector<T, 3>::DotProduct(v_, v_));
    T scale = (1 / length);
    s_ *= scale;
    v_ *= scale;
    return length;
  }

  inline Quaternion<T> Normalized() const {
    Quaternion<T> q(*this);
    q.Normalize();
    return q;
  }

  inline void ToAngleAxis(T* out_angle, Vector<T, 3>* out_axis) const {
    const Quaternion<T> q = (s_ > 0) ? *this : Quaternion<T>(-s_, -v_);
    q.ToAngleAxisFull(out_angle, out_axis);
  }

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

template <class T>
inline Quaternion<T> operator*(T s, const Quaternion<T>& q) {
  return q * s;
}
/// @}

}  // namespace zmath
#endif  // ZMATH_QUATERNION_H_
