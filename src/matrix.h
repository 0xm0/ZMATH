#ifndef ZMATH_MATRIX_H_
#define ZMATH_MATRIX_H_
#include "zmath/utilities.h"
#include "zmath/vector.h"
#include <cmath>
#include <assert.h>

#ifdef _MSC_VER
#pragma warning(push)
// The following disables warnings for ZMATH_MAT_OPERATION.
// The buffer overrun warning must be disabled as MSVC doesn't treat
// "Cols" as constant and therefore assumes that it's possible
// to overrun arrays indexed by "i".
// The conditional expression is constant warning is disabled since
// MSVC decides that "Cols" *is* constant when unrolling the operation
// loop.
#pragma warning(disable : 4127)  // conditional expression is constant
#pragma warning(disable : 4789)  // buffer overrun
#if _MSC_VER >= 1900             // MSVC 2015
#pragma warning(disable : 4456)  // allow shadowing in unrolled loops
#pragma warning(disable : 4723)  // suppress "potential divide by 0" warning
#endif                           // _MSC_VER >= 1900
#endif                           // _MSC_VER

#define ZMATH_VECTOR_STRIDE_FLOATS(vector) (sizeof(vector) / sizeof(float))

#define ZMATH_MAT_OPERATION(OP) ZMATH_UNROLLED_LOOP(i, Cols, OP)

#define ZMATH_MAT_OPERATOR(OP)                   \
  {                                               \
    Matrix<T, Rows, Cols> result;              \
    ZMATH_MAT_OPERATION(result.data_[i] = (OP)); \
    return result;                                \
  }

#define ZMATH_MAT_SELF_OPERATOR(OP) \
  {                                  \
    ZMATH_MAT_OPERATION(OP);        \
    return *this;                    \
  }

#define ZMATH_MATRIX_4X4_DOT(data1, data2, r)               \
  ((data1)[r] * (data2)[0] + (data1)[(r) + 4] * (data2)[1] + \
   (data1)[(r) + 8] * (data2)[2] + (data1)[(r) + 12] * (data2)[3])



#define ZMATH_MATRIX_3X3_DOT(data1, data2, r, size)              \
  ((data1)[r] * (data2)[0] + (data1)[(r) + (size)] * (data2)[1] + \
   (data1)[(r) + 2 * (size)] * (data2)[2])


namespace zmath {

template <class T, int Rows, int Cols = Rows>
class Matrix;
template <class T, int Rows, int Cols>
inline Matrix<T, Rows, Cols> IdentityHelper();
template <bool check_invertible, class T, int Rows, int Cols>
inline bool InverseHelper(
    const Matrix<T, Rows, Cols>& m, Matrix<T, Rows, Cols>* const inverse,
    T det_thresh);
template <class T, int size1, int size2, int size3>
inline void TimesHelper(const Matrix<T, size1, size2>& m1,
                        const Matrix<T, size2, size3>& m2,
                        Matrix<T, size1, size3>* out_m);
template <class T, int Rows, int Cols>
static inline Matrix<T, Rows, Cols> OuterProductHelper(
    const Vector<T, Rows>& v1, const Vector<T, Cols>& v2);
template <class T>
inline Matrix<T, 4, 4> PerspectiveHelper(T fovy, T aspect, T znear, T zfar,
                                         T handedness);
template <class T>
static inline Matrix<T, 4, 4> OrthoHelper(T left, T right, T bottom, T top,
                                          T znear, T zfar, T handedness);
template <class T>
static inline Matrix<T, 4, 4> LookAtHelper(const Vector<T, 3>& at,
                                           const Vector<T, 3>& eye,
                                           const Vector<T, 3>& up,
                                           T handedness);
template <class T>
static inline bool UnProjectHelper(const Vector<T, 3>& window_coord,
                                   const Matrix<T, 4, 4>& model_view,
                                   const Matrix<T, 4, 4>& projection,
                                   const float window_width,
                                   const float window_height,
                                   Vector<T, 3>& result);

template <typename T, int Rows, int Cols, typename CompatibleT>
static inline Matrix<T, Rows, Cols> FromTypeHelper(const CompatibleT& compatible);

template <typename T, int Rows, int Cols, typename CompatibleT>
static inline CompatibleT ToTypeHelper(const Matrix<T, Rows, Cols>& m);


template <class T>
class Constants {
 public:
  
  static T GetDeterminantThreshold() {
    assert(false);
    return 0;
  }
};

template <>
class Constants<float> {
 public:
  /// Effective uniform scale limit: ~(1/215)^3
  static float GetDeterminantThreshold() { return 1e-7f; }
};

/// Functions that return constants for <code>double</code> values.
template <>
class Constants<double> {
 public:
  /// Effective uniform scale limit: ~(1/100000)^3
  static double GetDeterminantThreshold() { return 1e-15; }
};

  /// Create a Matrix from the first row * column elements of an array.

  explicit inline Matrix(const T* const a) {
    ZMATH_MAT_OPERATION((data_[i] = Vector<T, Rows>(&a[i * Cols])));
  }

  /// Create a Matrix from an array of "Cols", "Rows" element packed vectors.

  explicit inline Matrix(const VectorPacked<T, Rows>* const vectors) {
    ZMATH_MAT_OPERATION((data_[i] = Vector<T, Rows>(vectors[i])));
  }

  inline const T& operator()(const int row, const int column) const {
    return data_[column][row];
  }

  inline T& operator()(const int row, const int column) {
    return data_[column][row];
  }

  inline const T& operator()(const int i) const { return operator[](i); }

  inline T& operator()(const int i) { return operator[](i); }

  inline const T& operator[](const int i) const {
    return const_cast<Matrix<T, Rows, Cols>*>(this)->operator[](i);
  }

  inline T& operator[](const int i) {
#if defined(ZMATH_COMPILE_WITH_PADDING)
    // In this case Vector<T, 3> is padded, so the element offset must be
    // accessed using the array operator.
    if (Rows == 3) {
      const int row = i % Rows;
      const int col = i / Rows;
      return data_[col][row];
    } else {
      return reinterpret_cast<T*>(data_)[i];
    }
#else
    return reinterpret_cast<T*>(data_)[i];
#endif  // defined(ZMATH_COMPILE_WITH_PADDING)
  }

 
  inline void Pack(VectorPacked<T, Rows>* const vector) const {
    ZMATH_MAT_OPERATION(GetColumn(i).Pack(&vector[i]));
  }

  inline Vector<T, Rows>& GetColumn(const int i) { return data_[i]; }

  inline const Vector<T, Rows>& GetColumn(const int i) const {
    return data_[i];
  }

  inline Matrix<T, Rows, Cols> operator-() const {
    ZMATH_MAT_OPERATOR(-data_[i]);
  }

  inline Matrix<T, Rows, Cols> operator+(
      const Matrix<T, Rows, Cols>& m) const {
    ZMATH_MAT_OPERATOR(data_[i] + m.data_[i]);
  }

  inline Matrix<T, Rows, Cols> operator-(
      const Matrix<T, Rows, Cols>& m) const {
    ZMATH_MAT_OPERATOR(data_[i] - m.data_[i]);
  }

  inline Matrix<T, Rows, Cols> operator+(T s) const {
    ZMATH_MAT_OPERATOR(data_[i] + s);
  }

  inline Matrix<T, Rows, Cols> operator-(T s) const {
    ZMATH_MAT_OPERATOR(data_[i] - s);
  }

  inline Matrix<T, Rows, Cols> operator*(T s) const {
    ZMATH_MAT_OPERATOR(data_[i] * s);
  }

  inline Matrix<T, Rows, Cols> operator/(T s) const {
    return (*this) * (1 / s);
  }

  inline Matrix<T, Rows, Cols> operator*(
      const Matrix<T, Rows, Cols>& m) const {
    Matrix<T, Rows, Cols> result;
    TimesHelper(*this, m, &result);
    return result;
  }

  inline Matrix<T, Rows, Cols>& operator+=(
      const Matrix<T, Rows, Cols>& m) {
    ZMATH_MAT_SELF_OPERATOR(data_[i] += m.data_[i]);
  }

  inline Matrix<T, Rows, Cols>& operator-=(
      const Matrix<T, Rows, Cols>& m) {
    ZMATH_MAT_SELF_OPERATOR(data_[i] -= m.data_[i]);
  }

 
  inline Matrix<T, Rows, Cols>& operator+=(T s) {
    ZMATH_MAT_SELF_OPERATOR(data_[i] += s);
  }

  inline Matrix<T, Rows, Cols>& operator-=(T s) {
    ZMATH_MAT_SELF_OPERATOR(data_[i] -= s);
  }

  inline Matrix<T, Rows, Cols>& operator*=(T s) {
    ZMATH_MAT_SELF_OPERATOR(data_[i] *= s);
  }

  inline Matrix<T, Rows, Cols>& operator/=(T s) {
    return (*this) *= (1 / s);
  }

  inline Matrix<T, Rows, Cols>& operator*=(
      const Matrix<T, Rows, Cols>& m) {
    const Matrix<T, Rows, Cols> copy_of_this(*this);
    TimesHelper(copy_of_this, m, this);
    return *this;
  }

  inline Matrix<T, Rows, Cols> Inverse() const {
    Matrix<T, Rows, Cols> inverse;
    InverseHelper<false>(*this, &inverse, static_cast<T>(0));
    return inverse;
  }

  inline bool InverseWithDeterminantCheck(
      Matrix<T, Rows, Cols>* const inverse,
      T det_thresh = Constants<T>::GetDeterminantThreshold()) const {
    return InverseHelper<true>(*this, inverse, det_thresh);
  }

  inline Matrix<T, Cols, Rows> Transpose() const {
    Matrix<T, Cols, Rows> transpose;
    ZMATH_UNROLLED_LOOP(
        i, Cols, ZMATH_UNROLLED_LOOP(
                        j, Rows, transpose.GetColumn(j)[i] = GetColumn(i)[j]))
    return transpose;
  }

  inline Vector<T, 2> TranslationVector2D() const {
    ZMATH_STATIC_ASSERT(Rows == 3 && Cols == 3);
    return Vector<T, 2>(data_[2][0], data_[2][1]);
  }

  inline Vector<T, 3> TranslationVector3D() const {
    ZMATH_STATIC_ASSERT(Rows == 4 && Cols == 4);
    return Vector<T, 3>(data_[3][0], data_[3][1], data_[3][2]);
  }

  inline Vector<T, 3> ScaleVector3D() const {
    ZMATH_STATIC_ASSERT(Rows >= 3 && Cols >= 3);
    return Vector<T, 3>(data_[0].xyz().Length(),
                        data_[1].xyz().Length(),
                        data_[2].xyz().Length());
  }

  template <typename CompatibleT>
  static inline Matrix<T, Rows, Cols> FromType(const CompatibleT& compatible) {
    return FromTypeHelper<T, Rows, Cols, CompatibleT>(compatible);
  }

  template <typename CompatibleT>
  static inline CompatibleT ToType(const Matrix<T, Rows, Cols>& m) {
    return ToTypeHelper<T, Rows, Cols, CompatibleT>(m);
  }

  static inline Matrix<T, Rows, Cols> OuterProduct(
      const Vector<T, Rows>& v1, const Vector<T, Cols>& v2) {
    return OuterProductHelper(v1, v2);
  }

  static inline Matrix<T, Rows, Cols> HadamardProduct(
      const Matrix<T, Rows, Cols>& m1, const Matrix<T, Rows, Cols>& m2) {
    ZMATH_MAT_OPERATOR(m1[i] * m2[i]);
  }

  static inline Matrix<T, Rows, Cols> Identity() {
    return IdentityHelper<T, Rows, Cols>();
  }

  static inline Matrix<T, 3> FromTranslationVector(const Vector<T, 2>& v) {
    return Matrix<T, 3>(1, 0, 0, 0, 1, 0, v[0], v[1], 1);
  }
  
  static inline Matrix<T, 4> FromTranslationVector(const Vector<T, 3>& v) {
    return Matrix<T, 4>(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, v[0], v[1], v[2],
                        1);
  }

  static inline Matrix<T, Rows> FromScaleVector(const Vector<T, Rows - 1>& v) {
 
    Matrix<T, Rows> return_matrix(Identity());
    for (int i = 0; i < Rows - 1; ++i) return_matrix(i, i) = v[i];
    return return_matrix;
  }

  static inline Matrix<T, 4> FromRotationMatrix(const Matrix<T, 3>& m) {
    return Matrix<T, 4>(m[0], m[1], m[2], 0, m[3], m[4], m[5], 0, m[6], m[7],
                        m[8], 0, 0, 0, 0, 1);
  }

  static inline Matrix<T, 3> ToRotationMatrix(const Matrix<T, 4>& m) {
    return Matrix<T, 3>(m[0], m[1], m[2], m[4], m[5], m[6], m[8], m[9],
                        m[10]);
  }

  static inline Matrix<T, 4> FromAffineTransform(
      const Matrix<T, 4, 3>& affine) {
    return Matrix<T, 4>(affine[0], affine[4], affine[8], static_cast<T>(0),
                        affine[1], affine[5], affine[9], static_cast<T>(0),
                        affine[2], affine[6], affine[10], static_cast<T>(0),
                        affine[3], affine[7], affine[11], static_cast<T>(1));
  }

  static inline Matrix<T, 4, 3> ToAffineTransform(const Matrix<T, 4>& m) {
    return Matrix<T, 4, 3>(m[0], m[4], m[8], m[12], m[1], m[5], m[9], m[13],
                           m[2], m[6], m[10], m[14]);
  }
  
  static inline Matrix<T, 3> RotationX(const Vector<T, 2>& v) {
    return Matrix<T, 3>(1, 0, 0, 0, v.x, v.y, 0, -v.y, v.x);
  }

  static inline Matrix<T, 3> RotationY(const Vector<T, 2>& v) {
    return Matrix<T, 3>(v.x, 0, -v.y, 0, 1, 0, v.y, 0, v.x);
  }

  static inline Matrix<T, 3> RotationZ(const Vector<T, 2>& v) {
    return Matrix<T, 3>(v.x, v.y, 0, -v.y, v.x, 0, 0, 0, 1);
  }
  
  static inline Matrix<T, 3> RotationX(T angle) {
    return RotationX(Vector<T, 2>(cosf(angle), sinf(angle)));
  }

  static inline Matrix<T, 3> RotationY(T angle) {
    return RotationY(Vector<T, 2>(cosf(angle), sinf(angle)));
  }

  static inline Matrix<T, 3> RotationZ(T angle) {
    return RotationZ(Vector<T, 2>(cosf(angle), sinf(angle)));
  }

  static inline Matrix<T, 4, 4> Perspective(T fovy, T aspect, T znear, T zfar,
                                            T handedness = 1) {
    return PerspectiveHelper(fovy, aspect, znear, zfar, handedness);
  }

  static inline Matrix<T, 4, 4> Ortho(T left, T right, T bottom, T top, T znear,
                                      T zfar, T handedness = 1) {
    return OrthoHelper(left, right, bottom, top, znear, zfar, handedness);
  }

  static inline Matrix<T, 4, 4> LookAt(const Vector<T, 3>& at,
                                       const Vector<T, 3>& eye,
                                       const Vector<T, 3>& up,
                                       T handedness = -1) {
    return LookAtHelper(at, eye, up, handedness);
  }

  static inline Matrix<T, 4, 4> Transform(const Vector<T, 3>& position,
                                          const Matrix<T, 3, 3>& rotation,
                                          const Vector<T, 3>& scale) {
    return TransformHelper(position, rotation, scale);
  }

  static inline Vector<T, 3> UnProject(const Vector<T, 3>& window_coord,
                                       const Matrix<T, 4, 4>& model_view,
                                       const Matrix<T, 4, 4>& projection,
                                       const float window_width,
                                       const float window_height) {
    Vector<T, 3> result;
    UnProjectHelper(window_coord, model_view, projection, window_width,
                    window_height, result);
    return result;
  }

  friend inline Vector<T, Cols> operator*(
      const Vector<T, Rows>& v, const Matrix<T, Rows, Cols>& m) {
    const int Dims = Cols;
    ZMATH_VECTOR_OPERATOR((Vector<T, Rows>::DotProduct(m.data_[i], v)));
  }

  static const int kRows = Rows;
  
  static const int kColumns = Cols;
  
  static const int kElements = Rows * Cols;

  ZMATH_DEFINE_CLASS_SIMD_AWARE_NEW_DELETE

  Vector<T, Rows> data_[Cols];
};

template <int index> struct zmathMatrixUnroller {
  template <class T, int Rows, int Cols>
  inline static bool NotEqual(const Matrix<T, Rows, Cols>& lhs,
                       const Matrix<T, Rows, Cols>& rhs) {
    return (lhs[index] != rhs[index]) |
           zmathMatrixUnroller<index-1>::NotEqual(lhs, rhs);
  }
};
template <> struct zmathMatrixUnroller<0> {
  template <class T, int Rows, int Cols>
  inline static bool NotEqual(const Matrix<T, Rows, Cols>& lhs,
                       const Matrix<T, Rows, Cols>& rhs) {
    return lhs[0] != rhs[0];
  }
};

template <class T, int Rows, int Cols>
inline bool operator!=(const Matrix<T, Rows, Cols>& lhs,
                       const Matrix<T, Rows, Cols>& rhs) {
  return zmathMatrixUnroller<Rows * Cols - 1>::NotEqual(lhs, rhs);
}

template <class T, int Rows, int Cols>
inline bool operator==(const Matrix<T, Rows, Cols>& lhs,
                       const Matrix<T, Rows, Cols>& rhs) {
  return !(lhs != rhs);
}

template <class T, int Rows, int Cols>
inline Matrix<T, Rows, Cols> operator*(T s, const Matrix<T, Cols, Rows>& m) {
  return m * s;
}

template <class T, int Rows, int Cols>
inline Vector<T, Rows> operator*(const Matrix<T, Rows, Cols>& m,
                                 const Vector<T, Cols>& v) {
  Vector<T, Rows> result(static_cast<T>(0));
  int offset = 0;
  for (int column = 0; column < Cols; column++) {
    for (int row = 0; row < Rows; row++) {
      result[row] += m[offset + row] * v[column];
    }
    offset += Rows;
  }
  return result;
}

template <class T>
inline Vector<T, 2> operator*(const Matrix<T, 2, 2>& m, const Vector<T, 2>& v) {
  return Vector<T, 2>(m[0] * v[0] + m[2] * v[1], m[1] * v[0] + m[3] * v[1]);
}

template <class T>
inline Vector<T, 3> operator*(const Matrix<T, 3, 3>& m, const Vector<T, 3>& v) {
  return Vector<T, 3>(ZMATH_MATRIX_3X3_DOT(&m.data_[0].data_[0], v, 0, 3),
                      ZMATH_MATRIX_3X3_DOT(&m.data_[0].data_[0], v, 1, 3),
                      ZMATH_MATRIX_3X3_DOT(&m.data_[0].data_[0], v, 2, 3));
}

template <>
inline Vector<float, 3> operator*(const Matrix<float, 3, 3>& m,
                                  const Vector<float, 3>& v) {
  return Vector<float, 3>(
      ZMATH_MATRIX_3X3_DOT(&m.data_[0].data_[0], v, 0,
                            ZMATH_VECTOR_STRIDE_FLOATS(v)),
      ZMATH_MATRIX_3X3_DOT(&m.data_[0].data_[0], v, 1,
                            ZMATH_VECTOR_STRIDE_FLOATS(v)),
      ZMATH_MATRIX_3X3_DOT(&m.data_[0].data_[0], v, 2,
                            ZMATH_VECTOR_STRIDE_FLOATS(v)));
}

template <class T>
inline Vector<T, 4> operator*(const Matrix<T, 4, 4>& m, const Vector<T, 4>& v) {
  return Vector<T, 4>(ZMATH_MATRIX_4X4_DOT(&m.data_[0].data_[0], v, 0),
                      ZMATH_MATRIX_4X4_DOT(&m.data_[0].data_[0], v, 1),
                      ZMATH_MATRIX_4X4_DOT(&m.data_[0].data_[0], v, 2),
                      ZMATH_MATRIX_4X4_DOT(&m.data_[0].data_[0], v, 3));
}

template <class T>
inline Vector<T, 3> operator*(const Matrix<T, 4, 4>& m, const Vector<T, 3>& v) {
  Vector<T, 4> v4(v[0], v[1], v[2], 1);
  v4 = m * v4;
  return Vector<T, 3>(v4[0] / v4[3], v4[1] / v4[3], v4[2] / v4[3]);
}

template <class T, int size1, int size2, int size3>


template <typename T, int Rows, int Cols, typename CompatibleT>
static inline Matrix<T, Rows, Cols> FromTypeHelper(const CompatibleT& compatible) {

#if __cplusplus >= 201103L
  
  union ConversionUnion {
    ConversionUnion() {}  
    CompatibleT compatible;
    VectorPacked<T, Rows> packed[Cols];
  } u;
  static_assert(sizeof(u.compatible) == sizeof(u.packed), "Conversion size mismatch.");
  u.compatible = compatible;
  return Matrix<T, Rows, Cols>(u.packed);
#else
  Matrix<T, Rows, Cols> m;
  assert(sizeof(m) == sizeof(compatible));
  memcpy(&m, &compatible, sizeof(m));
  return m;
#endif  
}
template <typename T, int Rows, int Cols, typename CompatibleT>
static inline CompatibleT ToTypeHelper(const Matrix<T, Rows, Cols>& m) {

#if __cplusplus >= 201103L
  union ConversionUnion {
    ConversionUnion() {}
    CompatibleT compatible;
    VectorPacked<T, Rows> packed[Cols];
  } u;
  static_assert(sizeof(u.compatible) == sizeof(u.packed), "Conversion size mismatch.");
  m.Pack(u.packed);
  return u.compatible;
#else
  CompatibleT compatible;
  assert(sizeof(m) == sizeof(compatible));
  memcpy(&compatible, &m, sizeof(compatible));
  return compatible;
#endif  
}
typedef Matrix<float, 4, 3> AffineTransform;
}  

#ifdef _MSC_VER
#pragma warning(pop)
#endif
#endif  