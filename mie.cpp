#include <cstdio>
#include <complex>
#include <assert.h>
#include <cmath>
#include <gsl/gsl_sf_bessel.h>

// be sure to check library functions, in case of change (i.e. bessel functions take double and return double)
#define FPVALUE double
#define VALUE std::complex<FPVALUE>
#define ASSERT(x) assert(x)
#define SQR(x) (x*x)
#define NORM(x) (sqrt(SQR(x.real ()) + SQR(x.imag())))
#define ACCURACY (0.000000001)
#define IS_FP_EXACT(x,y) (((x - y) < 0 ? y - x : x - y) < ACCURACY)
#define IS_COMPLEX_EXACT(x,y) (IS_FP_EXACT(x.real(), y.real()) && IS_FP_EXACT(x.imag(), y.imag()))

template <typename T>
class Vec3D
{
  T x;
  T y;
  T z;

public:

  Vec3D (T xx = T (0), T yy = T (0), T zz = T (0))
  : x (xx), y (yy), z (zz) {}

  T abs () const
  {
    return sqrt (SQR (x) + SQR (y) + SQR (z));
  }

  Vec3D operator + (const Vec3D &op2) const
  {
    return Vec3D (x + op2.x, y + op2.y, z + op2.z);
  }

  Vec3D operator - (const Vec3D &op2) const
  {
    return Vec3D (x - op2.x, y - op2.y, z - op2.z);
  }

  T getX () const
  {
    return x;
  }

  T getY () const
  {
    return y;
  }

  T getZ () const
  {
    return z;
  }

  void setX (T xx)
  {
    x = xx;
  }

  void setY (T yy)
  {
    y = yy;
  }

  void setZ (T zz)
  {
    z = zz;
  }
};

// calculate P_l^m (arg)
FPVALUE calc_associated_legendre_func (int m, int l, FPVALUE arg)
{
  ASSERT (m==1);

  if (l == 0)
  {
    return FPVALUE (0);
  }

  if (l == 1)
  {
    return sqrt (1 - SQR (arg));
  }

  ASSERT (l >= 2);

  FPVALUE prev = calc_associated_legendre_func (m, l - 1, arg);
  FPVALUE prevPrev = calc_associated_legendre_func (m, l - 2, arg);

  return (1.0/(l-1)) * ((2*l-1) * arg * prev - l * prevPrev);
}

// calculate P_l^m (cos(x)) / sin(x)
// also called Pi_n in Bohren/Huffman
FPVALUE calc_associated_legendre_func_special (int m, int l, FPVALUE arg)
{
  ASSERT (m==1);

  if (l == 0)
  {
    return FPVALUE (0);
  }

  if (l == 1)
  {
    return FPVALUE (1);
  }

  ASSERT (l >= 2);

  FPVALUE prev = calc_associated_legendre_func (m, l - 1, arg);
  FPVALUE prevPrev = calc_associated_legendre_func (m, l - 2, arg);

  return (1.0/(l-1)) * ((2*l-1) * arg * prev - l * prevPrev);
}

// calculate d(P_l^m(x))/dx
FPVALUE calc_derivative_associated_legendre_func (int m, int l, FPVALUE arg)
{
  ASSERT (m==1);

  if (l == 0)
  {
    return FPVALUE (0);
  }

  ASSERT (l >= 1);

  FPVALUE prev = calc_associated_legendre_func (m, l - 1, arg);
  FPVALUE cur = calc_associated_legendre_func (m, l, arg);

  return (1.0/(1 - SQR (arg))) * ((l + 1) * prev - l * arg * cur);
}

// calculate d(P_l^m(cos(x)))/d(cos(x)) * (-sin(x))
FPVALUE calc_derivative_associated_legendre_func_special (int m, int l, FPVALUE arg)
{
  ASSERT (m==1);

  if (l == 0)
  {
    return FPVALUE (0);
  }

  ASSERT (l >= 1);

  FPVALUE prev = calc_associated_legendre_func_special (m, l - 1, arg);
  FPVALUE cur = calc_associated_legendre_func_special (m, l, arg);

  return ((l + 1) * prev - l * arg * cur) * (-1);
}

// ================
// Bessel functions
// ================
VALUE calc_bessel_j (int l, FPVALUE arg)
{
  return VALUE (gsl_sf_bessel_jl (l, arg), 0);
}
VALUE calc_bessel_y (int l, FPVALUE arg)
{
  return VALUE (gsl_sf_bessel_yl (l, arg), 0);
}
VALUE calc_hankel_1 (int l, FPVALUE arg)
{
  return calc_bessel_j (l, arg) + VALUE (0, 1) * calc_bessel_y (l, arg);
}

// ===========================
// Bessel function derivatives
// ===========================
VALUE calc_derivative_bessel (int l, FPVALUE arg, VALUE prev, VALUE next)
{
  return 1.0/(2*l+1) * (FPVALUE(l) * prev - FPVALUE(l+1) * next);
}
VALUE calc_derivative_bessel_j (int l, FPVALUE arg)
{
  //return VALUE (0.5 * (gsl_sf_bessel_jl (l - 1, arg) - (1.0 / arg) * (gsl_sf_bessel_jl (l, arg) + arg * gsl_sf_bessel_jl (l + 1, arg))), 0);
  return calc_derivative_bessel (l, arg, calc_bessel_j (l - 1, arg), calc_bessel_j (l + 1, arg));
}
VALUE calc_derivative_bessel_y (int l, FPVALUE arg)
{
  return calc_derivative_bessel (l, arg, calc_bessel_y (l - 1, arg), calc_bessel_y (l + 1, arg));
}
VALUE calc_derivative_hankel_1 (int l, FPVALUE arg)
{
  return calc_derivative_bessel (l, arg, calc_hankel_1 (l - 1, arg), calc_hankel_1 (l + 1, arg));
}

// return vector in spherical coordinates
Vec3D<VALUE> calc_M (int m, int l, FPVALUE r, FPVALUE theta, FPVALUE phi, FPVALUE k, bool is_even, bool use_hankel)
{
  ASSERT (m == 1);

  VALUE z = VALUE (0, 0);
  if (use_hankel)
  {
    z = calc_hankel_1 (l, k * r);
  }
  else
  {
    z = calc_bessel_j (l, k * r);
  }

  VALUE e_r = VALUE (0, 0);
  VALUE e_theta = (is_even ? -1 : 1) * FPVALUE (m) * z * calc_associated_legendre_func_special (m, l, cos (theta)) * (is_even ? sin (m * phi) : cos (m * phi));
  VALUE e_phi = - z * calc_derivative_associated_legendre_func_special (m, l, cos (theta)) * (is_even ? cos (m * phi) : sin (m * phi));

  return Vec3D<VALUE> (e_r, e_theta, e_phi);
}

Vec3D<VALUE> calc_N (int m, int l, FPVALUE r, FPVALUE theta, FPVALUE phi, FPVALUE k, bool is_even, bool use_hankel)
{
  ASSERT (m == 1);

  VALUE z = VALUE (0, 0);
  VALUE z1 = VALUE (0, 0);
  if (use_hankel)
  {
    z = calc_hankel_1 (l, k * r);
    z1 = calc_derivative_hankel_1 (l, k * r);
  }
  else
  {
    z = calc_bessel_j (l, k * r);
    z1 = calc_derivative_bessel_j (l, k * r);
  }

  VALUE derivative = z + k * r * z1;

  VALUE e_r = (l * (l + 1) / (k * r)) * z * calc_associated_legendre_func (m, l, cos (theta)) * (is_even ? cos (m * phi) : sin (m * phi));
  VALUE e_theta = (1.0 / (k * r)) * (derivative) * calc_derivative_associated_legendre_func_special (m, l, cos (theta)) * (is_even ? cos (m * phi) : sin (m * phi));
  VALUE e_phi = (is_even ? -1 : 1) * (m / (k * r)) * (derivative) * calc_associated_legendre_func_special (m, l, cos (theta)) * (is_even ? sin (m * phi) : cos (m * phi));

  return Vec3D<VALUE> (e_r, e_theta, e_phi);
}

Vec3D<VALUE> calc_E_inc (int maxL, FPVALUE r, FPVALUE theta, FPVALUE phi, FPVALUE k)
{
  Vec3D<VALUE> res;

  VALUE imag (0, 1);
  VALUE imag_l = imag;

  for (int l = 1; l <= maxL; ++l)
  {
    FPVALUE multiplier = (2.0 * l + 1) / (l * (l + 1));

    Vec3D<VALUE> M = calc_M (1, l, r, theta, phi, k, false, false);
    Vec3D<VALUE> N = calc_N (1, l, r, theta, phi, k, true, false);

    VALUE e_r = imag_l * multiplier * (M.getX () - imag * N.getX ());
    VALUE e_theta = imag_l * multiplier * (M.getY () - imag * N.getY ());
    VALUE e_phi = imag_l * multiplier * (M.getZ () - imag * N.getZ ());

    res.setX (res.getX () + e_r);
    res.setY (res.getY () + e_theta);
    res.setZ (res.getZ () + e_phi);

    imag_l *= imag;
  }

  return res;
}

// pass angles of the polar coordinate system
Vec3D<VALUE> convert_polar_to_decart (const Vec3D<VALUE> &polar, FPVALUE theta, FPVALUE phi)
{
  VALUE x = polar.getX () * sin (theta) * cos (phi) + polar.getY () * cos (theta) * cos (phi) - polar.getZ () * sin (phi);
  VALUE y = polar.getX () * sin (theta) * sin (phi) + polar.getY () * cos (theta) * sin (phi) + polar.getZ () * cos (phi);
  VALUE z = polar.getX () * cos (theta) - polar.getY () * sin (theta);

  Vec3D<VALUE> res (x, y, z);

  return res;
}

void simple_test (FPVALUE theta, FPVALUE phi, Vec3D<VALUE> &correct_inc)
{
  FPVALUE lambda = 1.0;
  FPVALUE k = 2 * M_PI / lambda;

  int maxL = 30;

  Vec3D<VALUE> E_inc_polar = calc_E_inc (maxL, 1.0, theta, phi, k);
  Vec3D<VALUE> E_inc = convert_polar_to_decart (E_inc_polar, theta, phi);

  ASSERT (IS_COMPLEX_EXACT (E_inc.getX (), correct_inc.getX ()));
  ASSERT (IS_COMPLEX_EXACT (E_inc.getY (), correct_inc.getY ()));
  ASSERT (IS_COMPLEX_EXACT (E_inc.getZ (), correct_inc.getZ ()));
}

void test ()
{
  Vec3D<VALUE> test_zero (VALUE(1.0, 0.0), VALUE(0.0, 0.0), VALUE(0.0, 0.0));
  simple_test (M_PI / 2.0, 0.0, test_zero);
  simple_test (M_PI / 2.0, M_PI / 4.0, test_zero);
  simple_test (M_PI / 2.0, M_PI / 2.0, test_zero);

  Vec3D<VALUE> test_zero2 (VALUE(0.03799544386587661, 0.098683316029646362), VALUE(0.0, 0.0), VALUE(0.0, 0.0));
  simple_test (0, 0, test_zero2);
  // simple_test (0, M_PI / 4.0, test_zero2);
  // simple_test (0, M_PI / 2.0, test_zero2);

  // Vec3D<VALUE> test_zero3 (VALUE(0.03799544386587661, 0.098683316029646362), VALUE(0.0, 0.0), VALUE(0.0, 0.0));
  // simple_test (M_PI / 4.0, 0, test_zero3);
  // simple_test (M_PI / 4.0, M_PI / 2.0, test_zero3);
}

int main ()
{
  test ();

  // Vec3D<VALUE> E_inc_polar = calc_E_inc (maxL, 1.0, 0, 0, k);
  // Vec3D<VALUE> E_inc = convert_polar_to_decart (E_inc_polar, 0, 0);
  //
  // printf ("( {%f,%f}=|%f| , {%f,%f}=|%f| , {%f,%f}=|%f| )\n",
  //         E_inc_polar.getX ().real (), E_inc_polar.getX ().imag (), NORM (E_inc_polar.getX ()),
  //         E_inc_polar.getY ().real (), E_inc_polar.getY ().imag (), NORM (E_inc_polar.getY ()),
  //         E_inc_polar.getZ ().real (), E_inc_polar.getZ ().imag (), NORM (E_inc_polar.getZ ()));
  //
  // printf ("( {%f,%f}=|%f| , {%f,%f}=|%f| , {%f,%f}=|%f| )\n",
  //         E_inc.getX ().real (), E_inc.getX ().imag (), NORM (E_inc.getX ()),
  //         E_inc.getY ().real (), E_inc.getY ().imag (), NORM (E_inc.getY ()),
  //         E_inc.getZ ().real (), E_inc.getZ ().imag (), NORM (E_inc.getZ ()));

  return 0;
}
