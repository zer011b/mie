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

FPVALUE calc_bessel_j (int l, FPVALUE arg)
{
  return gsl_sf_bessel_jl (l, arg);
}

FPVALUE calc_derivative_bessel_j (int l, FPVALUE arg)
{
  // TODO
}

// return vector in spherical coordinates
Vec3D<FPVALUE> calc_M (int m, int l, FPVALUE r, FPVALUE theta, FPVALUE phi, FPVALUE k, bool is_even)
{
  ASSERT (m == 1);

  FPVALUE e_r = 0;
  FPVALUE e_theta = (is_even ? -1 : 1) * (m / sin (theta)) * calc_bessel_j (l, k * r) * calc_associated_legendre_func (m, l, cos (theta)) * (is_even ? sin (m * phi) : cos (m * phi));
  FPVALUE e_phi = calc_bessel_j (l, k * r) * calc_derivative_associated_legendre_func (m, l, cos (theta)) * (is_even ? cos (m * phi) : sin (m * phi));

  return Vec3D<FPVALUE> (e_r, e_theta, e_phi);
}

Vec3D<FPVALUE> calc_N (int m, int l, FPVALUE r, FPVALUE theta, FPVALUE phi, FPVALUE k, bool is_even)
{
  ASSERT (m == 1);

  FPVALUE derivative = calc_bessel_j (l, k * r) + r * calc_derivative_bessel_j (l, k * r);

  FPVALUE e_r = (l * (l + 1) / (k * r)) * calc_bessel_j (l, k * r) * calc_associated_legendre_func (m, l, cos (theta)) * (is_even ? cos (m * phi) : sin (m * phi));
  FPVALUE e_theta = (1.0 / (k * r)) * (derivative) * calc_associated_legendre_func (m, l, cos (theta)) * (is_even ? cos (m * phi) : sin (m * phi));
  FPVALUE e_phi = (is_even ? -1 : 1) * (m / (k * r * sin (theta))) * (derivative) * calc_associated_legendre_func (m, l, cos (theta)) * (is_even ? sin (m * phi) : cos (m * phi));

  return Vec3D<FPVALUE> (e_r, e_theta, e_phi);
}

Vec3D<VALUE> calc_E_inc (int maxL, FPVALUE r, FPVALUE theta, FPVALUE phi, FPVALUE k)
{
  Vec3D<VALUE> res;

  VALUE imag (0, 1);
  VALUE imag_l = imag;

  for (int l = 1; l <= maxL; ++l)
  {
    FPVALUE multiplier = (2.0 * l + 1) / (l * (l + 1));

    Vec3D<FPVALUE> M = calc_M (1, l, r, theta, phi, k, false);
    Vec3D<FPVALUE> N = calc_N (1, l, r, theta, phi, k, true);

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

int main ()
{
  //printf ("%f\n", calc_derivative_associated_legendre_func (1, 2, 0.5));
  FPVALUE lambda = 1.0;
  FPVALUE k = 2 * M_PI / lambda;

  int maxL = 40;

  Vec3D<VALUE> E_inc_polar = calc_E_inc (maxL, 1.0, M_PI / 2.0, 0.0, k);
  Vec3D<VALUE> E_inc = convert_polar_to_decart (E_inc_polar, M_PI / 2.0, 0.0);

  printf ("( {%f,%f} , {%f,%f} , {%f,%f} )\n",
          E_inc.getX ().real (), E_inc.getX ().imag (),
          E_inc.getY ().real (), E_inc.getY ().imag (),
          E_inc.getZ ().real (), E_inc.getZ ().imag ());

  return 0;
}
