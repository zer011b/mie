#include <cstdio>
#include <complex>
#include <assert.h>

#define FPVALUE double
#define VALUE std::complex<FPVALUE>
#define ASSERT(x) assert(x)
#define SQR(x) (x*x)

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

int main ()
{
  printf ("%f\n", calc_derivative_associated_legendre_func (1, 2, 0.5));

  return 0;
}
