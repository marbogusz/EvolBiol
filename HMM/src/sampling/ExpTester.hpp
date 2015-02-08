/*
 * ExpTester.hpp
 *
 *  Created on: Feb 3, 2015
 *      Author: marcin
 */

#ifndef SAMPLING_EXPTESTER_HPP_
#define SAMPLING_EXPTESTER_HPP_

// SSE2
#include <emmintrin.h>

#define SSE2_DOUBLE_CONST(Key, v) \
    static const  double _D_##Key[2] __attribute__((aligned(16))) = { v, v }
// repeat a const to be in packed form

#define SSE2_INT_CONST(Key, v) \
    static const  int _I_##Key[4] __attribute__((aligned(16))) = {v, v, v, v}

/*  SSE2 defs start */
SSE2_DOUBLE_CONST(1, 1.0);
SSE2_DOUBLE_CONST(0p5, 0.5);

SSE2_DOUBLE_CONST(cephes_LOG2, 1.44269504088896341);
SSE2_DOUBLE_CONST(cephes_C1, 0.693359375);
SSE2_DOUBLE_CONST(cephes_C2, -2.12194440e-4);

SSE2_DOUBLE_CONST(cephes_p0, 1.9875691500E-4);
SSE2_DOUBLE_CONST(cephes_p1, 1.3981999507E-3);
SSE2_DOUBLE_CONST(cephes_p2, 8.3334519073E-3);
SSE2_DOUBLE_CONST(cephes_p3, 4.1665795894E-2);
SSE2_DOUBLE_CONST(cephes_p4, 1.6666665459E-1);
SSE2_DOUBLE_CONST(cephes_p5, 5.0000001201E-1);

SSE2_DOUBLE_CONST(exp_upper,    88.3762626647949);
SSE2_DOUBLE_CONST(exp_lower,    -88.3762626647949);

SSE2_INT_CONST(0x7f, 0x7f);
/*  SSE2 defs end */
#define _PS_CONST(Name, Val)                                            \
  static const float _ps_##Name[4] __attribute__((aligned(16))) = { Val, Val, Val, Val }

#define _PI32_CONST(Name, Val)                                            \
  static const int _pi32_##Name[4] __attribute__((aligned(16))) = { Val, Val, Val, Val }

_PS_CONST(exp_hi,	88.3762626647949f);
_PS_CONST(exp_lo,	-88.3762626647949f);

_PS_CONST(cephes_LOG2EF, 1.44269504088896341);
_PS_CONST(cephes_exp_C1, 0.693359375);
_PS_CONST(cephes_exp_C2, -2.12194440e-4);

_PS_CONST(cephes_exp_p0, 1.9875691500E-4);
_PS_CONST(cephes_exp_p1, 1.3981999507E-3);
_PS_CONST(cephes_exp_p2, 8.3334519073E-3);
_PS_CONST(cephes_exp_p3, 4.1665795894E-2);
_PS_CONST(cephes_exp_p4, 1.6666665459E-1);
_PS_CONST(cephes_exp_p5, 5.0000001201E-1);
_PS_CONST(1  , 1.0f);
_PS_CONST(0p5, 0.5f);

_PI32_CONST(1, 1);
_PI32_CONST(inv1, ~1);
_PI32_CONST(2, 2);
_PI32_CONST(4, 4);
_PI32_CONST(0x7f, 0x7f);

typedef __m128 v4sf;  // vector of 4 float (sse1)
typedef __m128i v4si; // vector of 4 int (sse2)

typedef union {
  float f[4];
  int i[4];
  v4sf  v;
} __attribute__((aligned(16))) V4SF;

namespace EBC
{

class ExpTester {
private:
	double oneln2;
	double expConst;

public:
	ExpTester();
	virtual ~ExpTester();

	void exp_sse2(double in[2], double out[2]);

	v4sf exp_ps(v4sf x);

	double fexp1(double);

	float testExpFloat(unsigned int itr);

	double testExpDouble(unsigned int itr);
};

}
#endif /* SAMPLING_EXPTESTER_HPP_ */
