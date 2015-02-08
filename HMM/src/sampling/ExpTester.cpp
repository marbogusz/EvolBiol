/*
 * ExpTester.cpp
 *
 *  Created on: Feb 3, 2015
 *      Author: marcin
 */

#include <sampling/ExpTester.hpp>
#include <cmath>
#include <iostream>

namespace EBC
{

ExpTester::ExpTester() {
	// TODO Auto-generated constructor stub
	oneln2 = 1.0/log(2.0);
	expConst = exp(1.0);
}

ExpTester::~ExpTester() {
	// TODO Auto-generated destructor stub
}



void ExpTester::exp_sse2(double in[2], double out[2])
{
    __m128d x = _mm_load_pd(in);

    __m128d tmp = _mm_setzero_pd();

    __m128d fx;

    __m128i emm0;

    __m128d oneD = *(__m128d*) _D_1;

    x = _mm_min_pd(x, *(__m128d*)_D_exp_upper);
    x = _mm_max_pd(x, *(__m128d*)_D_exp_lower);

    // exp(x) = exp( z + n log(2) ) = exp(z) * n^2 , where 0<=z < 1

    fx = _mm_mul_pd(x, *(__m128d*)_D_cephes_LOG2);
    fx = _mm_add_pd(fx, *(__m128d*)_D_0p5);

    emm0 = _mm_cvttpd_epi32(fx); //overloaded to handle _m128d
    tmp = _mm_cvtepi32_pd(emm0);

    __m128d flag = _mm_cmpgt_pd(tmp, fx); // ? 0xffffffffffffffff : 0x0
    flag = _mm_and_pd(flag, oneD);
    fx = _mm_sub_pd(tmp, flag);

    tmp = _mm_mul_pd(fx, *(__m128d*)_D_cephes_C1);
    __m128d tmp2 = _mm_mul_pd(fx, *(__m128d*)_D_cephes_C2);
    x = _mm_sub_pd(x, tmp);

    x = _mm_mul_pd(x, tmp2);

    __m128d z = _mm_mul_pd(x, x);

    __m128d y = *(__m128d*)_D_cephes_p0;
    y = _mm_mul_pd(y, x);
    y = _mm_add_pd(y, *(__m128d*)_D_cephes_p1);
    y = _mm_mul_pd(y, x);
    y = _mm_add_pd(y, *(__m128d*)_D_cephes_p2);
    y = _mm_mul_pd(y, x);
    y = _mm_add_pd(y, *(__m128d*)_D_cephes_p3);
    y = _mm_mul_pd(y, x);
    y = _mm_add_pd(y, *(__m128d*)_D_cephes_p4);
    y = _mm_mul_pd(y, x);
    y = _mm_add_pd(y, *(__m128d*)_D_cephes_p5);
    y = _mm_mul_pd(y, z);
    y = _mm_add_pd(y, x);
    y = _mm_add_pd(y, oneD);

    /* build 2^n */

    emm0 = _mm_cvttpd_epi32(fx);
    emm0 = _mm_add_epi32(emm0, *(__m128i*)_I_0x7f);
    emm0 = _mm_slli_epi32(emm0, 23);
    __m128d pow2n = _mm_castsi128_pd(emm0);

  y = _mm_mul_pd(y, pow2n);
    _mm_store_pd(out, y);

}


double ExpTester::fexp1(double x)
{
	double res;
	double pos = x * -1.0;
	res = (1.0 + pos*(1.0 + (pos/2.0)*(1.0+(pos/3.0)*(1.0+(pos/4.0)))));

	return 1.0/res;

	//ef = 1 + f(1 + f/2(1 + f/3(1 + f/4)))

}

v4sf ExpTester::exp_ps(v4sf x) {
  v4sf tmp = _mm_setzero_ps(), fx;

  v4si emm0;

  v4sf one = *(v4sf*)_ps_1;

  x = _mm_min_ps(x, *(v4sf*)_ps_exp_hi);
  x = _mm_max_ps(x, *(v4sf*)_ps_exp_lo);

  /* express exp(x) as exp(g + n*log(2)) */
  fx = _mm_mul_ps(x, *(v4sf*)_ps_cephes_LOG2EF);
  fx = _mm_add_ps(fx, *(v4sf*)_ps_0p5);

  /* how to perform a floorf with SSE: just below */

  emm0 = _mm_cvttps_epi32(fx);
  tmp  = _mm_cvtepi32_ps(emm0);

  /* if greater, substract 1 */
  v4sf mask = _mm_cmpgt_ps(tmp, fx);
  mask = _mm_and_ps(mask, one);
  fx = _mm_sub_ps(tmp, mask);

  tmp = _mm_mul_ps(fx, *(v4sf*)_ps_cephes_exp_C1);
  v4sf z = _mm_mul_ps(fx, *(v4sf*)_ps_cephes_exp_C2);
  x = _mm_sub_ps(x, tmp);
  x = _mm_sub_ps(x, z);

  z = _mm_mul_ps(x,x);

  v4sf y = *(v4sf*)_ps_cephes_exp_p0;
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p1);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p2);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p3);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p4);
  y = _mm_mul_ps(y, x);
  y = _mm_add_ps(y, *(v4sf*)_ps_cephes_exp_p5);
  y = _mm_mul_ps(y, z);
  y = _mm_add_ps(y, x);
  y = _mm_add_ps(y, one);

  /* build 2^n */

  emm0 = _mm_cvttps_epi32(fx);
  emm0 = _mm_add_epi32(emm0, *(v4si*)_pi32_0x7f);
  emm0 = _mm_slli_epi32(emm0, 23);
  v4sf pow2n = _mm_castsi128_ps(emm0);

  y = _mm_mul_ps(y, pow2n);
  return y;
}

float ExpTester::testExpFloat(unsigned int itr) {

	unsigned int ctr = itr;
	float f = -0.001;
	float g = -1.001;
	float h = -10.001;
	float i = -100.001;

	float ex = 0.0;
	while(ctr > 0){
		ex += (exp(f) + exp(g) + exp(h)+ exp(i));
		ctr--;
	}
	return ex;
}

double ExpTester::testExpDouble(unsigned int itr) {
	unsigned int ctr = itr;

	V4SF vx, exp4;
	vx.f[0] = -0.001;
	vx.f[1] = -1.001;
	vx.f[2] = -10.001;
	vx.f[3] = -100.001;
	exp4.v = exp_ps(vx.v);
	std::cerr << exp4.f[0] << "\n";
	std::cerr << exp4.f[1] << "\n";
	std::cerr << exp4.f[2] << "\n";

	float ex = 0.0;
	while(ctr > 0){
		exp4.v = exp_ps(vx.v);
		ex += (exp4.f[0] + exp4.f[1] +exp4.f[2] +exp4.f[3]);
		ctr--;
	}
	return (double)ex;
}

}
