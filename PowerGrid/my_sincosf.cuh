    /*
      Copyright (c) 2016, Norbert Juffa
      All rights reserved.

      Redistribution and use in source and binary forms, with or without 
      modification, are permitted provided that the following conditions
      are met:

      1. Redistributions of source code must retain the above copyright 
         notice, this list of conditions and the following disclaimer.

      2. Redistributions in binary form must reproduce the above copyright
         notice, this list of conditions and the following disclaimer in the
         documentation and/or other materials provided with the distribution.

      THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 
      "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT 
      LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
      A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
      HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
      SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT 
      LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
      DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
      THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT 
      (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
      OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
    */
    /* 190 bits of 2/PI for Payne-Hanek style argument reduction. */

    #ifndef POWERGRID_my_sincosf_cuh
    #define POWERGRID_my_sincosf_cuh
   
    #include <math.h>
    #include "cudadevice.h"
    #define SANE_COMPILER 1
    extern "C" {
    static __constant__ unsigned int two_over_pi_f [] = 
    {
        0x28be60db,
        0x9391054a,
        0x7f09d5f4,
        0x7d4d3770,
        0x36d8a566,
        0x4f10e410
    };

    __forceinline__ __device__ float trig_red_slowpath_f (float a, int *quadrant)
    {
        unsigned long long int p;
        unsigned int ia, hi, mid, lo, tmp, i;
        int e, q;
        float r;

    #if PORTABLE
        ia = (unsigned int)(fabsf (frexpf (a, &e)) * 4.29496730e+9f); // 0x1.0p32f
    #else // !PORTABLE
        ia = __int_as_float ((__float_as_int (a) & 0x007fffff) + 0x4f000000);
        e = (((unsigned int)__float_as_int (a) >> 23) & 0xff) - 126;
    #endif // !PORTABLE
        
        /* compute product x * 2/pi in 2.62 fixed-point format */
        i = (unsigned int)e >> 5;
        e = (unsigned int)e & 31;

        hi  = i ? two_over_pi_f [i-1] : 0;
        mid = two_over_pi_f [i+0];
        lo  = two_over_pi_f [i+1];
        tmp = two_over_pi_f [i+2];
     
        if (e) {
            hi  = (hi  << e) | (mid >> (32 - e));
            mid = (mid << e) | (lo  >> (32 - e));
            lo  = (lo  << e) | (tmp >> (32 - e));
        }

        p = (unsigned long long int)ia * lo;
        p = (unsigned long long int)ia * mid + (p >> 32);
        p = ((unsigned long long int)(ia * hi) << 32) + p;

    #if SANE_COMPILER
        q = (int)(p >> 62);                // integral portion = quadrant<1:0>
        p = p & 0x3fffffffffffffffULL;     // fraction
        if (p & 0x2000000000000000ULL) {   // fraction >= 0.5
            p = p - 0x4000000000000000ULL; // fraction - 1.0
            q = q + 1;
        }
    #else // !SANE_COMPILER
        unsigned int phi, plo;
        asm ("mov.b64 {%0,%1}, %2;" : "=r"(plo), "=r"(phi) : "l"(p));
        q = phi >> 30;
        phi = phi & 0x3fffffff;
        if (phi & 0x20000000) {
            phi = phi - 0x40000000;
            q = q + 1;
        }
        asm ("mov.b64 %0, {%1,%2};" : "=l"(p) : "r"(plo), "r"(phi));
    #endif // !SANE_COMPILER

        /* compute remainder of x / (pi/2) */
        double d;

        d = (double)(long long int)p;
        d = d * 3.4061215800865545e-19; // 0x1.921fb54442d18p-62 // pi/2 * 0x1.p-62
        r = (float)d;

        if (a < 0.0f) {
            r = -r;
            q = -q;
        }

        *quadrant = q;
        return r;
    }

    /* Argument reduction for trigonometric functions that reduces the argument
       to the interval [-PI/4, +PI/4] and also returns the quadrant. It returns 
       -0.0f for an input of -0.0f 
    */
    __forceinline__ __device__ float trig_red_f (float a, float switch_over, int *q)
    {    
        float j, r;

        /* FMA-enhanced Cody-Waite style reduction. W. J. Cody and W. Waite, 
           "Software Manual for the Elementary Functions", Prentice-Hall 1980
        */
        j = fmaf (a, 0.636619747f, 12582912.0f); // 2/pi, 0x1.8p+23
        *q = __float_as_int (j);
        j = j - 12582912.0f; // 0x1.8p+23
        r = fmaf (j, -1.57079601e+00f, a); // -0x1.921fb0p+00 // pio2_high
        r = fmaf (j, -3.13916473e-07f, r); // -0x1.5110b4p-22 // pio2_mid
        r = fmaf (j, -5.39030253e-15f, r); // -0x1.846988p-48 // pio2_low
        if (fabsf (a) > switch_over) {
            /* Payne-Hanek style reduction. M. Payne and R. Hanek, Radian reduction
               for trigonometric functions. SIGNUM Newsletter, 18:19-24, 1983
            */
            r = trig_red_slowpath_f (a, q);
        }
        return r;
    }

    /* Approximate sine on [-PI/4,+PI/4]. Maximum ulp error = 0.64196
       Returns -0.0f for an argument of -0.0f
       Polynomial approximation based on T. Myklebust, "Computing accurate 
       Horner form approximations to special functions in finite precision
       arithmetic", http://arxiv.org/abs/1508.03211, retrieved on 8/29/2016
    */
    __forceinline__ __device__ float sinf_poly (float a, float s)
    {
        float r, t;
        r =              2.83960253e-6f;  //  0x1.7d2000p-19
        r = fmaf (r, s, -1.98562644e-4f); // -0x1.a06a82p-13
        r = fmaf (r, s,  8.33339617e-3f); //  0x1.111198p-07
        r = fmaf (r, s, -1.66666672e-1f); // -0x1.555556p-03
        t = fmaf (a, s, 0.0f); // ensure -0 is passed trough
        r = fmaf (r, t, a);
        return r;
    }

    /* Approximate cosine on [-PI/4,+PI/4]. Maximum ulp error = 0.87444 */
    __forceinline__ __device__ float cosf_poly (float s)
    {
        float r;
        r =              2.43484974e-5f;  //  0x1.988000p-16
        r = fmaf (r, s, -1.38865795e-3f); // -0x1.6c0742p-10
        r = fmaf (r, s,  4.16666307e-2f); //  0x1.555542p-05
        r = fmaf (r, s, -5.00000000e-1f); // -0x1.000000p-01
        r = fmaf (r, s,  1.00000000e+0f); //  0x1.000000p+00
        return r;
    }

    /* Compute sine and cosine simultaneously, based on quadrant */
    __forceinline__ __device__ void scf_core (float a, int i, float *sp, float *cp)
    {
        float c, s, t;

        s = a * a;
        c = cosf_poly (s);
        s = sinf_poly (a, s);
        if (i & 2) {
            s = 0.0f - s; // don't change "sign" of NaNs or create negative zeros
            c = 0.0f - c; // don't change "sign" of NaNs or create negative zeros
        }
        if (i & 1) {
            t = 0.0f - s; // don't change "sign" of NaNs or create negative zeros
            s = c;
            c = t;
        }
        *sp = s;
        *cp = c;
    }
    
    #pragma acc routine seq nohost
    /* maximum ulp error sin = 1.49241, maximum ulp error cos = 1.49510 */
    __device__ void my_sincosf (float a, float *sp, float *cp)
    {
        float r;
        int i;

        a = a * 0.0f + a; // inf -> NaN
        r = trig_red_f (a, 71476.0625f, &i);
        scf_core (r, i, sp, cp);
    }
    }
    #pragma acc routine seq nohost


    
    #endif