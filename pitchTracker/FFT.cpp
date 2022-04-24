/*
    Reference:
        https://github.com/JorenSix/TarsosDSP/blob/master/src/core/be/tarsos/dsp/util/fft/FloatFFT.java
*/

#include "FFT.h"
#include <cmath>
#include <bits/sysconf.h>
#include "androidDebug.h"

#define PI      3.14159265358979311599796346854418516f
#define TWO_PI  6.28318530717958623199592693708837032f

#define THREADS_BEGIN_N_1D_FFT_2THREADS     8192

#ifdef USE_CDFT_PTHREADS
#ifndef USE_CDFT_THREADS
#define USE_CDFT_THREADS
#endif
#ifndef CDFT_THREADS_BEGIN_N
#define CDFT_THREADS_BEGIN_N 8192
#endif
#ifndef CDFT_4THREADS_BEGIN_N
#define CDFT_4THREADS_BEGIN_N 65536
#endif
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#define cdft_thread_t pthread_t
#define cdft_thread_create(thp,func,argp) { \
    if(pthread_create(thp, NULL, func, (void *) argp) != 0) \
        LOGE("cdft thread error"); \
}
#define cdft_thread_wait(th) { \
    if (pthread_join(th, NULL) != 0) \
        LOGE("cdft thread error"); \
}
#endif /* USE_CDFT_PTHREADS */

/*
 * Creates new instance of floatFFT.
 * @param n size of data
 */
FFT::FFT( int n )
    : m_n(n), m_pW(nullptr), m_pIp(nullptr), m_nc(0), m_nBluestein(0)
    , m_pBk1(nullptr), m_pBk2(nullptr), m_pWtable(nullptr), m_pWtable_r(nullptr)
{
    m_vecFactors.push_back(4);
    m_vecFactors.push_back(2);
    m_vecFactors.push_back(3);
    m_vecFactors.push_back(5);

    if(n < 1)
        LOGE("n must be greater than 0");

    if(!IsPowerOf2(n))
    {
        if(getReminder(n) >= 211)
        {
            m_plan = Plans::BLUESTEIN;
            m_nBluestein = nextPow2(n*2-1);
            m_pBk1 = new float[2*m_nBluestein];
            m_pBk2 = new float[2*m_nBluestein];
            m_pIp = new int[2+(int)ceil(2+(1<<(int)(log(m_nBluestein + 0.5)/log(2))/2))];
            m_pW = new float[m_nBluestein];
            int twon = 2*m_nBluestein;
            m_nw = m_pIp[0];
            if (twon > (m_nw << 2))
            {
                m_nw = twon >> 2;
                makewt(m_nw, m_pIp, m_pW);
            }

            m_nc = m_pIp[1];
            if(m_nBluestein > (m_nc << 2))
            {
                m_nc = m_nBluestein >> 2;
                makect(m_nc, m_pIp, m_pW+m_nw);
            }

            bluesteini();
        }
        else
        {
            m_plan = Plans::MIXED_RADIX;
            m_pWtable = new float[4*n+15];
            m_pWtable_r = new float[2*n+15];
            cffti();
            rffti();
        }
    }
    else
    {
        m_plan = Plans::SPLIT_RADIX;
        m_pIp = new int[2+(int)ceil(2+(1<<(int)(log(n+0.5)/log(2))/2))];
        m_pW = new float[n];
        int twon = 2 * n;

        m_nw = m_pIp[0];
        if(twon > (m_nw << 2))
        {
            m_nw = twon >> 2;
            makewt(m_nw, m_pIp, m_pW);
        }

        m_nc = m_pIp[1];
        if (n > (m_nc << 2))
        {
            m_nc = n >> 2;
            makect(m_nc, m_pIp, m_pW+m_nw);
        }
    }
}

FFT::~FFT()
{
    if( m_pW != nullptr )
    {
        delete[] m_pW;
        m_pW = nullptr;
    }

    if( m_pIp != nullptr )
    {
        delete[] m_pIp;
        m_pIp = nullptr;
    }

    if( m_pBk1 != nullptr )
    {
        delete[] m_pBk1;
        m_pBk1 = nullptr;
    }

    if( m_pBk2 != nullptr )
    {
        delete[] m_pBk2;
        m_pBk2 = nullptr;
    }

    if( m_pWtable != nullptr )
    {
        delete[] m_pWtable;
        m_pWtable = nullptr;
    }

    if( m_pWtable_r != nullptr )
    {
        delete[] m_pWtable_r;
        m_pWtable_r = nullptr;
    }
}

/**
 * Computes 1D forward DFT of complex data leaving the result in
 * <code>a</code>. Complex number is stored as two float values in
 * sequence: the real and imaginary part, i.e. the size of the input array
 * must be greater or equal 2*n. The physical layout of the input data has
 * to be as follows:<br>
 *
 * <pre>
 * pArr[2*k] = Re[k],
 * pArr[2*k+1] = Im[k], 0&lt;=k&lt;n
 * </pre>
 *
 * @param pArr
 *            data to transform
 */
void FFT::complexForward(float* pArr)
{
    if(m_n == 1)
        return;

    switch(m_plan)
    {
        case SPLIT_RADIX:
            cftbsub(2*m_n, pArr, m_pIp, m_nw, m_pW);
            break;
        case MIXED_RADIX:
            cfftf(pArr, -1);
            break;
        case BLUESTEIN:
            bluestein_complex(pArr, -1);
            break;
    }
}

/**
 * Computes 1D inverse DFT of complex data leaving the result in
 * <code>a</code>. Complex number is stored as two float values in
 * sequence: the real and imaginary part, i.e. the size of the input array
 * must be greater or equal 2*n. The physical layout of the input data has
 * to be as follows:<br>
 *
 * <pre>
 * pArr[2*k] = Re[k],
 * pArr[2*k+1] = Im[k], 0&lt;=k&lt;n
 * </pre>
 *
 * @param pArr
 *            data to transform
 * @param scale
 *            if true then scaling is performed
 */
void FFT::complexInverse(float* pArr, bool bScale)
{
    if(m_n == 1)
        return;

    switch(m_plan)
    {
        case SPLIT_RADIX:
            cftfsub(2*m_n, pArr, m_pIp, m_nw, m_pW);
            break;
        case MIXED_RADIX:
            cfftf(pArr, +1);
            break;
        case BLUESTEIN:
            bluestein_complex(pArr, 1);
            break;
    }
    if(bScale)
        scale(m_n, pArr, true);
}

void FFT::cfftf(float* pArr, int isign)
{
    int idot;
    int l1, l2;
    int na, nf, ip, iw, ido, idl1;
    int nac = 0;
    int twon = 2 * m_n;

    int iw1, iw2;
    float* pCh = new float[twon];

    iw1 = twon;
    iw2 = 4 * m_n;
    nf = (int)m_pWtable[1+iw2];
    na = 0;
    l1 = 1;
    iw = iw1;
    for(int k1 = 2; k1 <= nf + 1; ++k1)
    {
        ip = (int)m_pWtable[k1 + iw2];
        l2 = ip*l1;
        ido = m_n/l2;
        idot = ido+ido;
        idl1 = idot*l1;
        switch(ip)
        {
            case 4:
                if(na == 0)
                    passf4(idot, l1, pArr, pCh, iw, isign);
                else
                    passf4(idot, l1, pCh, pArr, iw, isign);
                na = 1 - na;
                break;
            case 2:
                if(na == 0)
                    passf2(idot, l1, pArr, pCh, iw, isign);
                else
                    passf2(idot, l1, pCh, pArr, iw, isign);
                na = 1 - na;
                break;
            case 3:
                if (na == 0) {
                    passf3(idot, l1, pArr, pCh, iw, isign);
                } else {
                    passf3(idot, l1, pCh, pArr, iw, isign);
                }
                na = 1 - na;
                break;
            case 5:
                if(na == 0)
                    passf5(idot, l1, pArr, pCh, iw, isign);
                else
                    passf5(idot, l1, pCh, pArr, iw, isign);
                na = 1 - na;
                break;
            default:
                if(na == 0)
                    passfg(&nac, idot, ip, l1, idl1, pArr, pCh, iw, isign);
                else
                    passfg(&nac, idot, ip, l1, idl1, pCh, pArr, iw, isign);

                if(nac != 0)
                    na = 1 - na;
                break;
        }
        l1 = l2;
        iw += (ip - 1) * idot;
    }

    if(na != 0)
    {
        if( pCh != nullptr )
            memcpy( pArr, pCh, sizeof(float)*twon );
    }

    if( pCh != nullptr )
    {
        delete[] pCh;
        pCh = nullptr;
    }
}

/*----------------------------------------------------------------------
   passf2: Complex FFT's forward/backward processing of factor 2;
   isign is +1 for backward and -1 for forward transforms
  ----------------------------------------------------------------------*/

void FFT::passf2( int ido, int l1, float* pIn, float* pOut, int offset, int isign)
{
    float t1i, t1r;
    int iw1;
    iw1 = offset;
    int idx = ido * l1;
    if (ido <= 2)
    {
        for (int k = 0; k < l1; ++k)
        {
            int idx0 = k*ido;
            int iidx1 = 2*idx0;
            int iidx2 = iidx1+ido;
            float a1r = pIn[iidx1];
            float a1i = pIn[iidx1+1];
            float a2r = pIn[iidx2];
            float a2i = pIn[iidx2+1];

            int oidx1 = idx0+idx;
            pOut[idx0] = a1r+a2r;
            pOut[idx0+1] = a1i+a2i;
            pOut[oidx1] = a1r-a2r;
            pOut[oidx1+1] = a1i-a2i;
        }
    }
    else
    {
        for (int k = 0; k < l1; ++k)
        {
            for (int i = 0; i < ido - 1; i += 2)
            {
                int idx0 = k * ido;
                int iidx1 = i+2*idx0;
                int iidx2 = iidx1+ido;
                float i1r = pIn[iidx1];
                float i1i = pIn[iidx1+1];
                float i2r = pIn[iidx2];
                float i2i = pIn[iidx2+1];

                int widx1 = i + iw1;
                float w1r = m_pWtable[widx1];
                float w1i = isign * m_pWtable[widx1 + 1];

                t1r = i1r - i2r;
                t1i = i1i - i2i;

                int oidx1 = i+idx0;
                int oidx2 = oidx1+idx;
                pOut[oidx1] = i1r+i2r;
                pOut[oidx1+1] = i1i+i2i;
                pOut[oidx2] = w1r*t1r-w1i*t1i;
                pOut[oidx2+1] = w1r*t1i+w1i*t1r;
            }
        }
    }
}

/*----------------------------------------------------------------------
   passf3: Complex FFT's forward/backward processing of factor 3;
   isign is +1 for backward and -1 for forward transforms
  ----------------------------------------------------------------------*/
void FFT::passf3(int ido, int l1, float* pIn, float* pOut, int offset, int isign)
{
    float taur = -0.5f;
    float taui = 0.866025403784438707610604524234076962f;
    float ci2, ci3, di2, di3, cr2, cr3, dr2, dr3, ti2, tr2;
    int iw1, iw2;

    iw1 = offset;
    iw2 = iw1 + ido;

    int idxt = l1 * ido;

    if(ido == 2)
    {
        for(int k = 1; k <= l1; ++k)
        {
            int iidx1 = (3*k-2)*ido;
            int iidx2 = iidx1+ido;
            int iidx3 = iidx1-ido;
            float i1r = pIn[iidx1];
            float i1i = pIn[iidx1+1];
            float i2r = pIn[iidx2];
            float i2i = pIn[iidx2+1];
            float i3r = pIn[iidx3];
            float i3i = pIn[iidx3+1];

            tr2 = i1r+i2r;
            cr2 = i3r+taur*tr2;
            ti2 = i1i+i2i;
            ci2 = i3i+taur*ti2;
            cr3 = isign*taui*(i1r-i2r);
            ci3 = isign*taui*(i1i-i2i);

            int oidx1 = (k-1)*ido;
            int oidx2 = oidx1+idxt;
            int oidx3 = oidx2+idxt;
            pOut[oidx1] = pIn[iidx3]+tr2;
            pOut[oidx1+1] = i3i+ti2;
            pOut[oidx2] = cr2-ci3;
            pOut[oidx2+1] = ci2+cr3;
            pOut[oidx3] = cr2+ci3;
            pOut[oidx3+1] = ci2-cr3;
        }
    }
    else
    {
        for(int k = 1; k <= l1; ++k)
        {
            int idx1 = (3*k-2)*ido;
            int idx2 = (k-1)*ido;
            for(int idx = 0; idx < ido - 1; idx += 2)
            {
                int iidx1 = idx+idx1;
                int iidx2 = iidx1+ido;
                int iidx3 = iidx1-ido;
                float a1r = pIn[iidx1];
                float a1i = pIn[iidx1+1];
                float a2r = pIn[iidx2];
                float a2i = pIn[iidx2+1];
                float a3r = pIn[iidx3];
                float a3i = pIn[iidx3+1];

                tr2 = a1r+a2r;
                cr2 = a3r+taur*tr2;
                ti2 = a1i+a2i;
                ci2 = a3i+taur*ti2;
                cr3 = isign*taui*(a1r-a2r);
                ci3 = isign*taui*(a1i-a2i);
                dr2 = cr2-ci3;
                dr3 = cr2+ci3;
                di2 = ci2+cr3;
                di3 = ci2-cr3;

                int widx1 = idx+iw1;
                int widx2 = idx+iw2;
                float w1r = m_pWtable[widx1];
                float w1i = isign * m_pWtable[widx1 + 1];
                float w2r = m_pWtable[widx2];
                float w2i = isign * m_pWtable[widx2 + 1];

                int oidx1 = idx+idx2;
                int oidx2 = oidx1+idxt;
                int oidx3 = oidx2+idxt;
                pOut[oidx1] = a3r+tr2;
                pOut[oidx1+1] = a3i+ti2;
                pOut[oidx2] = w1r*dr2-w1i*di2;
                pOut[oidx2+1] = w1r*di2+w1i*dr2;
                pOut[oidx3] = w2r*dr3-w2i*di3;
                pOut[oidx3+1] = w2r*di3+w2i*dr3;
            }
        }
    }
}

/*----------------------------------------------------------------------
   passf4: Complex FFT's forward/backward processing of factor 4;
   isign is +1 for backward and -1 for forward transforms
  ----------------------------------------------------------------------*/
void FFT::passf4(int ido, int l1, float* pIn, float* pOut, int offset, int isign)
{
    float ci2, ci3, ci4, cr2, cr3, cr4, ti1, ti2, ti3, ti4, tr1, tr2, tr3, tr4;
    int iw1, iw2, iw3;
    iw1 = offset;
    iw2 = iw1 + ido;
    iw3 = iw2 + ido;

    int idx0 = l1 * ido;
    if(ido == 2)
    {
        for(int k = 0; k < l1; ++k)
        {
            int idxt1 = k*ido;
            int iidx1 = 4*idxt1+1;
            int iidx2 = iidx1+ido;
            int iidx3 = iidx2+ido;
            int iidx4 = iidx3+ido;

            float i1i = pIn[iidx1-1];
            float i1r = pIn[iidx1];
            float i2i = pIn[iidx2-1];
            float i2r = pIn[iidx2];
            float i3i = pIn[iidx3-1];
            float i3r = pIn[iidx3];
            float i4i = pIn[iidx4-1];
            float i4r = pIn[iidx4];

            ti1 = i1r-i3r;
            ti2 = i1r+i3r;
            tr4 = i4r-i2r;
            ti3 = i2r+i4r;
            tr1 = i1i-i3i;
            tr2 = i1i+i3i;
            ti4 = i2i-i4i;
            tr3 = i2i+i4i;

            int oidx1 = idxt1+idx0;
            int oidx2 = oidx1+idx0;
            int oidx3 = oidx2+idx0;
            pOut[idxt1] = tr2+tr3;
            pOut[idxt1+1] = ti2+ti3;
            pOut[oidx1] = tr1+isign*tr4;
            pOut[oidx1+1] = ti1+isign*ti4;
            pOut[oidx2] = tr2-tr3;
            pOut[oidx2+1] = ti2-ti3;
            pOut[oidx3] = tr1-isign*tr4;
            pOut[oidx3+1] = ti1-isign*ti4;
        }
    }
    else
    {
        for(int k = 0; k < l1; ++k)
        {
            int idx1 = k*ido;
            int idx2 = 1+4*idx1;
            for(int idx = 0; idx < ido - 1; idx += 2)
            {
                int iidx1 = idx+idx2;
                int iidx2 = iidx1+ido;
                int iidx3 = iidx2+ido;
                int iidx4 = iidx3+ido;
                float i1i = pIn[iidx1-1];
                float i1r = pIn[iidx1];
                float i2i = pIn[iidx2-1];
                float i2r = pIn[iidx2];
                float i3i = pIn[iidx3-1];
                float i3r = pIn[iidx3];
                float i4i = pIn[iidx4-1];
                float i4r = pIn[iidx4];

                ti1 = i1r-i3r;
                ti2 = i1r+i3r;
                ti3 = i2r+i4r;
                tr4 = i4r-i2r;
                tr1 = i1i-i3i;
                tr2 = i1i+i3i;
                ti4 = i2i-i4i;
                tr3 = i2i+i4i;
                cr3 = tr2-tr3;
                ci3 = ti2-ti3;
                cr2 = tr1+isign*tr4;
                cr4 = tr1-isign*tr4;
                ci2 = ti1+isign*ti4;
                ci4 = ti1-isign*ti4;

                int widx1 = idx+iw1;
                int widx2 = idx+iw2;
                int widx3 = idx+iw3;
                float w1r = m_pWtable[widx1];
                float w1i = isign * m_pWtable[widx1+1];
                float w2r = m_pWtable[widx2];
                float w2i = isign * m_pWtable[widx2+1];
                float w3r = m_pWtable[widx3];
                float w3i = isign * m_pWtable[widx3+1];

                int oidx1 = idx+idx1;
                int oidx2 = oidx1+idx0;
                int oidx3 = oidx2+idx0;
                int oidx4 = oidx3+idx0;
                pOut[oidx1] = tr2+tr3;
                pOut[oidx1+1] = ti2+ti3;
                pOut[oidx2] = w1r*cr2-w1i*ci2;
                pOut[oidx2+1] = w1r*ci2+w1i*cr2;
                pOut[oidx3] = w2r*cr3-w2i*ci3;
                pOut[oidx3+1] = w2r*ci3+w2i*cr3;
                pOut[oidx4] = w3r*cr4-w3i*ci4;
                pOut[oidx4+1] = w3r*ci4+w3i*cr4;
            }
        }
    }
}

/*----------------------------------------------------------------------
   passf5: Complex FFT's forward/backward processing of factor 5;
   isign is +1 for backward and -1 for forward transforms
  ----------------------------------------------------------------------*/
void FFT::passf5(int ido, int l1, float* pIn, float* pOut, int offset, int isign)
/* isign==-1 for forward transform and+1 for backward transform */
{
    float tr11 = 0.309016994374947451262869435595348477f;
    float ti11 = 0.951056516295153531181938433292089030f;
    float tr12 = -0.809016994374947340240566973079694435f;
    float ti12 = 0.587785252292473248125759255344746634f;
    float ci2, ci3, ci4, ci5, di3, di4, di5, di2, cr2, cr3, cr5, cr4, ti2, ti3, ti4, ti5, dr3, dr4, dr5, dr2, tr2, tr3, tr4, tr5;
    int iw1, iw2, iw3, iw4;

    iw1 = offset;
    iw2 = iw1 + ido;
    iw3 = iw2 + ido;
    iw4 = iw3 + ido;

    int idx0 = l1 * ido;

    if(ido == 2)
    {
        for (int k = 1; k <= l1; ++k)
        {
            int iidx1 = (5*k-4)*ido+1;
            int iidx2 = iidx1+ido;
            int iidx3 = iidx1-ido;
            int iidx4 = iidx2+ido;
            int iidx5 = iidx4+ido;

            float i1i = pIn[iidx1-1];
            float i1r = pIn[iidx1];
            float i2i = pIn[iidx2-1];
            float i2r = pIn[iidx2];
            float i3i = pIn[iidx3-1];
            float i3r = pIn[iidx3];
            float i4i = pIn[iidx4-1];
            float i4r = pIn[iidx4];
            float i5i = pIn[iidx5-1];
            float i5r = pIn[iidx5];

            ti5 = i1r-i5r;
            ti2 = i1r+i5r;
            ti4 = i2r-i4r;
            ti3 = i2r+i4r;
            tr5 = i1i-i5i;
            tr2 = i1i+i5i;
            tr4 = i2i-i4i;
            tr3 = i2i+i4i;
            cr2 = i3i+tr11*tr2+tr12*tr3;
            ci2 = i3r+tr11*ti2+tr12*ti3;
            cr3 = i3i+tr12*tr2+tr11*tr3;
            ci3 = i3r+tr12*ti2+tr11*ti3;
            cr5 = isign*(ti11*tr5+ti12*tr4);
            ci5 = isign*(ti11*ti5+ti12*ti4);
            cr4 = isign*(ti12*tr5-ti11*tr4);
            ci4 = isign*(ti12*ti5-ti11*ti4);

            int oidx1 = (k-1)*ido;
            int oidx2 = oidx1+idx0;
            int oidx3 = oidx2+idx0;
            int oidx4 = oidx3+idx0;
            int oidx5 = oidx4+idx0;
            pOut[oidx1] = i3i+tr2+tr3;
            pOut[oidx1+1] = i3r+ti2+ti3;
            pOut[oidx2] = cr2-ci5;
            pOut[oidx2+1] = ci2+cr5;
            pOut[oidx3] = cr3-ci4;
            pOut[oidx3+1] = ci3+cr4;
            pOut[oidx4] = cr3+ci4;
            pOut[oidx4+1] = ci3-cr4;
            pOut[oidx5] = cr2+ci5;
            pOut[oidx5+1] = ci2-cr5;
        }
    }
    else
    {
        for(int k = 1; k <= l1; ++k)
        {
            int idx1 = 1+(k*5-4)*ido;
            int idx2 = (k-1)*ido;
            for(int idx = 0; idx < ido - 1; idx += 2)
            {
                int iidx1 = idx+idx1;
                int iidx2 = iidx1+ido;
                int iidx3 = iidx1-ido;
                int iidx4 = iidx2+ido;
                int iidx5 = iidx4+ido;
                float i1i = pIn[iidx1-1];
                float i1r = pIn[iidx1];
                float i2i = pIn[iidx2-1];
                float i2r = pIn[iidx2];
                float i3i = pIn[iidx3-1];
                float i3r = pIn[iidx3];
                float i4i = pIn[iidx4-1];
                float i4r = pIn[iidx4];
                float i5i = pIn[iidx5-1];
                float i5r = pIn[iidx5];

                ti5 = i1r-i5r;
                ti2 = i1r+i5r;
                ti4 = i2r-i4r;
                ti3 = i2r+i4r;
                tr5 = i1i-i5i;
                tr2 = i1i+i5i;
                tr4 = i2i-i4i;
                tr3 = i2i+i4i;
                cr2 = i3i+tr11*tr2+tr12*tr3;
                ci2 = i3r+tr11*ti2+tr12*ti3;
                cr3 = i3i+tr12*tr2+tr11*tr3;
                ci3 = i3r+tr12*ti2+tr11*ti3;
                cr5 = isign*(ti11*tr5+ti12*tr4);
                ci5 = isign*(ti11*ti5+ti12*ti4);
                cr4 = isign*(ti12*tr5-ti11*tr4);
                ci4 = isign*(ti12*ti5-ti11*ti4);
                dr3 = cr3-ci4;
                dr4 = cr3+ci4;
                di3 = ci3+cr4;
                di4 = ci3-cr4;
                dr5 = cr2+ci5;
                dr2 = cr2-ci5;
                di5 = ci2-cr5;
                di2 = ci2+cr5;

                int widx1 = idx+iw1;
                int widx2 = idx+iw2;
                int widx3 = idx+iw3;
                int widx4 = idx+iw4;
                float w1r = m_pWtable[widx1];
                float w1i = isign * m_pWtable[widx1+1];
                float w2r = m_pWtable[widx2];
                float w2i = isign * m_pWtable[widx2+1];
                float w3r = m_pWtable[widx3];
                float w3i = isign * m_pWtable[widx3+1];
                float w4r = m_pWtable[widx4];
                float w4i = isign * m_pWtable[widx4+1];

                int oidx1 = idx+idx2;
                int oidx2 = oidx1+idx0;
                int oidx3 = oidx2+idx0;
                int oidx4 = oidx3+idx0;
                int oidx5 = oidx4+idx0;
                pOut[oidx1] = i3i+tr2+tr3;
                pOut[oidx1+1] = i3r+ti2+ti3;
                pOut[oidx2] = w1r*dr2-w1i*di2;
                pOut[oidx2+1] = w1r*di2+w1i*dr2;
                pOut[oidx3] = w2r*dr3-w2i*di3;
                pOut[oidx3+1] = w2r*di3+w2i*dr3;
                pOut[oidx4] = w3r*dr4-w3i*di4;
                pOut[oidx4+1] = w3r*di4+w3i*dr4;
                pOut[oidx5] = w4r*dr5-w4i*di5;
                pOut[oidx5+1] = w4r*di5+w4i*dr5;
            }
        }
    }
}

/*----------------------------------------------------------------------
   passfg: Complex FFT's forward/backward processing of general factor;
   isign is +1 for backward and -1 for forward transforms
  ----------------------------------------------------------------------*/
void FFT::passfg(int* pNac, int ido, int ip, int l1, int idl1, float* pIn
                , float* pOut, int offset, int isign)
{
    int idij, idlj, idot, ipph, l, jc, lc, idj, idl, inc, idp;
    float w1r, w1i, w2i, w2r;
    int iw1;

    iw1 = offset;
    idot = ido/2;
    ipph = (ip+1)/2;
    idp = ip*ido;

    if(ido >= l1)
    {
        for(int j = 1; j < ipph; ++j)
        {
            jc = ip-j;
            int idx1 = j*ido;
            int idx2 = jc*ido;
            for(int k = 0; k < l1; ++k)
            {
                int idx3 = k*ido;
                int idx4 = idx3+idx1*l1;
                int idx5 = idx3+idx2*l1;
                int idx6 = idx3*ip;
                for(int idx = 0; idx < ido; ++idx)
                {
                    float i1r = pIn[idx+idx1+idx6];
                    float i2r = pIn[idx+idx2+idx6];
                    pOut[idx+idx4] = i1r+i2r;
                    pOut[idx+idx5] = i1r-i2r;
                }
            }
        }
        for(int k = 0; k < l1; ++k)
        {
            int idxt1 = k*ido;
            int idxt2 = idxt1*ip;
            for(int idx = 0; idx < ido; ++idx)
                pOut[idx+idxt1] = pIn[idx+idxt2];
        }
    }
    else
    {
        for(int j = 1; j < ipph; ++j)
        {
            jc = ip-j;
            int idxt1 = j*l1*ido;
            int idxt2 = jc*l1*ido;
            int idxt3 = j*ido;
            int idxt4 = jc*ido;
            for(int idx = 0; idx < ido; ++idx)
            {
                for(int k = 0; k < l1; ++k)
                {
                    int idx1 = k*ido;
                    int idx2 = idx1*ip;
                    float i1r = pIn[idx+idxt3+idx2];
                    float i2r = pIn[idx+idxt4+idx2];
                    pOut[idx+idx1+idxt1] = i1r+i2r;
                    pOut[idx+idx1+idxt2] = i1r-i2r;
                }
            }
        }
        for(int idx = 0; idx < ido; ++idx)
        {
            for(int k = 0; k < l1; ++k)
            {
                int idx1 = k*ido;
                pOut[idx+idx1] = pIn[idx+idx1*ip];
            }
        }
    }

    idl = 2-ido;
    inc = 0;
    int idxt0 = (ip-1)*idl1;
    for(l = 1; l < ipph; ++l)
    {
        lc = ip-l;
        idl += ido;
        int idxt1 = l*idl1;
        int idxt2 = lc*idl1;
        int idxt3 = idl+iw1;
        w1r = m_pWtable[idxt3-2];
        w1i = isign * m_pWtable[idxt3-1];
        for(int ik = 0; ik < idl1; ++ik)
        {
            pIn[ik+idxt1] = pOut[ik]+w1r*pOut[ik+idl1];
            pIn[ik+idxt2] = w1i*pOut[ik+idxt0];
        }
        idlj = idl;
        inc += ido;
        for(int j = 2; j < ipph; ++j)
        {
            jc = ip-j;
            idlj += inc;
            if (idlj > idp)
                idlj -= idp;
            int idxt4 = idlj+iw1;
            w2r = m_pWtable[idxt4-2];
            w2i = isign*m_pWtable[idxt4-1];
            int idxt5 = j*idl1;
            int idxt6 = jc*idl1;
            for(int ik = 0; ik < idl1; ++ik)
            {
                pIn[ik+idxt1] += w2r*pOut[ik+idxt5];
                pIn[ik+idxt2] += w2i*pOut[ik+idxt6];
            }
        }
    }

    for(int j = 1; j < ipph; ++j)
    {
        int idxt1 = j*idl1;
        for(int ik = 0; ik < idl1; ++ik)
        {
            pOut[ik] += pOut[ik+idxt1];
        }
    }

    for(int j = 1; j < ipph; ++j)
    {
        jc = ip-j;
        int idx1 = j*idl1;
        int idx2 = jc*idl1;
        for(int ik = 1; ik < idl1; ik += 2)
        {
            int iidx1 = ik+idx1;
            int iidx2 = ik+idx2;
            float i1i = pIn[iidx1-1];
            float i1r = pIn[iidx1];
            float i2i = pIn[iidx2-1];
            float i2r = pIn[iidx2];

            int oidx1 = ik+idx1;
            int oidx2 = ik+idx2;
            pOut[oidx1-1] = i1i-i2r;
            pOut[oidx2-1] = i1i+i2r;
            pOut[oidx1] = i1r+i2i;
            pOut[oidx2] = i1r-i2i;
        }
    }
    pNac[0] = 1;
    if (ido == 2)
        return;
    pNac[0] = 0;
    memcpy( pIn, pOut, sizeof(float)*idl1 );
    int idx0 = l1*ido;
    for(int j = 1; j < ip; ++j)
    {
        int idx1 = j*idx0;
        for(int k = 0; k < l1; ++k)
        {
            int idx2 = k*ido;
            int idx3 = idx2+idx1;
            pIn[idx3] = pOut[idx3];
            pIn[idx3+1] = pOut[idx3+1];
        }
    }
    if(idot <= l1)
    {
        idij = 0;
        for(int j = 1; j < ip; ++j)
        {
            idij += 2;
            int idx1 = j*l1*ido;
            for(int idx = 3; idx < ido; idx += 2)
            {
                idij += 2;
                int idx2 = idij+iw1-1;
                w1r = m_pWtable[idx2-1];
                w1i = isign*m_pWtable[idx2];
                for(int k = 0; k < l1; ++k)
                {
                    int idx3 = idx+(k*ido+idx1);
                    float o1i = pOut[idx3-1];
                    float o1r = pOut[idx3];
                    pIn[idx3-1] = w1r*o1i-w1i*o1r;
                    pIn[idx3] = w1r*o1r+w1i*o1i;
                }
            }
        }
    }
    else
    {
        idj = 2-ido;
        for(int j = 1; j < ip; ++j)
        {
            idj += ido;
            int idx1 = j*l1*ido;
            for(int k = 0; k < l1; ++k)
            {
                idij = idj;
                int idx3 = k*ido+idx1;
                for(int idx = 3; idx < ido; idx += 2)
                {
                    idij += 2;
                    int idx2 = idij-1+iw1;
                    w1r = m_pWtable[idx2-1];
                    w1i = isign * m_pWtable[idx2];
                    int idx4 = idx+idx3;
                    float o1i = pOut[idx4-1];
                    float o1r = pOut[idx4];
                    pIn[idx4-1] = w1r*o1i-w1i*o1r;
                    pIn[idx4] = w1r*o1r+w1i*o1i;
                }
            }
        }
    }
}

bool FFT::IsPowerOf2( int x )
{
    if( x <= 0 )
        return false;

    return ((x & (x-1)) == 0);
}

int FFT::prevPow2( int x )
{
    if( x < 1 )
    {
        LOGE( "x must be greater or equal 1" );
        return 0;
    }

    return (int)pow(2, floor(log(x)/log(2)));
}

int FFT::nextPow2( int x )
{
    if( x < 1 )
    {
        LOGE( "x must be greater or equal 1" );
        return 0;
    }

    if( (x & (x-1)) == 0 )
    {
        return x; // x is already a power-of-two number
    }
    x |= (x >> 1);
    x |= (x >> 2);
    x |= (x >> 4);
    x |= (x >> 8);
    x |= (x >> 16);
    x |= (int)((int64_t)x >> 32);

    return x+1;
}

int FFT::getReminder(int n)
{
    int reminder = n;

    if(n <= 0)
        LOGE("n must be positive integer");

    for(int idx = 0; idx < m_vecFactors.size() && reminder != 1; ++idx)
    {
        int factor = m_vecFactors[idx];
        while ((reminder % factor) == 0)
            reminder/=factor;
    }

    return reminder;
}

void FFT::cffti(int n, int offw)
{
    if (n == 1)
        return;

    int twon = 2*n;
    int fourn = 4*n;
    float argh;
    int idot, ntry = 0, i = 0;
    float argld;
    int i1, k1, l1, l2, ib;
    float fi;
    int ld, ii, ip, nq, nr;
    float arg;
    int ido, ipm;

    int nl = n;
    int nf = 0;
    int j = 0;

    while(true)
    {
        ++j;
        if (j <= 4)
            ntry = m_vecFactors[j-1];
        else
            ntry += 2;

        bool bContinue = false;
        do
        {
            nq = nl/ntry;
            nr = nl-ntry*nq;

            if (nr != 0)
            {
                bContinue = true;
                break;
            }
            nf++;
            m_pWtable[offw+nf+1+fourn] = ntry;
            nl = nq;

            if(ntry == 2 && nf != 1)
            {
                for(i = 2; i <= nf; ++i)
                {
                    ib = nf-i+2;
                    int idx = ib+fourn;
                    m_pWtable[offw+idx+1] = m_pWtable[offw+idx];
                }
                m_pWtable[offw+2+fourn] = 2;
            }
        } while (nl != 1);

        if( bContinue == false )
            break;
    }

    m_pWtable[offw+fourn] = n;
    m_pWtable[offw+1+fourn] = nf;
    argh = TWO_PI/(float)n;
    i = 1;
    l1 = 1;

    for (k1 = 1; k1 <= nf; ++k1)
    {
        ip = (int)m_pWtable[offw+k1+1+fourn];
        ld = 0;
        l2 = l1*ip;
        ido = n/l2;
        idot = ido + ido + 2;
        ipm = ip - 1;

        for(j = 1; j <= ipm; ++j)
        {
            i1 = i;
            m_pWtable[offw+i-1+twon] = 1;
            m_pWtable[offw+i+twon] = 0;
            ld += l1;
            fi = 0;
            argld = ld * argh;
            for (ii = 4; ii <= idot; ii += 2)
            {
                i += 2;
                fi += 1;
                arg = fi*argld;
                int idx = i+twon;
                m_pWtable[offw+idx-1] = (float)cos(arg);
                m_pWtable[offw+idx] = (float)sin(arg);
            }

            if (ip > 5)
            {
                int idx1 = i1 + twon;
                int idx2 = i + twon;
                m_pWtable[offw+idx1-1] = m_pWtable[offw+idx2-1];
                m_pWtable[offw+idx1] = m_pWtable[offw+idx2];
            }
        }

        l1 = l2;
    }

}

void FFT::cffti()
{
    if (m_n == 1)
        return;

    int twon = 2*m_n;
    int fourn = 4*m_n;
    float argh;
    int idot, ntry = 0;
    float argld;
    int i1, k1, l1, l2, ib;
    float fi;
    int ld, ii, ip, nq, nr;
    float arg;
    int ido, ipm;

    int nl = m_n;
    int nf = 0;
    int j = 0;

    while (true)
    {
        ++j;
        if (j <= 4)
            ntry = m_vecFactors[j-1];
        else
            ntry += 2;

        bool bContinue = false;
        do
        {
            nq = nl/ntry;
            nr = nl- ntry*nq;
            if(nr != 0)
            {
                bContinue = true;
                break;
            }

            ++nf;
            m_pWtable[nf+1+fourn] = ntry;
            nl = nq;

            if (ntry == 2 && nf != 1)
            {
                for(int i = 2; i <= nf; ++i)
                {
                    ib = nf - i+2;
                    int idx = ib+fourn;
                    m_pWtable[idx+1] = m_pWtable[idx];
                }
                m_pWtable[2+fourn] = 2;
            }
        } while (nl != 1);

        if( bContinue == false )
            break;
    }

    m_pWtable[fourn] = m_n;
    m_pWtable[1+fourn] = nf;
    argh = TWO_PI/(float) m_n;
    int i = 1;
    l1 = 1;

    for (k1 = 1; k1 <= nf; ++k1)
    {
        ip = (int)m_pWtable[k1+1+fourn];
        ld = 0;
        l2 = l1*ip;
        ido = m_n/l2;
        idot = ido+ido+2;
        ipm = ip-1;

        for(j = 1; j <= ipm; j++)
        {
            i1 = i;
            m_pWtable[i-1+twon] = 1;
            m_pWtable[i+twon] = 0;
            ld += l1;
            fi = 0;
            argld = ld * argh;
            for(ii = 4; ii <= idot; ii += 2)
            {
                i += 2;
                fi += 1;
                arg = fi*argld;
                int idx = i+twon;
                m_pWtable[idx-1] = (float)cos(arg);
                m_pWtable[idx] = (float)sin(arg);
            }

            if(ip > 5)
            {
                int idx1 = i1+twon;
                int idx2 = i+twon;
                m_pWtable[idx1-1] = m_pWtable[idx2-1];
                m_pWtable[idx1] = m_pWtable[idx2];
            }
        }

        l1 = l2;
    }

}

void FFT::rffti()
{
    if(m_n == 1)
        return;

    int twon = 2*m_n;
    float argh;
    int ntry = 0, i = 0;
    float argld;
    int k1, l1, l2, ib;
    float fi;
    int ld, ii, ip, is, nq, nr;
    float arg;
    int ido, ipm;
    int nfm1;

    int nl = m_n;
    int nf = 0;
    int j = 0;

    while(true)
    {
        ++j;
        if (j <= 4)
            ntry = m_vecFactors[j-1];
        else
            ntry += 2;

        bool bContinue = false;

        do
        {
            nq = nl/ntry;
            nr = nl-ntry*nq;
            if (nr != 0)
            {
                bContinue = true;
                break;
            }

            ++nf;
            m_pWtable_r[nf+1+twon] = ntry;

            nl = nq;

            if(ntry == 2 && nf != 1)
            {
                for(i = 2; i <= nf; ++i)
                {
                    ib = nf-i+2;
                    int idx = ib+twon;
                    m_pWtable_r[idx+1] = m_pWtable_r[idx];
                }
                m_pWtable_r[2+twon] = 2;
            }
        } while (nl != 1);

        if( bContinue == false )
            break;

    }

    m_pWtable_r[twon] = m_n;
    m_pWtable_r[1+twon] = nf;
    argh = TWO_PI/(float)(m_n);
    is = 0;
    nfm1 = nf - 1;
    l1 = 1;

    if (nfm1 == 0)
        return;

    for (k1 = 1; k1 <= nfm1; ++k1)
    {
        ip = (int) m_pWtable_r[k1+1+twon];
        ld = 0;
        l2 = l1*ip;
        ido = m_n/l2;
        ipm = ip-1;

        for(j = 1; j <= ipm; ++j)
        {
            ld += l1;
            i = is;
            argld = (float) ld * argh;

            fi = 0;
            for (ii = 3; ii <= ido; ii += 2)
            {
                i += 2;
                fi += 1;
                arg = fi * argld;
                int idx = i + m_n;
                m_pWtable_r[idx-2] = (float)cos(arg);
                m_pWtable_r[idx-1] = (float)sin(arg);
            }
            is += ido;
        }
        l1 = l2;
    }
}

void FFT::bluesteini()
{
    int k = 0;
    float arg;
    float pi_n = PI/m_n;
    m_pBk1[0] = 1;
    m_pBk1[1] = 0;

    for (int idx = 1; idx < m_n; ++idx)
    {
        k += 2*idx-1;
        if (k >= 2*m_n)
            k -= 2*m_n;
        arg = pi_n*k;
        m_pBk1[2*idx] = (float)cos(arg);
        m_pBk1[2*idx+1] = (float)sin(arg);
    }

    float scale = (float)(1.0/m_nBluestein);
    m_pBk2[0] = m_pBk1[0]*scale;
    m_pBk2[1] = m_pBk1[1]*scale;

    for(int idx = 2; idx < 2 * m_n; idx += 2)
    {
        m_pBk2[idx] = m_pBk1[idx]*scale;
        m_pBk2[idx+1] = m_pBk1[idx+1]*scale;
        m_pBk2[2*m_nBluestein-idx] = m_pBk2[idx];
        m_pBk2[2*m_nBluestein-idx+1] = m_pBk2[idx+1];
    }

    cftbsub(2*m_nBluestein, m_pBk2, m_pIp, m_nw, m_pW);
}

void FFT::cdft(int n, int isgn, float* pArr, int* pIp, float* pW)
{
    int nw = pIp[0];

    if(n > (nw << 2))
    {
        nw = n >> 2;
        makewt(nw, pIp, pW);
    }

    if(isgn >= 0)
        cftfsub(n, pArr, pIp, nw, pW);
    else
        cftbsub(n, pArr, pIp, nw, pW);
}

void FFT::rdft( int n, int isgn, float* pArr, int* pIp, float* pW )
{
    int nw = pIp[0];

    if(n > (nw << 2))
    {
        nw = n >> 2;
        makewt(nw, pIp, pW);
    }

    int nc = pIp[1];

    if(n > (nc << 2))
    {
        nc = n >> 2;
        makect(nc, pIp, pW + nw);
    }

    if(isgn >= 0)
    {
        if (n > 4)
        {
            cftfsub(n, pArr, pIp, nw, pW);
            rftfsub(n, pArr, nc, pW + nw);
        }
        else if(n == 4)
            cftfsub(n, pArr, pIp, nw, pW);

        float xi = pArr[0] - pArr[1];
        pArr[0] += pArr[1];
        pArr[1] = xi;
    }
    else
    {
        pArr[1] = 0.5 * (pArr[0] - pArr[1]);
        pArr[0] -= pArr[1];

        if(n > 4)
        {
            rftbsub( n, pArr, nc, pW+nw );
            cftbsub( n, pArr, pIp, nw, pW );
        }
        else if(n == 4)
            cftbsub(n, pArr, pIp, nw, pW);
    }
}

void FFT::ddct( int n, int isgn, float* pArr, int* pIp, float* pW )
{
    int nw = pIp[0];

    if(n > (nw << 2))
    {
        nw = n >> 2;
        makewt(nw, pIp, pW);
    }

    int nc = pIp[1];

    if(n > nc)
    {
        nc = n;
        makect(nc, pIp, pW + nw);
    }

    if(isgn < 0)
    {
        float xr = pArr[n-1];

        for(int j = n-2; j >= 2; j -= 2)
        {
            pArr[j+1] = pArr[j] - pArr[j-1];
            pArr[j] += pArr[j-1];
        }

        pArr[1] = pArr[0] - xr;
        pArr[0] += xr;

        if(n > 4)
        {
            rftbsub(n, pArr, nc, pW+nw);
            cftbsub(n, pArr, pIp, nw, pW);
        }
        else if(n == 4)
            cftbsub(n, pArr, pIp, nw, pW);
    }

    dctsub(n, pArr, nc, pW + nw);

    if(isgn >= 0)
    {
        if(n > 4)
        {
            cftfsub(n, pArr, pIp, nw, pW);
            rftfsub(n, pArr, nc, pW+nw);
        }
        else if(n == 4)
        {
            cftfsub(n, pArr, pIp, nw, pW);
        }

        float xr = pArr[0] - pArr[1];
        pArr[0] += pArr[1];

        for(int j = 2; j < n; j += 2)
        {
            pArr[j-1] = pArr[j] - pArr[j+1];
            pArr[j] += pArr[j+1];
        }

        pArr[n-1] = xr;
    }
}

void FFT::ddst( int n, int isgn, float* pArr, int* pIp, float* pW )
{
    int nw = pIp[0];
    if (n > (nw << 2))
    {
        nw = n >> 2;
        makewt(nw, pIp, pW);
    }

    int nc = pIp[1];
    if(n > nc)
    {
        nc = n;
        makect(nc, pIp, pW + nw);
    }

    if(isgn < 0)
    {
        float xr = pArr[n-1];
        for(int j = n - 2; j >= 2; j -= 2)
        {
            pArr[j+1] = -pArr[j]-pArr[j-1];
            pArr[j] -= pArr[j-1];
        }

        pArr[1] = pArr[0]+xr;
        pArr[0] -= xr;

        if(n > 4)
        {
            rftbsub(n, pArr, nc, pW + nw);
            cftbsub(n, pArr, pIp, nw, pW);
        }
        else if(n == 4)
            cftbsub(n, pArr, pIp, nw, pW);
    }

    dstsub(n, pArr, nc, pW + nw);

    if (isgn >= 0)
    {
        if (n > 4)
        {
            cftfsub(n, pArr, pIp, nw, pW);
            rftfsub(n, pArr, nc, pW + nw);
        }
        else if (n == 4)
            cftfsub(n, pArr, pIp, nw, pW);

        float xr = pArr[0] - pArr[1];
        pArr[0] += pArr[1];

        for(int j = 2; j < n; j += 2)
        {
            pArr[j-1] = -pArr[j]-pArr[j+1];
            pArr[j] -= pArr[j+1];
        }

        pArr[n-1] = -xr;
    }
}

void FFT::dfct( int n, float* pArr, float* pT, int* pIp, float* pW )
{
    int nw = pIp[0];

    if (n > (nw << 3))
    {
        nw = n >> 3;
        makewt(nw, pIp, pW);
    }

    int nc = pIp[1];

    if (n > (nc << 1))
    {
        nc = n >> 1;
        makect(nc, pIp, pW + nw);
    }

    int m = n >> 1;
    float yi = pArr[m];
    float xi = pArr[0]+pArr[n];

    pArr[0] -= pArr[n];
    pT[0] = xi-yi;
    pT[m] = xi+yi;

    if(n > 2)
    {
        int mh = m >> 1;

        for(int j = 1; j < mh; ++j)
        {
            int k = m-j;
            float xr = pArr[j]-pArr[n-j];
            xi = pArr[j] + pArr[n-j];
            float yr = pArr[k]-pArr[n-k];
            yi = pArr[k]+pArr[n-k];
            pArr[j] = xr;
            pArr[k] = yr;
            pT[j] = xi-yi;
            pT[k] = xi+yi;
        }

        pT[mh] = pArr[mh]+pArr[n-mh];
        pArr[mh] -= pArr[n-mh];
        dctsub(m, pArr, nc, pW+nw);

        if(m > 4)
        {
            cftfsub(m, pArr, pIp, nw, pW);
            rftfsub(m, pArr, nc, pW+nw);
        }
        else if (m == 4)
            cftfsub(m, pArr, pIp, nw, pW);

        pArr[n-1] = pArr[0]-pArr[1];
        pArr[1] = pArr[0]+pArr[1];

        for(int j = m - 2; j >= 2; j -= 2)
        {
            pArr[2*j+1] = pArr[j]+pArr[j+1];
            pArr[2*j-1] = pArr[j]-pArr[j+1];
        }

        int l = 2;
        m = mh;

        while(m >= 2)
        {
            dctsub(m, pT, nc, pW + nw);
            if(m > 4)
            {
                cftfsub(m, pT, pIp, nw, pW);
                rftfsub(m, pT, nc, pW+nw);
            }
            else if(m == 4)
                cftfsub(m, pT, pIp, nw, pW);

            pArr[n-l] = pT[0]-pT[1];
            pArr[l] = pT[0]+pT[1];
            int k = 0;

            for(int j = 2; j < m; j += 2)
            {
                k += l << 2;
                pArr[k-l] = pT[j]-pT[j+1];
                pArr[k+l] = pT[j]+pT[j+1];
            }

            l <<= 1;
            mh = m >> 1;

            for(int j = 0; j < mh; ++j)
            {
                k = m-j;
                pT[j] = pT[m+k]-pT[m+j];
                pT[k] = pT[m+k]+pT[m+j];
            }

            pT[mh] = pT[m+mh];
            m = mh;
        }

        pArr[l] = pT[0];
        pArr[n] = pT[2]-pT[1];
        pArr[0] = pT[2]+pT[1];
    }
    else
    {
        pArr[1] = pArr[0];
        pArr[2] = pT[0];
        pArr[0] = pT[1];
    }
}

void FFT::dfst( int n, float* pArr, float* pT, int* pIp, float* pW )
{
    int nw = pIp[0];

    if(n > (nw << 3))
    {
        nw = n >> 3;
        makewt(nw, pIp, pW);
    }

    int nc = pIp[1];

    if(n > (nc << 1))
    {
        nc = n >> 1;
        makect(nc, pIp, pW+nw);
    }

    if(n > 2)
    {
        int m = n >> 1;
        int mh = m >> 1;

        for(int j = 1; j < mh; ++j)
        {
            int k = m-j;
            float xr = pArr[j]+pArr[n-j];
            float xi = pArr[j]-pArr[n-j];
            float yr = pArr[k]+pArr[n-k];
            float yi = pArr[k]-pArr[n-k];

            pArr[j] = xr;
            pArr[k] = yr;
            pT[j] = xi+yi;
            pT[k] = xi-yi;
        }

        pT[0] = pArr[mh] - pArr[n-mh];
        pArr[mh] += pArr[n-mh];
        pArr[0] = pArr[m];
        dstsub(m, pArr, nc, pW+nw);

        if(m > 4)
        {
            cftfsub(m, pArr, pIp, nw, pW);
            rftfsub(m, pArr, nc, pW+nw);
        }
        else if(m == 4)
            cftfsub(m, pArr, pIp, nw, pW);

        pArr[n-1] = pArr[1]-pArr[0];
        pArr[1] = pArr[0]+pArr[1];

        for(int j = m - 2; j >= 2; j -= 2)
        {
            pArr[2*j+1] = pArr[j]-pArr[j+1];
            pArr[2*j-1] = -pArr[j]-pArr[j+1];
        }

        int l = 2;
        m = mh;

        while (m >= 2)
        {
            dstsub(m, pT, nc, pW + nw);

            if (m > 4)
            {
                cftfsub(m, pT, pIp, nw, pW);
                rftfsub(m, pT, nc, pW+nw);
            }
            else if (m == 4)
                cftfsub(m, pT, pIp, nw, pW);

            pArr[n-l] = pT[1]-pT[0];
            pArr[l] = pT[0]+pT[1];
            int k = 0;

            for(int j = 2; j < m; j += 2)
            {
                k += l << 2;
                pArr[k-l] = -pT[j]-pT[j+1];
                pArr[k+l] = pT[j]-pT[j+1];
            }

            l <<= 1;
            mh = m >> 1;

            for(int j = 1; j < mh; ++j)
            {
                k = m - j;
                pT[j] = pT[m+k]+pT[m+j];
                pT[k] = pT[m+k]-pT[m+j];
            }

            pT[0] = pT[m+mh];
            m = mh;
        }
        pArr[l] = pT[0];
    }
    pArr[0] = 0;
}

void FFT::makewt(int nw, int* pIp, float* pW )
{
    pIp[0] = nw;
    pIp[1] = 1;

    if(nw > 2)
    {
        int nwh = nw >> 1;
        float delta = atan(1.0)/nwh;
        float wn4r = cos(delta*nwh);

        pW[0] = 1;
        pW[1] = wn4r;

        if(nwh == 4)
        {
            pW[2] = cos(delta*2);
            pW[3] = sin(delta*2);
        }
        else if(nwh > 4)
        {
            makeipt(nw, pIp);
            pW[2] = 0.5/cos(delta*2);
            pW[3] = 0.5/cos(delta*6);
            for(int j = 4; j < nwh; j += 4)
            {
                pW[j] = cos(delta*j);
                pW[j+1] = sin(delta*j);
                pW[j+2] = cos(3*delta*j);
                pW[j+3] = -sin(3*delta*j);
            }
        }

        int nw0 = 0;

        while(nwh > 2)
        {
            int nw1 = nw0+nwh;
            nwh >>= 1;
            pW[nw1] = 1;
            pW[nw1+1] = wn4r;

            float wk1r = 0.f;
            float wk1i = 0.f;

            if(nwh == 4)
            {
                wk1r = pW[nw0+4];
                wk1i = pW[nw0+5];
                pW[nw1+2] = wk1r;
                pW[nw1+3] = wk1i;
            }
            else if(nwh > 4)
            {
                wk1r = pW[nw0+4];
                float wk3r = pW[nw0+6];
                float wk3i = 0.f;

                pW[nw1+2] = 0.5/wk1r;
                pW[nw1+3] = 0.5/wk3r;

                for(int j = 4; j < nwh; j += 4)
                {
                    int idx1 = nw0+2*j;
                    int idx2 = nw1+j;

                    wk1r = pW[idx1];
                    wk1i = pW[idx1+1];
                    wk3r = pW[idx1+2];
                    wk3i = pW[idx1+3];

                    pW[idx2] = wk1r;
                    pW[idx2+1] = wk1i;
                    pW[idx2+2] = wk3r;
                    pW[idx2+3] = wk3i;
                }
            }
            nw0 = nw1;
        }
    }
}

void FFT::makeipt(int nw, int* pIp)
{
    pIp[2] = 0;
    pIp[3] = 16;
    int m = 2;

    for(int l = nw; l > 32; l >>= 2)
    {
        int m2 = m << 1;
        int q = m2 << 3;

        for(int j = m; j < m2; ++j)
        {
            int p = pIp[j] << 2;
            pIp[m+j] = p;
            pIp[m2+j] = p+q;
        }

        m = m2;
    }
}

void FFT::makect( int nc, int* pIp, float* pC )
{
    pIp[1] = nc;

    if(nc > 1)
    {
        int nch = nc >> 1;
        float delta = atan(1.0)/nch;
        pC[0] = cos(delta*nch);
        pC[nch] = 0.5*pC[0];

        for(int j = 1; j < nch; ++j)
        {
            pC[j] = 0.5 * cos(delta*j);
            pC[nc-j] = 0.5 * sin(delta*j);
        }
    }
}

void FFT::bluestein_complex(float* pArr, int isign)
{
    float* pAk = new float[2*m_nBluestein];

#ifdef USE_CDFT_THREADS
    if( m_n > CDFT_4THREADS_BEGIN_N )
    {
        int nthread = 4;

        {
            cdft_thread_t th[4];
            bluestein_complex_arg1_st ag[4];

            int k = m_n/nthread;
            for( int idx = 0; idx < nthread; ++idx )
            {
                int firstIdx = idx*k;
                int lastIdx = (idx == (nthread-1)) ? m_n : firstIdx+k;

                ag[idx].maxIdx = lastIdx-firstIdx;
                ag[idx].pAk = pAk[firstIdx];
                ag[idx].pArr = pArr[firstIdx];
                ag[idx].pBk1 = m_pBk1[firstIdx];
                ag[idx].isign = isign;

                cdft_thread_create( &th[idx], bluestein_complex1_th, &ag[idx] );
            }

            for( int idx = 0; idx < nthread; ++idx )
                cdft_thread_wait( th[idx] );
        }

        cftbsub(2*m_nBluestein, pAk, m_pIp, m_nw, m_pW);

        {
            cdft_thread_t th[4];
            bluestein_complex_arg2_st ag[4];

            int k = m_nBluestein/nthread;
            for( int idx = 0; idx < nthread; ++idx )
            {
                int firstIdx = idx*k;
                int lastIdx = (idx == (nthread-1)) ? m_nBluestein : firstIdx+k;

                ag[idx].maxIdx = lastIdx-firstIdx;
                ag[idx].pAk = pAk[firstIdx];
                ag[idx].pBk2 = m_pBk2[firstIdx];
                ag[idx].isign = isign;

                cdft_thread_create( &th[idx], bluestein_complex2_th, &ag[idx] );
            }

            for( int idx = 0; idx < nthread; ++idx )
                cdft_thread_wait( th[idx] );
        }

        cftfsub(2*m_nBluestein, pAk, m_pIp, m_nw, m_pW);

        {
            cdft_thread_t th[4];
            bluestein_complex_arg1_st ag[4];

            int k = m_n/nthread;
            for(int idx = 0; idx < nthread; ++idx)
            {
                int firstIdx = idx*k;
                int lastIdx = (idx == (nthread-1)) ? m_n : firstIdx+k;

                ag[idx].maxIdx = lastIdx-firstIdx;
                ag[idx].pAk = pAk[firstIdx];
                ag[idx].pArr = pArr[firstIdx];
                ag[idx].pBk1 = m_pBk1[firstIdx];
                ag[idx].isign = isign;

                cdft_thread_create( &th[idx], bluestein_complex3_th, &ag[idx] );
            }

            for( int idx = 0; idx < nthread; ++idx )
                cdft_thread_wait( th[idx] );
        }
    }
    else
#endif /* USE_CDFT_THREADS */
    {
        if(isign > 0)
        {
            for(int idx = 0; idx < m_n; ++idx)
            {
                int idx1 = 2*idx;
                int idx2 = idx1+1;
                int idx3 = idx1;
                int idx4 = idx2;
                pAk[idx1] = pArr[idx3]*m_pBk1[idx1]-pArr[idx4]*m_pBk1[idx2];
                pAk[idx2] = pArr[idx3]*m_pBk1[idx2]+pArr[idx4]*m_pBk1[idx1];
            }
        }
        else
        {
            for(int idx = 0; idx < m_n; ++idx)
            {
                int idx1 = 2*idx;
                int idx2 = idx1+1;
                pAk[idx1] = pArr[idx1]*m_pBk1[idx1]+pArr[idx2]*m_pBk1[idx2];
                pAk[idx2] = -pArr[idx1]*m_pBk1[idx2]+pArr[idx2]*m_pBk1[idx1];
            }
        }

        cftbsub(2*m_nBluestein, pAk, m_pIp, m_nw, m_pW);

        if(isign > 0)
        {
            for(int idx = 0; idx < m_nBluestein; ++idx)
            {
                int idx1 = 2*idx;
                int idx2 = idx1+1;
                float im = -pAk[idx1]*m_pBk2[idx2]+pAk[idx2]*m_pBk2[idx1];
                pAk[idx1] = pAk[idx1]*m_pBk2[idx1]+pAk[idx2]*m_pBk2[idx2];
                pAk[idx2] = im;
            }
        }
        else
        {
            for(int idx = 0; idx < m_nBluestein; ++idx)
            {
                int idx1 = 2*idx;
                int idx2 = idx1+1;
                float im = pAk[idx1]*m_pBk2[idx2]+pAk[idx2]*m_pBk2[idx1];
                pAk[idx1] = pAk[idx1]*m_pBk2[idx1]-pAk[idx2]*m_pBk2[idx2];
                pAk[idx2] = im;
            }
        }

        cftfsub(2*m_nBluestein, pAk, m_pIp, m_nw, m_pW);

        if(isign > 0)
        {
            for(int idx = 0; idx < m_n; ++idx)
            {
                int idx1 = 2*idx;
                int idx2 = idx1+1;
                pArr[idx1] = m_pBk1[idx1]*pAk[idx1]-m_pBk1[idx2]*pAk[idx2];
                pArr[idx2] = m_pBk1[idx2]*pAk[idx1]+m_pBk1[idx1]*pAk[idx2];
            }
        }
        else
        {
            for(int idx = 0; idx < m_n; ++idx)
            {
                int idx1 = 2*idx;
                int idx2 = idx1+1;
                pArr[idx1] = m_pBk1[idx1]*pAk[idx1]+m_pBk1[idx2]*pAk[idx2];
                pArr[idx2] = -m_pBk1[idx2]*pAk[idx1]+m_pBk1[idx1]*pAk[idx2];
            }
        }
    }

    if( pAk != nullptr )
    {
        delete[] pAk;
        pAk = nullptr;
    }
}

#ifdef USE_CDFT_THREADS
void* FFT::bluestein_complex1_th(void* p)
{
    int maxIdx = ((bluestein_complex_arg1_t*) p)->maxIdx;
    float* pAk = ((bluestein_complex_arg1_t*) p)->pAk;
    float* pArr = ((bluestein_complex_arg1_t*) p)->pArr;
    float* pBk1 = ((bluestein_complex_arg1_t*) p)->pBk1;
    int isign = ((bluestein_complex_arg1_t*) p)->isign;

    if(isign > 0)
    {
        for(int idx = 0; idx < lastIdx; ++idx)
        {
            int idx1 = 2 * idx;
            int idx2 = idx1 + 1;
            pAk[idx1] = pArr[idx1] * pBk1[idx1] - pArr[idx2] * pBk1[idx2];
            pAk[idx2] = pArr[idx1] * pBk1[idx2] + pArr[idx2] * pBk1[idx1];
        }
    }
    else
    {
        for(int idx = 0; idx < lastIdx; ++idx)
        {
            int idx1 = 2 * idx;
            int idx2 = idx1 + 1;
            pAk[idx1] = pArr[idx1] * pBk1[idx1] + pArr[idx2] * pBk1[idx2];
            pAk[idx2] = -pArr[idx1] * pBk1[idx2] + pArr[idx2] * pBk1[idx1];
        }
    }

    return (void*)0;
}

void* FFT::bluestein_complex2_th(void* p)
{
    int maxIdx = ((bluestein_complex_arg2_t*) p)->maxIdx;
    float* pAk = ((bluestein_complex_arg2_t*) p)->pAk;
    float* pBk2 = ((bluestein_complex_arg2_t*) p)->pBk2;
    int isign = ((bluestein_complex_arg2_t*) p)->isign;

    if(isign > 0)
    {
        for(int idx = 0; idx < lastIdx; ++idx)
        {
            int idx1 = 2*idx;
            int idx2 = idx1+1;
            float im = -pAk[idx1]*pBk2[idx2]+pAk[idx2]*pBk2[idx1];
            pAk[idx1] = pAk[idx1]*pBk2[idx1]+pAk[idx2]*pBk2[idx2];
            pAk[idx2] = im;
        }
    }
    else
    {
        for(int idx = 0; idx < lastIdx; ++idx)
        {
            int idx1 = 2*idx;
            int idx2 = idx1+1;
            float im = pAk[idx1]*pBk2[idx2]+pAk[idx2]*pBk2[idx1];
            pAk[idx1] = pAk[idx1]*pBk2[idx1]-pAk[idx2]*pBk2[idx2];
            pAk[idx2] = im;
        }
    }

    return (void*)0;
}

void* FFT::bluestein_complex3_th(void* p)
{
    int maxIdx = ((bluestein_complex_arg1_t*) p)->maxIdx;
    float* pAk = ((bluestein_complex_arg1_t*) p)->pAk;
    float* pArr = ((bluestein_complex_arg1_t*) p)->pArr;
    float* pBk1 = ((bluestein_complex_arg1_t*) p)->pBk1;
    int isign = ((bluestein_complex_arg1_t*) p)->isign;

    if(isign > 0)
    {
        for(int idx = 0; idx < lastIdx; ++idx)
        {
            int idx1 = 2*idx;
            int idx2 = idx1+1;
            pArr[idx1] = pBk1[idx1]*pAk[idx1]-pBk1[idx2]*pAk[idx2];
            pArr[idx2] = pBk1[idx2]*pAk[idx1]+pBk1[idx1]*pAk[idx2];
        }
    }
    else
    {
        for(int idx = 0; idx < lastIdx; ++idx)
        {
            int idx1 = 2*idx;
            int idx2 = idx1+1;
            pArr[idx1] = pBk1[idx1]*pAk[idx1]+pBk1[idx2]*pAk[idx2];
            pArr[idx2] = -pBk1[idx2]*pAk[idx1]+pBk1[idx1]*pAk[idx2];
        }
    }

    return (void*)0;
}
#endif

void FFT::cftfsub( int n, float* pArr, int* pIp, int nw, float* pW )
{
    if(n > 8)
    {
        if(n > 32)
        {
            cftf1st(n, pArr, &pW[nw-(n>>2)]);
#ifdef USE_CDFT_THREADS
            if(n > CDFT_THREADS_BEGIN_N)
                cftrec4_th(n, pArr, nw, pW);
            else
#endif /* USE_CDFT_THREADS */
            if( n > 512 )
                cftrec4( n, pArr, nw, pW );
            else if( n > 128 )
                cftleaf( n, 1, pArr, nw, pW );
            else
                cftfx41( n, pArr, nw, pW );
            bitrv2(n, pIp, pArr);
        }
        else if( n == 32 )
        {
            cftf161(pArr, &pW[nw-8]);
            bitrv216(pArr);
        }
        else
        {
            cftf081(pArr, pW);
            bitrv208(pArr);
        }
    }
    else if(n == 8)
        cftf040(pArr);
    else if(n == 4)
        cftx020(pArr);
}

void FFT::cftbsub( int n, float* pArr, int*pIp, int nw, float* pW )
{
    if( n > 8 )
    {
        if( n > 32 )
        {
            cftb1st(n, pArr, &pW[nw-(n>>2)]);
#ifdef USE_CDFT_THREADS
            if( n > CDFT_THREADS_BEGIN_N )
                cftrec4_th(n, a, nw, w);
            else
#endif /* USE_CDFT_THREADS */
            if( n > 512 )
                cftrec4( n, pArr, nw, pW );
            else if( n > 128 )
                cftleaf( n, 1, pArr, nw, pW );
            else
                cftfx41(n, pArr, nw, pW);

            bitrv2conj(n, pIp, pArr);
        }
        else if( n == 32 )
        {
            cftf161(pArr, &pW[nw-8]);
            bitrv216neg(pArr);
        }
        else
        {
            cftf081(pArr, pW);
            bitrv208neg(pArr);
        }
    }
    else if( n == 8 )
        cftb040(pArr);
    else if( n == 4 )
        cftx020(pArr);
}

void FFT::bitrv2( int n, int* pIp, float* pArr )
{
    int m = 1;
    int l = n >> 2;
    for( ; l > 8; l >>= 2)
        m <<= 1;

    int nh = n >> 1;
    int nm = 4 * m;
    int j1 = 0, k1 = 0;
    float xr = 0.f, xi = 0.f, yr = 0.f, yi = 0.f;

    if(l == 8)
    {
        for(int k = 0; k < m; ++k)
        {
            for(int j = 0; j < k; ++j)
            {
                j1 = 4*j+2*pIp[m+k];
                k1 = 4*k+2*pIp[m+j];
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += 2*nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 -= nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += 2 * nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nh;
                k1 += 2;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= 2 * nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 += nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= 2 * nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += 2;
                k1 += nh;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += 2 * nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 -= nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += 2 * nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nh;
                k1 -= 2;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= 2*nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 += nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= 2 * nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
            }
            k1 = 4*k+2*pIp[m+k];
            j1 = k1+2;
            k1 += nh;
            xr = pArr[j1];
            xi = pArr[j1+1];
            yr = pArr[k1];
            yi = pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 += nm;
            k1 += 2 * nm;
            xr = pArr[j1];
            xi = pArr[j1+1];
            yr = pArr[k1];
            yi = pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 += nm;
            k1 -= nm;
            xr = pArr[j1];
            xi = pArr[j1+1];
            yr = pArr[k1];
            yi = pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 -= 2;
            k1 -= nh;
            xr = pArr[j1];
            xi = pArr[j1+1];
            yr = pArr[k1];
            yi = pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 += nh + 2;
            k1 += nh + 2;
            xr = pArr[j1];
            xi = pArr[j1+1];
            yr = pArr[k1];
            yi = pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 -= nh-nm;
            k1 += 2*nm-2;
            xr = pArr[j1];
            xi = pArr[j1+1];
            yr = pArr[k1];
            yi = pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
        }
    }
    else
    {
        for(int k = 0; k < m; ++k)
        {
            for(int j = 0; j < k; ++j)
            {
                j1 = 4*j+pIp[m+k];
                k1 = 4*k+pIp[m+j];
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nh;
                k1 += 2;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += 2;
                k1 += nh;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nh;
                k1 -= 2;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= nm;
                xr = pArr[j1];
                xi = pArr[j1+1];
                yr = pArr[k1];
                yi = pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
            }

            k1 = 4*k+pIp[m+k];
            j1 = k1+2;
            k1 += nh;
            xr = pArr[j1];
            xi = pArr[j1+1];
            yr = pArr[k1];
            yi = pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 += nm;
            k1 += nm;
            xr = pArr[j1];
            xi = pArr[j1+1];
            yr = pArr[k1];
            yi = pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
        }
    }
}

void FFT::bitrv2conj( int n, int* pIp, float* pArr )
{
    int m = 1;

    int l = n >> 2;
    for(; l > 8; l >>= 2)
        m <<= 1;

    int nh = n >> 1;
    int nm = 4*m;
    int j1 = 0, k1 = 0;
    float xr = 0.f, xi = 0.f, yr = 0.f, yi = 0.f;

    if(l == 8)
    {
        for(int k = 0; k < m; ++k)
        {
            for(int j = 0; j < k; ++j)
            {
                j1 = 4*j+2*pIp[m+k];
                k1 = 4*k+2*pIp[m+j];
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += 2 * nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 -= nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += 2 * nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nh;
                k1 += 2;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= 2 * nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 += nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= 2 * nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += 2;
                k1 += nh;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += 2 * nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 -= nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += 2 * nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1 + 1] = yi;
                pArr[k1] = xr;
                pArr[k1 + 1] = xi;
                j1 -= nh;
                k1 -= 2;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= 2 * nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 += nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= 2 * nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
            }

            k1 = 4*k+2*pIp[m+k];
            j1 = k1+2;
            k1 += nh;
            pArr[j1-1] = -pArr[j1-1];
            xr = pArr[j1];
            xi = -pArr[j1+1];
            yr = pArr[k1];
            yi = -pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            pArr[k1+3] = -pArr[k1+3];
            j1 += nm;
            k1 += 2 * nm;
            xr = pArr[j1];
            xi = -pArr[j1+1];
            yr = pArr[k1];
            yi = -pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 += nm;
            k1 -= nm;
            xr = pArr[j1];
            xi = -pArr[j1+1];
            yr = pArr[k1];
            yi = -pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 -= 2;
            k1 -= nh;
            xr = pArr[j1];
            xi = -pArr[j1+1];
            yr = pArr[k1];
            yi = -pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 += nh + 2;
            k1 += nh + 2;
            xr = pArr[j1];
            xi = -pArr[j1 + 1];
            yr = pArr[k1];
            yi = -pArr[k1 + 1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            j1 -= nh - nm;
            k1 += 2 * nm - 2;
            pArr[j1-1] = -pArr[j1-1];
            xr = pArr[j1];
            xi = -pArr[j1+1];
            yr = pArr[k1];
            yi = -pArr[k1+1];
            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            pArr[k1+3] = -pArr[k1+3];
        }
    }
    else
    {
        for(int k = 0; k < m; ++k)
        {
            for(int j = 0; j < k; ++j)
            {
                j1 = 4*j+pIp[m+k];
                k1 = 4*k+pIp[m+j];
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nh;
                k1 += 2;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += 2;
                k1 += nh;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 += nm;
                k1 += nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nh;
                k1 -= 2;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
                j1 -= nm;
                k1 -= nm;
                xr = pArr[j1];
                xi = -pArr[j1+1];
                yr = pArr[k1];
                yi = -pArr[k1+1];
                pArr[j1] = yr;
                pArr[j1+1] = yi;
                pArr[k1] = xr;
                pArr[k1+1] = xi;
            }

            k1 = 4*k+pIp[m+k];

            j1 = k1+2;
            k1 += nh;

            pArr[j1-1] = -pArr[j1-1];

            xr = pArr[j1];
            xi = -pArr[j1+1];
            yr = pArr[k1];
            yi = -pArr[k1+1];

            pArr[j1] = yr;
            pArr[j1 + 1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            pArr[k1+3] = -pArr[k1+3];

            j1 += nm;
            k1 += nm;

            pArr[j1-1] = -pArr[j1-1];

            xr = pArr[j1];
            xi = -pArr[j1+1];
            yr = pArr[k1];
            yi = -pArr[k1+1];

            pArr[j1] = yr;
            pArr[j1+1] = yi;
            pArr[k1] = xr;
            pArr[k1+1] = xi;
            pArr[k1+3] = -pArr[k1+3];
        }
    }
}

void FFT::bitrv216(float* pArr)
{
    float x1r = pArr[2];
    float x1i = pArr[3];
    float x2r = pArr[4];
    float x2i = pArr[5];
    float x3r = pArr[6];
    float x3i = pArr[7];
    float x4r = pArr[8];
    float x4i = pArr[9];
    float x5r = pArr[10];
    float x5i = pArr[11];
    float x7r = pArr[14];
    float x7i = pArr[15];
    float x8r = pArr[16];
    float x8i = pArr[17];
    float x10r = pArr[20];
    float x10i = pArr[21];
    float x11r = pArr[22];
    float x11i = pArr[23];
    float x12r = pArr[24];
    float x12i = pArr[25];
    float x13r = pArr[26];
    float x13i = pArr[27];
    float x14r = pArr[28];
    float x14i = pArr[29];

    pArr[2] = x8r;
    pArr[3] = x8i;
    pArr[4] = x4r;
    pArr[5] = x4i;
    pArr[6] = x12r;
    pArr[7] = x12i;
    pArr[8] = x2r;
    pArr[9] = x2i;
    pArr[10] = x10r;
    pArr[11] = x10i;
    pArr[14] = x14r;
    pArr[15] = x14i;
    pArr[16] = x1r;
    pArr[17] = x1i;
    pArr[20] = x5r;
    pArr[21] = x5i;
    pArr[22] = x13r;
    pArr[23] = x13i;
    pArr[24] = x3r;
    pArr[25] = x3i;
    pArr[26] = x11r;
    pArr[27] = x11i;
    pArr[28] = x7r;
    pArr[29] = x7i;
}

void FFT::bitrv216neg( float*pArr )
{
    float x1r = pArr[2];
    float x1i = pArr[3];
    float x2r = pArr[4];
    float x2i = pArr[5];
    float x3r = pArr[6];
    float x3i = pArr[7];
    float x4r = pArr[8];
    float x4i = pArr[9];
    float x5r = pArr[10];
    float x5i = pArr[11];
    float x6r = pArr[12];
    float x6i = pArr[13];
    float x7r = pArr[14];
    float x7i = pArr[15];
    float x8r = pArr[16];
    float x8i = pArr[17];
    float x9r = pArr[18];
    float x9i = pArr[19];
    float x10r = pArr[20];
    float x10i = pArr[21];
    float x11r = pArr[22];
    float x11i = pArr[23];
    float x12r = pArr[24];
    float x12i = pArr[25];
    float x13r = pArr[26];
    float x13i = pArr[27];
    float x14r = pArr[28];
    float x14i = pArr[29];
    float x15r = pArr[30];
    float x15i = pArr[31];

    pArr[2] = x15r;
    pArr[3] = x15i;
    pArr[4] = x7r;
    pArr[5] = x7i;
    pArr[6] = x11r;
    pArr[7] = x11i;
    pArr[8] = x3r;
    pArr[9] = x3i;
    pArr[10] = x13r;
    pArr[11] = x13i;
    pArr[12] = x5r;
    pArr[13] = x5i;
    pArr[14] = x9r;
    pArr[15] = x9i;
    pArr[16] = x1r;
    pArr[17] = x1i;
    pArr[18] = x14r;
    pArr[19] = x14i;
    pArr[20] = x6r;
    pArr[21] = x6i;
    pArr[22] = x10r;
    pArr[23] = x10i;
    pArr[24] = x2r;
    pArr[25] = x2i;
    pArr[26] = x12r;
    pArr[27] = x12i;
    pArr[28] = x4r;
    pArr[29] = x4i;
    pArr[30] = x8r;
    pArr[31] = x8i;
}

void FFT::bitrv208(float* pArr)
{
    float x1r = pArr[2];
    float x1i = pArr[3];
    float x3r = pArr[6];
    float x3i = pArr[7];
    float x4r = pArr[8];
    float x4i = pArr[9];
    float x6r = pArr[12];
    float x6i = pArr[13];

    pArr[2] = x4r;
    pArr[3] = x4i;
    pArr[6] = x6r;
    pArr[7] = x6i;
    pArr[8] = x1r;
    pArr[9] = x1i;
    pArr[12] = x3r;
    pArr[13] = x3i;
}

void FFT::bitrv208neg(float* pArr)
{
    float x1r = pArr[2];
    float x1i = pArr[3];
    float x2r = pArr[4];
    float x2i = pArr[5];
    float x3r = pArr[6];
    float x3i = pArr[7];
    float x4r = pArr[8];
    float x4i = pArr[9];
    float x5r = pArr[10];
    float x5i = pArr[11];
    float x6r = pArr[12];
    float x6i = pArr[13];
    float x7r = pArr[14];
    float x7i = pArr[15];

    pArr[2] = x7r;
    pArr[3] = x7i;
    pArr[4] = x3r;
    pArr[5] = x3i;
    pArr[6] = x5r;
    pArr[7] = x5i;
    pArr[8] = x1r;
    pArr[9] = x1i;
    pArr[10] = x6r;
    pArr[11] = x6i;
    pArr[12] = x2r;
    pArr[13] = x2i;
    pArr[14] = x4r;
    pArr[15] = x4i;
}

void FFT::cftf1st(int n, float* pArr, float* pW)
{
    int mh = n >> 3;
    int m = 2 * mh;
    int j0 = 0;
    int j1 = m;
    int j2 = j1 + m;
    int j3 = j2 + m;

    float x0r = pArr[0]+pArr[j2];
    float x0i = pArr[1]+pArr[j2+1];
    float x1r = pArr[0]-pArr[j2];
    float x1i = pArr[1]-pArr[j2+1];
    float x2r = pArr[j1]+pArr[j3];
    float x2i = pArr[j1+1]+pArr[j3+1];
    float x3r = pArr[j1]-pArr[j3];
    float x3i = pArr[j1+1]-pArr[j3+1];

    pArr[0] = x0r+x2r;
    pArr[1] = x0i+x2i;
    pArr[j1] = x0r-x2r;
    pArr[j1+1] = x0i-x2i;
    pArr[j2] = x1r-x3i;
    pArr[j2+1] = x1i+x3r;
    pArr[j3] = x1r+x3i;
    pArr[j3+1] = x1i-x3r;

    float wn4r = pW[1];
    float csc1 = pW[2];
    float csc3 = pW[3];
    float wd1r = 1;
    float wd1i = 0;
    float wd3r = 1;
    float wd3i = 0;

    float wk1r = 0.f;
    float wk1i = 0.f;
    float wk3r = 0.f;
    float wk3i = 0.f;

    int k = 0;
    for(int j = 2; j < mh-2; j += 4)
    {
        k += 4;

        wk1r = csc1*(wd1r+pW[k]);
        wk1i = csc1*(wd1i+pW[k+1]);
        wk3r = csc3*(wd3r+pW[k+2]);
        wk3i = csc3*(wd3i+pW[k+3]);

        wd1r = pW[k];
        wd1i = pW[k+1];
        wd3r = pW[k+2];
        wd3i = pW[k+3];

        j1 = j + m;
        j2 = j1 + m;
        j3 = j2 + m;

        x0r = pArr[j]+pArr[j2];
        x0i = pArr[j+1]+pArr[j2+1];
        x1r = pArr[j]-pArr[j2];
        x1i = pArr[j+1]-pArr[j2+1];
        x2r = pArr[j1]+pArr[j3];
        x2i = pArr[j1+1]+pArr[j3+1];
        x3r = pArr[j1]-pArr[j3];
        x3i = pArr[j1+1]-pArr[j3+1];

        float y0r = pArr[j+2]+pArr[j2+2];
        float y0i = pArr[j+3]+pArr[j2+3];
        float y1r = pArr[j+2]-pArr[j2+2];
        float y1i = pArr[j+3]-pArr[j2+3];
        float y2r = pArr[j1+2]+pArr[j3+2];
        float y2i = pArr[j1+3]+pArr[j3+3];
        float y3r = pArr[j1+2]-pArr[j3+2];
        float y3i = pArr[j1+3]-pArr[j3+3];

        pArr[j] = x0r+x2r;
        pArr[j+1] = x0i+x2i;
        pArr[j+2] = y0r+y2r;
        pArr[j+3] = y0i+y2i;
        pArr[j1] = x0r-x2r;
        pArr[j1+1] = x0i-x2i;
        pArr[j1+2] = y0r-y2r;
        pArr[j1+3] = y0i-y2i;

        x0r = x1r - x3i;
        x0i = x1i + x3r;

        pArr[j2] = wk1r*x0r-wk1i*x0i;
        pArr[j2+1] = wk1r*x0i+wk1i*x0r;

        x0r = y1r-y3i;
        x0i = y1i+y3r;

        pArr[j2+2] = wd1r*x0r-wd1i*x0i;
        pArr[j2+3] = wd1r*x0i+wd1i*x0r;

        x0r = x1r+x3i;
        x0i = x1i-x3r;

        pArr[j3] = wk3r*x0r+wk3i*x0i;
        pArr[j3+1] = wk3r*x0i-wk3i*x0r;

        x0r = y1r+y3i;
        x0i = y1i-y3r;

        pArr[j3+2] = wd3r*x0r+wd3i*x0i;
        pArr[j3+3] = wd3r*x0i-wd3i*x0r;

        j0 = m-j;
        j1 = j0+m;
        j2 = j1+m;
        j3 = j2+m;

        x0r = pArr[j0]+pArr[j2];
        x0i = pArr[j0+1]+pArr[j2+1];
        x1r = pArr[j0]-pArr[j2];
        x1i = pArr[j0+1]-pArr[j2+1];
        y0r = pArr[j0-2]+pArr[j2-2];
        y0i = pArr[j0-1]+pArr[j2-1];
        y1r = pArr[j0-2]-pArr[j2-2];
        y1i = pArr[j0-1]-pArr[j2-1];
        x2r = pArr[j1]+pArr[j3];
        x2i = pArr[j1+1]+pArr[j3+1];
        x3r = pArr[j1]-pArr[j3];
        x3i = pArr[j1+1]-pArr[j3+1];
        y2r = pArr[j1-2]+pArr[j3-2];
        y2i = pArr[j1-1]+pArr[j3-1];
        y3r = pArr[j1-2]-pArr[j3-2];
        y3i = pArr[j1-1]-pArr[j3-1];

        pArr[j0] = x0r+x2r;
        pArr[j0+1] = x0i+x2i;
        pArr[j0-2] = y0r+y2r;
        pArr[j0-1] = y0i+y2i;
        pArr[j1] = x0r-x2r;
        pArr[j1+1] = x0i-x2i;
        pArr[j1-2] = y0r-y2r;
        pArr[j1-1] = y0i-y2i;
        x0r = x1r-x3i;
        x0i = x1i+x3r;
        pArr[j2] = wk1i*x0r-wk1r*x0i;
        pArr[j2+1] = wk1i*x0i+wk1r*x0r;
        x0r = y1r - y3i;
        x0i = y1i + y3r;
        pArr[j2-2] = wd1i*x0r-wd1r*x0i;
        pArr[j2-1] = wd1i*x0i+wd1r*x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        pArr[j3] = wk3i*x0r+wk3r*x0i;
        pArr[j3+1] = wk3i*x0i-wk3r*x0r;
        x0r = y1r+y3i;
        x0i = y1i-y3r;
        pArr[j3-2] = wd3i*x0r+wd3r*x0i;
        pArr[j3-1] = wd3i*x0i-wd3r*x0r;
    }

    wk1r = csc1*(wd1r+wn4r);
    wk1i = csc1*(wd1i+wn4r);
    wk3r = csc3*(wd3r-wn4r);
    wk3i = csc3*(wd3i-wn4r);

    j0 = mh;
    j1 = j0+m;
    j2 = j1+m;
    j3 = j2+m;

    x0r = pArr[j0-2]+pArr[j2-2];
    x0i = pArr[j0-1]+pArr[j2-1];
    x1r = pArr[j0-2]-pArr[j2-2];
    x1i = pArr[j0-1]-pArr[j2-1];
    x2r = pArr[j1-2]+pArr[j3-2];
    x2i = pArr[j1-1]+pArr[j3-1];
    x3r = pArr[j1-2]-pArr[j3-2];
    x3i = pArr[j1-1]-pArr[j3-1];

    pArr[j0-2] = x0r+x2r;
    pArr[j0-1] = x0i+x2i;
    pArr[j1-2] = x0r-x2r;
    pArr[j1-1] = x0i-x2i;

    x0r = x1r-x3i;
    x0i = x1i+x3r;

    pArr[j2-2] = wk1r*x0r-wk1i*x0i;
    pArr[j2-1] = wk1r*x0i+wk1i*x0r;

    x0r = x1r+x3i;
    x0i = x1i-x3r;

    pArr[j3-2] = wk3r*x0r+wk3i*x0i;
    pArr[j3-1] = wk3r*x0i-wk3i*x0r;

    x0r = pArr[j0]+pArr[j2];
    x0i = pArr[j0+1]+pArr[j2+1];
    x1r = pArr[j0]-pArr[j2];
    x1i = pArr[j0+1]-pArr[j2+1];
    x2r = pArr[j1]+pArr[j3];
    x2i = pArr[j1+1]+pArr[j3+1];
    x3r = pArr[j1]-pArr[j3];
    x3i = pArr[j1+1]-pArr[j3+1];

    pArr[j0] = x0r+x2r;
    pArr[j0+1] = x0i+x2i;
    pArr[j1] = x0r-x2r;
    pArr[j1+1] = x0i-x2i;

    x0r = x1r-x3i;
    x0i = x1i+x3r;

    pArr[j2] = wn4r*(x0r-x0i);
    pArr[j2+1] = wn4r*(x0i+x0r);

    x0r = x1r+x3i;
    x0i = x1i-x3r;

    pArr[j3] = -wn4r*(x0r+x0i);
    pArr[j3+1] = -wn4r*(x0i-x0r);

    x0r = pArr[j0+2]+pArr[j2+2];
    x0i = pArr[j0+3]+pArr[j2+3];
    x1r = pArr[j0+2]-pArr[j2+2];
    x1i = pArr[j0+3]-pArr[j2+3];
    x2r = pArr[j1+2]+pArr[j3+2];
    x2i = pArr[j1+3]+pArr[j3+3];
    x3r = pArr[j1+2]-pArr[j3+2];
    x3i = pArr[j1+3]-pArr[j3+3];

    pArr[j0+2] = x0r+x2r;
    pArr[j0+3] = x0i+x2i;
    pArr[j1+2] = x0r-x2r;
    pArr[j1+3] = x0i-x2i;

    x0r = x1r-x3i;
    x0i = x1i+x3r;

    pArr[j2+2] = wk1i*x0r-wk1r*x0i;
    pArr[j2+3] = wk1i*x0i+wk1r*x0r;

    x0r = x1r+x3i;
    x0i = x1i-x3r;

    pArr[j3+2] = wk3i*x0r+wk3r*x0i;
    pArr[j3+3] = wk3i*x0i-wk3r*x0r;
}

void FFT::cftb1st(int n, float* pArr, float* pW)
{
    int mh = n >> 3;
    int m = 2 * mh;

    int j0 = 0;
    int j1 = m;
    int j2 = j1+m;
    int j3 = j2+m;

    float x0r = pArr[0]+pArr[j2];
    float x0i = -pArr[1]-pArr[j2+1];
    float x1r = pArr[0]-pArr[j2];
    float x1i = -pArr[1]+pArr[j2+1];
    float x2r = pArr[j1]+pArr[j3];
    float x2i = pArr[j1+1]+pArr[j3+1];
    float x3r = pArr[j1]-pArr[j3];
    float x3i = pArr[j1+1]-pArr[j3+1];

    pArr[0] = x0r+x2r;
    pArr[1] = x0i-x2i;
    pArr[j1] = x0r-x2r;
    pArr[j1+1] = x0i+x2i;
    pArr[j2] = x1r+x3i;
    pArr[j2+1] = x1i+x3r;
    pArr[j3] = x1r-x3i;
    pArr[j3+1] = x1i-x3r;

    float wn4r = pW[1];
    float csc1 = pW[2];
    float csc3 = pW[3];
    float wd1r = 1;
    float wd1i = 0;
    float wd3r = 1;
    float wd3i = 0;
    float wk1r = 0.f;
    float wk1i = 0.f;
    float wk3r = 0.f;
    float wk3i = 0.f;

    int k = 0;
    for(int j = 2; j < mh - 2; j += 4)
    {
        k += 4;
        wk1r = csc1*(wd1r+pW[k]);
        wk1i = csc1*(wd1i+pW[k+1]);
        wk3r = csc3*(wd3r+pW[k+2]);
        wk3i = csc3*(wd3i+pW[k+3]);
        wd1r = pW[k];
        wd1i = pW[k+1];
        wd3r = pW[k+2];
        wd3i = pW[k+3];

        j1 = j+m;
        j2 = j1+m;
        j3 = j2+m;

        x0r = pArr[j]+pArr[j2];
        x0i = -pArr[j+1]-pArr[j2+1];
        x1r = pArr[j]-pArr[j2];
        x1i = -pArr[j+1]+pArr[j2+1];

        float y0r = pArr[j+2]+pArr[j2+2];
        float y0i = -pArr[j+3]-pArr[j2+3];
        float y1r = pArr[j+2]-pArr[j2+2];
        float y1i = -pArr[j+3]+pArr[j2+3];

        x2r = pArr[j1]+pArr[j3];
        x2i = pArr[j1+1]+pArr[j3+1];
        x3r = pArr[j1]- pArr[j3];
        x3i = pArr[j1+1]-pArr[j3+1];

        float y2r = pArr[j1+2]+pArr[j3+2];
        float y2i = pArr[j1+3]+pArr[j3+3];
        float y3r = pArr[j1+2]-pArr[j3+2];
        float y3i = pArr[j1+3]-pArr[j3+3];

        pArr[j] = x0r+x2r;
        pArr[j+1] = x0i-x2i;
        pArr[j+2] = y0r+y2r;
        pArr[j+3] = y0i-y2i;
        pArr[j1] = x0r-x2r;
        pArr[j1+1] = x0i+x2i;
        pArr[j1+2] = y0r-y2r;
        pArr[j1+3] = y0i+y2i;
        x0r = x1r+x3i;
        x0i = x1i+x3r;
        pArr[j2] = wk1r*x0r-wk1i*x0i;
        pArr[j2+1] = wk1r*x0i+wk1i*x0r;
        x0r = y1r+y3i;
        x0i = y1i+y3r;
        pArr[j2+2] = wd1r*x0r-wd1i*x0i;
        pArr[j2+3] = wd1r*x0i+wd1i*x0r;
        x0r = x1r-x3i;
        x0i = x1i-x3r;
        pArr[j3] = wk3r*x0r+wk3i*x0i;
        pArr[j3+1] = wk3r*x0i-wk3i*x0r;
        x0r = y1r-y3i;
        x0i = y1i-y3r;
        pArr[j3+2] = wd3r*x0r+wd3i*x0i;
        pArr[j3+3] = wd3r*x0i-wd3i*x0r;
        j0 = m-j;
        j1 = j0+m;
        j2 = j1+m;
        j3 = j2+m;
        x0r = pArr[j0]+pArr[j2];
        x0i = -pArr[j0+1]-pArr[j2+1];
        x1r = pArr[j0]-pArr[j2];
        x1i = -pArr[j0+1]+pArr[j2+1];
        y0r = pArr[j0-2]+pArr[j2-2];
        y0i = -pArr[j0-1]-pArr[j2-1];
        y1r = pArr[j0-2]-pArr[j2-2];
        y1i = -pArr[j0-1]+pArr[j2-1];
        x2r = pArr[j1]+pArr[j3];
        x2i = pArr[j1+1]+pArr[j3+1];
        x3r = pArr[j1]-pArr[j3];
        x3i = pArr[j1+1]-pArr[j3+1];
        y2r = pArr[j1-2]+pArr[j3-2];
        y2i = pArr[j1-1]+pArr[j3-1];
        y3r = pArr[j1-2]-pArr[j3-2];
        y3i = pArr[j1-1]-pArr[j3-1];

        pArr[j0] = x0r+x2r;
        pArr[j0+1] = x0i-x2i;
        pArr[j0-2] = y0r+y2r;
        pArr[j0-1] = y0i-y2i;
        pArr[j1] = x0r-x2r;
        pArr[j1+1] = x0i+x2i;
        pArr[j1-2] = y0r-y2r;
        pArr[j1-1] = y0i+y2i;

        x0r = x1r+x3i;
        x0i = x1i+x3r;
        pArr[j2] = wk1i*x0r-wk1r*x0i;
        pArr[j2+1] = wk1i*x0i+wk1r*x0r;
        x0r = y1r+y3i;
        x0i = y1i+y3r;
        pArr[j2-2] = wd1i*x0r-wd1r*x0i;
        pArr[j2-1] = wd1i*x0i+wd1r*x0r;
        x0r = x1r-x3i;
        x0i = x1i-x3r;
        pArr[j3] = wk3i*x0r+wk3r*x0i;
        pArr[j3+1] = wk3i*x0i-wk3r*x0r;
        x0r = y1r-y3i;
        x0i = y1i-y3r;
        pArr[j3-2] = wd3i*x0r+wd3r*x0i;
        pArr[j3-1] = wd3i*x0i-wd3r*x0r;
    }

    wk1r = csc1*(wd1r+wn4r);
    wk1i = csc1*(wd1i+wn4r);
    wk3r = csc3*(wd3r-wn4r);
    wk3i = csc3*(wd3i-wn4r);

    j0 = mh;
    j1 = j0+m;
    j2 = j1+m;
    j3 = j2+m;

    x0r = pArr[j0-2]+pArr[j2-2];
    x0i = -pArr[j0-1]-pArr[j2-1];
    x1r = pArr[j0-2]-pArr[j2-2];
    x1i = -pArr[j0-1]+pArr[j2-1];
    x2r = pArr[j1-2]+pArr[j3-2];
    x2i = pArr[j1-1]+pArr[j3-1];
    x3r = pArr[j1-2]-pArr[j3-2];
    x3i = pArr[j1-1]-pArr[j3-1];
    pArr[j0-2] = x0r+x2r;
    pArr[j0-1] = x0i-x2i;
    pArr[j1-2] = x0r-x2r;
    pArr[j1-1] = x0i+x2i;
    x0r = x1r + x3i;
    x0i = x1i + x3r;
    pArr[j2-2] = wk1r*x0r-wk1i*x0i;
    pArr[j2-1] = wk1r*x0i+wk1i*x0r;
    x0r = x1r-x3i;
    x0i = x1i-x3r;
    pArr[j3-2] = wk3r*x0r+wk3i*x0i;
    pArr[j3-1] = wk3r*x0i-wk3i*x0r;
    x0r = pArr[j0]+pArr[j2];
    x0i = -pArr[j0+1]-pArr[j2+1];
    x1r = pArr[j0]-pArr[j2];
    x1i = -pArr[j0+1]+pArr[j2+1];
    x2r = pArr[j1]+pArr[j3];
    x2i = pArr[j1+1]+pArr[j3+1];
    x3r = pArr[j1]-pArr[j3];
    x3i = pArr[j1+1]-pArr[j3+1];
    pArr[j0] = x0r+x2r;
    pArr[j0+1] = x0i-x2i;
    pArr[j1] = x0r-x2r;
    pArr[j1+1] = x0i+x2i;
    x0r = x1r+x3i;
    x0i = x1i+x3r;
    pArr[j2] = wn4r*(x0r-x0i);
    pArr[j2+1] = wn4r*(x0i+x0r);
    x0r = x1r-x3i;
    x0i = x1i-x3r;
    pArr[j3] = -wn4r*(x0r+x0i);
    pArr[j3+1] = -wn4r*(x0i-x0r);
    x0r = pArr[j0+2]+pArr[j2+2];
    x0i = -pArr[j0+3]-pArr[j2+3];
    x1r = pArr[j0+2]-pArr[j2+2];
    x1i = -pArr[j0+3]+pArr[j2+3];
    x2r = pArr[j1+2]+pArr[j3+2];
    x2i = pArr[j1+3]+pArr[j3+3];
    x3r = pArr[j1+2]-pArr[j3+2];
    x3i = pArr[j1+3]-pArr[j3+3];
    pArr[j0+2] = x0r+x2r;
    pArr[j0+3] = x0i-x2i;
    pArr[j1+2] = x0r-x2r;
    pArr[j1+3] = x0i+x2i;
    x0r = x1r+x3i;
    x0i = x1i+x3r;
    pArr[j2+2] = wk1i*x0r-wk1r*x0i;
    pArr[j2+3] = wk1i*x0i+wk1r*x0r;
    x0r = x1r-x3i;
    x0i = x1i-x3r;
    pArr[j3+2] = wk3i*x0r+wk3r*x0i;
    pArr[j3+3] = wk3i*x0i-wk3r*x0r;
}

#ifdef USE_CDFT_THREADS
void FFT::cftrec4_th(int n, float* pArr, int nw, float* pW)
{
    int idiv4, m, nthread;
    cdft_thread_t th[4];
    cdft_arg_t ag[4];

    nthread = 2;
    idiv4 = 0;
    m = n >> 1;

    if(n > CDFT_4THREADS_BEGIN_N)
    {
        nthread = 4;
        idiv4 = 1;
        m >>= 1;
    }

    for(int idx = 0; idx < nthread; ++idx)
    {
        ag[idx].n0 = n;
        ag[idx].n = m;
        ag[idx].pArr = &pArr[idx*m];
        ag[idx].nw = nw;
        ag[idx].pW = pW;

        if(idx != idiv4)
            cdft_thread_create(&th[idx], cftrec1_th, &ag[idx]);
        else
            cdft_thread_create(&th[idx], cftrec2_th, &ag[idx]);
    }

    for(int idx = 0; idx < nthread; ++idx)
        cdft_thread_wait(th[idx]);
}

void* FFT::cftrec1_th(void* p)
{
    int n0 = ((cdft_arg_t *) p)->n0;
    int n = ((cdft_arg_t *) p)->n;
    float* pArr = ((cdft_arg_t *) p)->pArr;
    int nw = ((cdft_arg_t *) p)->nw;
    float* pW = ((cdft_arg_t *) p)->pW;

    int m = n0;
    while (m > 512)
    {
        m >>= 2;
        cftmdl1(m, &pArr[n-m], &pW[nw-(m>>1)]);
    }
    cftleaf(m, 1, &pArr[n-m], nw, pW);

    int k = 0;
    for(int j = n - m; j > 0; j -= m)
    {
        k++;
        int isplt = cfttree(m, j, k, pArr, nw, pW);
        cftleaf(m, isplt, &pArr[j-m], nw, pW);
    }
    return (void *) 0;
}

void* FFT::cftrec2_th(void* p)
{
    int n0 = ((cdft_arg_t *) p)->n0;
    int n = ((cdft_arg_t *) p)->n;
    float* pArr = ((cdft_arg_t *) p)->pArr;
    int nw = ((cdft_arg_t *) p)->nw;
    float* pW = ((cdft_arg_t *) p)->pW;
    int k = 1;
    int m = n0;

    while(m > 512)
    {
        m >>= 2;
        k <<= 2;
        cftmdl2(m, &pArr[n-m], &pW[nw-m]);
    }

    cftleaf(m, 0, &pArr[n-m], nw, pW);

    k >>= 1;

    for(int j = n - m; j > 0; j -= m)
    {
        ++k;
        int isplt = cfttree(m, j, k, pArr, nw, pW);
        cftleaf(m, isplt, &pArr[j-m], nw, pW);
    }

    return (void *)0;
}
#endif /* USE_CDFT_THREADS */

void FFT::cftrec4( int n, float* pArr, int nw, float* pW )
{
    int m = n;

    while(m > 512)
    {
        m >>= 2;
        cftmdl1(m, &pArr[n-m], &pW[nw-(m>>1)]);
    }

    cftleaf(m, 1, &pArr[n-m], nw, pW);

    int k = 0;
    for(int j = n - m; j > 0; j -= m)
    {
        k++;
        int isplt = cfttree(m, j, k, pArr, nw, pW);
        cftleaf(m, isplt, &pArr[j-m], nw, pW);
    }
}

int FFT::cfttree(int n, int j, int k, float* pArr, int nw, float* pW)
{
    int isplt = 0;

    if((k & 3) != 0)
    {
        isplt = k&1;
        if (isplt != 0)
            cftmdl1(n, &pArr[j-n], &pW[nw-(n>>1)]);
        else
            cftmdl2(n, &pArr[j-n], &pW[nw-n]);
    }
    else
    {
        int m = n;
        int idx = k;
        for( ; (idx&3) == 0; idx>>=2)
            m <<= 2;

        isplt = idx&1;
        if(isplt != 0)
        {
            while(m > 128)
            {
                cftmdl1(m, &pArr[j-m], &pW[nw-(m>>1)]);
                m >>= 2;
            }
        }
        else
        {
            while(m > 128)
            {
                cftmdl2(m, &pArr[j-m], &pW[nw-m]);
                m >>= 2;
            }
        }
    }
    return isplt;
}

void FFT::cftleaf( int n, int isplt, float* pArr, int nw, float* pW )
{
    if(n == 512)
    {
        cftmdl1(128, pArr, &pW[nw-64]);
        cftf161(pArr, &pW[nw - 8]);
        cftf162(&pArr[32], &pW[nw-32]);
        cftf161(&pArr[64], &pW[nw-8]);
        cftf161(&pArr[96], &pW[nw-8]);
        cftmdl2(128, &pArr[128], &pW[nw-128]);
        cftf161(&pArr[128], &pW[nw-8]);
        cftf162(&pArr[160], &pW[nw-32]);
        cftf161(&pArr[192], &pW[nw-8]);
        cftf162(&pArr[224], &pW[nw-32]);
        cftmdl1(128, &pArr[256], &pW[nw-64]);
        cftf161(&pArr[256], &pW[nw-8]);
        cftf162(&pArr[288], &pW[nw-32]);
        cftf161(&pArr[320], &pW[nw-8]);
        cftf161(&pArr[352], &pW[nw-8]);

        if (isplt != 0)
        {
            cftmdl1(128, &pArr[384], &pW[nw-64]);
            cftf161(&pArr[480], &pW[nw-8]);
        }
        else
        {
            cftmdl2(128, &pArr[384], &pW[nw-128]);
            cftf162(&pArr[480], &pW[nw-32]);
        }

        cftf161(&pArr[384], &pW[nw-8]);
        cftf162(&pArr[416], &pW[nw-32]);
        cftf161(&pArr[448], &pW[nw-8]);
    }
    else
    {
        cftmdl1(64, pArr, &pW[nw-32]);
        cftf081(pArr, &pW[nw-8]);
        cftf082(&pArr[16], &pW[nw-8]);
        cftf081(&pArr[32], &pW[nw-8]);
        cftf081(&pArr[48], &pW[nw-8]);
        cftmdl2(64, &pArr[64], &pW[nw-64]);
        cftf081(&pArr[64], &pW[nw-8]);
        cftf082(&pArr[80], &pW[nw-8]);
        cftf081(&pArr[96], &pW[nw-8]);
        cftf082(&pArr[112], &pW[nw-8]);
        cftmdl1(64, &pArr[128], &pW[nw-32]);
        cftf081(&pArr[128], &pW[nw-8]);
        cftf082(&pArr[144], &pW[nw-8]);
        cftf081(&pArr[160], &pW[nw-8]);
        cftf081(&pArr[176], &pW[nw-8]);

        if (isplt != 0)
        {
            cftmdl1(64, &pArr[192], &pW[nw-32]);
            cftf081(&pArr[240], &pW[nw-8]);
        }
        else
        {
            cftmdl2(64, &pArr[192], &pW[nw-64]);
            cftf082(&pArr[240], &pW[nw-8]);
        }

        cftf081(&pArr[192], &pW[nw-8]);
        cftf082(&pArr[208], &pW[nw-8]);
        cftf081(&pArr[224], &pW[nw-8]);
    }
}

void FFT::cftmdl1(int n, float* pArr, float* pW)
{
    float wn4r, wk1r, wk1i, wk3r, wk3i;
    float x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    int mh = n >> 3;
    int m = 2 * mh;

    int j0 = 0;
    int j1 = m;
    int j2 = j1+m;
    int j3 = j2+m;

    x0r = pArr[0]+pArr[j2];
    x0i = pArr[1]+pArr[j2+1];
    x1r = pArr[0]-pArr[j2];
    x1i = pArr[1]-pArr[j2+1];
    x2r = pArr[j1]+pArr[j3];
    x2i = pArr[j1+1]+pArr[j3+1];
    x3r = pArr[j1]-pArr[j3];
    x3i = pArr[j1+1]-pArr[j3+1];
    pArr[0] = x0r+x2r;
    pArr[1] = x0i+x2i;
    pArr[j1] = x0r-x2r;
    pArr[j1+1] = x0i-x2i;
    pArr[j2] = x1r-x3i;
    pArr[j2+1] = x1i+x3r;
    pArr[j3] = x1r+x3i;
    pArr[j3+1] = x1i-x3r;
    wn4r = pW[1];

    int k = 0;

    for(int j = 2; j < mh; j += 2)
    {
        k += 4;
        wk1r = pW[k];
        wk1i = pW[k+1];
        wk3r = pW[k+2];
        wk3i = pW[k+3];
        j1 = j+m;
        j2 = j1+m;
        j3 = j2+m;

        x0r = pArr[j]+pArr[j2];
        x0i = pArr[j+1]+pArr[j2+1];
        x1r = pArr[j]-pArr[j2];
        x1i = pArr[j+1]-pArr[j2+1];
        x2r = pArr[j1]+pArr[j3];
        x2i = pArr[j1+1]+pArr[j3+1];
        x3r = pArr[j1]-pArr[j3];
        x3i = pArr[j1+1]-pArr[j3+1];
        pArr[j] = x0r+x2r;
        pArr[j+1] = x0i+x2i;
        pArr[j1] = x0r-x2r;
        pArr[j1+1] = x0i-x2i;
        x0r = x1r-x3i;
        x0i = x1i+x3r;
        pArr[j2] = wk1r*x0r-wk1i*x0i;
        pArr[j2+1] = wk1r*x0i+wk1i*x0r;
        x0r = x1r+x3i;
        x0i = x1i-x3r;
        pArr[j3] = wk3r*x0r+wk3i*x0i;
        pArr[j3+1] = wk3r*x0i-wk3i*x0r;

        j0 = m- j;
        j1 = j0+m;
        j2 = j1+m;
        j3 = j2+m;

        x0r = pArr[j0]+pArr[j2];
        x0i = pArr[j0+1]+pArr[j2+1];
        x1r = pArr[j0]-pArr[j2];
        x1i = pArr[j0+1]-pArr[j2+1];
        x2r = pArr[j1]+pArr[j3];
        x2i = pArr[j1+1]+pArr[j3+1];
        x3r = pArr[j1]-pArr[j3];
        x3i = pArr[j1+1]-pArr[j3+1];
        pArr[j0] = x0r+x2r;
        pArr[j0+1] = x0i+x2i;
        pArr[j1] = x0r-x2r;
        pArr[j1+1] = x0i-x2i;
        x0r = x1r-x3i;
        x0i = x1i+x3r;
        pArr[j2] = wk1i*x0r-wk1r*x0i;
        pArr[j2+1] = wk1i*x0i+wk1r*x0r;
        x0r = x1r+x3i;
        x0i = x1i-x3r;
        pArr[j3] = wk3i*x0r+wk3r*x0i;
        pArr[j3+1] = wk3i*x0i-wk3r*x0r;
    }

    j0 = mh;
    j1 = j0+m;
    j2 = j1+m;
    j3 = j2+m;

    x0r = pArr[j0]+pArr[j2];
    x0i = pArr[j0+1]+pArr[j2+1];
    x1r = pArr[j0]-pArr[j2];
    x1i = pArr[j0+1]-pArr[j2+1];
    x2r = pArr[j1]+pArr[j3];
    x2i = pArr[j1+1]+pArr[j3+1];
    x3r = pArr[j1]-pArr[j3];
    x3i = pArr[j1+1]-pArr[j3+1];
    pArr[j0] = x0r+x2r;
    pArr[j0+1] = x0i+x2i;
    pArr[j1] = x0r-x2r;
    pArr[j1+1] = x0i-x2i;
    x0r = x1r-x3i;
    x0i = x1i+x3r;
    pArr[j2] = wn4r*(x0r-x0i);
    pArr[j2+1] = wn4r*(x0i+x0r);
    x0r = x1r+x3i;
    x0i = x1i-x3r;
    pArr[j3] = -wn4r*(x0r+x0i);
    pArr[j3+1] = -wn4r*(x0i-x0r);
}

void FFT::cftmdl2(int n, float* pArr, float* pW)
{
    int mh = n >> 3;
    int m = 2 * mh;
    float wn4r = pW[1];

    int j0 = 0;
    int j1 = m;
    int j2 = j1+m;
    int j3 = j2+m;

    float x0r = pArr[0]-pArr[j2 + 1];
    float x0i = pArr[1]+pArr[j2];
    float x1r = pArr[0]+pArr[j2 + 1];
    float x1i = pArr[1]-pArr[j2];
    float x2r = pArr[j1]-pArr[j3 + 1];
    float x2i = pArr[j1+1]+pArr[j3];
    float x3r = pArr[j1]+pArr[j3 + 1];
    float x3i = pArr[j1+1]-pArr[j3];
    float y0r = wn4r*(x2r-x2i);
    float y0i = wn4r*(x2i+x2r);
    float y2r = 0.f;
    float y2i = 0.f;

    float wk1r = 0.f, wk1i = 0.f, wk3r = 0.f, wk3i = 0.f;

    pArr[0] = x0r+y0r;
    pArr[1] = x0i+y0i;
    pArr[j1] = x0r-y0r;
    pArr[j1+1] = x0i-y0i;
    y0r = wn4r * (x3r-x3i);
    y0i = wn4r * (x3i+x3r);
    pArr[j2] = x1r-y0i;
    pArr[j2+1] = x1i+y0r;
    pArr[j3] = x1r+y0i;
    pArr[j3+1] = x1i-y0r;

    int k = 0;
    int kr = 2 * m;
    for(int j = 2; j < mh; j += 2)
    {
        k += 4;
        wk1r = pW[k];
        wk1i = pW[k+1];
        wk3r = pW[k+2];
        wk3i = pW[k+3];
        kr -= 4;

        float wd1i = pW[kr];
        float wd1r = pW[kr+1];
        float wd3i = pW[kr+2];
        float wd3r = pW[kr+3];

        j1 = j+m;
        j2 = j1+m;
        j3 = j2+m;

        x0r = pArr[j]-pArr[j2+1];
        x0i = pArr[j+1]+pArr[j2];
        x1r = pArr[j]+pArr[j2+1];
        x1i = pArr[j+1]-pArr[j2];
        x2r = pArr[j1]-pArr[j3+1];
        x2i = pArr[j1+1]+pArr[j3];
        x3r = pArr[j1]+pArr[j3+1];
        x3i = pArr[j1+1]-pArr[j3];
        y0r = wk1r*x0r-wk1i*x0i;
        y0i = wk1r*x0i+wk1i*x0r;
        y2r = wd1r*x2r-wd1i*x2i;
        y2i = wd1r*x2i+wd1i*x2r;
        pArr[j] = y0r+y2r;
        pArr[j+1] = y0i+y2i;
        pArr[j1] = y0r-y2r;
        pArr[j1+1] = y0i-y2i;
        y0r = wk3r*x1r+wk3i*x1i;
        y0i = wk3r*x1i-wk3i*x1r;
        y2r = wd3r*x3r+wd3i*x3i;
        y2i = wd3r*x3i-wd3i*x3r;
        pArr[j2] = y0r+y2r;
        pArr[j2+1] = y0i+y2i;
        pArr[j3] = y0r-y2r;
        pArr[j3+1] = y0i-y2i;

        j0 = m-j;
        j1 = j0+m;
        j2 = j1+m;
        j3 = j2+m;

        x0r = pArr[j0]-pArr[j2+1];
        x0i = pArr[j0+1]+pArr[j2];
        x1r = pArr[j0]+pArr[j2+1];
        x1i = pArr[j0+1]-pArr[j2];
        x2r = pArr[j1]-pArr[j3+1];
        x2i = pArr[j1+1]+pArr[j3];
        x3r = pArr[j1]+pArr[j3+1];
        x3i = pArr[j1+1]-pArr[j3];
        y0r = wd1i*x0r-wd1r*x0i;
        y0i = wd1i*x0i+wd1r*x0r;
        y2r = wk1i*x2r-wk1r*x2i;
        y2i = wk1i*x2i+wk1r*x2r;

        pArr[j0] = y0r+y2r;
        pArr[j0+1] = y0i+y2i;
        pArr[j1] = y0r-y2r;
        pArr[j1+1] = y0i-y2i;
        y0r = wd3i*x1r+wd3r*x1i;
        y0i = wd3i*x1i-wd3r*x1r;
        y2r = wk3i*x3r+wk3r*x3i;
        y2i = wk3i*x3i-wk3r*x3r;
        pArr[j2] = y0r+y2r;
        pArr[j2+1] = y0i+y2i;
        pArr[j3] = y0r-y2r;
        pArr[j3+1] = y0i-y2i;
    }

    wk1r = pW[m];
    wk1i = pW[m+1];

    j0 = mh;
    j1 = j0+m;
    j2 = j1+m;
    j3 = j2+m;

    x0r = pArr[j0]-pArr[j2+1];
    x0i = pArr[j0+1]+pArr[j2];
    x1r = pArr[j0]+pArr[j2+1];
    x1i = pArr[j0+1]-pArr[j2];
    x2r = pArr[j1]-pArr[j3+1];
    x2i = pArr[j1+1]+pArr[j3];
    x3r = pArr[j1]+pArr[j3+1];
    x3i = pArr[j1+1]-pArr[j3];
    y0r = wk1r*x0r-wk1i*x0i;
    y0i = wk1r*x0i+wk1i*x0r;
    y2r = wk1i*x2r-wk1r*x2i;
    y2i = wk1i*x2i+wk1r*x2r;
    pArr[j0] = y0r+y2r;
    pArr[j0+1] = y0i+y2i;
    pArr[j1] = y0r-y2r;
    pArr[j1+1] = y0i-y2i;
    y0r = wk1i*x1r-wk1r*x1i;
    y0i = wk1i*x1i+wk1r*x1r;
    y2r = wk1r*x3r-wk1i*x3i;
    y2i = wk1r*x3i+wk1i*x3r;
    pArr[j2] = y0r-y2r;
    pArr[j2+1] = y0i-y2i;
    pArr[j3] = y0r+y2r;
    pArr[j3+1] = y0i+y2i;
}

void FFT::cftfx41(int n, float* pArr, int nw, float* pW)
{
    if(n == 128)
    {
        cftf161(pArr, &pW[nw-8]);
        cftf162(&pArr[32], &pW[nw-32]);
        cftf161(&pArr[64], &pW[nw-8]);
        cftf161(&pArr[96], &pW[nw-8]);
    }
    else
    {
        cftf081(pArr, &pW[nw-8]);
        cftf082(&pArr[16], &pW[nw-8]);
        cftf081(&pArr[32], &pW[nw-8]);
        cftf081(&pArr[48], &pW[nw-8]);
    }
}


void FFT::cftf161(float* pArr, float* pW)
{
    float wn4r = pW[1];
    float wk1r = pW[2];
    float wk1i = pW[3];

    float x0r = pArr[0]+pArr[16];
    float x0i = pArr[1]+pArr[17];
    float x1r = pArr[0]-pArr[16];
    float x1i = pArr[1]-pArr[17];
    float x2r = pArr[8]+pArr[24];
    float x2i = pArr[9]+pArr[25];
    float x3r = pArr[8]-pArr[24];
    float x3i = pArr[9]-pArr[25];

    float y0r = x0r+x2r;
    float y0i = x0i+x2i;
    float y4r = x0r-x2r;
    float y4i = x0i-x2i;
    float y8r = x1r-x3i;
    float y8i = x1i+x3r;
    float y12r = x1r+x3i;
    float y12i = x1i-x3r;

    x0r = pArr[2]+pArr[18];
    x0i = pArr[3]+pArr[19];
    x1r = pArr[2]-pArr[18];
    x1i = pArr[3]-pArr[19];
    x2r = pArr[10]+pArr[26];
    x2i = pArr[11]+pArr[27];
    x3r = pArr[10]-pArr[26];
    x3i = pArr[11]-pArr[27];

    float y1r = x0r+x2r;
    float y1i = x0i+x2i;
    float y5r = x0r-x2r;
    float y5i = x0i-x2i;

    x0r = x1r-x3i;
    x0i = x1i+x3r;

    float y9r = wk1r*x0r-wk1i*x0i;
    float y9i = wk1r*x0i+wk1i*x0r;

    x0r = x1r+x3i;
    x0i = x1i-x3r;

    float y13r = wk1i*x0r-wk1r*x0i;
    float y13i = wk1i*x0i+wk1r*x0r;

    x0r = pArr[4]+pArr[20];
    x0i = pArr[5]+pArr[21];
    x1r = pArr[4]-pArr[20];
    x1i = pArr[5]-pArr[21];
    x2r = pArr[12]+pArr[28];
    x2i = pArr[13]+pArr[29];
    x3r = pArr[12]-pArr[28];
    x3i = pArr[13]-pArr[29];

    float y2r = x0r+x2r;
    float y2i = x0i+x2i;
    float y6r = x0r-x2r;
    float y6i = x0i-x2i;

    x0r = x1r-x3i;
    x0i = x1i+x3r;

    float y10r = wn4r*(x0r-x0i);
    float y10i = wn4r*(x0i+x0r);

    x0r = x1r+x3i;
    x0i = x1i-x3r;

    float y14r = wn4r*(x0r+x0i);
    float y14i = wn4r*(x0i-x0r);

    x0r = pArr[6]+pArr[22];
    x0i = pArr[7]+pArr[23];
    x1r = pArr[6]-pArr[22];
    x1i = pArr[7]-pArr[23];
    x2r = pArr[14]+pArr[30];
    x2i = pArr[15]+pArr[31];
    x3r = pArr[14]-pArr[30];
    x3i = pArr[15]-pArr[31];

    float y3r = x0r+x2r;
    float y3i = x0i+x2i;
    float y7r = x0r-x2r;
    float y7i = x0i-x2i;

    x0r = x1r-x3i;
    x0i = x1i+x3r;

    float y11r = wk1i*x0r-wk1r*x0i;
    float y11i = wk1i*x0i+wk1r*x0r;
    x0r = x1r+x3i;
    x0i = x1i-x3r;

    float y15r = wk1r*x0r-wk1i*x0i;
    float y15i = wk1r*x0i+wk1i*x0r;

    x0r = y12r-y14r;
    x0i = y12i-y14i;
    x1r = y12r+y14r;
    x1i = y12i+y14i;
    x2r = y13r-y15r;
    x2i = y13i-y15i;
    x3r = y13r+y15r;
    x3i = y13i+y15i;
    pArr[24] = x0r+x2r;
    pArr[25] = x0i+x2i;
    pArr[26] = x0r-x2r;
    pArr[27] = x0i-x2i;
    pArr[28] = x1r-x3i;
    pArr[29] = x1i+x3r;
    pArr[30] = x1r+x3i;
    pArr[31] = x1i-x3r;
    x0r = y8r+y10r;
    x0i = y8i+y10i;
    x1r = y8r-y10r;
    x1i = y8i-y10i;
    x2r = y9r+y11r;
    x2i = y9i+y11i;
    x3r = y9r-y11r;
    x3i = y9i-y11i;
    pArr[16] = x0r+x2r;
    pArr[17] = x0i+x2i;
    pArr[18] = x0r-x2r;
    pArr[19] = x0i-x2i;
    pArr[20] = x1r-x3i;
    pArr[21] = x1i+x3r;
    pArr[22] = x1r+x3i;
    pArr[23] = x1i-x3r;
    x0r = y5r-y7i;
    x0i = y5i+y7r;
    x2r = wn4r*(x0r-x0i);
    x2i = wn4r*(x0i+x0r);
    x0r = y5r+y7i;
    x0i = y5i-y7r;
    x3r = wn4r*(x0r-x0i);
    x3i = wn4r*(x0i+x0r);
    x0r = y4r-y6i;
    x0i = y4i+y6r;
    x1r = y4r+y6i;
    x1i = y4i-y6r;
    pArr[8] = x0r+x2r;
    pArr[9] = x0i+x2i;
    pArr[10] = x0r-x2r;
    pArr[11] = x0i-x2i;
    pArr[12] = x1r-x3i;
    pArr[13] = x1i+x3r;
    pArr[14] = x1r+x3i;
    pArr[15] = x1i-x3r;
    x0r = y0r+y2r;
    x0i = y0i+y2i;
    x1r = y0r-y2r;
    x1i = y0i-y2i;
    x2r = y1r+y3r;
    x2i = y1i+y3i;
    x3r = y1r-y3r;
    x3i = y1i-y3i;
    pArr[0] = x0r+x2r;
    pArr[1] = x0i+x2i;
    pArr[2] = x0r-x2r;
    pArr[3] = x0i-x2i;
    pArr[4] = x1r-x3i;
    pArr[5] = x1i+x3r;
    pArr[6] = x1r+x3i;
    pArr[7] = x1i-x3r;
}


void FFT::cftf162(float* pArr, float* pW)
{
    float wn4r = pW[1];
    float wk1r = pW[4];
    float wk1i = pW[5];
    float wk3r = pW[6];
    float wk3i = -pW[7];
    float wk2r = pW[8];
    float wk2i = pW[9];
    float x1r = pArr[0]-pArr[17];
    float x1i = pArr[1]+pArr[16];
    float x0r = pArr[8]-pArr[25];
    float x0i = pArr[9]+pArr[24];
    float x2r = wn4r*(x0r-x0i);
    float x2i = wn4r*(x0i+x0r);
    float y0r = x1r+x2r;
    float y0i = x1i+x2i;
    float y4r = x1r-x2r;
    float y4i = x1i-x2i;

    x1r = pArr[0]+pArr[17];
    x1i = pArr[1]-pArr[16];
    x0r = pArr[8]+pArr[25];
    x0i = pArr[9]-pArr[24];
    x2r = wn4r*(x0r-x0i);
    x2i = wn4r*(x0i+x0r);

    float y8r = x1r-x2i;
    float y8i = x1i+x2r;
    float y12r = x1r+x2i;
    float y12i = x1i-x2r;

    x0r = pArr[2]-pArr[19];
    x0i = pArr[3]+pArr[18];
    x1r = wk1r*x0r-wk1i*x0i;
    x1i = wk1r*x0i+wk1i*x0r;
    x0r = pArr[10]-pArr[27];
    x0i = pArr[11]+pArr[26];
    x2r = wk3i*x0r-wk3r*x0i;
    x2i = wk3i*x0i+wk3r*x0r;

    float y1r = x1r+x2r;
    float y1i = x1i+x2i;
    float y5r = x1r-x2r;
    float y5i = x1i-x2i;

    x0r = pArr[2]+pArr[19];
    x0i = pArr[3]-pArr[18];
    x1r = wk3r*x0r-wk3i*x0i;
    x1i = wk3r*x0i+wk3i*x0r;
    x0r = pArr[10]+pArr[27];
    x0i = pArr[11]-pArr[26];
    x2r = wk1r*x0r+wk1i*x0i;
    x2i = wk1r*x0i-wk1i*x0r;

    float y9r = x1r-x2r;
    float y9i = x1i-x2i;
    float y13r = x1r+x2r;
    float y13i = x1i+x2i;

    x0r = pArr[4]-pArr[21];
    x0i = pArr[5]+pArr[20];
    x1r = wk2r*x0r-wk2i*x0i;
    x1i = wk2r*x0i+wk2i*x0r;
    x0r = pArr[12]-pArr[29];
    x0i = pArr[13]+pArr[28];
    x2r = wk2i*x0r-wk2r*x0i;
    x2i = wk2i*x0i+wk2r*x0r;
    float y2r = x1r+x2r;
    float y2i = x1i+x2i;
    float y6r = x1r-x2r;
    float y6i = x1i-x2i;

    x0r = pArr[4]+pArr[21];
    x0i = pArr[5]-pArr[20];
    x1r = wk2i*x0r-wk2r*x0i;
    x1i = wk2i*x0i+wk2r*x0r;
    x0r = pArr[12]+pArr[29];
    x0i = pArr[13]-pArr[28];
    x2r = wk2r*x0r-wk2i*x0i;
    x2i = wk2r*x0i+wk2i*x0r;

    float y10r = x1r-x2r;
    float y10i = x1i-x2i;
    float y14r = x1r+x2r;
    float y14i = x1i+x2i;

    x0r = pArr[6]-pArr[23];
    x0i = pArr[7]+pArr[22];
    x1r = wk3r*x0r-wk3i*x0i;
    x1i = wk3r*x0i+wk3i*x0r;
    x0r = pArr[14]-pArr[31];
    x0i = pArr[15]+pArr[30];
    x2r = wk1i*x0r-wk1r*x0i;
    x2i = wk1i*x0i+wk1r*x0r;

    float y3r = x1r+x2r;
    float y3i = x1i+x2i;
    float y7r = x1r-x2r;
    float y7i = x1i-x2i;

    x0r = pArr[6]+pArr[23];
    x0i = pArr[7]-pArr[22];
    x1r = wk1i*x0r+wk1r*x0i;
    x1i = wk1i*x0i-wk1r*x0r;
    x0r = pArr[14]+pArr[31];
    x0i = pArr[15]-pArr[30];
    x2r = wk3i*x0r-wk3r*x0i;
    x2i = wk3i*x0i+wk3r*x0r;

    float y11r = x1r+x2r;
    float y11i = x1i+x2i;
    float y15r = x1r-x2r;
    float y15i = x1i-x2i;

    x1r = y0r+y2r;
    x1i = y0i+y2i;
    x2r = y1r+y3r;
    x2i = y1i+y3i;
    pArr[0] = x1r+x2r;
    pArr[1] = x1i+x2i;
    pArr[2] = x1r-x2r;
    pArr[3] = x1i-x2i;
    x1r = y0r-y2r;
    x1i = y0i-y2i;
    x2r = y1r-y3r;
    x2i = y1i-y3i;
    pArr[4] = x1r-x2i;
    pArr[5] = x1i+x2r;
    pArr[6] = x1r+x2i;
    pArr[7] = x1i-x2r;
    x1r = y4r-y6i;
    x1i = y4i+y6r;
    x0r = y5r-y7i;
    x0i = y5i+y7r;
    x2r = wn4r*(x0r-x0i);
    x2i = wn4r*(x0i+x0r);
    pArr[8] = x1r+x2r;
    pArr[9] = x1i+x2i;
    pArr[10] = x1r-x2r;
    pArr[11] = x1i-x2i;
    x1r = y4r+y6i;
    x1i = y4i-y6r;
    x0r = y5r+y7i;
    x0i = y5i-y7r;
    x2r = wn4r*(x0r-x0i);
    x2i = wn4r*(x0i+x0r);
    pArr[12] = x1r-x2i;
    pArr[13] = x1i+x2r;
    pArr[14] = x1r+x2i;
    pArr[15] = x1i-x2r;
    x1r = y8r+y10r;
    x1i = y8i+y10i;
    x2r = y9r-y11r;
    x2i = y9i-y11i;
    pArr[16] = x1r+x2r;
    pArr[17] = x1i+x2i;
    pArr[18] = x1r-x2r;
    pArr[19] = x1i-x2i;
    x1r = y8r-y10r;
    x1i = y8i-y10i;
    x2r = y9r+y11r;
    x2i = y9i+y11i;
    pArr[20] = x1r-x2i;
    pArr[21] = x1i+x2r;
    pArr[22] = x1r+x2i;
    pArr[23] = x1i-x2r;
    x1r = y12r-y14i;
    x1i = y12i+y14r;
    x0r = y13r+y15i;
    x0i = y13i-y15r;
    x2r = wn4r*(x0r-x0i);
    x2i = wn4r*(x0i+x0r);
    pArr[24] = x1r+x2r;
    pArr[25] = x1i+x2i;
    pArr[26] = x1r-x2r;
    pArr[27] = x1i-x2i;
    x1r = y12r+y14i;
    x1i = y12i-y14r;
    x0r = y13r-y15i;
    x0i = y13i+y15r;
    x2r = wn4r*(x0r-x0i);
    x2i = wn4r*(x0i+x0r);
    pArr[28] = x1r-x2i;
    pArr[29] = x1i+x2r;
    pArr[30] = x1r+x2i;
    pArr[31] = x1i-x2r;
}


void FFT::cftf081(float* pArr, float* pW)
{
    float wn4r = pW[1];
    float x0r = pArr[0]+pArr[8];
    float x0i = pArr[1]+pArr[9];
    float x1r = pArr[0]-pArr[8];
    float x1i = pArr[1]-pArr[9];
    float x2r = pArr[4]+pArr[12];
    float x2i = pArr[5]+pArr[13];
    float x3r = pArr[4]-pArr[12];
    float x3i = pArr[5]-pArr[13];
    float y0r = x0r+x2r;
    float y0i = x0i+x2i;
    float y2r = x0r-x2r;
    float y2i = x0i-x2i;
    float y1r = x1r-x3i;
    float y1i = x1i+x3r;
    float y3r = x1r+x3i;
    float y3i = x1i-x3r;

    x0r = pArr[2]+pArr[10];
    x0i = pArr[3]+pArr[11];
    x1r = pArr[2]-pArr[10];
    x1i = pArr[3]-pArr[11];
    x2r = pArr[6]+pArr[14];
    x2i = pArr[7]+pArr[15];
    x3r = pArr[6]-pArr[14];
    x3i = pArr[7]-pArr[15];

    float y4r = x0r+x2r;
    float y4i = x0i+x2i;
    float y6r = x0r-x2r;
    float y6i = x0i-x2i;

    x0r = x1r-x3i;
    x0i = x1i+x3r;
    x2r = x1r+x3i;
    x2i = x1i-x3r;

    float y5r = wn4r*(x0r-x0i);
    float y5i = wn4r*(x0r+x0i);
    float y7r = wn4r*(x2r-x2i);
    float y7i = wn4r*(x2r+x2i);

    pArr[8] = y1r+y5r;
    pArr[9] = y1i+y5i;
    pArr[10] = y1r-y5r;
    pArr[11] = y1i-y5i;
    pArr[12] = y3r-y7i;
    pArr[13] = y3i+y7r;
    pArr[14] = y3r+y7i;
    pArr[15] = y3i-y7r;
    pArr[0] = y0r+y4r;
    pArr[1] = y0i+y4i;
    pArr[2] = y0r-y4r;
    pArr[3] = y0i-y4i;
    pArr[4] = y2r-y6i;
    pArr[5] = y2i+y6r;
    pArr[6] = y2r+y6i;
    pArr[7] = y2i-y6r;
}

void FFT::cftf082(float* pArr, float* pW)
{
    float wn4r = pW[1];
    float wk1r = pW[2];
    float wk1i = pW[3];
    float y0r = pArr[0]-pArr[9];
    float y0i = pArr[1]+pArr[8];
    float y1r = pArr[0]+pArr[9];
    float y1i = pArr[1]-pArr[8];
    float x0r = pArr[4]-pArr[13];
    float x0i = pArr[5]+pArr[12];
    float y2r = wn4r*(x0r-x0i);
    float y2i = wn4r*(x0i+x0r);

    x0r = pArr[4]+pArr[13];
    x0i = pArr[5]-pArr[12];

    float y3r = wn4r*(x0r-x0i);
    float y3i = wn4r*(x0i+x0r);

    x0r = pArr[2]-pArr[11];
    x0i = pArr[3]+pArr[10];

    float y4r = wk1r*x0r-wk1i*x0i;
    float y4i = wk1r*x0i+wk1i*x0r;

    x0r = pArr[2]+pArr[11];
    x0i = pArr[3]-pArr[10];

    float y5r = wk1i*x0r-wk1r*x0i;
    float y5i = wk1i*x0i+wk1r*x0r;

    x0r = pArr[6]-pArr[15];
    x0i = pArr[7]+pArr[14];

    float y6r = wk1i*x0r-wk1r*x0i;
    float y6i = wk1i*x0i+wk1r*x0r;

    x0r = pArr[6]+pArr[15];
    x0i = pArr[7]-pArr[14];

    float y7r = wk1r*x0r-wk1i*x0i;
    float y7i = wk1r*x0i+wk1i*x0r;

    x0r = y0r+y2r;
    x0i = y0i+y2i;

    float x1r = y4r+y6r;
    float x1i = y4i+y6i;

    pArr[0] = x0r+x1r;
    pArr[1] = x0i+x1i;
    pArr[2] = x0r-x1r;
    pArr[3] = x0i-x1i;
    x0r = y0r-y2r;
    x0i = y0i-y2i;
    x1r = y4r-y6r;
    x1i = y4i-y6i;
    pArr[4] = x0r-x1i;
    pArr[5] = x0i+x1r;
    pArr[6] = x0r+x1i;
    pArr[7] = x0i-x1r;
    x0r = y1r-y3i;
    x0i = y1i+y3r;
    x1r = y5r-y7r;
    x1i = y5i-y7i;
    pArr[8] = x0r+x1r;
    pArr[9] = x0i+x1i;
    pArr[10] = x0r-x1r;
    pArr[11] = x0i-x1i;
    x0r = y1r+y3i;
    x0i = y1i-y3r;
    x1r = y5r+y7r;
    x1i = y5i+y7i;
    pArr[12] = x0r-x1i;
    pArr[13] = x0i+x1r;
    pArr[14] = x0r+x1i;
    pArr[15] = x0i-x1r;
}

void FFT::cftf040(float* pArr)
{
    float x0r = pArr[0]+pArr[4];
    float x0i = pArr[1]+pArr[5];
    float x1r = pArr[0]-pArr[4];
    float x1i = pArr[1]-pArr[5];
    float x2r = pArr[2]+pArr[6];
    float x2i = pArr[3]+pArr[7];
    float x3r = pArr[2]-pArr[6];
    float x3i = pArr[3]-pArr[7];

    pArr[0] = x0r+x2r;
    pArr[1] = x0i+x2i;
    pArr[2] = x1r-x3i;
    pArr[3] = x1i+x3r;
    pArr[4] = x0r-x2r;
    pArr[5] = x0i-x2i;
    pArr[6] = x1r+x3i;
    pArr[7] = x1i-x3r;
}

void FFT::cftb040(float* pArr)
{
    float x0r = pArr[0]+pArr[4];
    float x0i = pArr[1]+pArr[5];
    float x1r = pArr[0]-pArr[4];
    float x1i = pArr[1]-pArr[5];
    float x2r = pArr[2]+pArr[6];
    float x2i = pArr[3]+pArr[7];
    float x3r = pArr[2]-pArr[6];
    float x3i = pArr[3]-pArr[7];

    pArr[0] = x0r+x2r;
    pArr[1] = x0i+x2i;
    pArr[2] = x1r+x3i;
    pArr[3] = x1i-x3r;
    pArr[4] = x0r-x2r;
    pArr[5] = x0i-x2i;
    pArr[6] = x1r-x3i;
    pArr[7] = x1i+x3r;
}

void FFT::cftx020(float* pArr)
{
    float x0r = pArr[0]-pArr[2];
    float x0i = pArr[1]-pArr[3];
    pArr[0] += pArr[2];
    pArr[1] += pArr[3];
    pArr[2] = x0r;
    pArr[3] = x0i;
}


void FFT::rftfsub(int n, float* pArr, int nc, float* pC)
{
    int m = n >> 1;
    int ks = 2*nc/m;
    int kk = 0;

    for(int j = 2; j < m; j += 2)
    {
        kk += ks;

        int k = n-j;
        float wkr = 0.5 - pC[nc-kk];
        float wki = pC[kk];
        float xr = pArr[j]-pArr[k];
        float xi = pArr[j+1]+pArr[k+1];
        float yr = wkr*xr-wki*xi;
        float yi = wkr*xi+wki*xr;
        pArr[j] -= yr;
        pArr[j+1] -= yi;
        pArr[k] += yr;
        pArr[k+1] -= yi;
    }
}


void FFT::rftbsub(int n, float* pArr, int nc, float* pC)
{
    int m = n >> 1;
    int ks = 2*nc/m;
    int kk = 0;

    for(int j = 2; j < m; j += 2)
    {
        int k = n-j;

        kk += ks;
        float wkr = 0.5-pC[nc-kk];
        float wki = pC[kk];
        float xr = pArr[j]-pArr[k];
        float xi = pArr[j+1]+pArr[k+1];
        float yr = wkr*xr+wki*xi;
        float yi = wkr*xi-wki*xr;

        pArr[j] -= yr;
        pArr[j+1] -= yi;
        pArr[k] += yr;
        pArr[k+1] -= yi;
    }
}

void FFT::dctsub(int n, float* pArr, int nc, float* pC)
{
    int m = n >> 1;
    int ks = nc/n;
    int kk = 0;

    for(int j = 1; j < m; ++j)
    {
        kk += ks;

        int k = n-j;
        float wkr = pC[kk]-pC[nc-kk];
        float wki = pC[kk]+pC[nc-kk];
        float xr = wki*pArr[j]-wkr*pArr[k];

        pArr[j] = wkr * pArr[j]+wki*pArr[k];
        pArr[k] = xr;
    }
    pArr[m] *= pC[0];
}


void FFT::dstsub(int n, float* pArr, int nc, float* pC)
{
    int m = n >> 1;
    int ks = nc / n;
    int kk = 0;

    for(int j = 1; j < m; ++j)
    {
        kk += ks;

        int k = n-j;
        float wkr = pC[kk]-pC[nc-kk];
        float wki = pC[kk]+pC[nc-kk];
        float xr = wki*pArr[k]-wkr*pArr[j];

        pArr[k] = wkr*pArr[k]+wki*pArr[j];
        pArr[j] = xr;
    }
    pArr[m] *= pC[0];
}

void FFT::scale(float m, float* pArr, bool bComplex)
{
    float norm = (float)(1.0/m);
    int n2;
    if(bComplex)
        n2 = 2*m_n;
    else
        n2 = m_n;

#ifdef USE_CDFT_THREADS
    if(n2 >= CDFT_THREADS_BEGIN_N)
    {
        cdft_thread_t th[4];
        scale_arg_st ag[4];
        int nthread = 4;

        for(int idx = 0; idx < nthread; ++idx)
        {
            int firstIdx = idx * k;
            int lastIdx = (idx == (nthread - 1)) ? n2 : firstIdx + k;

            ag[idx].pArr = pArr[firstIdx];
            ag[idx].maxIdx = lastIdx-firstIdx;
            ag[idx].norm = norm;

            cdft_thread_create(&th[idx], scale_th, &ag[idx]);
        }

        for(int idx = 0; idx < nthread; ++idx)
            cdft_thread_wait(th[idx]);
    }
    else
#endif /* USE_CDFT_THREADS */
    {
        for(int idx = 0; idx < n2; ++idx)
            pArr[idx] *= norm;

    }
}

#ifdef USE_CDFT_THREADS
void* FFT::scale_th(void* p)
{
    int maxIdx = ((scale_arg_t*) p)->maxIdx;
    float norm = ((scale_arg_t*) p)->norm;
    float* pArr = ((scale_arg_t*) p)->pArr;

    for (int idx = 0; idx < maxIdx; ++idx)
        pArr[idx] *= norm;

    return (void*)0;
}
#endif /* USE_CDFT_THREADS */