/*
    Reference:
        https://github.com/JorenSix/TarsosDSP/blob/master/src/core/be/tarsos/dsp/util/fft/FloatFFT.java
*/

#ifndef NATIVE_FFT_H_
#define NATIVE_FFT_H_

#include <complex>
#include <vector>

#ifdef USE_CDFT_PTHREADS
#ifndef USE_CDFT_THREADS
#define USE_CDFT_THREADS
#endif
#endif


#ifdef USE_CDFT_THREADS
typedef struct cdft_arg_st {
    int n0;
    int n;
    float* pArr;
    int nw;
    float* pW;

    cdft_arg_t(): n0(0), n(0), pArr(nullptr), nw(0), pW(nullptr){}
} cdft_arg_t;

typedef struct scale_arg_st {
    int maxIdx;
    float* pArr;
    float norm;

    scale_arg_t(): firstIdx(0), maxIdx(0), pArr(nullptr), norm(0.f){}
} scale_arg_t;

typedef struct bluestein_complex_arg1_st {
    int maxIdx;

    float* pAk;
    float* pArr;
    float* pBk1;

    int isign;

    bluestein_complex_arg1_t(): maxIdx(0), pAk(nullptr), pArr(nullptr)
                                , pBk1(nullptr), isign(0){}
} bluestein_complex_arg1_t;

typedef struct bluestein_complex_arg2_st {
    int maxIdx;

    float* pAk;
    float* pBk2;

    int isign;

    bluestein_complex_arg2_t(): maxIdx(0), pAk(nullptr), pBk2(nullptr), isign(0){}
} bluestein_complex_arg2_t;
#endif /* USE_CDFT_THREADS */

class FFT
{
public:
    enum Plans {
        SPLIT_RADIX, MIXED_RADIX, BLUESTEIN
    };
    typedef std::complex<float> Complex;
    /* Initializes FFT. n must be a power of 2. */

public:
    FFT( int n );
    ~FFT();

    void    complexForward(float* pArr);

    void    complexInverse(float* pArr, bool scale);
private:
    bool    IsPowerOf2(int x);
    int     prevPow2(int x);
    int     nextPow2(int x);
    int     getReminder(int n);
    void    cffti(int n, int offw);
    void    cffti();
    void    rffti();
    void    bluesteini();

    void    cfftf(float* pArr, int isign);
    void    passf2(int ido, int l1, float* pIn, float* pOut, int offset, int isign);
    void    passf3(int ido, int l1, float* pIn, float* pOut, int offset, int isign);
    void    passf4(int ido, int l1, float* pIn, float* pOut, int offset, int isign);
    void    passf5(int ido, int l1, float* pIn, float* pOut, int offset, int isign);
    void    passfg(int* pNac, int ido, int ip, int l1, int idl1, float* pIn
                    , float* pOut, int offset, int isign);

    void    cdft(int n, int isgn, float* pArr, int* pIp, float* pW);      // Complex Discrete Fourier Transform
    void    rdft(int n, int isgn, float* pArr, int* pIp, float* pW);      // Real Discrete Fourier Transform
    void    ddct(int n, int isgn, float* pArr, int* pIp, float* pW);      // Discrete Cosine Transform
    void    ddst(int n, int isgn, float* pArr, int* pIp, float* pW);      // Discrete Sine Transform
    void    dfct(int n, float* pArr, float* pT, int* pIp, float* pW);    // Cosine Transform of RDFT (Real Symmetric DFT)
    void    dfst(int n, float* pArr, float* pT, int* pIp, float* pW);    // Sine Transform of RDFT (Real Anti-symmetric DFT)

    void    makewt(int nw, int* pIp, float* pW );
    void    makeipt(int nw, int* pIp);
    void    makect(int nc, int* pIp, float* pC);

    void    bluestein_complex(float* pArr, int isign);

    void    cftfsub(int n, float* pArr, int* pIp, int nw, float* pW);
    void    cftbsub(int n, float* pArr, int* pIp, int nw, float* pW);
    void    bitrv2(int n, int* pIp, float* pArr);
    void    bitrv2conj(int n, int* pIp, float* pArr);
    void    bitrv216(float* pArr);
    void    bitrv216neg(float* pArr);
    void    bitrv208(float* pArr);
    void    bitrv208neg(float* pArr);
    void    cftf1st(int n, float* pArr, float* pW);
    void    cftb1st(int n, float* pArr, float* pW);

#ifdef USE_CDFT_THREADS
    void    cftrec4_th(int n, float* pArr, int nw, float* pW);
    void*   cftrec1_th(void* p);
    void*   cftrec2_th(void* p);

    void*   scale_th(void* p);

    void*   bluestein_complex1_th(void* p);
    void*   bluestein_complex2_th(void* p);
    void*   bluestein_complex3_th(void* p);
#endif /* USE_CDFT_THREADS */

    void    cftrec4(int n, float* pArr, int nw, float* pW);
    int     cfttree(int n, int j, int k, float* pArr, int nw, float* pW);
    void    cftleaf(int n, int isplt, float* pArr, int nw, float* pW);
    void    cftmdl1(int n, float* pArr, float* pW);
    void    cftmdl2(int n, float* pArr, float* pW);
    void    cftfx41(int n, float* pArr, int nw, float* pW);
    void    cftf161(float* pArr, float* pW);
    void    cftf162(float* pArr, float* pW);
    void    cftf081(float* pArr, float* pW);
    void    cftf082(float* pArr, float* pW);
    void    cftf040(float* pArr);
    void    cftb040(float* pArr);
    void    cftx020(float* pArr);
    void    rftfsub(int n, float* pArr, int nc, float* pC);
    void    rftbsub(int n, float* pArr, int nc, float* pC);
    void    dctsub(int n, float* pArr, int nc, float* pC);
    void    dstsub(int n, float* pArr, int nc, float* pC);
    void    scale(float m, float* pArr, bool bComplex);

private:
    int                 m_n;
    float*              m_pW;
    int*                m_pIp;
    int                 m_nw;
    int                 m_nc;

    std::vector<int>    m_vecFactors;
    Plans               m_plan;
    int                 m_nBluestein;

    float*              m_pBk1;
    float*              m_pBk2;
    float*              m_pWtable;
    float*              m_pWtable_r;
};

#endif //NATIVE_FFT_H_