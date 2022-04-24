/*
    Reference:
        https://github.com/JorenSix/TarsosDSP/blob/master/src/core/be/tarsos/dsp/pitch/FastYin.java
*/

#include "pitchTracker.h"
#include <cmath>
#include "FFT.h"

PitchTracker::PitchTracker()
        : m_nSampleSize(0), m_nYinBufSize(0), m_pBuf(nullptr), m_nCurBufSize(0), m_pYinBuf(nullptr)
        , m_fThreshold(0.f), m_maxValue(0), m_frequency(0.f), m_linearVolume(0.f)
        , m_fProbability(0.f), m_note(-1), m_dB(DECIBEL_MIN), m_bReadyToCompute(false)
        , m_bComputed(false), m_bIsSet(false)
{

}

PitchTracker::~PitchTracker()
{
    if( m_pBuf != nullptr )
    {
        delete[] m_pBuf;
        m_pBuf = nullptr;
    }

    if( m_pYinBuf != nullptr )
    {
        delete[] m_pYinBuf;
        m_pYinBuf = nullptr;
    }
}

bool PitchTracker::InitPitchTracker( SampleFormat* sampleFormat )
{
    if( m_bIsSet )
    {
        LOGI( "PeachTracker already set" );
        return true;
    }

    bool isSuccess = false;
    if( sampleFormat == nullptr )
        return isSuccess;

    m_sampleInfo = (*sampleFormat);

    m_nSampleSize = (uint32_t)round(((float)(m_sampleInfo.pcmFormat>>3)*((float)m_sampleInfo.sampleRate/1000.0)/2.0/*buf_per_ms*/)*COMPUTE_BUFFER_MS);
    m_nYinBufSize = m_nSampleSize/2;

    m_pBuf = new int16_t[m_nSampleSize];
    if( m_pBuf == nullptr )
    {
        LOGE( "Failed to create wavelet sample array. nSampleSize: %d", m_nSampleSize );
        return isSuccess;
    }

    m_pYinBuf = new float[m_nYinBufSize];
    if( m_pYinBuf == nullptr )
    {
        LOGE( "Failed to create calc wavelet sample array. nSampleSize: %d", m_nSampleSize );
        return isSuccess;
    }

    memset( m_pBuf, 0, sizeof(int16_t)*m_nSampleSize );
    memset( m_pYinBuf, 0, sizeof(float)*m_nYinBufSize );

    m_bIsSet = true;

    isSuccess = true;

    return isSuccess;
}

int PitchTracker::AddBuffer( int16_t* pSampleBuf, uint32_t nBufSize, uint32_t& nWritten )
{
    nWritten = 0;

    if( m_bIsSet == false )
    {
        LOGE( "Pitch tracker not set yet" );
        return ERROR_PITCH_TRACKER_NOT_SET;
    }

    if( m_bReadyToCompute )
        return ERROR_BUFFER_FULL;

    if( pSampleBuf == nullptr || nBufSize <= 0 )
    {
        LOGE( "Wrong sample buffer. pSampleBuf is null(%s) or nBufSize(%u) is less than or equal to 0"
        , (pSampleBuf==nullptr?"T":"F"), nBufSize );
        return ERROR_INVALID_BUFFER;
    }

    uint32_t idxCurBuf = m_nCurBufSize;
    uint32_t idxCurSampleBuf = 0;

    if( m_sampleInfo.channels == 2 )
    {
        while( idxCurBuf < m_nSampleSize && idxCurSampleBuf < nBufSize )
        {
            int32_t mix = (((int32_t)pSampleBuf[idxCurSampleBuf]+(int32_t)pSampleBuf[idxCurSampleBuf+1])/(int32_t)2);
            if( mix < SHRT_MIN )
                mix = SHRT_MIN;
            else if( mix > SHRT_MAX )
                mix = SHRT_MAX;

            m_pBuf[idxCurBuf] = (int16_t)mix;

            int32_t localMax = abs( m_pBuf[idxCurBuf] );
            if( localMax > m_maxValue )
                m_maxValue = localMax;

            idxCurSampleBuf += 2;
            ++idxCurBuf;
            ++m_nCurBufSize;
        }
    }
    else
    {
        while( idxCurBuf < m_nSampleSize && idxCurSampleBuf < nBufSize )
        {
            m_pBuf[idxCurBuf] = pSampleBuf[idxCurSampleBuf++];

            int32_t localMax = abs( m_pBuf[idxCurBuf] );
            if( localMax > m_maxValue )
                m_maxValue = localMax;

            ++idxCurBuf;
            ++m_nCurBufSize;
        }
    }

    m_linearVolume = (float)m_maxValue/(float)SHRT_MAX;
    m_bComputed = false;
    nWritten = idxCurSampleBuf;

    if( m_nCurBufSize == m_nSampleSize )
    {
        m_bReadyToCompute = true;
        m_dB = 20 * (log(m_linearVolume)/log(10));
    }
    else
        m_bReadyToCompute = false;

    if( idxCurSampleBuf < nBufSize )
        return ERROR_BUFFER_FULL;

    return ERROR_NONE;
}

int PitchTracker::AddBuffer( float* pSampleBuf, uint32_t nBufSize, uint32_t& nWritten )
{
    nWritten = 0;

    if( m_bIsSet == false )
    {
        LOGE( "Pitch tracker not set yet" );
        return ERROR_PITCH_TRACKER_NOT_SET;
    }

    if( m_bReadyToCompute )
        return ERROR_BUFFER_FULL;

    if( pSampleBuf == nullptr || nBufSize <= 0 )
    {
        LOGE( "Wrong sample buffer. pSampleBuf is null(%s) or nBufSize(%u) is less than or equal to 0"
        , (pSampleBuf==nullptr?"T":"F"), nBufSize );
        return ERROR_INVALID_BUFFER;
    }

    uint32_t idxCurBuf = m_nCurBufSize;
    uint32_t idxCurSampleBuf = 0;

    if( m_sampleInfo.channels == 2 )
    {
        while( idxCurBuf < m_nSampleSize && idxCurSampleBuf < nBufSize )
        {
            float mix = ((pSampleBuf[idxCurSampleBuf]+pSampleBuf[idxCurSampleBuf+1])/2)*32768.0;
            if( mix < SHRT_MIN )
                mix = SHRT_MIN;
            else if( mix > SHRT_MAX )
                mix = SHRT_MAX;

            m_pBuf[idxCurBuf] = (int16_t)mix;

            int32_t localMax = abs( m_pBuf[idxCurBuf] );
            if( localMax > m_maxValue )
                m_maxValue = localMax;

            idxCurSampleBuf += 2;
            ++idxCurBuf;
            ++m_nCurBufSize;
        }
    }
    else
    {
        while( idxCurBuf < m_nSampleSize && idxCurSampleBuf < nBufSize )
        {
            float fSample = pSampleBuf[idxCurSampleBuf++]*32768.0;
            if( fSample < SHRT_MIN )
                fSample = SHRT_MIN;
            else if( fSample > SHRT_MAX )
                fSample = SHRT_MAX;

            m_pBuf[idxCurBuf] = (int16_t)fSample;

            int32_t localMax = abs( m_pBuf[idxCurBuf] );
            if( localMax > m_maxValue )
                m_maxValue = localMax;

            ++idxCurBuf;
            ++m_nCurBufSize;
        }
    }

    m_linearVolume = (float)m_maxValue/(float)SHRT_MAX;
    m_bComputed = false;
    nWritten = idxCurSampleBuf;

    if( m_nCurBufSize == m_nSampleSize )
    {
        m_bReadyToCompute = true;
        m_dB = 20 * (log(m_linearVolume)/log(10));
    }
    else
        m_bReadyToCompute = false;

    if( idxCurSampleBuf < nBufSize )
        return ERROR_BUFFER_FULL;

    return ERROR_NONE;
}

uint32_t PitchTracker::CurrentBufSize()
{
    return m_nCurBufSize;
}

bool PitchTracker::ReadyToComputePitch()
{
    if( m_bIsSet == false )
        return false;

    if( m_bComputed )
        return false;

    return m_bReadyToCompute;
}

bool PitchTracker::PitchComputed()
{
    return m_bComputed;
}

bool PitchTracker::ComputePitch( float fThreshold )
{
    if( m_bIsSet == false || m_bReadyToCompute == false )
        return false;

    if( m_bComputed )
        return true;

    m_fThreshold = fThreshold;

    bool isSuccess = false;

    /* Step 1: Calculates the squared difference of the signal with a shifted version of itself. */
    difference();

    /* Step 2: Calculate the cumulative mean on the normalised difference calculated in step 1 */
    cumulativeMeanNormalizedDifference();

    /* Step 3: Search through the normalised cumulative mean array and find values that are over the threshold */
    int32_t nTauEstimate = absoluteThreshold();

    /* Step 5: Interpolate the shift value (tau) to improve the pitch estimate. */
    if( nTauEstimate <= -1 )
        return isSuccess;

    m_frequency = ((float)m_sampleInfo.sampleRate)/parabolicInterpolation(nTauEstimate);

    isSuccess = measurePitch();

    return isSuccess;
}

void PitchTracker::difference()
{
/*
    // Calculate the difference for difference shift values (tau) for the half of the samples
    for( int32_t tau = 0 ; tau < m_nSampleSize; ++tau )
    {
        // Take the difference of the signal with a shifted version of itself, then square it.
        // (This is the Yin algorithm's tweak on autocorellation)
        float entry = 0.f;
        for( int32_t idx = 0; idx < m_nSampleSize-tau; ++idx )
        {
            float delta = m_pBuf[idx] - m_pBuf[idx+tau];
            entry += delta*delta;
        }
        m_pYinBuf[tau] = entry;
    }
*/
    uint32_t nCalcBufSize = m_nSampleSize*2;

    float* pPowerTerms = new float[m_nYinBufSize];

    float* pAudioBufferFFT = new float[nCalcBufSize];
    float* pKernel = new float[nCalcBufSize];
    float* pYinStyleACF = new float[nCalcBufSize];

    FFT fft( m_nSampleSize );

    memset( pPowerTerms, 0, sizeof(float)*m_nYinBufSize );

    memset( pAudioBufferFFT, 0, sizeof(float)*nCalcBufSize );
    memset( pKernel, 0, sizeof(float)*nCalcBufSize );
    memset( pYinStyleACF, 0, sizeof(float)*nCalcBufSize );

    for( uint32_t idx = 0; idx < m_nYinBufSize; ++idx )
        pPowerTerms[0] += m_pBuf[idx]*m_pBuf[idx];

    for( uint32_t tau = 1; tau < m_nYinBufSize; ++tau )
        pPowerTerms[tau] = pPowerTerms[tau-1]-m_pBuf[tau-1]*m_pBuf[tau-1]+m_pBuf[tau+m_nYinBufSize]*m_pBuf[tau+m_nYinBufSize];

    // YIN-STYLE AUTOCORRELATION via FFT
    // 1. data
    for(uint32_t idx = 0; idx < m_nSampleSize; ++idx)
    {
        pAudioBufferFFT[2*idx] = m_pBuf[idx];
        pAudioBufferFFT[2*idx+1] = 0;
    }
    fft.complexForward( pAudioBufferFFT );

    // 2. half of the data, disguised as a convolution kernel
    for(uint32_t idx = 0; idx < m_nYinBufSize; ++idx)
    {
        pKernel[2*idx] = m_pBuf[(m_nYinBufSize-1)-idx];
        pKernel[2*idx+1] = 0;
        pKernel[2*idx+m_nSampleSize] = 0;
        pKernel[2*idx+m_nSampleSize+1] = 0;
    }
    fft.complexForward( pKernel );

    // 3. convolution via complex multiplication
    for(uint32_t idx = 0; idx < m_nSampleSize; ++idx)
    {
        pYinStyleACF[2*idx] = pAudioBufferFFT[2*idx]*pKernel[2*idx]-pAudioBufferFFT[2*idx+1]*pKernel[2*idx+1]; // real
        pYinStyleACF[2*idx+1] = pAudioBufferFFT[2*idx+1]*pKernel[2*idx]+pAudioBufferFFT[2*idx]*pKernel[2*idx+1]; // imaginary
    }

    fft.complexInverse( pYinStyleACF, true );

    // CALCULATION OF difference function
    // ... according to (7) in the Yin paper.
    for(uint32_t idx = 0; idx < m_nYinBufSize; ++idx)
    {
        // taking only the real part
        m_pYinBuf[idx] = pPowerTerms[0]+pPowerTerms[idx]-2*pYinStyleACF[2*(m_nYinBufSize-1+idx)];
    }

    if( pPowerTerms != nullptr )
        delete[] pPowerTerms;
    if( pAudioBufferFFT != nullptr )
        delete[] pAudioBufferFFT;
    if( pKernel != nullptr )
        delete[] pKernel;
    if( pYinStyleACF != nullptr )
        delete[] pYinStyleACF;
}

void PitchTracker::cumulativeMeanNormalizedDifference()
{
    float runningSum = 0.f;
    m_pYinBuf[0] = 1; //if tau == 0, d't(tau) = 1

    // otherwise
    /* Sum all the values in the autocorellation buffer and nomalise the result, replacing
     * the value in the autocorellation buffer with a cumulative mean of the normalised difference */
    for(int32_t tau = 1; tau < m_nYinBufSize; ++tau)
    {
        runningSum += m_pYinBuf[tau];
        m_pYinBuf[tau] *= tau / runningSum;
    }
}

int32_t PitchTracker::absoluteThreshold()
{
    int32_t tau = 0;

    /* Search through the array of cumulative mean values, and look for ones that are over the threshold
     * The first two positions in yinBuffer are always so start at the third (index 2) */
    for( tau = 2; tau < m_nYinBufSize; ++tau )
    {
        if( m_pYinBuf[tau] < m_fThreshold )
        {
            while( tau+1 < m_nYinBufSize && m_pYinBuf[tau+1] < m_pYinBuf[tau] )
                ++tau;

            /* found tau, exit loop and return
             * store the probability
             * From the YIN paper: The yin->threshold determines the list of
             * candidates admitted to the set, and can be interpreted as the
             * proportion of aperiodic power tolerated
             * within a periodic signal.
             *
             * Since we want the periodicity and and not aperiodicity:
             * periodicity = 1 - aperiodicity */

            m_fProbability = 1-m_pYinBuf[tau];
            break;
        }
    }

    /* if no pitch found, tau => -1 */
    if( tau == m_nYinBufSize || m_pYinBuf[tau] >= m_fThreshold || m_fProbability > 1.f )
    {
        tau = -1;
        m_fProbability = 0;
    }

    return tau;
}

float PitchTracker::parabolicInterpolation( int32_t nTauEstimate )
{
    float fBetterTau = 0.f;
    int32_t x0 = 0;
    int32_t x2 = 0;

    if( nTauEstimate < 1 )
        x0 = nTauEstimate;
    else
        x0 = nTauEstimate-1;

    if( nTauEstimate +1 < m_nYinBufSize )
        x2 = nTauEstimate+1;
    else
        x2 = nTauEstimate;

    if( x0 == nTauEstimate )
    {
        if( m_pYinBuf[nTauEstimate] <= m_pYinBuf[x2] )
            fBetterTau = (float)nTauEstimate;
        else
            fBetterTau = (float)x2;
    }
    else if( x2 == nTauEstimate )
    {
        if( m_pYinBuf[nTauEstimate] <= m_pYinBuf[x0] )
            fBetterTau = (float)nTauEstimate;
        else
            fBetterTau = (float)x0;
    }
    else
    {
        float s0 = m_pYinBuf[x0];
        float s1 = m_pYinBuf[nTauEstimate];
        float s2 = m_pYinBuf[x2];
        /* fixed AUBIO implementation, thanks to Karl Helgason:
           (2.0f * s1 - s2 - s0) was incorrectly multiplied with -1 */
        fBetterTau = ((float)nTauEstimate+((s2-s0)/(2*(2*s1-s2-s0))));
    }

    return fBetterTau;
}

bool PitchTracker::measurePitch()
{
    double kPitchA = 440.0; // Hz.= midi 69
    m_note = (int)round(12.0*log((double)m_frequency/kPitchA)/log(2)+69.0);
    if(m_note >= 0 && m_note <= SCHAR_MAX )
        return true;

    m_note = -1;
    return false;
}

void PitchTracker::ClearCache()
{
    m_bReadyToCompute = false;
    m_bComputed = false;
    m_maxValue = 0;
    m_frequency = 0.f;
    m_linearVolume = 0.f;
    m_dB = DECIBEL_MIN;
    m_nCurBufSize = 0;
    m_fThreshold = 0.f;
    m_fProbability = 0.f;

    memset( m_pBuf, 0, sizeof(int16_t)*m_nSampleSize );
    memset( m_pYinBuf, 0, sizeof(float)*m_nYinBufSize );
}

float PitchTracker::GetFrequency()
{
    return m_frequency;
}

float PitchTracker::GetLinearVolume()
{
    return m_linearVolume;
}

int8_t PitchTracker::GetNote()
{
    return m_note;
}

short PitchTracker::GetDecibel()
{
    return m_dB;
}