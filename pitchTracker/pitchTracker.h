/*
    Reference:
        https://github.com/JorenSix/TarsosDSP/blob/master/src/core/be/tarsos/dsp/pitch/FastYin.java
*/

#ifndef NATIVE_AUDIO_PITCH_TRACKER_H
#define NATIVE_AUDIO_PITCH_TRACKER_H

#include <cstdint>
#include "audioCommon.h"

#define DECIBEL_MIN                     -90
#define DECIBEL_VOICE                   -24
#define DECIBEL_MAX                     0

#define ERROR_NONE                      0
#define ERROR_BUFFER_FULL               1
#define ERROR_INVALID_BUFFER            2
#define ERROR_PITCH_TRACKER_NOT_SET     3

#define COMPUTE_BUFFER_MS               50  // 50 ms

#include <stdint.h>

#define YIN_DEFAULT_THRESHOLD 0.20

class PitchTracker
{
public:
    explicit PitchTracker();
    ~PitchTracker();

public:
    bool                InitPitchTracker( SampleFormat* sampleFormat );
    int                 AddBuffer( int16_t* pSampleBuf, uint32_t nBufSize, uint32_t& nWritten );
    int                 AddBuffer( float* pSampleBuf, uint32_t nBufSize, uint32_t& nWritten );
    uint32_t            CurrentBufSize();

    bool                ReadyToComputePitch();
    bool                PitchComputed();
    bool                ComputePitch( float fThreshold );

    void                ClearCache();

    float               GetFrequency();
    float               GetLinearVolume();
    int8_t              GetNote();
    short               GetDecibel();

private:
    void                difference();
    void                cumulativeMeanNormalizedDifference();
    int32_t             absoluteThreshold();
    float               parabolicInterpolation(int32_t nTauEstimate);
    bool                measurePitch();

private:
    SampleFormat        m_sampleInfo;
    uint32_t            m_nSampleSize;
    uint32_t            m_nYinBufSize;

    int16_t*            m_pBuf;
    uint32_t            m_nCurBufSize;
    float*              m_pYinBuf;
    float               m_fThreshold;

    int32_t             m_maxValue;
    float               m_frequency;
    float               m_linearVolume;
    float               m_fProbability;
    int8_t              m_note;
    short               m_dB;

    bool                m_bReadyToCompute;
    bool                m_bComputed;

    bool                m_bIsSet;
};

#endif  // NATIVE_AUDIO_PITCH_TRACKER_H
