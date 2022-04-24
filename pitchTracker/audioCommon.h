/*
 * Copyright 2015 The Android Open Source Project
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

#ifndef NATIVE_AUDIO_AUDIO_COMMON_H
#define NATIVE_AUDIO_AUDIO_COMMON_H

#include <SLES/OpenSLES.h>
#include <SLES/OpenSLES_Android.h>

#include "androidDebug.h"
#include "debugUtils.h"
#include "bufManager.h"

struct SampleFormat
{
    uint32_t sampleRate;
    uint32_t framesPerBuf;
    uint16_t channels;
    uint16_t pcmFormat;  // 8 bit, 16 bit, 24 bit ...
    uint32_t representation;  // android extensions

    SampleFormat():sampleRate(0), framesPerBuf(0), channels(0)
                    , pcmFormat(0), representation(0){}
};

extern void ConvertToSLSampleFormat(SLAndroidDataFormat_PCM_EX* pFormat
                                    , SampleFormat* format);

/*
 * GetSystemTicks(void):  return the time in micro sec
 */
__inline__ uint64_t GetSystemTicks(void) {
  struct timeval Time;
  gettimeofday(&Time, NULL);

  return (static_cast<uint64_t>(1000000) * Time.tv_sec + Time.tv_usec);
}

#define SLASSERT(x)                   \
  do {                                \
    assert(SL_RESULT_SUCCESS == (x)); \
    (void)(x);                        \
  } while (0)


#endif  // NATIVE_AUDIO_AUDIO_COMMON_H
