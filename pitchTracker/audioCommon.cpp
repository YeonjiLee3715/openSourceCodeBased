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

#include "audioCommon.h"

void ConvertToSLSampleFormat(SLAndroidDataFormat_PCM_EX* pFormat,
                             SampleFormat* pSampleInfo_) {
  if( pFormat == nullptr )
      return;

  memset(pFormat, 0, sizeof(*pFormat));

  pFormat->formatType = SL_ANDROID_DATAFORMAT_PCM_EX;
  // Only support 2 channels
  // For channelMask, refer to wilhelm/src/android/channels.c for details
  if (pSampleInfo_->channels <= 1) {
    pFormat->numChannels = 1;
    pFormat->channelMask = SL_SPEAKER_FRONT_LEFT;
  } else {
    pFormat->numChannels = 2;
    pFormat->channelMask = SL_SPEAKER_FRONT_LEFT | SL_SPEAKER_FRONT_RIGHT;
  }
  pFormat->sampleRate = pSampleInfo_->sampleRate*1000;

  /*
   * fixup for android extended representations...
   */
  pFormat->endianness = SL_BYTEORDER_LITTLEENDIAN;

  pFormat->representation = pSampleInfo_->representation;
  switch (pFormat->representation) {
    case SL_ANDROID_PCM_REPRESENTATION_UNSIGNED_INT:
      pFormat->bitsPerSample = SL_PCMSAMPLEFORMAT_FIXED_8;
          pFormat->containerSize = SL_PCMSAMPLEFORMAT_FIXED_8;
          break;
    case SL_ANDROID_PCM_REPRESENTATION_SIGNED_INT:
      pFormat->bitsPerSample =
              SL_PCMSAMPLEFORMAT_FIXED_16;  // supports 16, 24, and 32
          pFormat->containerSize = SL_PCMSAMPLEFORMAT_FIXED_16;
          break;
    case SL_ANDROID_PCM_REPRESENTATION_FLOAT:
      pFormat->bitsPerSample = SL_PCMSAMPLEFORMAT_FIXED_32;
          pFormat->containerSize = SL_PCMSAMPLEFORMAT_FIXED_32;
          break;
    case 0:
          pFormat->bitsPerSample = pSampleInfo_->pcmFormat;
          pFormat->containerSize = pSampleInfo_->pcmFormat;
          if( pSampleInfo_->pcmFormat == SL_PCMSAMPLEFORMAT_FIXED_8 )
            pFormat->representation = SL_ANDROID_PCM_REPRESENTATION_UNSIGNED_INT;
          else if( pSampleInfo_->pcmFormat == SL_PCMSAMPLEFORMAT_FIXED_16 )
            pFormat->representation = SL_ANDROID_PCM_REPRESENTATION_SIGNED_INT;
          else if( pSampleInfo_->pcmFormat == SL_PCMSAMPLEFORMAT_FIXED_32 )
            pFormat->representation = SL_ANDROID_PCM_REPRESENTATION_FLOAT;
      break;
    default:
      assert(0);
  }
}
