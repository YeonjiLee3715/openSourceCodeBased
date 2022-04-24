# openSourceCodeBased

open source projects based project

- pitchTracker:
  - c++ based pitch tracker
  - Applied FastYIN theory based on YIN pitch detection theory - Using FFT for Autocorrelation and Difference to speed up the calculation.

- fluidsynth:
  - Version: 1.1.11
  - Can work on Android
  - Gain adjustment function for each channel
  - Master volume, reverb volume, chorus volume and preset volume control
  - Note Seeking
  - Roland exclusive message handling
  - midi channel mute and solo
  - Minor errors fixed

# Reference:

- pitchTracker:
  - Based project: 
        1. [TarsosDSP](https://github.com/JorenSix/TarsosDSP)
        2. [sevagh/pitch-detection](https://github.com/sevagh/pitch-detection/tree/master/misc/yin)
  - [Yin](http://audition.ens.fr/adc/pdf/2002_JASA_YIN.pdf)
  - [PYin](https://www.eecs.qmul.ac.uk/~simond/pub/2014/MauchDixon-PYIN-ICASSP2014.pdf)
  - [MPM](https://www.cs.otago.ac.nz/students/postgrads/tartini/papers/A_Smarter_Way_to_Find_Pitch.pdf)
  - [Wavelets](https://courses.physics.illinois.edu/phys406/sp2017/NSF_REU_Reports/2005_reu/Real-Time_Time-Domain_Pitch_Tracking_Using_Wavelets.pdf)
  - [FFT](https://towardsdatascience.com/fast-fourier-transform-937926e591cb)

- fluidsynth:
  - Based project: [fluidsynth](https://github.com/FluidSynth/fluidsynth)
  - [SoundFont® Technical Specification 2.04](http://www.synthfont.com/sfspec24.pdf)
  - [Roland System Exclusive Implementation](http://www.chromakinetics.com/handsonic/rolSysEx.htm)
  - [MIDI Implementation](http://cdn.roland.com/assets/media/pdf/F-20_MIDI_Imple_e01_W.pdf)
  - [PPQN](http://midi.teragonaudio.com/tech/midifile/ppqn.htm)
  - [MIDI system exclusive message](https://www.recordingblogs.com/wiki/midi-system-exclusive-message)
  - [MIDI 1.0 Universal System Exclusive Messages](https://www.midi.org/specifications-old/item/table-4-universal-system-exclusive-messages)