/* FluidSynth - A Software Synthesizer
 *
 * Copyright (C) 2003  Peter Hanappe and others.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * as published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this library; if not, write to the Free
 * Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 * 02110-1301, USA
 */

#ifndef _FLUID_CHAN_H
#define _FLUID_CHAN_H

#include "fluidsynth_priv.h"
#include "fluid_midi.h"
#include "fluid_tuning.h"

/*
 * fluid_channel_t
 *
 * Mutual exclusion notes (as of 1.1.2):
 * None - everything should have been synchronized by the synth.
 */
struct _fluid_channel_t
{
  fluid_mutex_t mutex;                  /* Lock for thread sensitive parameters */

  fluid_synth_t* synth;                 /**< Parent synthesizer instance */
  int channum;                          /**< MIDI channel number */

  /* Poly Mono variables see macro access description */
  int mode;								/**< Poly Mono mode */
  int mode_val;							/**< number of channel in basic channel group */

  int sfont_bank_prog;                  /**< SoundFont ID (bit 21-31), bank (bit 7-20), program (bit 0-6) */
  fluid_preset_t* preset;               /**< Selected preset */

  int channel_pressure;                 /**< MIDI channel pressure */
  int pitch_bend;                       /**< Current pitch bend value */
  int pitch_wheel_sensitivity;          /**< Current pitch wheel sensitivity */

  fluid_real_t chan_gain;               /**< default 1.0 */

  int cc[128];                          /**< MIDI controller values */
  unsigned char key_pressure[128];      /**< MIDI polyphonic key pressure from [0;127] */

  /* Sostenuto order id gives the order of SostenutoOn event.
     This value is useful to known when the sostenuto pedal is depressed
     (before or after a key note). We need to compare SostenutoOrderId with voice id.
   */
  unsigned int  sostenuto_orderid;
  int interp_method;                    /**< Interpolation method (enum fluid_interp) */
  fluid_tuning_t* tuning;               /**< Micro tuning */
  int tuning_bank;                      /**< Current tuning bank number */
  int tuning_prog;                      /**< Current tuning program number */

  /* NRPN system */
  int nrpn_select;      /* Generator ID of SoundFont NRPN message */
  int nrpn_active;      /* 1 if data entry CCs are for NRPN, 0 if RPN */

  /* The values of the generators, set by NRPN messages, or by
   * fluid_synth_set_gen(), are cached in the channel so they can be
   * applied to future notes. They are copied to a voice's generators
   * in fluid_voice_init(), which calls fluid_gen_init().  */
  fluid_real_t gen[GEN_LAST];

  /* By default, the NRPN values are relative to the values of the
   * generators set in the SoundFont. For example, if the NRPN
   * specifies an attack of 100 msec then 100 msec will be added to the
   * combined attack time of the sound font and the modulators.
   *
   * However, it is useful to be able to specify the generator value
   * absolutely, completely ignoring the generators of the SoundFont
   * and the values of modulators. The gen_abs field, is a boolean
   * flag indicating whether the NRPN value is absolute or not.
   */
  char gen_abs[GEN_LAST];

  /* Drum channel flag, CHANNEL_TYPE_MELODIC, or CHANNEL_TYPE_DRUM. */
  int channel_type;

  int mute;

};

fluid_channel_t* new_fluid_channel(fluid_synth_t* synth, int num);
void fluid_channel_init_ctrl(fluid_channel_t* chan, int is_all_ctrl_off);
int delete_fluid_channel(fluid_channel_t* chan);
void fluid_channel_reset(fluid_channel_t* chan);
int fluid_channel_set_preset(fluid_channel_t* chan, fluid_preset_t* preset);
fluid_preset_t* fluid_channel_get_preset(fluid_channel_t* chan);
void fluid_channel_set_sfont_bank_prog(fluid_channel_t* chan, int sfont,
                                       int bank, int prog);
void fluid_channel_set_bank_lsb(fluid_channel_t* chan, int banklsb);
void fluid_channel_set_bank_msb(fluid_channel_t* chan, int bankmsb);
void fluid_channel_get_sfont_bank_prog(fluid_channel_t* chan, int *sfont,
                                       int *bank, int *prog);
int fluid_channel_get_num(fluid_channel_t* chan);
void fluid_channel_set_interp_method(fluid_channel_t* chan, int new_method);
int fluid_channel_get_interp_method(fluid_channel_t* chan);

#define fluid_channel_get_preset(chan)          ((chan)->preset)
#define fluid_channel_set_cc(chan, num, val) \
  ((chan)->cc[num] = (val))
#define fluid_channel_get_cc(chan, num) \
  ((chan)->cc[num])
#define fluid_channel_get_key_pressure(chan, key) \
  ((chan)->key_pressure[key])
#define fluid_channel_set_key_pressure(chan, key, val) \
  ((chan)->key_pressure[key] = (val))
#define fluid_channel_get_channel_pressure(chan) \
  ((chan)->channel_pressure)
#define fluid_channel_set_channel_pressure(chan, val) \
  ((chan)->channel_pressure = (val))
#define fluid_channel_get_pitch_bend(chan) \
  ((chan)->pitch_bend)
#define fluid_channel_set_pitch_bend(chan, val) \
  ((chan)->pitch_bend = (val))
#define fluid_channel_get_pitch_wheel_sensitivity(chan) \
  ((chan)->pitch_wheel_sensitivity)
#define fluid_channel_set_pitch_wheel_sensitivity(chan, val) \
  ((chan)->pitch_wheel_sensitivity = (val))
#define fluid_channel_get_num(chan)             ((chan)->channum)
#define fluid_channel_set_interp_method(chan, new_method) \
  ((chan)->interp_method = (new_method))
#define fluid_channel_get_interp_method(chan) \
  ((chan)->interp_method);
#define fluid_channel_set_tuning(_c, _t)        { (_c)->tuning = _t; }
#define fluid_channel_has_tuning(_c)            ((_c)->tuning != NULL)
#define fluid_channel_get_tuning(_c)            ((_c)->tuning)
#define fluid_channel_get_tuning_bank(chan)     \
  ((chan)->tuning_bank)
#define fluid_channel_set_tuning_bank(chan, bank) \
  ((chan)->tuning_bank = (bank))
#define fluid_channel_get_tuning_prog(chan)     \
  ((chan)->tuning_prog)
#define fluid_channel_set_tuning_prog(chan, prog) \
  ((chan)->tuning_prog = (prog))
#define fluid_channel_sustained(_c)             ((_c)->cc[SUSTAIN_SWITCH] >= 64)
#define fluid_channel_sostenuto(_c)             ((_c)->cc[SOSTENUTO_SWITCH] >= 64)
#define fluid_channel_set_gen(_c, _n, _v, _a)   { (_c)->gen[_n] = _v; (_c)->gen_abs[_n] = _a; }
#define fluid_channel_get_gen(_c, _n)           ((_c)->gen[_n])
#define fluid_channel_get_gen_abs(_c, _n)       ((_c)->gen_abs[_n])
#define fluid_channel_get_min_note_length_ticks(chan) \
  ((chan)->synth->min_note_length_ticks)

/* Interface to poly/mono mode variables */
enum fluid_channel_mode_flags_internal
{
    FLUID_CHANNEL_BASIC = 0x04,    /**< if flag set the corresponding midi channel is a basic channel */
    FLUID_CHANNEL_ENABLED = 0x08,  /**< if flag set the corresponding midi channel is enabled, else disabled, i.e. channel ignores any MIDI messages */

    /*
      FLUID_CHANNEL_LEGATO_PLAYING bit of channel mode keeps trace of the legato /staccato
      state playing.
      FLUID_CHANNEL_LEGATO_PLAYING bit is updated on noteOn/noteOff mono by the legato detector:
      - On noteOn, before inserting a new note into the monolist.
      - On noteOff, after removing a note out of the monolist.

      - On noteOn, this state is used by fluid_synth_noteon_mono_LOCAL()
      to play the current  note legato or staccato.
      - On noteOff, this state is used by fluid_synth_noteoff_mono_LOCAL()
      to play the current noteOff legato with the most recent note.
    */
    /* bit7, 1: means legato playing , 0: means staccato playing */
            FLUID_CHANNEL_LEGATO_PLAYING = 0x80
};

/* End of interface to monophonic list variables */

#endif /* _FLUID_CHAN_H */
