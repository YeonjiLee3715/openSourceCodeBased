# FluidSynth - A Software Synthesizer
#
# Copyright (C) 2003-2010 Peter Hanappe and others.
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public License
# as published by the Free Software Foundation; either version 2.1 of
# the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
# 02111-1307, USA

# CMake based build system. Pedro Lopez-Cabanillas <plcl@users.sf.net>

# FluidSynth package name
set ( PACKAGE "fluidsynth" )

# FluidSynth package version
set ( FLUIDSYNTH_VERSION_MAJOR 1 )
set ( FLUIDSYNTH_VERSION_MINOR 1 )
set ( FLUIDSYNTH_VERSION_MICRO 11 )
set ( VERSION "${FLUIDSYNTH_VERSION_MAJOR}.${FLUIDSYNTH_VERSION_MINOR}.${FLUIDSYNTH_VERSION_MICRO}" )
set ( FLUIDSYNTH_VERSION "\"${VERSION}\"" )

# *** NOTICE ***
# Update library version upon each release (follow these steps in order)
# if any source code changes: REVISION++
# if any interfaces added/removed/changed: REVISION=0
# if any interfaces removed/changed (compatibility broken): CURRENT++
# if any interfaces have been added: AGE++
# if any interfaces have been removed/changed (compatibility broken): AGE=0
# This is not exactly the same algorithm as the libtool one, but the results are the same.
set ( LIB_VERSION_CURRENT 1 )
set ( LIB_VERSION_AGE 7 )
set ( LIB_VERSION_REVISION 2 )
set ( LIB_VERSION_INFO
        "${LIB_VERSION_CURRENT}.${LIB_VERSION_AGE}.${LIB_VERSION_REVISION}" )

#set(OPENSLES_SUPPORT 1)
set(HAVE_DLFCN_H 1)
set(HAVE_ERRNO_H 1)
set(HAVE_FCNTL_H 1)
set(HAVE_INTTYPES_H 1)
set(HAVE_LIBPTHREAD TRUE)
set(HAVE_LIMITS_H 1)
set(HAVE_MATH_H 1)
set(HAVE_MEMORY_H 1)
set(HAVE_NETINET_IN_H 1)
set(HAVE_NETINET_TCP_H 1)
set(HAVE_PTHREAD_H 1)
set(HAVE_SIGNAL_H 1)
set(HAVE_STDARG_H 1)
set(HAVE_STDINT_H 1)
set(HAVE_STDIO_H 1)
set(HAVE_STDLIB_H 1)
set(HAVE_STRINGS_H 1)
set(HAVE_STRING_H 1)
set(HAVE_SYS_IOCTL_H 1)
set(HAVE_SYS_MMAN_H 1)
set(HAVE_SYS_SOCKET_H 1)
set(HAVE_NETINET_IN_H 1)
set(HAVE_NETINET_TCP_H 1)
set(HAVE_ARPA_INET_H 1)
set(HAVE_SYS_STAT_H 1)
set(HAVE_SYS_TIME_H 1)
set(HAVE_SYS_TYPES_H 1)
set(HAVE_UNISTD_H 1)
set(HAVE_GETOPT_H 1)
set(STDC_HEADERS 1)
set(SUPPORTS_VLA 1)

set( FLUID_INCLUDE_DIRECTORIES
    ${CMAKE_BINARY_DIR}
    ${MIDIPLAYER_SOURCE_DIR}/synth
    ${MIDIPLAYER_SOURCE_DIR}/rvoice
    ${MIDIPLAYER_SOURCE_DIR}/midi
    ${MIDIPLAYER_SOURCE_DIR}/utils
    ${MIDIPLAYER_SOURCE_DIR}/sfloader
    ${MIDIPLAYER_SOURCE_DIR}/bindings
    ${MIDIPLAYER_SOURCE_DIR}/../include
    ${CMAKE_BINARY_DIR}/include
    ${PTHREADS_INCLUDE_DIR}
)

#include_directories (
#    SYSTEM
#    ${GLIB_INCLUDEDIR}
#    ${GLIB_INCLUDE_DIRS}
#)

# ************ library ************

if ( DBUS_SUPPORT )
  set ( fluid_dbus_SOURCES ${MIDIPLAYER_SOURCE_DIR}/bindings/fluid_rtkit.c ${MIDIPLAYER_SOURCE_DIR}/bindings/fluid_rtkit.h )
  set ( FLUID_INCLUDE_DIRECTORIES ${FLUID_INCLUDE_DIRECTORIES} ${DBUS_INCLUDEDIR} ${DBUS_INCLUDE_DIRS} )
endif ( DBUS_SUPPORT )

if ( LASH_SUPPORT OR LADCCA_SUPPORT )
  set ( fluid_lash_SOURCES ${MIDIPLAYER_SOURCE_DIR}/bindings/fluid_lash.c ${MIDIPLAYER_SOURCE_DIR}/bindings/fluid_lash.h )
  set ( FLUID_INCLUDE_DIRECTORIES ${FLUID_INCLUDE_DIRECTORIES} ${LASH_INCLUDEDIR} ${LADCCA_INCLUDEDIR}
                        ${LASH_INCLUDE_DIRS} ${LADCCA_INCLUDE_DIRS} )
endif ( LASH_SUPPORT OR LADCCA_SUPPORT )

if ( LIBSNDFILE_SUPPORT )
    set ( FLUID_INCLUDE_DIRECTORIES ${FLUID_INCLUDE_DIRECTORIES} ${LIBSNDFILE_INCLUDEDIR} ${LIBSNDFILE_INCLUDE_DIRS} )
endif ( LIBSNDFILE_SUPPORT )

if ( LADSPA_SUPPORT )
  set ( fluid_ladspa_SOURCES ${MIDIPLAYER_SOURCE_DIR}/bindings/fluid_ladspa.c ${MIDIPLAYER_SOURCE_DIR}/bindings/fluid_ladspa.h )
endif ( LADSPA_SUPPORT )

set ( public_main_HEADER ${CMAKE_BINARY_DIR}/include/fluidsynth.h )
configure_file ( ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth.cmake ${public_main_HEADER} )

set ( config_SOURCES ${CMAKE_BINARY_DIR}/config.h )
configure_file ( ${MIDIPLAYER_SOURCE_DIR}/config.cmake ${CMAKE_BINARY_DIR}/config.h )

configure_file ( ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/version.h.in
        ${CMAKE_BINARY_DIR}/include/fluidsynth/version.h )

set( FLUID_COMPILE_OPTIONS -DHAVE_CONFIG_H -DENABLE_MIXER_THREADS)

set ( fluidsynth_SOURCES
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_conv.c
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_conv.h
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_hash.c
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_hash.h
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_list.c
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_list.h
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_ringbuffer.c
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_ringbuffer.h
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_settings.c
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_settings.h
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluidsynth_priv.h
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_sys.c
    ${MIDIPLAYER_SOURCE_DIR}/utils/fluid_sys.h
    ${MIDIPLAYER_SOURCE_DIR}/sfloader/fluid_defsfont.c
    ${MIDIPLAYER_SOURCE_DIR}/sfloader/fluid_defsfont.h
    ${MIDIPLAYER_SOURCE_DIR}/sfloader/fluid_ramsfont.c
    ${MIDIPLAYER_SOURCE_DIR}/sfloader/fluid_ramsfont.h
    ${MIDIPLAYER_SOURCE_DIR}/sfloader/fluid_sfont.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_adsr_env.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_adsr_env.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_chorus.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_chorus.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_iir_filter.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_iir_filter.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_lfo.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_lfo.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rvoice.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rvoice.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rvoice_dsp.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rvoice_event.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rvoice_event.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rvoice_mixer.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rvoice_mixer.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_phase.h
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rev.c
    ${MIDIPLAYER_SOURCE_DIR}/rvoice/fluid_rev.h
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_chan.c
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_chan.h
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_event.c
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_event_priv.h
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_gen.c
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_gen.h
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_mod.c
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_mod.h
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_synth.c
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_synth.h
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_tuning.c
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_tuning.h
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_voice.c
    ${MIDIPLAYER_SOURCE_DIR}/synth/fluid_voice.h
    ${MIDIPLAYER_SOURCE_DIR}/midi/fluid_midi.c
    ${MIDIPLAYER_SOURCE_DIR}/midi/fluid_midi.h
    ${MIDIPLAYER_SOURCE_DIR}/midi/fluid_midi_router.c
    ${MIDIPLAYER_SOURCE_DIR}/midi/fluid_midi_router.h
    ${MIDIPLAYER_SOURCE_DIR}/midi/fluid_seqbind.c
    ${MIDIPLAYER_SOURCE_DIR}/midi/fluid_seq.c
    ${MIDIPLAYER_SOURCE_DIR}/bindings/fluid_filerenderer.c
)

set ( public_HEADERS
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/audio.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/event.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/gen.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/log.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/midi.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/misc.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/mod.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/ramsfont.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/seq.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/seqbind.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/settings.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/sfont.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/shell.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/synth.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/types.h
    ${MIDIPLAYER_SOURCE_DIR}/../include/fluidsynth/voice.h
    ${CMAKE_BINARY_DIR}/include/fluidsynth/version.h
)

set( FLUID_LINK_DIRECTORIES
    ${LASH_LIBDIR}
    ${LASH_LIBRARY_DIRS}
    ${LADCCA_LIBDIR}
    ${LADCCA_LIBRARY_DIRS}
    ${LIBSNDFILE_LIBDIR}
    ${LIBSNDFILE_LIBRARY_DIRS}
    ${DBUS_LIBDIR}
    ${DBUS_LIBRARY_DIRS}
)

# ************ CLI program ************

if ( FLUID_CPPFLAGS )
  set_target_properties ( fluidsynth
    PROPERTIES COMPILE_FLAGS ${FLUID_CPPFLAGS} )
endif ( FLUID_CPPFLAGS )

set (FLUIDSYNTH_INTERFACE_INCLUDE_DIRECTORIES
        ${MIDIPLAYER_SOURCE_DIR}/../include
        ${CMAKE_BINARY_DIR}/include
        )

set (FLUID_LINK_LIBRARIES
        ${LASH_LIBRARIES}
        ${LADCCA_LIBRARIES}
        ${LIBSNDFILE_LIBRARIES}
        ${DBUS_LIBRARIES}
)

