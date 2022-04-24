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

#include "fluid_midi.h"
#include "fluid_sys.h"
#include "fluid_synth.h"
#include "fluid_settings.h"


static int fluid_midi_event_length(unsigned char event);

/* Read the entire contents of a file into memory, allocating enough memory
 * for the file, and returning the length and the buffer.
 * Note: This rewinds the file to the start before reading.
 * Returns NULL if there was an error reading or allocating memory.
 */
static char* fluid_file_read_full(fluid_file fp, size_t* length);
static void fluid_midi_event_set_sysex_LOCAL(fluid_midi_event_t *evt, int type, void *data, int size, int dynamic);
static void fluid_midi_event_get_sysex_LOCAL(fluid_midi_event_t *evt, void **data, int *size);
#define READ_FULL_INITIAL_BUFLEN 1024


/***************************************************************
 *
 *                      MIDIFILE
 */

/**
 * Return a new MIDI file handle for parsing an already-loaded MIDI file.
 * @internal
 * @param buffer Pointer to full contents of MIDI file (borrows the pointer).
 *  The caller must not free buffer until after the fluid_midi_file is deleted.
 * @param length Size of the buffer in bytes.
 * @return New MIDI file handle or NULL on error.
 */
fluid_midi_file *
new_fluid_midi_file(const char* buffer, size_t length)
{
    fluid_midi_file *mf;

    mf = FLUID_NEW(fluid_midi_file);
    if (mf == NULL) {
        FLUID_LOG(FLUID_ERR, "Out of memory");
        return NULL;
    }
    FLUID_MEMSET(mf, 0, sizeof(fluid_midi_file));

    mf->c = -1;
    mf->running_status = -1;

    mf->buffer = buffer;
    mf->buf_len = length;
    mf->buf_pos = 0;
    mf->eof = FALSE;

    if (fluid_midi_file_read_mthd(mf) != FLUID_OK) {
        FLUID_FREE(mf);
        return NULL;
    }
    return mf;
}

static char*
fluid_file_read_full(fluid_file fp, size_t* length)
{
    size_t buflen;
    char* buffer;
    size_t n;
    /* Work out the length of the file in advance */
    if (FLUID_FSEEK(fp, 0, SEEK_END) != 0)
    {
        FLUID_LOG(FLUID_ERR, "File load: Could not seek within file");
        return NULL;
    }
    buflen = ftell(fp);
    if (FLUID_FSEEK(fp, 0, SEEK_SET) != 0)
    {
        FLUID_LOG(FLUID_ERR, "File load: Could not seek within file");
        return NULL;
    }
    FLUID_LOG(FLUID_DBG, "File load: Allocating %d bytes", buflen);
    buffer = FLUID_MALLOC(buflen);
    if (buffer == NULL) {
        FLUID_LOG(FLUID_PANIC, "Out of memory");
        return NULL;
    }
    FLUID_MEMSET( buffer, 0, sizeof(char)*buflen );

    n = FLUID_FREAD(buffer, 1, buflen, fp);
    if (n != buflen) {
        FLUID_LOG(FLUID_ERR, "Only read %d bytes; expected %d", n,
                  buflen);
        FLUID_FREE(buffer);
        return NULL;
    };
    *length = n;
    return buffer;
}

/**
 * Delete a MIDI file handle.
 * @internal
 * @param mf MIDI file handle to close and free.
 */
void
delete_fluid_midi_file (fluid_midi_file *mf)
{
    if (mf == NULL) {
        return;
    }
    FLUID_FREE(mf);
    return;
}

/*
 * Gets the next byte in a MIDI file, taking into account previous running status.
 *
 * returns FLUID_FAILED if EOF or read error
 */
int
fluid_midi_file_getc (fluid_midi_file *mf)
{
    unsigned char c;
    if (mf->c >= 0) {
        c = mf->c;
        mf->c = -1;
    } else {
        if (mf->buf_pos >= mf->buf_len) {
            mf->eof = TRUE;
            return FLUID_FAILED;
        }
        c = mf->buffer[mf->buf_pos++];
        mf->trackpos++;
    }
    return (int) c;
}

/*
 * Saves a byte to be returned the next time fluid_midi_file_getc() is called,
 * when it is necessary according to running status.
 */
int
fluid_midi_file_push(fluid_midi_file *mf, int c)
{
    mf->c = c;
    return FLUID_OK;
}

/*
 * fluid_midi_file_read
 */
int
fluid_midi_file_read(fluid_midi_file *mf, void *buf, int len)
{
    int num = len < mf->buf_len - mf->buf_pos
        ? len : mf->buf_len - mf->buf_pos;
    if (num != len) {
        mf->eof = TRUE;
    }
    if (num < 0) {
        num = 0;
    }
    /* Note: Read bytes, even if there aren't enough, but only increment
     * trackpos if successful (emulates old behaviour of fluid_midi_file_read)
     */
    FLUID_MEMCPY(buf, mf->buffer+mf->buf_pos, num);
    mf->buf_pos += num;
    if (num == len)
        mf->trackpos += num;
#if DEBUG
    else
        FLUID_LOG(FLUID_DBG, "Could not read the requested number of bytes");
#endif
    return (num != len) ? FLUID_FAILED : FLUID_OK;
}

/*
 * fluid_midi_file_skip
 */
int
fluid_midi_file_skip (fluid_midi_file *mf, int skip)
{
    int new_pos = mf->buf_pos + skip;
    /* Mimic the behaviour of fseek: Error to seek past the start of file, but
     * OK to seek past end (this just puts it into the EOF state). */
    if (new_pos < 0) {
        FLUID_LOG(FLUID_ERR, "Failed to seek position in file");
        return FLUID_FAILED;
    }
    /* Clear the EOF flag, even if moved past the end of the file (this is
     * consistent with the behaviour of fseek). */
    mf->eof = FALSE;
    mf->buf_pos = new_pos;
    return FLUID_OK;
}

/*
 * fluid_midi_file_eof
 */
int fluid_midi_file_eof(fluid_midi_file* mf)
{
	/* Note: This does not simply test whether the file read pointer is past
	 * the end of the file. It mimics the behaviour of feof by actually
	 * testing the stateful EOF condition, which is set to TRUE if getc or
	 * fread have attempted to read past the end (but not if they have
	 * precisely reached the end), but reset to FALSE upon a successful seek.
	 */
	return mf->eof;
}

/*
 * fluid_midi_file_read_mthd
 */
int
fluid_midi_file_read_mthd(fluid_midi_file *mf)
{
    char mthd[14];
    if (fluid_midi_file_read(mf, mthd, sizeof(mthd)) != FLUID_OK) {
        return FLUID_FAILED;
    }
    if ((FLUID_STRNCMP(mthd, "MThd", 4) != 0) || (mthd[7] != 6)
            || (mthd[9] > 2)) {
        FLUID_LOG(FLUID_ERR,
                "Doesn't look like a MIDI file: invalid MThd header");
        return FLUID_FAILED;
    }
    mf->type = mthd[9];
    mf->ntracks = (unsigned) mthd[11];
    mf->ntracks += (unsigned int) (mthd[10]) << 16;
    if ((signed char)mthd[12] < 0) {
        mf->uses_smpte = 1;
        mf->smpte_fps = -(signed char)mthd[12];
        mf->smpte_res = (unsigned) mthd[13];
        FLUID_LOG(FLUID_ERR, "File uses SMPTE timing -- Not implemented yet");
        return FLUID_FAILED;
    } else {
        mf->uses_smpte = 0;
        mf->division = ((unsigned)mthd[12] << 8) | ((unsigned)mthd[13] & 0xff);
        FLUID_LOG(FLUID_DBG, "Division=%d", mf->division);
    }
    return FLUID_OK;
}

/*
 * fluid_midi_file_load_tracks
 */
int
fluid_midi_file_load_tracks(fluid_midi_file *mf, fluid_player_t *player)
{
    int i;
    for (i = 0; i < mf->ntracks; i++) {
        if (fluid_midi_file_read_track(mf, player, i) != FLUID_OK) {
            return FLUID_FAILED;
        }
    }
    return FLUID_OK;
}

/*
 * fluid_isasciistring
 */
int
fluid_isasciistring(char *s)
{
    int i;
    int len = (int) FLUID_STRLEN(s);
    for (i = 0; i < len; i++) {
        if (!fluid_isascii(s[i])) {
            return 0;
        }
    }
    return 1;
}

/*
 * fluid_getlength
 */
long
fluid_getlength(unsigned char *s)
{
    long i = 0;
    i = s[3] | (s[2] << 8) | (s[1] << 16) | (s[0] << 24);
    return i;
}

/*
 * fluid_midi_file_read_tracklen
 */
int
fluid_midi_file_read_tracklen(fluid_midi_file *mf)
{
    unsigned char length[5];
    if (fluid_midi_file_read(mf, length, 4) != FLUID_OK) {
        return FLUID_FAILED;
    }
    mf->tracklen = fluid_getlength(length);
    mf->trackpos = 0;
    mf->eot = 0;
    return FLUID_OK;
}

/*
 * fluid_midi_file_eot
 */
int
fluid_midi_file_eot(fluid_midi_file *mf)
{
#if DEBUG
    if (mf->trackpos > mf->tracklen) {
        printf("track overrun: %d > %d\n", mf->trackpos, mf->tracklen);
    }
#endif
    return mf->eot || (mf->trackpos >= mf->tracklen);
}

static int fluid_tempo_compare( const fluid_tempo_t* a, const fluid_tempo_t* b )
{
    if( a == NULL && b == NULL )
        return 0;

    if( a == NULL )
        return 1;

    if( b == NULL )
        return -1;

    //ticks -> track -> channel
    if( a->ticks > b->ticks )
        return 1;

    if( a->ticks < b->ticks )
        return -1;

    if( a->track > b->track )
        return 1;

    if( a->track < b->track )
        return -1;

    if( a->channel > b->channel )
        return 1;

    if( a->channel < b->channel )
        return -1;

    return 0;
}

/*
 * fluid_midi_file_read_track
 */
int
fluid_midi_file_read_track(fluid_midi_file *mf, fluid_player_t *player, int num)
{
    fluid_track_t *track;
    unsigned char id[5], length[5];
    int found_track = 0;
    int skip;
    int trackTicks = 0;
    fluid_list_t* tempo_list = NULL;
    fluid_list_t* on_note_list = NULL;
    fluid_evt_with_ticks_t* tmp_evt = NULL;
    int isError = 0;

    if (fluid_midi_file_read(mf, id, 4) != FLUID_OK) {
        return FLUID_FAILED;
    }
    id[4] = '\0';
    mf->dtime = 0;

    while (!found_track) {

        if (fluid_isasciistring((char *) id) == 0) {
            FLUID_LOG(FLUID_ERR,
                    "An non-ascii track header found, corrupt file");
            return FLUID_FAILED;

        } else if (strcmp((char *) id, "MTrk") == 0) {

            found_track = 1;

            if (fluid_midi_file_read_tracklen(mf) != FLUID_OK) {
                return FLUID_FAILED;
            }

            track = new_fluid_track(num);
            if (track == NULL) {
                FLUID_LOG(FLUID_ERR, "Out of memory");
                return FLUID_FAILED;
            }

            trackTicks = 0;
            tempo_list = NULL;
            on_note_list = NULL;

            while (!fluid_midi_file_eot(mf)) {
                if (fluid_midi_file_read_event(mf, track, &tempo_list, &trackTicks, &on_note_list) != FLUID_OK) {
                    isError = TRUE;
                    break;
                }
            }

            if( on_note_list != NULL )
            {
                FLUID_LOG( FLUID_WARN, "track %d, on_note_list is not empty. Incorrect note on/note off event", num );
                while( on_note_list != NULL )
                {
                    tmp_evt = (fluid_evt_with_ticks_t*)fluid_list_get(on_note_list);
                    FLUID_LOG( FLUID_WARN, "ticks: %u, channel: %u key: %u, vel: %u"
                        , tmp_evt->ticks, tmp_evt->evt->channel, tmp_evt->evt->param1, tmp_evt->evt->param2 );
                    on_note_list = fluid_list_remove_link( on_note_list, on_note_list );
                    free( tmp_evt );
                    tmp_evt = NULL;
                }
            }

            if( isError )
            {
                delete_fluid_track(track);
                delete_fluid_tempo_list(tempo_list);
                return FLUID_FAILED;
            }

            /* Skip remaining track data, if any */
            if (mf->trackpos < mf->tracklen) {
                if (fluid_midi_file_skip(mf, mf->tracklen - mf->trackpos) != FLUID_OK) {
                    delete_fluid_track(track);
                    delete_fluid_tempo_list(tempo_list);
                    return FLUID_FAILED;
                }
            }

            if (fluid_player_add_track(player, track) != FLUID_OK) {
                delete_fluid_track(track);
                delete_fluid_tempo_list(tempo_list);
                return FLUID_FAILED;
            }

            if( tempo_list != NULL )
            {
                if (fluid_player_add_temp_list(player, tempo_list) != FLUID_OK) {
                    delete_fluid_tempo_list(tempo_list);
                    return FLUID_FAILED;
                }
            }

            if( player->total_ticks < trackTicks )
                player->total_ticks = trackTicks;

        } else {
            found_track = 0;
            if (fluid_midi_file_read(mf, length, 4) != FLUID_OK) {
                return FLUID_FAILED;
            }
            skip = fluid_getlength(length);
            /* fseek(mf->fp, skip, SEEK_CUR); */
            if (fluid_midi_file_skip(mf, skip) != FLUID_OK) {
                return FLUID_FAILED;
            }
        }
    }
    if (fluid_midi_file_eof(mf)) {
        FLUID_LOG(FLUID_ERR, "Unexpected end of file");
        return FLUID_FAILED;
    }

    return FLUID_OK;
}

/*
 * fluid_midi_file_read_varlen
 */
int
fluid_midi_file_read_varlen(fluid_midi_file *mf)
{
    int i;
    int c;
    mf->varlen = 0;
    for (i = 0;; i++) {
        if (i == 4) {
            FLUID_LOG(FLUID_ERR, "Invalid variable length number");
            return FLUID_FAILED;
        }
        c = fluid_midi_file_getc(mf);
        if (c < 0) {
            FLUID_LOG(FLUID_ERR, "Unexpected end of file");
            return FLUID_FAILED;
        }
        if (c & 0x80) {
            mf->varlen |= (int) (c & 0x7F);
            mf->varlen <<= 7;
        } else {
            mf->varlen += c;
            break;
        }
    }
    return FLUID_OK;
}

/*
 * fluid_midi_file_read_event
 */
int
fluid_midi_file_read_event(fluid_midi_file *mf, fluid_track_t *track, fluid_list_t** tempo_list, int* trackTicks, fluid_list_t** on_note_list)
{
    int status;
    int type;

    unsigned char *metadata = NULL;
    unsigned char *dyn_buf = NULL;
    unsigned char static_buf[256];
    int nominator, denominator, clocks, notes;
    fluid_midi_event_t *evt;
    int channel = 0;
    int param1 = 0;
    int param2 = 0;
    int size;
    fluid_tempo_t* tempo = NULL;
    fluid_list_t* list = NULL;
    fluid_evt_with_ticks_t* tmp_evt = NULL;

    /* read the delta-time of the event */
    if(fluid_midi_file_read_varlen(mf) != FLUID_OK)
    {
        return FLUID_FAILED;
    }

    mf->dtime += mf->varlen;

    /* read the status byte */
    status = fluid_midi_file_getc(mf);
    (*trackTicks) += mf->varlen;

    if(status < 0)
    {
        FLUID_LOG(FLUID_ERR, "Unexpected end of file");
        return FLUID_FAILED;
    }

    /* not a valid status byte: use the running status instead */
    if((status & 0x80) == 0)
    {
        if((mf->running_status & 0x80) == 0)
        {
            FLUID_LOG(FLUID_ERR, "Undefined status and invalid running status");
            return FLUID_FAILED;
        }

        fluid_midi_file_push(mf, status);
        status = mf->running_status;
    }

    /* check what message we have */

    mf->running_status = status;

    if(status == MIDI_SYSEX)    /* system exclusif */
    {
        /* read the length of the message */
        if(fluid_midi_file_read_varlen(mf) != FLUID_OK)
        {
            return FLUID_FAILED;
        }

        if(mf->varlen)
        {
            FLUID_LOG(FLUID_DBG, "%s: %d: alloc metadata, len = %d", __FILE__,
                      __LINE__, mf->varlen);
            metadata = FLUID_MALLOC(mf->varlen + 1);

            if(metadata == NULL)
            {
                FLUID_LOG(FLUID_PANIC, "Out of memory");
                return FLUID_FAILED;
            }
            FLUID_MEMSET( metadata, 0, mf->varlen + 1 );

            /* read the data of the message */
            if(fluid_midi_file_read(mf, metadata, mf->varlen) != FLUID_OK)
            {
                FLUID_FREE(metadata);
                return FLUID_FAILED;
            }

            evt = new_fluid_midi_event();

            if(evt == NULL)
            {
                FLUID_LOG(FLUID_ERR, "Out of memory");
                FLUID_FREE(metadata);
                return FLUID_FAILED;
            }

            evt->dtime = mf->dtime;
            size = mf->varlen;

            if(metadata[mf->varlen - 1] == MIDI_EOX)
            {
                size--;
            }

            /* Add SYSEX event and indicate that its dynamically allocated and should be freed with event */
            fluid_midi_event_set_sysex(evt, metadata, size, TRUE);
            fluid_track_add_event(track, evt);
            mf->dtime = 0;
        }

        return FLUID_OK;

    }
    else if(status == MIDI_META_EVENT)      /* meta events */
    {

        int result = FLUID_OK;

        /* get the type of the meta message */
        type = fluid_midi_file_getc(mf);

        if(type < 0)
        {
            FLUID_LOG(FLUID_ERR, "Unexpected end of file");
            return FLUID_FAILED;
        }

        /* get the length of the data part */
        if(fluid_midi_file_read_varlen(mf) != FLUID_OK)
        {
            return FLUID_FAILED;
        }

        if(mf->varlen < 255)
        {
            metadata = &static_buf[0];
        }
        else
        {
            FLUID_LOG(FLUID_DBG, "%s: %d: alloc metadata, len = %d", __FILE__,
                      __LINE__, mf->varlen);
            dyn_buf = FLUID_MALLOC(mf->varlen + 1);

            if(dyn_buf == NULL)
            {
                FLUID_LOG(FLUID_PANIC, "Out of memory");
                return FLUID_FAILED;
            }
            FLUID_MEMSET( dyn_buf, 0, mf->varlen + 1 );
            metadata = dyn_buf;
        }

        /* read the data */
        if(mf->varlen)
        {
            if(fluid_midi_file_read(mf, metadata, mf->varlen) != FLUID_OK)
            {
                if(dyn_buf)
                {
                    FLUID_FREE(dyn_buf);
                }

                return FLUID_FAILED;
            }
        }

        /* handle meta data */
        switch(type)
        {

            case MIDI_COPYRIGHT:
                metadata[mf->varlen] = 0;
                break;

            case MIDI_TRACK_NAME:
                metadata[mf->varlen] = 0;
                fluid_track_set_name(track, (char *) metadata);
                break;

            case MIDI_INST_NAME:
                metadata[mf->varlen] = 0;
                break;

            case MIDI_LYRIC:
            case MIDI_TEXT:
            {
                void *tmp;
                int size = mf->varlen + 1;

                /* NULL terminate strings for safety */
                metadata[size - 1] = '\0';

                evt = new_fluid_midi_event();

                if(evt == NULL)
                {
                    FLUID_LOG(FLUID_ERR, "Out of memory");
                    result = FLUID_FAILED;
                    break;
                }

                evt->dtime = mf->dtime;

                tmp = FLUID_MALLOC(size);

                if(tmp == NULL)
                {
                    FLUID_LOG(FLUID_PANIC, "Out of memory");
                    delete_fluid_midi_event(evt);
                    evt = NULL;
                    result = FLUID_FAILED;
                    break;
                }

                FLUID_MEMCPY(tmp, metadata, size);

                fluid_midi_event_set_sysex_LOCAL(evt, type, tmp, size, TRUE);
                fluid_track_add_event(track, evt);
                mf->dtime = 0;
            }
                break;

            case MIDI_MARKER:
                break;

            case MIDI_CUE_POINT:
                break; /* don't care much for text events */

            case MIDI_EOT:
                if(mf->varlen != 0)
                {
                    FLUID_LOG(FLUID_ERR, "Invalid length for EndOfTrack event");
                    result = FLUID_FAILED;
                    break;
                }

                mf->eot = 1;
                evt = new_fluid_midi_event();

                if(evt == NULL)
                {
                    FLUID_LOG(FLUID_ERR, "Out of memory");
                    result = FLUID_FAILED;
                    break;
                }

                evt->dtime = mf->dtime;
                evt->type = MIDI_EOT;
                fluid_track_add_event(track, evt);
                mf->dtime = 0;
                break;

            case MIDI_SET_TEMPO:
                if (mf->varlen != 3) {
                    FLUID_LOG(FLUID_ERR,
                              "Invalid length for SetTempo meta event");
                    result = FLUID_FAILED;
                    break;
                }
                tempo = new_fluid_tempo();
                if (tempo == NULL) {
                    FLUID_LOG(FLUID_ERR, "Out of memory");
                    result = FLUID_FAILED;
                    break;
                }

                tempo->tempo = (metadata[0] << 16) + (metadata[1] << 8) + metadata[2];
                tempo->channel = 0;
                tempo->track = track->num;
                tempo->ticks = (*trackTicks);

                (*tempo_list) = fluid_list_append( (*tempo_list), tempo );

                evt = new_fluid_midi_event();
                if (evt == NULL) {
                    FLUID_LOG(FLUID_ERR, "Out of memory");
                    result = FLUID_FAILED;
                    break;
                }
                evt->dtime = mf->dtime;
                evt->type = MIDI_SET_TEMPO;
                evt->channel = tempo->channel;
                evt->param1 = tempo->tempo;
                evt->param2 = 0;
                evt->param3 = 0;
                fluid_track_add_event(track, evt);
                mf->dtime = 0;

                break;

            case MIDI_SMPTE_OFFSET:
                if(mf->varlen != 5)
                {
                    FLUID_LOG(FLUID_ERR,
                              "Invalid length for SMPTE Offset meta event");
                    result = FLUID_FAILED;
                    break;
                }

                break; /* we don't use smtp */

            case MIDI_TIME_SIGNATURE:
                if(mf->varlen != 4)
                {
                    FLUID_LOG(FLUID_ERR,
                              "Invalid length for TimeSignature meta event");
                    result = FLUID_FAILED;
                    break;
                }

                nominator = metadata[0];
                denominator = pow(2.0, (double) metadata[1]);
                clocks = metadata[2];
                notes = metadata[3];

                FLUID_LOG(FLUID_DBG,
                          "signature=%d/%d, metronome=%d, 32nd-notes=%d",
                          nominator, denominator, clocks, notes);

                break;

            case MIDI_KEY_SIGNATURE:
                if(mf->varlen != 2)
                {
                    FLUID_LOG(FLUID_ERR,
                              "Invalid length for KeySignature meta event");
                    result = FLUID_FAILED;
                    break;
                }

                /* We don't care about key signatures anyway */
                /* sf = metadata[0];
                mi = metadata[1]; */
                break;

            case MIDI_SEQUENCER_EVENT:
                break;

            default:
                break;
        }

        if(dyn_buf)
        {
            FLUID_LOG(FLUID_DBG, "%s: %d: free metadata", __FILE__, __LINE__);
            FLUID_FREE(dyn_buf);
        }

        return result;

    }
    else     /* channel messages */
    {

        type = status & 0xf0;
        channel = status & 0x0f;

        /* all channel message have at least 1 byte of associated data */
        if((param1 = fluid_midi_file_getc(mf)) < 0)
        {
            FLUID_LOG(FLUID_ERR, "Unexpected end of file");
            return FLUID_FAILED;
        }

        switch(type)
        {

            case NOTE_ON:
                if((param2 = fluid_midi_file_getc(mf)) < 0)
                {
                    FLUID_LOG(FLUID_ERR, "Unexpected end of file");
                    return FLUID_FAILED;
                }

                break;

            case NOTE_OFF:
                if((param2 = fluid_midi_file_getc(mf)) < 0)
                {
                    FLUID_LOG(FLUID_ERR, "Unexpected end of file");
                    return FLUID_FAILED;
                }

                break;

            case KEY_PRESSURE:
                if((param2 = fluid_midi_file_getc(mf)) < 0)
                {
                    FLUID_LOG(FLUID_ERR, "Unexpected end of file");
                    return FLUID_FAILED;
                }

                break;

            case CONTROL_CHANGE:
                if((param2 = fluid_midi_file_getc(mf)) < 0)
                {
                    FLUID_LOG(FLUID_ERR, "Unexpected end of file");
                    return FLUID_FAILED;
                }

                break;

            case PROGRAM_CHANGE:
                break;

            case CHANNEL_PRESSURE:
                break;

            case PITCH_BEND:
                if((param2 = fluid_midi_file_getc(mf)) < 0)
                {
                    FLUID_LOG(FLUID_ERR, "Unexpected end of file");
                    return FLUID_FAILED;
                }

                param1 = ((param2 & 0x7f) << 7) | (param1 & 0x7f);
                param2 = 0;
                break;

            default:
                /* Can't possibly happen !? */
                FLUID_LOG(FLUID_ERR, "Unrecognized MIDI event");
                return FLUID_FAILED;
        }

        evt = new_fluid_midi_event();

        if(evt == NULL)
        {
            FLUID_LOG(FLUID_ERR, "Out of memory");
            return FLUID_FAILED;
        }

        evt->dtime = mf->dtime;
        evt->type = type;
        evt->channel = channel;
        evt->param1 = param1;
        evt->param2 = param2;
        evt->param3 = 0;
        fluid_track_add_event(track, evt);

        if( type == NOTE_ON )
        {
            if( param2 > 0 )
            {
                tmp_evt = FLUID_NEW(fluid_evt_with_ticks_t);
                if( tmp_evt != NULL )
                {
                    FLUID_MEMSET( tmp_evt, 0, sizeof( fluid_evt_with_ticks_t ) );
                    tmp_evt->ticks = ( *trackTicks );
                    tmp_evt->evt = evt;
                    ( *on_note_list ) = fluid_list_append( ( *on_note_list ), tmp_evt );
                    tmp_evt = NULL;
                }
            }
            else
            {
                for( list = (*on_note_list); list; list = fluid_list_next(list))
                {
                    tmp_evt = (fluid_evt_with_ticks_t*)fluid_list_get(list);
                    if( tmp_evt->evt->channel == channel && tmp_evt->evt->param1 == param1 )
                    {
                        tmp_evt->evt->param3 = (*trackTicks)-tmp_evt->ticks;
                        (*on_note_list) = fluid_list_remove_link( (*on_note_list), list );
                        free( tmp_evt );
                        tmp_evt = NULL;
                        break;
                    }
                    tmp_evt = NULL;
                }
            }
        }
        else if( evt->type == NOTE_OFF )
        {
            for( list = (*on_note_list); list; list = fluid_list_next(list))
            {
                tmp_evt = (fluid_evt_with_ticks_t*)fluid_list_get(list);
                if( tmp_evt->evt->channel == channel && tmp_evt->evt->param1 == param1 )
                {
                    tmp_evt->evt->param3 = (*trackTicks)-tmp_evt->ticks;
                    (*on_note_list) = fluid_list_remove_link( (*on_note_list), list );
                    free( tmp_evt );
                    tmp_evt = NULL;
                    break;
                }
                tmp_evt = NULL;
            }
        }

        mf->dtime = 0;
    }

    return FLUID_OK;
}

/*
 * fluid_midi_file_get_division
 */
int
fluid_midi_file_get_division(fluid_midi_file *midifile)
{
    return midifile->division;
}

/******************************************************
 *
 *     fluid_track_t
 */

/**
 * Create a MIDI event structure.
 * @return New MIDI event structure or NULL when out of memory.
 */
fluid_midi_event_t *
new_fluid_midi_event ()
{
    fluid_midi_event_t* evt;
    evt = FLUID_NEW(fluid_midi_event_t);
    if (evt == NULL) {
        FLUID_LOG(FLUID_ERR, "Out of memory");
        return NULL;
    }
    FLUID_MEMSET(evt, 0, sizeof(fluid_midi_event_t));

    evt->dtime = 0;
    evt->type = 0;
    evt->channel = 0;
    evt->param1 = 0;
    evt->param2 = 0;
    evt->param3 = 0;
    evt->next = NULL;
    evt->paramptr = NULL;
    return evt;
}

/**
 * Delete MIDI event structure.
 * @param evt MIDI event structure
 * @return Always returns #FLUID_OK
 */
int
delete_fluid_midi_event(fluid_midi_event_t *evt)
{
    fluid_midi_event_t *temp;

    while (evt) {
        temp = evt->next;

        /* Dynamic SYSEX event? - free (param2 indicates if dynamic) */
        if((evt->type == MIDI_SYSEX || (evt-> type == MIDI_TEXT) || (evt->type == MIDI_LYRIC)) &&
           evt->paramptr && evt->param2)
        {
            FLUID_FREE(evt->paramptr);
        }

        FLUID_FREE(evt);
        evt = temp;
    }
    return FLUID_OK;
}

/**
 * Get the event type field of a MIDI event structure.
 * @param evt MIDI event structure
 * @return Event type field (MIDI status byte without channel)
 */
int
fluid_midi_event_get_type(fluid_midi_event_t *evt)
{
    return evt->type;
}

/**
 * Set the event type field of a MIDI event structure.
 * @param evt MIDI event structure
 * @param type Event type field (MIDI status byte without channel)
 * @return Always returns #FLUID_OK
 */
int
fluid_midi_event_set_type(fluid_midi_event_t *evt, int type)
{
    evt->type = type;
    return FLUID_OK;
}

/**
 * Get the channel field of a MIDI event structure.
 * @param evt MIDI event structure
 * @return Channel field
 */
int
fluid_midi_event_get_channel(fluid_midi_event_t *evt)
{
    return evt->channel;
}

/**
 * Set the channel field of a MIDI event structure.
 * @param evt MIDI event structure
 * @param chan MIDI channel field
 * @return Always returns #FLUID_OK
 */
int
fluid_midi_event_set_channel(fluid_midi_event_t *evt, int chan)
{
    evt->channel = chan;
    return FLUID_OK;
}

/**
 * Get the key field of a MIDI event structure.
 * @param evt MIDI event structure
 * @return MIDI note number (0-127)
 */
int
fluid_midi_event_get_key(fluid_midi_event_t *evt)
{
    return evt->param1;
}

/**
 * Set the key field of a MIDI event structure.
 * @param evt MIDI event structure
 * @param v MIDI note number (0-127)
 * @return Always returns #FLUID_OK
 */
int
fluid_midi_event_set_key(fluid_midi_event_t *evt, int v)
{
    evt->param1 = v;
    return FLUID_OK;
}

/**
 * Get the velocity field of a MIDI event structure.
 * @param evt MIDI event structure
 * @return MIDI velocity number (0-127)
 */
int
fluid_midi_event_get_velocity(fluid_midi_event_t *evt)
{
    return evt->param2;
}

/**
 * Set the velocity field of a MIDI event structure.
 * @param evt MIDI event structure
 * @param v MIDI velocity value
 * @return Always returns #FLUID_OK
 */
int
fluid_midi_event_set_velocity(fluid_midi_event_t *evt, int v)
{
    evt->param2 = v;
    return FLUID_OK;
}

int
fluid_midi_event_get_duration(fluid_midi_event_t *evt)
{
    return evt->param3;
}

int
fluid_midi_event_set_duration(fluid_midi_event_t *evt, int d)
{
    evt->param3 = d;
    return FLUID_OK;
}

/**
 * Get the control number of a MIDI event structure.
 * @param evt MIDI event structure
 * @return MIDI control number
 */
int
fluid_midi_event_get_control(fluid_midi_event_t *evt)
{
    return evt->param1;
}

/**
 * Set the control field of a MIDI event structure.
 * @param evt MIDI event structure
 * @param v MIDI control number
 * @return Always returns #FLUID_OK
 */
int
fluid_midi_event_set_control(fluid_midi_event_t *evt, int v)
{
    evt->param1 = v;
    return FLUID_OK;
}

/**
 * Get the value field from a MIDI event structure.
 * @param evt MIDI event structure
 * @return Value field
 */
int
fluid_midi_event_get_value(fluid_midi_event_t *evt)
{
    return evt->param2;
}

/**
 * Set the value field of a MIDI event structure.
 * @param evt MIDI event structure
 * @param v Value to assign
 * @return Always returns #FLUID_OK
 */
int
fluid_midi_event_set_value(fluid_midi_event_t *evt, int v)
{
    evt->param2 = v;
    return FLUID_OK;
}

/**
 * Get the program field of a MIDI event structure.
 * @param evt MIDI event structure
 * @return MIDI program number (0-127)
 */
int
fluid_midi_event_get_program(fluid_midi_event_t *evt)
{
    return evt->param1;
}

/**
 * Set the program field of a MIDI event structure.
 * @param evt MIDI event structure
 * @param val MIDI program number (0-127)
 * @return Always returns #FLUID_OK
 */
int
fluid_midi_event_set_program(fluid_midi_event_t *evt, int val)
{
    evt->param1 = val;
    return FLUID_OK;
}

/**
 * Get the pitch field of a MIDI event structure.
 * @param evt MIDI event structure
 * @return Pitch value (14 bit value, 0-16383, 8192 is center)
 */
int
fluid_midi_event_get_pitch(fluid_midi_event_t *evt)
{
    return evt->param1;
}

/**
 * Set the pitch field of a MIDI event structure.
 * @param evt MIDI event structure
 * @param val Pitch value (14 bit value, 0-16383, 8192 is center)
 * @return Always returns FLUID_OK
 */
int
fluid_midi_event_set_pitch(fluid_midi_event_t *evt, int val)
{
    evt->param1 = val;
    return FLUID_OK;
}

/**
 * Assign sysex data to a MIDI event structure.
 * @param evt MIDI event structure
 * @param data Pointer to SYSEX data
 * @param size Size of SYSEX data in bytes
 * @param dynamic TRUE if the SYSEX data has been dynamically allocated and
 *   should be freed when the event is freed (only applies if event gets destroyed
 *   with delete_fluid_midi_event())
 * @return Always returns #FLUID_OK
 */
int
fluid_midi_event_set_sysex(fluid_midi_event_t *evt, void *data, int size, int dynamic)
{
    fluid_midi_event_set_sysex_LOCAL(evt, MIDI_SYSEX, data, size, dynamic);
    return FLUID_OK;
}

/**
 * Assign text data to a MIDI event structure.
 * @param evt MIDI event structure
 * @param data Pointer to text data
 * @param size Size of text data in bytes
 * @param dynamic TRUE if the data has been dynamically allocated and
 *   should be freed when the event is freed via delete_fluid_midi_event()
 * @return Always returns #FLUID_OK
 *
 * @since 2.0.0
 */
int
fluid_midi_event_set_text(fluid_midi_event_t *evt, void *data, int size, int dynamic)
{
    fluid_midi_event_set_sysex_LOCAL(evt, MIDI_TEXT, data, size, dynamic);
    return FLUID_OK;
}

/**
 * Get the text of a MIDI event structure.
 * @param evt MIDI event structure
 * @param data Pointer to return text data on.
 * @param size Pointer to return text size on.
 * @return Returns #FLUID_OK if \p data and \p size previously set by
 * fluid_midi_event_set_text() have been successfully retrieved.
 * Else #FLUID_FAILED is returned and \p data and \p size are not changed.
 * @since 2.0.3
 */
int fluid_midi_event_get_text(fluid_midi_event_t *evt, void **data, int *size)
{
    fluid_return_val_if_fail(evt != NULL, FLUID_FAILED);
    fluid_return_val_if_fail(evt->type == MIDI_TEXT, FLUID_FAILED);

    fluid_midi_event_get_sysex_LOCAL(evt, data, size);
    return FLUID_OK;
}

/**
 * Assign lyric data to a MIDI event structure.
 * @param evt MIDI event structure
 * @param data Pointer to lyric data
 * @param size Size of lyric data in bytes
 * @param dynamic TRUE if the data has been dynamically allocated and
 *   should be freed when the event is freed via delete_fluid_midi_event()
 * @return Always returns #FLUID_OK
 *
 * @since 2.0.0
 */
int
fluid_midi_event_set_lyrics(fluid_midi_event_t *evt, void *data, int size, int dynamic)
{
    fluid_midi_event_set_sysex_LOCAL(evt, MIDI_LYRIC, data, size, dynamic);
    return FLUID_OK;
}

/**
 * Get the lyric of a MIDI event structure.
 * @param evt MIDI event structure
 * @param data Pointer to return lyric data on.
 * @param size Pointer to return lyric size on.
 * @return Returns #FLUID_OK if \p data and \p size previously set by
 * fluid_midi_event_set_lyrics() have been successfully retrieved.
 * Else #FLUID_FAILED is returned and \p data and \p size are not changed.
 * @since 2.0.3
 */
int fluid_midi_event_get_lyrics(fluid_midi_event_t *evt, void **data, int *size)
{
    fluid_return_val_if_fail(evt != NULL, FLUID_FAILED);
    fluid_return_val_if_fail(evt->type == MIDI_LYRIC, FLUID_FAILED);

    fluid_midi_event_get_sysex_LOCAL(evt, data, size);
    return FLUID_OK;
}

static void fluid_midi_event_set_sysex_LOCAL(fluid_midi_event_t *evt, int type, void *data, int size, int dynamic)
{
    evt->type = type;
    evt->paramptr = data;
    evt->param1 = size;
    evt->param2 = dynamic;
}

static void fluid_midi_event_get_sysex_LOCAL(fluid_midi_event_t *evt, void **data, int *size)
{
    if(data)
    {
        *data = evt->paramptr;
    }

    if(size)
    {
        *size = evt->param1;
    }
}

/******************************************************
 *
 *     fluid_track_t
 */

/*
 * new_fluid_track
 */
fluid_track_t *
new_fluid_track(int num)
{
    fluid_track_t *track;
    track = FLUID_NEW(fluid_track_t);
    if (track == NULL) {
        return NULL;
    }
    FLUID_MEMSET(track, 0, sizeof(fluid_track_t));

    track->name = NULL;
    track->num = num;
    track->first = NULL;
    track->cur = NULL;
    track->last = NULL;
    track->ticks = 0;
    return track;
}

/*
 * delete_fluid_track
 */
int
delete_fluid_track(fluid_track_t *track)
{
    if (track->name != NULL) {
        FLUID_FREE(track->name);
    }
    if (track->first != NULL) {
        delete_fluid_midi_event(track->first);
    }
    FLUID_FREE(track);
    return FLUID_OK;
}

/*
 * fluid_track_set_name
 */
int
fluid_track_set_name(fluid_track_t *track, char *name)
{
    int len;
    if (track->name != NULL) {
        FLUID_FREE(track->name);
    }
    if (name == NULL) {
        track->name = NULL;
        return FLUID_OK;
    }
    len = FLUID_STRLEN(name);
    track->name = FLUID_MALLOC(len + 1);
    if (track->name == NULL) {
        FLUID_LOG(FLUID_ERR, "Out of memory");
        return FLUID_FAILED;
    }
    FLUID_MEMSET( track->name, 0, sizeof(char)*(len+1) );
    FLUID_STRCPY(track->name, name);
    return FLUID_OK;
}

/*
 * fluid_track_get_name
 */
char *
fluid_track_get_name(fluid_track_t *track)
{
    return track->name;
}

/*
 * fluid_track_get_duration
 */
int
fluid_track_get_duration(fluid_track_t *track)
{
    int time = 0;
    fluid_midi_event_t *evt = track->first;
    while (evt != NULL) {
        time += evt->dtime;
        evt = evt->next;
    }
    return time;
}

/*
 * fluid_track_count_events
 */
int
fluid_track_count_events(fluid_track_t *track, int *on, int *off)
{
    fluid_midi_event_t *evt = track->first;
    while (evt != NULL) {
        if (evt->type == NOTE_ON) {
            (*on)++;
        } else if (evt->type == NOTE_OFF) {
            (*off)++;
        }
        evt = evt->next;
    }
    return FLUID_OK;
}

/*
 * fluid_track_add_event
 */
int
fluid_track_add_event(fluid_track_t *track, fluid_midi_event_t *evt)
{
    evt->next = NULL;
    if (track->first == NULL) {
        track->first = evt;
        track->cur = evt;
        track->last = evt;
    } else {
        track->last->next = evt;
        track->last = evt;
    }
    return FLUID_OK;
}

/*
 * fluid_track_first_event
 */
fluid_midi_event_t *
fluid_track_first_event(fluid_track_t *track)
{
    track->cur = track->first;
    return track->cur;
}

/*
 * fluid_track_next_event
 */
fluid_midi_event_t *
fluid_track_next_event(fluid_track_t *track)
{
    if (track->cur != NULL) {
        track->cur = track->cur->next;
    }
    return track->cur;
}

/*
 * fluid_track_reset
 */
int
fluid_track_reset(fluid_track_t *track)
{
    track->ticks = 0;
    track->cur = track->first;
    return FLUID_OK;
}

/*
 * fluid_track_send_events
 */
int
fluid_track_send_events(fluid_track_t *track,
			fluid_synth_t *synth,
			fluid_player_t *player,
			int is_seeking,
			unsigned int ticks)
{
    int status = FLUID_OK;
    fluid_midi_event_t *event;
    int seek_forward = FALSE;

    if(is_seeking)
    {
        if(track->ticks > ticks)
        {
            fluid_track_reset(track);    /* reset track if seeking backwards */
            seek_forward = TRUE;
        }
    }

    while (1) {

        event = track->cur;
        if (event == NULL) {
            return status;
        }

        /* 		printf("track=%02d\tticks=%05u\ttrack=%05u\tdtime=%05u\tnext=%05u\n", */
        /* 		       track->num, */
        /* 		       ticks, */
        /* 		       track->ticks, */
        /* 		       event->dtime, */
        /* 		       track->ticks + event->dtime); */

        if (track->ticks + event->dtime > ticks) {
            return status;
        }

        track->ticks += event->dtime;

        if (!player || event->type == MIDI_EOT) {
        }
        else if(is_seeking && (event->type == NOTE_ON || event->type == NOTE_OFF)) {
            if( seek_forward )
            {
                /* skip on/off messages */
                if( track->ticks+event->param3 > ticks )
                {
                    if (player->seeking_callback)
                        player->seeking_callback(player->playback_userdata, event, ticks-track->ticks);
                }
            }
            else
            {
                if( event->type == NOTE_ON && event->param3 > 0 )
                {
                    if( track->ticks+event->param3 > ticks )
                    {
                        if (player->seeking_callback)
                            player->seeking_callback(player->playback_userdata, event, ticks-track->ticks);
                    }
                }
                else if( event->type == NOTE_ON && event->param3 == 0 || event->type == NOTE_OFF )
                {
                    if (player->playback_callback)
                        player->playback_callback(player->playback_userdata, event);
                }
            }
        }
        else {
            if (player->playback_callback)
                player->playback_callback(player->playback_userdata, event);
        }

        fluid_track_next_event(track);

    }
    return status;
}

fluid_tempo_t* new_fluid_tempo()
{
    fluid_tempo_t* tempo;
    tempo = FLUID_NEW(fluid_tempo_t);
    if (tempo == NULL) {
        return NULL;
    }
    FLUID_MEMSET(tempo, 0, sizeof(fluid_tempo_t));

    tempo->track = -1;
    tempo->ticks = 0;
    tempo->channel = 0;
    tempo->tempo = 500000;

    return tempo;
}

void delete_fluid_tempo( fluid_tempo_t* tempo )
{
    if( tempo != NULL )
        free( tempo );
}

void delete_fluid_tempo_list( fluid_list_t* tempo_list )
{
    fluid_list_t* list = tempo_list;
    fluid_list_t* list_tmp = NULL;

    while( list != NULL )
    {
        list_tmp = list;
        list = list->next;

        delete_fluid_tempo( (fluid_tempo_t*)list_tmp->data );
        delete1_fluid_list( list_tmp );
        list_tmp = NULL;
    }
}

/******************************************************
 *
 *     fluid_player
 */
static void
fluid_player_handle_reset_synth(void *data, int value)
{
    fluid_player_t *player = data;
    fluid_return_if_fail(player != NULL);

    player->reset_synth_between_songs = value;
}

/**
 * Create a new MIDI player.
 * @param synth Fluid synthesizer instance to create player for
 * @return New MIDI player instance or NULL on error (out of memory)
 */
fluid_player_t *
new_fluid_player(fluid_synth_t *synth)
{
    int i, play_next_auto = TRUE;
    fluid_player_t *player;
    player = FLUID_NEW(fluid_player_t);
    if (player == NULL) {
        FLUID_LOG(FLUID_ERR, "Out of memory");
        return NULL;
    }
    FLUID_MEMSET(player, 0, sizeof(fluid_player_t));

    player->status = FLUID_PLAYER_STOPPED;
    player->loop = 1;
    player->ntracks = 0;
    for (i = 0; i < MAX_NUMBER_OF_TRACKS; i++) {
        player->track[i] = NULL;
    }
    player->synth = synth;
    player->system_timer = NULL;
    player->sample_timer = NULL;
    player->playlist = NULL;
    player->currentfile = NULL;
    player->division = 0;
    player->send_program_change = 1;
    player->miditempo = 500000;
    player->playtempo = 500000;
    player->deltatime = 4.0;
    player->midi_duration_msec = 0.0f;
    player->play_duration_msec = 0.0f;
    player->duration_ticks = 0;
    player->tick_duration_msec = 0.0f;
    player->last_sync_msec = 0;
    player->playdeltatime = 4.0;
    fluid_settings_getint(synth->settings, "player.play-next-auto", &play_next_auto);
    player->play_next_auto = play_next_auto;
    fluid_player_set_playback_callback(player, fluid_synth_handle_midi_event, fluid_synth_handle_midi_seek_event, synth);
    player->tempo = NULL;
    player->cur_tempo = NULL;
    player->tempo_offset_bpm = 0;
    player->tempo_offset_timeCoeff = 1.0f;
    player->total_ticks = 0;
    player->total_msec = 0;

    player->is_seeking = FALSE;
    fluid_mutex_init(player->player_mutex);

    player->use_system_timer = fluid_settings_str_equal(synth->settings,
            "player.timing-source", "system");

    fluid_settings_getint(synth->settings, "player.reset-synth", &i);
    fluid_player_handle_reset_synth(player, i);

    return player;
}

/**
 * Delete a MIDI player instance.
 * @param player MIDI player instance
 * @return Always returns #FLUID_OK
 */
int
delete_fluid_player(fluid_player_t *player)
{
    fluid_list_t *q;
    fluid_playlist_item* pi;

    if (player == NULL) {
        return FLUID_OK;
    }
    fluid_player_stop(player);
    fluid_player_reset(player);

    while (player->playlist != NULL) {
        q = player->playlist->next;
        pi = (fluid_playlist_item*) player->playlist->data;
        FLUID_FREE(pi->filename);
        FLUID_FREE(pi->buffer);
        FLUID_FREE(pi);
        delete1_fluid_list(player->playlist);
        player->playlist = q;
    }

    fluid_mutex_destroy(player->player_mutex);

    FLUID_FREE(player);
    return FLUID_OK;
}

/**
 * Registers settings related to the MIDI player
 */
void
fluid_player_settings(fluid_settings_t *settings)
{
    /* player.timing-source can be either "system" (use system timer)
     or "sample" (use timer based on number of written samples) */
    fluid_settings_register_str(settings, "player.timing-source", "sample", 0,
            NULL, NULL);
    fluid_settings_add_option(settings, "player.timing-source", "sample");
    fluid_settings_add_option(settings, "player.timing-source", "system");

    /* Selects whether the player should reset the synth between songs, or not. */
    fluid_settings_register_int(settings, "player.reset-synth", 1, 0, 1,
            FLUID_HINT_TOGGLED, NULL, NULL);
    fluid_settings_register_int(settings, "player.play-next-auto", 1, 0, 1,
                                FLUID_HINT_TOGGLED, NULL, NULL);
}


int
fluid_player_reset(fluid_player_t *player)
{
    int i;

    for (i = 0; i < MAX_NUMBER_OF_TRACKS; i++) {
        if (player->track[i] != NULL) {
            delete_fluid_track(player->track[i]);
            player->track[i] = NULL;
        }
    }
    /*	player->current_file = NULL; */
    /*	player->status = FLUID_PLAYER_READY; */
    /*	player->loop = 1; */
    player->ntracks = 0;
    player->division = 0;
    player->send_program_change = 1;
    player->miditempo = 500000;
    player->playtempo = 500000;
    player->deltatime = 4.0;
    player->playdeltatime = 4.0;
    player->cur_tempo = NULL;
    delete_fluid_tempo_list( player->tempo );
    player->tempo = NULL;
    player->tempo_offset_bpm = 0;
    player->tempo_offset_timeCoeff = 1.0f;
    player->total_ticks = 0;
    player->total_msec = 0;
    player->midi_duration_msec = 0.0f;
    player->play_duration_msec = 0.0f;
    player->duration_ticks = 0;
    player->tick_duration_msec = 0.0f;
    player->last_sync_msec = 0;

    return 0;
}

/*
 * fluid_player_add_track
 */
int
fluid_player_add_track(fluid_player_t *player, fluid_track_t *track)
{
    if (player->ntracks < MAX_NUMBER_OF_TRACKS) {
        player->track[player->ntracks++] = track;
        return FLUID_OK;
    } else {
        return FLUID_FAILED;
    }
}

int
fluid_player_add_temp_list(fluid_player_t *player, fluid_list_t *tempo_list)
{
    fluid_list_t* last = NULL;

    if (player->tempo != NULL) {
        last = fluid_list_last( player->tempo );
        last->next = tempo_list;
    } else {
        player->tempo = tempo_list;
    }

    return FLUID_OK;
}

/*
 * fluid_player_count_tracks
 */
int
fluid_player_count_tracks(fluid_player_t *player)
{
    return player->ntracks;
}

/*
 * fluid_player_get_track
 */
fluid_track_t *
fluid_player_get_track(fluid_player_t *player, int i)
{
    if ((i >= 0) && (i < MAX_NUMBER_OF_TRACKS)) {
        return player->track[i];
    } else {
        return NULL;
    }
}

/**
 * Change the MIDI callback function. This is usually set to 
 * fluid_synth_handle_midi_event, but can optionally be changed
 * to a user-defined function instead, for intercepting all MIDI
 * messages sent to the synth. You can also use a midi router as 
 * the callback function to modify the MIDI messages before sending
 * them to the synth. 
 * @param player MIDI player instance
 * @param handler Pointer to callback function
 * @param handler_data Parameter sent to the callback function
 * @returns FLUID_OK
 * @since 1.1.4
 */
int 
fluid_player_set_playback_callback(fluid_player_t* player, 
    handle_midi_event_func_t handler, handle_midi_seek_event_func_t seek_handler, void* handler_data)
{
    player->playback_callback = handler;
    player->seeking_callback = seek_handler;
    player->playback_userdata = handler_data;
    return FLUID_OK;
}

/**
 * Add a MIDI file to a player queue.
 * @param player MIDI player instance
 * @param midifile File name of the MIDI file to add
 * @return #FLUID_OK or #FLUID_FAILED
 */
int
fluid_player_add(fluid_player_t *player, const char *midifile)
{
    fluid_playlist_item *pi = FLUID_MALLOC(sizeof(fluid_playlist_item));
    char* f = FLUID_STRDUP(midifile);
    if (!pi || !f) {
        FLUID_FREE(pi);
        FLUID_FREE(f);
        FLUID_LOG(FLUID_PANIC, "Out of memory");
        return FLUID_FAILED;
    }
    FLUID_MEMSET( pi, 0, sizeof(fluid_playlist_item) );

    pi->filename = f;
    pi->buffer = NULL;
    pi->buffer_len = 0;
    player->playlist = fluid_list_append(player->playlist, pi);
    player->status = FLUID_PLAYER_READY;
    return FLUID_OK;
}

/**
 * Add a MIDI file to a player queue, from a buffer in memory.
 * @param player MIDI player instance
 * @param buffer Pointer to memory containing the bytes of a complete MIDI
 *   file. The data is copied, so the caller may free or modify it immediately
 *   without affecting the playlist.
 * @param len Length of the buffer, in bytes.
 * @return #FLUID_OK or #FLUID_FAILED
 */
int
fluid_player_add_mem(fluid_player_t* player, const void *buffer, size_t len)
{
    /* Take a copy of the buffer, so the caller can free immediately. */
    fluid_playlist_item *pi = FLUID_MALLOC(sizeof(fluid_playlist_item));
    void *buf_copy = FLUID_MALLOC(len);
    if (!pi || !buf_copy) {
        if( pi != NULL ) {
            FLUID_FREE(pi);
            pi = NULL;
        }
        if( buf_copy != NULL ) {
            FLUID_FREE(buf_copy);
            buf_copy = NULL;
        }
        FLUID_LOG(FLUID_PANIC, "Out of memory");
        return FLUID_FAILED;
    }

    FLUID_MEMSET( pi, 0, sizeof(fluid_playlist_item) );
    FLUID_MEMSET( buf_copy, 0, len);

    FLUID_MEMCPY(buf_copy, buffer, len);
    pi->filename = NULL;
    pi->buffer = buf_copy;
    pi->buffer_len = len;
    player->playlist = fluid_list_append(player->playlist, pi);
    player->status = FLUID_PLAYER_READY;
    return FLUID_OK;
}

/*
 * fluid_player_load
 */
int
fluid_player_load(fluid_player_t *player, fluid_playlist_item *item)
{
    fluid_midi_file *midifile;
    char* buffer;
    size_t buffer_length;
    int buffer_owned;

    if (item->filename != NULL)
    {
        fluid_file fp;
        /* This file is specified by filename; load the file from disk */
        FLUID_LOG(FLUID_DBG, "%s: %d: Loading midifile %s", __FILE__, __LINE__,
                item->filename);
        /* Read the entire contents of the file into the buffer */
        fp = FLUID_FOPEN(item->filename, "rb");
        if (fp == NULL) {
            FLUID_LOG(FLUID_ERR, "Couldn't open the MIDI file");
            return FLUID_FAILED;
        }
        buffer = fluid_file_read_full(fp, &buffer_length);
        FLUID_FCLOSE(fp);
        if (buffer == NULL)
        {
            FLUID_FCLOSE(fp);
            return FLUID_FAILED;
        }
        buffer_owned = 1;
    }
    else
    {
        /* This file is specified by a pre-loaded buffer; load from memory */
        FLUID_LOG(FLUID_DBG, "%s: %d: Loading midifile from memory (%p)",
                __FILE__, __LINE__, item->buffer);
        buffer = (char *) item->buffer;
        buffer_length = item->buffer_len;
        /* Do not free the buffer (it is owned by the playlist) */
        buffer_owned = 0;
    }

    midifile = new_fluid_midi_file(buffer, buffer_length);
    if (midifile == NULL) {
        if (buffer_owned) {
            if( buffer != NULL ) {
                FLUID_FREE(buffer);
                buffer = NULL;
            }
        }
        return FLUID_FAILED;
    }
    player->division = fluid_midi_file_get_division(midifile);
    fluid_player_set_midi_tempo(player, player->miditempo); // Update deltatime
    /*FLUID_LOG(FLUID_DBG, "quarter note division=%d\n", player->division); */

    if (fluid_midi_file_load_tracks(midifile, player) != FLUID_OK) {
        if (buffer_owned) {
            if( buffer != NULL ) {
                FLUID_FREE(buffer);
                buffer = NULL;
            }
        }
        delete_fluid_midi_file(midifile);
        return FLUID_FAILED;
    }
    player->total_msec = fluid_player_convert_ticks_to_msec(player, player->total_ticks);
    if( player->tempo != NULL )
    {
        player->tempo = fluid_list_sort( player->tempo, fluid_tempo_compare );
        player->cur_tempo = player->tempo;
    }

    delete_fluid_midi_file(midifile);
    if (buffer_owned) {
        if( buffer != NULL ) {
            FLUID_FREE(buffer);
            buffer = NULL;
        }
    }
    return FLUID_OK;
}

void
fluid_player_advancefile(fluid_player_t *player)
{
    if (player->playlist == NULL) {
        return; /* No files to play */
    }
    if (player->currentfile != NULL) {
        player->currentfile = fluid_list_next(player->currentfile);
    }
    if (player->currentfile == NULL) {
        if (player->loop == 0) {
            return; /* We're done playing */
        }
        if (player->loop > 0) {
            player->loop--;
        }
        player->currentfile = player->playlist;
    }
}

void
fluid_player_playlist_load(fluid_player_t *player, unsigned int msec)
{
    fluid_playlist_item* current_playitem;
    int i;

    do {
        fluid_player_advancefile(player);
        if (player->currentfile == NULL) {
            /* Failed to find next song, probably since we're finished */
            player->status = FLUID_PLAYER_DONE;
            return;
        }

        fluid_player_reset(player);
        current_playitem = (fluid_playlist_item *) player->currentfile->data;
    } while (fluid_player_load(player, current_playitem) != FLUID_OK);

    /* Successfully loaded midi file */
    player->status = FLUID_PLAYER_READY;

    if( player->tempo != NULL )
        player->cur_tempo = player->tempo;
    player->miditempo = 500000;
    player->playtempo = rint(60000000.0/((60000000.0/(double)player->miditempo)+(double)player->tempo_offset_bpm));
    player->deltatime = (double)player->miditempo / (double)player->division / 1000.0; /* in milliseconds */
    player->playdeltatime = ((double)player->playtempo / (double)player->division / 1000.0)*player->tempo_offset_timeCoeff; /* in milliseconds */
    player->midi_duration_msec = 0.0f;
    player->play_duration_msec = 0.0f;
    player->duration_ticks = 0;
    player->tick_duration_msec = 0.0f;
    player->last_sync_msec = msec;

    if (player->reset_synth_between_songs) {
        fluid_synth_system_reset(player->synth);
    }

    for (i = 0; i < player->ntracks; i++) {
        if (player->track[i] != NULL) {
            fluid_track_reset(player->track[i]);
        }
    }
}

int fluid_player_update_duration( fluid_player_t* player, int msec )
{
    int is_updated = FALSE;

    int pre_duration_msec = 0;
    int next_tempo_delta = 0;
    int remain_ticks = 0;
    fluid_real_t remain_delta_msec = 0.0f;
    fluid_tempo_t* tempo = NULL;
    fluid_real_t real_remain_ticks = 0.0f;

    pre_duration_msec = player->midi_duration_msec;

    remain_delta_msec = (fluid_real_t)(msec-player->last_sync_msec);
    if( player->midi_duration_msec - player->tick_duration_msec > 0 )
        remain_delta_msec += (((player->midi_duration_msec-player->tick_duration_msec)/player->deltatime)*player->playdeltatime);

    do
    {
        while( player->cur_tempo != NULL )
        {
            tempo = (fluid_tempo_t*)player->cur_tempo->data;
            if( tempo == NULL )
            {
                player->cur_tempo = fluid_list_next(player->cur_tempo);
                continue;
            }

            if( tempo->ticks > player->duration_ticks )
            {
                next_tempo_delta = (tempo->ticks-player->duration_ticks);
                break;
            }

            fluid_player_set_midi_tempo( player, tempo->tempo );

            player->cur_tempo = fluid_list_next(player->cur_tempo);
        }

        if( remain_delta_msec > 0 )
            remain_ticks = (int)floor(remain_delta_msec/player->playdeltatime);

        if( remain_ticks <= 0 )
            break;

        if( next_tempo_delta > 0 && remain_ticks > next_tempo_delta )
        {
            player->tick_duration_msec += (player->deltatime*(fluid_real_t)next_tempo_delta);
            player->duration_ticks += next_tempo_delta;
            remain_delta_msec -= (player->playdeltatime*(fluid_real_t)next_tempo_delta);
        }
        else
        {
            player->tick_duration_msec += (player->deltatime*(fluid_real_t)remain_ticks);
            player->duration_ticks += remain_ticks;
            remain_delta_msec -= (player->playdeltatime*(fluid_real_t)remain_ticks);
        }

        next_tempo_delta = 0;
        remain_ticks = 0;
    }while( remain_delta_msec > 0 );

    if( remain_delta_msec > 0 )
    {
        real_remain_ticks = remain_delta_msec/player->playdeltatime;
        player->midi_duration_msec = player->tick_duration_msec+(real_remain_ticks*player->deltatime);
    }

    player->play_duration_msec += (msec-player->last_sync_msec);
    player->last_sync_msec = msec;

    if( pre_duration_msec < (int)player->midi_duration_msec )
        is_updated = true;

    return is_updated;
}

/*
 * fluid_player_callback
 */
int
fluid_player_callback(void *data, unsigned int msec)
{
    int i;
    int loadnextfile;
    int status = FLUID_PLAYER_DONE;
    fluid_player_t *player;
    fluid_synth_t *synth;

    player = (fluid_player_t *) data;
    synth = player->synth;

    if( player->is_seeking )
        return 1;

    fluid_mutex_lock( player->player_mutex );

    loadnextfile = player->currentfile == NULL ? 1 : 0;

    if(player->status == FLUID_PLAYER_STOPPED || player->status == FLUID_PLAYER_DONE)
    {
        fluid_synth_all_notes_off(synth, -1);
        fluid_mutex_unlock( player->player_mutex );
        return 1;
    }

    do {
        if (loadnextfile) {
            loadnextfile = 0;
            fluid_player_playlist_load(player, msec);
            if (player->currentfile == NULL) {
                fluid_mutex_unlock( player->player_mutex );
                return 0;
            }
        }

        fluid_player_update_duration( player, msec );

        for (i = 0; i < player->ntracks; i++) {
            if (fluid_track_eot(player->track[i]) == FALSE) {
                status = FLUID_PLAYER_PLAYING;
                if (fluid_track_send_events(player->track[i], synth, player, FALSE,
                        player->duration_ticks) != FLUID_OK) {
                    /* */
                }
            }
        }

        if (status == FLUID_PLAYER_DONE) {
            FLUID_LOG(FLUID_INFO, "%s: %d: Duration=%.3f sec, Real time duration:%d sec", __FILE__,
                    __LINE__, player->midi_duration_msec, player->play_duration_msec);
            if( player->play_next_auto )
                loadnextfile = 1;
        }


    } while (loadnextfile);

    player->status = status;
    fluid_mutex_unlock( player->player_mutex );

    return 1;
}

/**
 * Activates play mode for a MIDI player if not already playing.
 * @param player MIDI player instance
 * @return #FLUID_OK on success, #FLUID_FAILED otherwise
 */
int
fluid_player_play(fluid_player_t *player)
{
    if (player->status == FLUID_PLAYER_PLAYING) {
        return FLUID_OK;
    }

    if (player->status != FLUID_PLAYER_READY) {
        FLUID_LOG( FLUID_ERR, "player not yet ready" );
        return FLUID_FAILED;
    }

    if (player->playlist == NULL) {
        return FLUID_OK;
    }

    player->status = FLUID_PLAYER_PLAYING;

    if (player->use_system_timer) {
        player->system_timer = new_fluid_timer((int) player->playdeltatime,
                fluid_player_callback, (void *) player, TRUE, FALSE, TRUE);
        if (player->system_timer == NULL) {
            return FLUID_FAILED;
        }
    } else {
        player->sample_timer = new_fluid_sample_timer(player->synth,
                fluid_player_callback, (void *) player);

        if (player->sample_timer == NULL) {
            return FLUID_FAILED;
        }
    }
    return FLUID_OK;
}

/**
 * Stops a MIDI player.
 * @param player MIDI player instance
 * @return Always returns #FLUID_OK
 */
int
fluid_player_stop(fluid_player_t *player)
{
    if (player->system_timer != NULL) {
        delete_fluid_timer(player->system_timer);
        player->system_timer = NULL;
    }
    if (player->sample_timer != NULL) {
        delete_fluid_sample_timer(player->synth, player->sample_timer);
        player->sample_timer = NULL;
    }
    player->status = FLUID_PLAYER_STOPPED;
    fluid_player_seek_ticks(player, fluid_player_get_duration_ticks(player));
    return FLUID_OK;
}

/**
 * Get MIDI player status.
 * @param player MIDI player instance
 * @return Player status (#fluid_player_status)
 * @since 1.1.0
 */
int
fluid_player_get_status(fluid_player_t *player)
{
    return player->status;
}

/**
 * Seek in the currently playing file.
 * @param player MIDI player instance
 * @param ticks the position to seek to in the current file
 * @return #FLUID_FAILED if ticks is negative or after the latest tick of the file,
 *   #FLUID_OK otherwise
 * @since 2.0.0
 *
 * The actual seek is performed during the player_callback.
 */
int fluid_player_seek_ticks(fluid_player_t* player, int ticks )
{
    int i = 0;
    int status = FLUID_PLAYER_DONE;
    int cur_delta_ticks = 0;
    int next_tempo_delta = 0;
    fluid_tempo_t* tempo = NULL;
    fluid_real_t delta_msec = 0.0f;

    if(ticks < 0 || ticks > fluid_player_get_total_ticks(player))
    {
        return FLUID_FAILED;
    }

    if( ticks == player->duration_ticks )
    {
        return FLUID_OK;
    }

    player->is_seeking = TRUE;
    fluid_mutex_lock( player->player_mutex );

    if( ticks < player->duration_ticks )
    {
        if( player->tempo != NULL )
            player->cur_tempo = player->tempo;
        player->miditempo = 500000;
        player->playtempo = rint(60000000.0/((60000000.0/(double)player->miditempo)+(double)player->tempo_offset_bpm));
        player->deltatime = (double)player->miditempo / (double)player->division / 1000.0; /* in milliseconds */
        player->playdeltatime = ((double)player->playtempo / (double)player->division / 1000.0)*player->tempo_offset_timeCoeff; /* in milliseconds */
        player->midi_duration_msec = 0.0f;
        player->duration_ticks = 0;
        player->tick_duration_msec = 0.0f;

        fluid_synth_all_sounds_off(player->synth, -1); /* avoid hanging notes */
    }

    if( player->midi_duration_msec - player->tick_duration_msec > 0 )
        delta_msec = (player->midi_duration_msec-player->tick_duration_msec);

    cur_delta_ticks = ticks-player->duration_ticks;

    do
    {
        while( player->cur_tempo != NULL )
        {
            tempo = (fluid_tempo_t*)player->cur_tempo->data;
            if( tempo == NULL )
            {
                player->cur_tempo = fluid_list_next(player->cur_tempo);
                continue;
            }

            if( tempo->ticks > player->duration_ticks )
            {
                next_tempo_delta = (tempo->ticks-player->duration_ticks);
                break;
            }

            fluid_player_set_midi_tempo( player, tempo->tempo );
            player->cur_tempo = fluid_list_next(player->cur_tempo);
        }

        if( cur_delta_ticks <= 0 )
            break;

        if( next_tempo_delta > 0 && cur_delta_ticks > next_tempo_delta )
        {
            player->tick_duration_msec += (player->deltatime*(fluid_real_t)next_tempo_delta);
            player->duration_ticks += next_tempo_delta;
            player->tick_duration_msec += (player->deltatime*(fluid_real_t)next_tempo_delta);
            cur_delta_ticks -= next_tempo_delta;
        }
        else
        {
            player->tick_duration_msec += (player->deltatime*(fluid_real_t)cur_delta_ticks);
            player->duration_ticks += cur_delta_ticks;
            cur_delta_ticks = 0;
        }

        next_tempo_delta = 0;
    }while(1);

    player->midi_duration_msec = player->tick_duration_msec+delta_msec;

    for (i = 0; i < player->ntracks; i++) {
        status = FLUID_PLAYER_PLAYING;
        if (fluid_track_send_events(player->track[i], player->synth, player, TRUE,
                                    player->duration_ticks) != FLUID_OK) {
            /* */
        }
    }

    player->status = status;

    if( player->play_next_auto == FALSE && player->status == FLUID_PLAYER_DONE )
        player->status = FLUID_PLAYER_READY;

    player->is_seeking = FALSE;
    fluid_mutex_unlock( player->player_mutex );

    return FLUID_OK;
}

int fluid_player_seek_msec(fluid_player_t* player, int msec )
{
    int i = 0;
    int status = FLUID_PLAYER_DONE;
    int next_tempo_delta = 0;
    int remain_ticks = 0;
    fluid_real_t remain_delta_msec = 0.0f;
    fluid_tempo_t* tempo = NULL;

    if(msec < 0 || msec > fluid_player_get_total_msec(player))
    {
        return FLUID_FAILED;
    }

    if(msec == fluid_player_get_duration_msec(player))
    {
        return FLUID_OK;
    }

    player->is_seeking = TRUE;
    fluid_mutex_lock( player->player_mutex );

    if( msec < player->midi_duration_msec )
    {
        if( player->tempo != NULL )
            player->cur_tempo = player->tempo;
        player->miditempo = 500000;
        player->playtempo = rint(60000000.0/((60000000.0/(double)player->miditempo)+(double)player->tempo_offset_bpm));
        player->deltatime = (double)player->miditempo / (double)player->division / 1000.0; /* in milliseconds */
        player->playdeltatime = ((double)player->playtempo / (double)player->division / 1000.0)*player->tempo_offset_timeCoeff; /* in milliseconds */
        player->midi_duration_msec = 0.0f;
        player->duration_ticks = 0;
        player->tick_duration_msec = 0.0f;

        fluid_synth_all_sounds_off(player->synth, -1); /* avoid hanging notes */
    }

    remain_delta_msec = ((fluid_real_t)msec-player->midi_duration_msec);

    if( player->midi_duration_msec - player->tick_duration_msec > 0 )
        remain_delta_msec += (player->midi_duration_msec-player->tick_duration_msec);

    do
    {
        while( player->cur_tempo != NULL )
        {
            tempo = (fluid_tempo_t*)player->cur_tempo->data;
            if( tempo == NULL )
            {
                player->cur_tempo = fluid_list_next(player->cur_tempo);
                continue;
            }

            if( tempo->ticks > player->duration_ticks )
            {
                next_tempo_delta = (tempo->ticks-player->duration_ticks);
                break;
            }

            fluid_player_set_midi_tempo( player, tempo->tempo );

            player->cur_tempo = fluid_list_next(player->cur_tempo);
        }

        if( remain_delta_msec > 0 )
            remain_ticks = (int)floor(remain_delta_msec/player->deltatime);

        if( remain_ticks <= 0 )
            break;

        if( next_tempo_delta > 0 && remain_ticks > next_tempo_delta )
        {
            player->tick_duration_msec += (player->deltatime*(fluid_real_t)next_tempo_delta);
            player->duration_ticks += next_tempo_delta;
            remain_delta_msec -= (player->deltatime*(fluid_real_t)next_tempo_delta);
        }
        else
        {
            player->tick_duration_msec += (player->deltatime*(fluid_real_t)remain_ticks);
            player->duration_ticks += remain_ticks;
            remain_delta_msec -= (player->deltatime*(fluid_real_t)remain_ticks);
        }

        next_tempo_delta = 0;
        remain_ticks = 0;
    }while( remain_delta_msec > 0 );

    player->midi_duration_msec = msec;

    for (i = 0; i < player->ntracks; i++) {
        status = FLUID_PLAYER_PLAYING;
        if (fluid_track_send_events(player->track[i], player->synth, player, TRUE,
                                    player->duration_ticks) != FLUID_OK) {
            /* */
        }
    }

    player->status = status;

    if( player->play_next_auto == FALSE && player->status == FLUID_PLAYER_DONE )
        player->status = FLUID_PLAYER_READY;

    player->is_seeking = FALSE;
    fluid_mutex_unlock( player->player_mutex );

    return FLUID_OK;
}

int fluid_player_get_channel_index_from_track_index(fluid_player_t* player, int ntrack)
{
    int channum = -1;
    fluid_midi_event_t* evt = NULL;

    if( ntrack < 0 || ntrack >= player->ntracks )
        return channum;

    evt = player->track[ntrack]->first;
    if( evt != NULL )
        channum = evt->channel;

    return channum;
}

/**
 * Enable looping of a MIDI player 
 * @param player MIDI player instance
 * @param loop Times left to loop the playlist. -1 means loop infinitely.
 * @return Always returns #FLUID_OK
 * @since 1.1.0
 *
 * For example, if you want to loop the playlist twice, set loop to 2 
 * and call this function before you start the player.
 */
int fluid_player_set_loop(fluid_player_t *player, int loop)
{
    player->loop = loop;
    return FLUID_OK;
}

/**
 * Set the tempo of a MIDI player.
 * @param player MIDI player instance
 * @param tempo Tempo to set playback speed to (in microseconds per quarter note, as per MIDI file spec)
 * @return Always returns #FLUID_OK
 */
int fluid_player_set_midi_tempo(fluid_player_t *player, int tempo)
{
    player->miditempo = tempo;
    player->playtempo = rint(60000000.0/((60000000.0/(double)player->miditempo)+(double)player->tempo_offset_bpm));
    player->deltatime = (double)player->miditempo / (double)player->division / 1000.0; /* in milliseconds */
    player->playdeltatime = ((double)player->playtempo / (double)player->division / 1000.0)*player->tempo_offset_timeCoeff; /* in milliseconds */

    if( player->use_system_timer )
    {
        if( player->system_timer )
            fluid_timer_change_interval( player->system_timer, player->playdeltatime );
    }

    FLUID_LOG(FLUID_INFO,
            "tempo=%d, tempo_offset_bpm=%d, tick time=%f msec, player->realtime_duration_msec, midi_duration_msec=%d msec, duration_ticks=%d",
            tempo, player->tempo_offset_bpm, player->deltatime, player->play_duration_msec, player->midi_duration_msec, player->duration_ticks);

    return FLUID_OK;
}

/**
 * Set the tempo of a MIDI player in beats per minute.
 * @param player MIDI player instance
 * @param bpm Tempo in beats per minute
 * @return Always returns #FLUID_OK
 */
int fluid_player_set_bpm(fluid_player_t *player, int bpm)
{
    return fluid_player_set_midi_tempo(player, (int) ((double) 60 * 1e6 / bpm));
}

int fluid_player_set_tempo_offset_bpm(fluid_player_t* player, int bpm)
{
    player->tempo_offset_timeCoeff = 1.0f;
    player->tempo_offset_bpm = bpm;

    return fluid_player_set_midi_tempo(player, player->miditempo);
}

int fluid_player_set_tempo_offset_timeCoeff(fluid_player_t* player, float coeff)
{
    player->tempo_offset_bpm = 0;
    player->tempo_offset_timeCoeff = coeff;

    return fluid_player_set_midi_tempo(player, player->miditempo);
}

/**
 * Wait for a MIDI player to terminate (when done playing).
 * @param player MIDI player instance
 * @return #FLUID_OK on success, #FLUID_FAILED otherwise
 */
int
fluid_player_join(fluid_player_t *player)
{
    if (player->system_timer) {
        return fluid_timer_join(player->system_timer);
    } else if (player->sample_timer) {
        /* Busy-wait loop, since there's no thread to wait for... */
        while (player->status != FLUID_PLAYER_STOPPED && player->status != FLUID_PLAYER_DONE) {
#if defined(WIN32)
            Sleep(10);
#else
            usleep(10000);
#endif
        }
    }
    return FLUID_OK;
}

/**
 * Get the number of tempo ticks passed.
 * @param player MIDI player instance
 * @return The number of tempo ticks passed
 * @since 1.1.7
 */
int fluid_player_get_duration_ticks(fluid_player_t * player)
{
    return player->duration_ticks;
}

int fluid_player_get_duration_msec(fluid_player_t * player)
{
    return (int)floor(player->midi_duration_msec);
}

int fluid_player_get_play_duration_msec(fluid_player_t * player)
{
    return player->play_duration_msec;
}

/**
 * Looks through all available MIDI tracks and gets the absolute tick of the very last event to play. 
 * @param player MIDI player instance
 * @return Total tick count of the sequence
 * @since 1.1.7
 */
int fluid_player_get_total_ticks(fluid_player_t * player)
{
    return player->total_ticks;
}

/**
 * Looks through all available MIDI tracks and gets the absolute tick of the very last event to play.
 * @param player MIDI player instance
 * @return Total tick count of the sequence
 * @since 1.1.7
 */
int fluid_player_get_total_msec(fluid_player_t * player)
{
    return player->total_msec;
}


int fluid_player_convert_msec_to_ticks(fluid_player_t * player, unsigned int msec)
{
    int ticks = 0;
    fluid_real_t duration = 0.0f;
    fluid_real_t duration_delta = 0.0f;
    fluid_list_t* cur_tempo = NULL;
    fluid_tempo_t* defaultTempo = new_fluid_tempo();
    fluid_tempo_t* tempo = NULL;
    fluid_real_t pre_deltatime = (fluid_real_t)defaultTempo->tempo/player->division / 1000.0;
    fluid_real_t deltatime = 0;

    if( msec > 0 )
    {
        if( player->tempo )
        {
            for (cur_tempo = player->tempo; cur_tempo; cur_tempo = fluid_list_next(cur_tempo))
            {
                if( cur_tempo->data == NULL )
                    continue;

                tempo = (fluid_tempo_t*)cur_tempo->data;
                deltatime = (fluid_real_t)tempo->tempo/player->division / 1000.0;
                duration_delta = ((fluid_real_t)((int)tempo->ticks-ticks)*pre_deltatime);

                if( (unsigned int)(rint(duration+duration_delta)) > msec )
                    break;

                duration += duration_delta;
                ticks = tempo->ticks;
                pre_deltatime = deltatime;
                duration_delta = 0;
                tempo = NULL;

                if( duration == msec )
                    break;
            }

            if( (unsigned int)(rint(duration)) < msec )
                ticks += (int)(((fluid_real_t)msec-duration) / pre_deltatime + 0.5);
        }
        else
        {
            ticks = (int)((fluid_real_t) msec / pre_deltatime + 0.5); /* 0.5 to average overall error when casting */
        }
    }

    delete_fluid_tempo( defaultTempo );

    return ticks;
}

unsigned int fluid_player_convert_ticks_to_msec(fluid_player_t * player, int ticks)
{
    int pre_ticks = 0;
    fluid_real_t duration = 0;
    fluid_real_t delta_duration = 0;
    fluid_list_t* cur_tempo = NULL;
    fluid_tempo_t* defaultTempo = new_fluid_tempo();
    fluid_tempo_t* tempo = NULL;
    fluid_real_t pre_deltatime = (fluid_real_t)defaultTempo->tempo/player->division / 1000.0;
    fluid_real_t deltatime = 0;

    if( ticks > 0 )
    {
        if( player->tempo )
        {
            for (cur_tempo = player->tempo; cur_tempo; cur_tempo = fluid_list_next(cur_tempo))
            {
                if( cur_tempo->data == NULL )
                    continue;

                tempo = (fluid_tempo_t*)cur_tempo->data;

                deltatime = (fluid_real_t)tempo->tempo/player->division / 1000.0;

                delta_duration = ((fluid_real_t)(tempo->ticks-pre_ticks)*pre_deltatime);
                if( tempo->ticks > ticks )
                    break;

                duration += delta_duration;
                pre_ticks = tempo->ticks;
                pre_deltatime = deltatime;
                delta_duration = 0;
                tempo = NULL;

                if( pre_ticks == ticks )
                    break;
            }

            if( pre_ticks < ticks )
            {
                delta_duration = ((fluid_real_t)(ticks-pre_ticks)*pre_deltatime );
                duration += delta_duration;
            }
        }
        else
        {
            duration = ((fluid_real_t)(ticks)*pre_deltatime);
        }
    }

    delete_fluid_tempo( defaultTempo );

    return (unsigned int)rint(duration);
}

/**
 * Get the tempo of a MIDI player in beats per minute.
 * @param player MIDI player instance
 * @return MIDI player tempo in BPM
 * @since 1.1.7
 */
int fluid_player_get_bpm(fluid_player_t * player)
{
    return (int)(60000000.0 / (double)player->miditempo);
}

/**
 * Get the tempo of a MIDI player.
 * @param player MIDI player instance
 * @return Tempo of the MIDI player (in microseconds per quarter note, as per MIDI file spec)
 * @since 1.1.7
 */
int fluid_player_get_midi_tempo(fluid_player_t * player)
{
    return player->miditempo;
}

int fluid_player_get_tempo_offset_bpm(fluid_player_t* player)
{
    return player->tempo_offset_bpm;
}

int fluid_player_get_play_tempo(fluid_player_t* player)
{
    return player->playtempo;
}

int fluid_player_get_player_bpm(fluid_player_t* player)
{
    return rint((60000000.0/(double)player->miditempo)+(double)player->tempo_offset_bpm);
}

double fluid_player_get_player_deltatime(fluid_player_t* player)
{
    return player->deltatime;
}

double fluid_player_get_player_playdeltatime(fluid_player_t* player)
{
    return player->playdeltatime;
}

/************************************************************************
 *       MIDI PARSER
 *
 */

/*
 * new_fluid_midi_parser
 */
fluid_midi_parser_t *
new_fluid_midi_parser ()
{
    fluid_midi_parser_t *parser;
    parser = FLUID_NEW(fluid_midi_parser_t);
    if (parser == NULL) {
        FLUID_LOG(FLUID_ERR, "Out of memory");
        return NULL;
    }
    FLUID_MEMSET (parser, 0, sizeof (fluid_midi_parser_t));

    parser->status = 0; /* As long as the status is 0, the parser won't do anything -> no need to initialize all the fields. */
    return parser;
}

/*
 * delete_fluid_midi_parser
 */
int
delete_fluid_midi_parser(fluid_midi_parser_t *parser)
{
    FLUID_FREE(parser);
    return FLUID_OK;
}

/**
 * Parse a MIDI stream one character at a time.
 * @param parser Parser instance
 * @param c Next character in MIDI stream
 * @return A parsed MIDI event or NULL if none.  Event is internal and should
 *   not be modified or freed and is only valid until next call to this function.
 */
fluid_midi_event_t *
fluid_midi_parser_parse(fluid_midi_parser_t *parser, unsigned char c)
{
    fluid_midi_event_t *event;

    /* Real-time messages (0xF8-0xFF) can occur anywhere, even in the middle
     * of another message. */
    if (c >= 0xF8) {
        if (c == MIDI_SYSTEM_RESET) {
            parser->event.type = c;
            parser->status = 0; /* clear the status */
            return &parser->event;
        }

        return NULL;
    }

    /* Status byte? - If previous message not yet complete, it is discarded (re-sync). */
    if (c & 0x80) {
        /* Any status byte terminates SYSEX messages (not just 0xF7) */
        if (parser->status == MIDI_SYSEX && parser->nr_bytes > 0) {
            event = &parser->event;
            fluid_midi_event_set_sysex(event, parser->data, parser->nr_bytes,
                    FALSE);
        } else
            event = NULL;

        if (c < 0xF0) /* Voice category message? */
        {
            parser->channel = c & 0x0F;
            parser->status = c & 0xF0;

            /* The event consumes x bytes of data... (subtract 1 for the status byte) */
            parser->nr_bytes_total = fluid_midi_event_length(parser->status)
                    - 1;

            parser->nr_bytes = 0; /* 0  bytes read so far */
        } else if (c == MIDI_SYSEX) {
            parser->status = MIDI_SYSEX;
            parser->nr_bytes = 0;
        } else
            parser->status = 0; /* Discard other system messages (0xF1-0xF7) */

        return event; /* Return SYSEX event or NULL */
    }

    /* Data/parameter byte */

    /* Discard data bytes for events we don't care about */
    if (parser->status == 0)
        return NULL;

    /* Max data size exceeded? (SYSEX messages only really) */
    if (parser->nr_bytes == FLUID_MIDI_PARSER_MAX_DATA_SIZE) {
        parser->status = 0; /* Discard the rest of the message */
        return NULL;
    }

    /* Store next byte */
    parser->data[parser->nr_bytes++] = c;

    /* Do we still need more data to get this event complete? */
    if (parser->status == MIDI_SYSEX || parser->nr_bytes < parser->nr_bytes_total)
        return NULL;

    /* Event is complete, return it.
     * Running status byte MIDI feature is also handled here. */
    parser->event.type = parser->status;
    parser->event.channel = parser->channel;
    parser->nr_bytes = 0; /* Reset data size, in case there are additional running status messages */

    switch (parser->status) {
        case NOTE_OFF:
        case NOTE_ON:
        case KEY_PRESSURE:
        case CONTROL_CHANGE:
        case PROGRAM_CHANGE:
        case CHANNEL_PRESSURE:
            parser->event.param1 = parser->data[0]; /* For example key number */
            parser->event.param2 = parser->data[1]; /* For example velocity */
            parser->event.param3 = 0;
            break;
        case PITCH_BEND:
            /* Pitch-bend is transmitted with 14-bit precision. */
            parser->event.param1 = (parser->data[1] << 7) | parser->data[0];
            break;
        default: /* Unlikely */
            return NULL;
    }

    return &parser->event;
}

/* Purpose:
 * Returns the length of a MIDI message. */
static int
fluid_midi_event_length(unsigned char event)
{
    switch (event & 0xF0) {
        case NOTE_OFF:
        case NOTE_ON:
        case KEY_PRESSURE:
        case CONTROL_CHANGE:
        case PITCH_BEND:
            return 3;
        case PROGRAM_CHANGE:
        case CHANNEL_PRESSURE:
            return 2;
    }
    switch (event) {
        case MIDI_TIME_CODE:
        case MIDI_SONG_SELECT:
        case 0xF4:
        case 0xF5:
            return 2;
        case MIDI_TUNE_REQUEST:
            return 1;
        case MIDI_SONG_POSITION:
            return 3;
    }
    return 1;
}
