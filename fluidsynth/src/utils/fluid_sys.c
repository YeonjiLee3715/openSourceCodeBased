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

#include "fluid_sys.h"


#if WITH_READLINE
#include <readline/readline.h>
#include <readline/history.h>
#endif

#ifdef DBUS_SUPPORT
#include "fluid_rtkit.h"
#endif

#ifdef ANDROID
#include <android/log.h>
#endif

/* WIN32 HACK - Flag used to differentiate between a file descriptor and a socket.
 * Should work, so long as no SOCKET or file descriptor ends up with this bit set. - JG */
#define WIN32_SOCKET_FLAG       0x40000000

/* SCHED_FIFO priority for high priority timer threads */
#define FLUID_SYS_TIMER_HIGH_PRIO_LEVEL         10


typedef struct
{
  fluid_thread_func_t func;
  void *data;
  int prio_level;
} fluid_thread_info_t;

struct _fluid_timer_t
{
  long last_timeout_msec;
  int timeout_count;
  long duration_msec;
  long interval_msec;
  fluid_timer_callback_t callback;
  void *data;
  fluid_thread_t *thread;
  int cont;
  int auto_destroy;
};

struct _fluid_server_socket_t
{
  fluid_socket_t socket;
  fluid_thread_t *thread;
  int cont;
  fluid_server_func_t func;
  void *data;
};


static int fluid_istream_gets(fluid_istream_t in, char* buf, int len);


static char fluid_errbuf[512];  /* buffer for error message */

static fluid_log_function_t fluid_log_function[LAST_LOG_LEVEL];
static void* fluid_log_user_data[LAST_LOG_LEVEL];
static int fluid_log_initialized = 0;

static const char fluid_libname[] = "fluidsynth";


void fluid_sys_config()
{
  fluid_log_config();
}


unsigned int fluid_debug_flags = 0;

#if DEBUG
/*
 * fluid_debug
 */
int fluid_debug(int level, char * fmt, ...)
{
  if (fluid_debug_flags & level) {
    fluid_log_function_t fun;
    va_list args;

    va_start (args, fmt);
    vsnprintf(fluid_errbuf, sizeof (fluid_errbuf), fmt, args);
    va_end (args);

    fun = fluid_log_function[FLUID_DBG];
    if (fun != NULL) {
      (*fun)(level, fluid_errbuf, fluid_log_user_data[FLUID_DBG]);
    }
  }
  return 0;
}
#endif

/**
 * Installs a new log function for a specified log level.
 * @param level Log level to install handler for.
 * @param fun Callback function handler to call for logged messages
 * @param data User supplied data pointer to pass to log function
 * @return The previously installed function.
 */
fluid_log_function_t
fluid_set_log_function(int level, fluid_log_function_t fun, void* data)
{
  fluid_log_function_t old = NULL;

  if ((level >= 0) && (level < LAST_LOG_LEVEL)) {
    old = fluid_log_function[level];
    fluid_log_function[level] = fun;
    fluid_log_user_data[level] = data;
  }
  return old;
}

/**
 * Default log function which prints to the stderr.
 * @param level Log level
 * @param message Log message
 * @param data User supplied data (not used)
 */
void
fluid_default_log_function(int level, const char* message, void* data)
{
  FILE* out;

#if defined(WIN32)
  out = stdout;
#else
  out = stderr;
#endif

  if (fluid_log_initialized == 0) {
    fluid_log_config();
  }

  switch (level) {
    case FLUID_PANIC:
#ifdef ANDROID
      __android_log_write(ANDROID_LOG_FATAL, fluid_libname, message);
#else
          FLUID_FPRINTF(out, "%s: panic: %s\n", fluid_libname, message);
#endif
          break;
    case FLUID_ERR:
#ifdef ANDROID
      __android_log_write(ANDROID_LOG_ERROR, fluid_libname, message);
#else
          FLUID_FPRINTF(out, "%s: error: %s\n", fluid_libname, message);
#endif
          break;
    case FLUID_WARN:
#ifdef ANDROID
      __android_log_write(ANDROID_LOG_WARN, fluid_libname, message);
#else
          FLUID_FPRINTF(out, "%s: warning: %s\n", fluid_libname, message);
#endif
          break;
    case FLUID_INFO:
#ifdef ANDROID
      __android_log_write(ANDROID_LOG_INFO, fluid_libname, message);
#else
          FLUID_FPRINTF(out, "%s: %s\n", fluid_libname, message);
#endif
          break;
    case FLUID_DBG:
#if DEBUG
      #ifdef ANDROID
          __android_log_write(ANDROID_LOG_DEBUG, fluid_libname, message);
        #else
          FLUID_FPRINTF(out, "%s: debug: %s\n", fluid_libname, message);
        #endif
#endif
      break;
    default:
#ifdef ANDROID
      __android_log_write(ANDROID_LOG_VERBOSE, fluid_libname, message);
#else
          FLUID_FPRINTF(out, "%s: %s\n", fluid_libname, message);
#endif
          break;
  }
  fflush(out);
}

/*
 * fluid_init_log
 */
void
fluid_log_config(void)
{
  if (fluid_log_initialized == 0) {

    fluid_log_initialized = 1;

    if (fluid_log_function[FLUID_PANIC] == NULL) {
      fluid_set_log_function(FLUID_PANIC, fluid_default_log_function, NULL);
    }

    if (fluid_log_function[FLUID_ERR] == NULL) {
      fluid_set_log_function(FLUID_ERR, fluid_default_log_function, NULL);
    }

    if (fluid_log_function[FLUID_WARN] == NULL) {
      fluid_set_log_function(FLUID_WARN, fluid_default_log_function, NULL);
    }

    if (fluid_log_function[FLUID_INFO] == NULL) {
      fluid_set_log_function(FLUID_INFO, fluid_default_log_function, NULL);
    }

    if (fluid_log_function[FLUID_DBG] == NULL) {
      fluid_set_log_function(FLUID_DBG, fluid_default_log_function, NULL);
    }
  }
}

/**
 * Print a message to the log.
 * @param level Log level (#fluid_log_level).
 * @param fmt Printf style format string for log message
 * @param ... Arguments for printf 'fmt' message string
 * @return Always returns #FLUID_FAILED
 */
int
fluid_log(int level, const char* fmt, ...)
{
  fluid_log_function_t fun = NULL;

  va_list args;
  va_start (args, fmt);
  vsnprintf(fluid_errbuf, sizeof (fluid_errbuf), fmt, args);
  va_end (args);

  if ((level >= 0) && (level < LAST_LOG_LEVEL)) {
    fun = fluid_log_function[level];
    if (fun != NULL) {
      (*fun)(level, fluid_errbuf, fluid_log_user_data[level]);
    }
  }
  return FLUID_FAILED;
}

/**
 * An improved strtok, still trashes the input string, but is portable and
 * thread safe.  Also skips token chars at beginning of token string and never
 * returns an empty token (will return NULL if source ends in token chars though).
 * NOTE: NOT part of public API
 * @internal
 * @param str Pointer to a string pointer of source to tokenize.  Pointer gets
 *   updated on each invocation to point to beginning of next token.  Note that
 *   token char get's overwritten with a 0 byte.  String pointer is set to NULL
 *   when final token is returned.
 * @param delim String of delimiter chars.
 * @return Pointer to the next token or NULL if no more tokens.
 */
char *fluid_strtok (char **str, char *delim)
{
  char *s, *d, *token;
  char c;

  if (str == NULL || delim == NULL || !*delim)
  {
    FLUID_LOG(FLUID_ERR, "Null pointer");
    return NULL;
  }

  s = *str;
  if (!s) return NULL;  /* str points to a NULL pointer? (tokenize already ended) */

  /* skip delimiter chars at beginning of token */
  do
  {
    c = *s;
    if (!c)     /* end of source string? */
    {
      *str = NULL;
      return NULL;
    }

    for (d = delim; *d; d++)    /* is source char a token char? */
    {
      if (c == *d)      /* token char match? */
      {
        s++;            /* advance to next source char */
        break;
      }
    }
  } while (*d);         /* while token char match */

  token = s;            /* start of token found */

  /* search for next token char or end of source string */
  for (s = s+1; *s; s++)
  {
    c = *s;

    for (d = delim; *d; d++)    /* is source char a token char? */
    {
      if (c == *d)      /* token char match? */
      {
        *s = '\0';      /* overwrite token char with zero byte to terminate token */
        *str = s+1;     /* update str to point to beginning of next token */
        return token;
      }
    }
  }

  /* we get here only if source string ended */
  *str = NULL;
  return token;
}

/*
 * fluid_error
 */
char*
fluid_error()
{
  return fluid_errbuf;
}

/**
 * Suspend the execution of the current thread for the specified amount of time.
 * @param milliseconds to wait.
 */
void fluid_msleep(unsigned int msecs)
{
  usleep(msecs * 1000);
}

/**
 * Check if a file is a MIDI file.
 * @param filename Path to the file to check
 * @return TRUE if it could be a MIDI file, FALSE otherwise
 *
 * The current implementation only checks for the "MThd" header in the file.
 * It is useful only to distinguish between SoundFont and MIDI files.
 */
int
fluid_is_midifile(const char *filename)
{
  FILE* fp = fopen(filename, "rb");
  char id[4];

  if (fp == NULL) {
    return 0;
  }
  if (fread((void*) id, 1, 4, fp) != 4) {
    fclose(fp);
    return 0;
  }
  fclose(fp);

  return strncmp(id, "MThd", 4) == 0;
}

/**
 * Check if a file is a SoundFont file.
 * @param filename Path to the file to check
 * @return TRUE if it could be a SoundFont, FALSE otherwise
 *
 * @note The current implementation only checks for the "RIFF" header in the file.
 * Version >1.1.9 also checks for "sfbk" header.
 * It is useful to distinguish between SoundFont and other (e.g. MIDI) files.
 */
int
fluid_is_soundfont(const char *filename)
{
  FILE* fp = fopen(filename, "rb");
  char riff_id[4], sfbk_id[4];

  if (fp == NULL) {
    return 0;
  }
  if((fread((void*) riff_id, 1, sizeof(riff_id), fp) != sizeof(riff_id)) ||
     (fseek(fp, 4, SEEK_CUR) != 0) || 
     (fread((void*) sfbk_id, 1, sizeof(sfbk_id), fp) != sizeof(sfbk_id)))
  {
      goto error_rec;
  }
  
  fclose(fp);
  return (FLUID_STRNCMP(riff_id, "RIFF", sizeof(riff_id)) == 0) &&
         (FLUID_STRNCMP(sfbk_id, "sfbk", sizeof(sfbk_id)) == 0);
  
error_rec:
    fclose(fp);
    return 0;
}

/**
 * Get time in milliseconds to be used in relative timing operations.
 * @return Unix time in milliseconds.
 */
unsigned int fluid_curtime(void)
{
  static long initial_seconds = 0;
  struct timespec timeval;

  if (initial_seconds == 0) {
    clock_gettime(CLOCK_REALTIME, &timeval);
    initial_seconds = timeval.tv_sec;
  }

  clock_gettime(CLOCK_REALTIME, &timeval);

  return (unsigned int)((timeval.tv_sec - initial_seconds) * 1000.0 + timeval.tv_nsec / 1000000.0);
}

/**
 * Get time in microseconds to be used in relative timing operations.
 * @return Unix time in microseconds.
 */
double
fluid_utime (void)
{
  struct timespec timeval;

  clock_gettime(CLOCK_REALTIME, &timeval);

  return (timeval.tv_sec * 1000000.0 + timeval.tv_nsec / 1000.0);
}


void
fluid_thread_self_set_prio (int prio_level)
{
  struct sched_param priority;

  if (prio_level > 0)
  {

    memset(&priority, 0, sizeof(priority));
    priority.sched_priority = prio_level;

    if (pthread_setschedparam (pthread_self (), SCHED_FIFO, &priority) == 0) {
      return;
    }
    FLUID_LOG(FLUID_WARN, "Failed to set thread to high priority");
  }
}

#ifdef FPE_CHECK

/***************************************************************
 *
 *               Floating point exceptions
 *
 *  The floating point exception functions were taken from Ircam's
 *  jMax source code. http://www.ircam.fr/jmax
 *
 *  FIXME: check in config for i386 machine
 *
 *  Currently not used. I leave the code here in case we want to pick
 *  this up again some time later.
 */

/* Exception flags */
#define _FPU_STATUS_IE    0x001  /* Invalid Operation */
#define _FPU_STATUS_DE    0x002  /* Denormalized Operand */
#define _FPU_STATUS_ZE    0x004  /* Zero Divide */
#define _FPU_STATUS_OE    0x008  /* Overflow */
#define _FPU_STATUS_UE    0x010  /* Underflow */
#define _FPU_STATUS_PE    0x020  /* Precision */
#define _FPU_STATUS_SF    0x040  /* Stack Fault */
#define _FPU_STATUS_ES    0x080  /* Error Summary Status */

/* Macros for accessing the FPU status word.  */

/* get the FPU status */
#define _FPU_GET_SW(sw) __asm__ ("fnstsw %0" : "=m" (*&sw))

/* clear the FPU status */
#define _FPU_CLR_SW() __asm__ ("fnclex" : : )

/* Purpose:
 * Checks, if the floating point unit has produced an exception, print a message
 * if so and clear the exception.
 */
unsigned int fluid_check_fpe_i386(char* explanation)
{
  unsigned int s;

  _FPU_GET_SW(s);
  _FPU_CLR_SW();

  s &= _FPU_STATUS_IE | _FPU_STATUS_DE | _FPU_STATUS_ZE | _FPU_STATUS_OE | _FPU_STATUS_UE;

  if (s)
  {
      FLUID_LOG(FLUID_WARN, "FPE exception (before or in %s): %s%s%s%s%s", explanation,
               (s & _FPU_STATUS_IE) ? "Invalid operation " : "",
               (s & _FPU_STATUS_DE) ? "Denormal number " : "",
               (s & _FPU_STATUS_ZE) ? "Zero divide " : "",
               (s & _FPU_STATUS_OE) ? "Overflow " : "",
               (s & _FPU_STATUS_UE) ? "Underflow " : "");
  }

  return s;
}

/* Purpose:
 * Clear floating point exception.
 */
void fluid_clear_fpe_i386 (void)
{
  _FPU_CLR_SW();
}

#endif  // ifdef FPE_CHECK

/***************************************************************
 *
 *               Profiling (Linux, i586 only)
 *
 */

#if WITH_PROFILING

fluid_profile_data_t fluid_profile_data[] =
{
  { FLUID_PROF_WRITE,            "fluid_synth_write_*             ", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_ONE_BLOCK,        "fluid_synth_one_block           ", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_ONE_BLOCK_CLEAR,  "fluid_synth_one_block:clear     ", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_ONE_BLOCK_VOICE,  "fluid_synth_one_block:one voice ", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_ONE_BLOCK_VOICES, "fluid_synth_one_block:all voices", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_ONE_BLOCK_REVERB, "fluid_synth_one_block:reverb    ", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_ONE_BLOCK_CHORUS, "fluid_synth_one_block:chorus    ", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_VOICE_NOTE,       "fluid_voice:note                ", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_VOICE_RELEASE,    "fluid_voice:release             ", 1e10, 0.0, 0.0, 0},
  { FLUID_PROF_LAST, "last", 1e100, 0.0, 0.0, 0}
};


void fluid_profiling_print(void)
{
  int i;

  printf("fluid_profiling_print\n");

  FLUID_LOG(FLUID_INFO, "Estimated times: min/avg/max (micro seconds)");

  for (i = 0; i < FLUID_PROF_LAST; i++) {
    if (fluid_profile_data[i].count > 0) {
      FLUID_LOG(FLUID_INFO, "%s: %.3f/%.3f/%.3f",
               fluid_profile_data[i].description,
               fluid_profile_data[i].min,
               fluid_profile_data[i].total / fluid_profile_data[i].count,
               fluid_profile_data[i].max);
    } else {
      FLUID_LOG(FLUID_DBG, "%s: no profiling available", fluid_profile_data[i].description);
    }
  }
}


#endif /* WITH_PROFILING */



/***************************************************************
 *
 *               Threads
 *
 */

static fluid_thread_return_t
fluid_thread_high_prio (void *data)
{
  fluid_thread_info_t *info = data;

  fluid_thread_self_set_prio (info->prio_level);

  info->func (info->data);
  if( info != NULL )
  {
    FLUID_FREE (info);
    info = NULL;
  }

  return FLUID_THREAD_RETURN_VALUE;
}

/**
 * Create a new thread.
 * @param func Function to execute in new thread context
 * @param data User defined data to pass to func
 * @param prio_level Priority level.  If greater than 0 then high priority scheduling will
 *   be used, with the given priority level (used by pthreads only).  0 uses normal scheduling.
 * @param detach If TRUE, 'join' does not work and the thread destroys itself when finished.
 * @return New thread pointer or NULL on error
 */
fluid_thread_t *
new_fluid_thread (const char *name, fluid_thread_func_t func, void *data, int prio_level, int detach)
{
  fluid_thread_t *thread;
  fluid_thread_info_t *info;

  fluid_return_val_if_fail (func != NULL, NULL);

  thread = FLUID_NEW(fluid_thread_t);
  if (prio_level > 0) {
    info = FLUID_NEW (fluid_thread_info_t);
    if (!info) {
      FLUID_LOG(FLUID_ERR, "Out of memory");
      return NULL;
    }
    FLUID_MEMSET (info, 0, sizeof (fluid_thread_info_t));

    info->func = func;
    info->data = data;
    info->prio_level = prio_level;
    data = info;
    func = fluid_thread_high_prio;
  }

  pthread_create(thread, NULL, func, data);

  if (!thread) {
    FLUID_LOG(FLUID_ERR, "Failed to create the thread");
    return NULL;
  }

  if (detach) {
    pthread_detach(*thread);
  }

  return thread;
}

/**
 * Frees data associated with a thread (does not actually stop thread).
 * @param thread Thread to free
 */
void
delete_fluid_thread(fluid_thread_t* thread)
{
  /* Threads free themselves when they quit, nothing to do */
  if( thread != NULL )
    free( thread );
}

/**
 * Join a thread (wait for it to terminate).
 * @param thread Thread to join
 * @return FLUID_OK
 */
int
fluid_thread_join(fluid_thread_t* thread)
{
  pthread_join(*thread, NULL);

  return FLUID_OK;
}


fluid_thread_return_t
fluid_timer_run (void *data)
{
  fluid_timer_t *timer;
  int cont;
  long delay;
  long cur_timeout_msec = 0;
  long delta_msec = 0;
  long cur_msec = 0;

  timer = (fluid_timer_t *)data;

  /* keep track of the start time for absolute positioning */
  timer->last_timeout_msec = fluid_curtime();

  while (timer->cont)
  {
    cur_msec = fluid_curtime();
    cur_timeout_msec = cur_msec;
    delta_msec = (long)rint((double)(cur_msec-timer->last_timeout_msec));
    if( delta_msec > 0 )
    {
      timer->duration_msec += delta_msec;
      timer->last_timeout_msec = cur_msec;
      cont = (*timer->callback)(timer->data, timer->duration_msec);
    }

    if (!cont) break;

    /* to avoid incremental time errors, calculate the delay between
       two callbacks bringing in the "absolute" time (count *
       timer->msec) */

    cur_msec = fluid_curtime();
    delta_msec = (long)rint((double)(cur_msec-timer->last_timeout_msec));
    if( delta_msec > 0 )
    {
      timer->duration_msec += delta_msec;
      timer->last_timeout_msec = cur_msec;
    }
    delay = timer->interval_msec-(cur_msec-cur_timeout_msec);
    if (delay > 0) fluid_msleep(delay);
  }

  FLUID_LOG (FLUID_DBG, "Timer thread finished");

  if (timer->auto_destroy)
    if( timer != NULL )
    {
      FLUID_FREE (timer);
      timer = NULL;
    }

  return FLUID_THREAD_RETURN_VALUE;
}

fluid_timer_t*
new_fluid_timer (int msec, fluid_timer_callback_t callback, void* data,
                 int new_thread, int auto_destroy, int high_priority)
{
  fluid_timer_t *timer = NULL;

  timer = FLUID_NEW (fluid_timer_t);
  if (timer == NULL)
  {
    FLUID_LOG (FLUID_ERR, "Out of memory");
    return NULL;
  }
  FLUID_MEMSET (timer, 0, sizeof (fluid_timer_t));

  timer->last_timeout_msec = 0;
  timer->duration_msec = 0;
  timer->interval_msec = msec;
  timer->callback = callback;
  timer->data = data;
  timer->cont = TRUE ;
  timer->thread = NULL;
  timer->auto_destroy = auto_destroy;

  if (new_thread)
  {
    timer->thread = new_fluid_thread ("timer", fluid_timer_run, timer, high_priority
                                      ? FLUID_SYS_TIMER_HIGH_PRIO_LEVEL : 0, FALSE);
    if (!timer->thread)
    {
      if( timer != NULL )
      {
        FLUID_FREE (timer);
        timer = NULL;
      }
      return NULL;
    }
  }
  else
  {
      fluid_timer_run (timer);  /* Run directly, instead of as a separate thread */
      if(auto_destroy)
      {
          /* do NOT return freed memory */
          return NULL;
      }
  }

  return timer;
}

int fluid_timer_change_interval(fluid_timer_t* timer, int msec)
{
  if( timer == NULL )
    return FLUID_FAILED;

  timer->interval_msec = msec;

  return FLUID_OK;
}

int fluid_timer_seek(fluid_timer_t* timer, int msec)
{
  if( timer == NULL )
    return FLUID_FAILED;

  timer->duration_msec = msec;
  timer->last_timeout_msec = fluid_curtime();

  return FLUID_OK;
}

int
delete_fluid_timer (fluid_timer_t *timer)
{
  int auto_destroy = timer->auto_destroy;

  timer->cont = 0;
  fluid_timer_join (timer);

  /* Shouldn't access timer now if auto_destroy enabled, since it has been destroyed */

  if (!auto_destroy){
    if( timer != NULL )
    {
      FLUID_FREE (timer);
      timer = NULL;
    }
  }

  return FLUID_OK;
}

int
fluid_timer_join (fluid_timer_t *timer)
{
  int auto_destroy;

  if (timer->thread)
  {
    auto_destroy = timer->auto_destroy;
    fluid_thread_join (timer->thread);

    if (!auto_destroy) timer->thread = NULL;
  }

  return FLUID_OK;
}


/***************************************************************
 *
 *               Sockets and I/O
 *
 */

/**
 * Get standard in stream handle.
 * @return Standard in stream.
 */
fluid_istream_t
fluid_get_stdin (void)
{
  return STDIN_FILENO;
}

/**
 * Get standard output stream handle.
 * @return Standard out stream.
 */
fluid_ostream_t
fluid_get_stdout (void)
{
  return STDOUT_FILENO;
}

/**
 * Read a line from an input stream.
 * @return 0 if end-of-stream, -1 if error, non zero otherwise
 */
int
fluid_istream_readline (fluid_istream_t in, fluid_ostream_t out, char* prompt,
                        char* buf, int len)
{
#if WITH_READLINE
  if (in == fluid_get_stdin ())
  {
    char *line;

    line = readline (prompt);

    if (line == NULL)
      return -1;

    snprintf(buf, len, "%s", line);
    buf[len - 1] = 0;

    free(line);
    return 1;
  }
  else
#endif
  {
    fluid_ostream_printf (out, "%s", prompt);
    return fluid_istream_gets (in, buf, len);
  }
}

/**
 * Reads a line from an input stream (socket).
 * @param in The input socket
 * @param buf Buffer to store data to
 * @param len Maximum length to store to buf
 * @return 1 if a line was read, 0 on end of stream, -1 on error
 */
static int
fluid_istream_gets (fluid_istream_t in, char* buf, int len)
{
  char c;
  int n;

  buf[len - 1] = 0;

  while (--len > 0)
  {
#ifndef WIN32
    n = read(in, &c, 1);
    if (n == -1) return -1;
#else
    /* Handle read differently depending on if its a socket or file descriptor */
    if (!(in & WIN32_SOCKET_FLAG))
    {
      n = read(in, &c, 1);
      if (n == -1) return -1;
    }
    else
    {
      n = recv(in & ~WIN32_SOCKET_FLAG, &c, 1, 0);
      if (n == SOCKET_ERROR) return -1;
    }
#endif

    if (n == 0)
    {
      *buf++ = 0;
      return 0;
    }

    if (c == '\n')
    {
      *buf++ = 0;
      return 1;
    }

    /* Store all characters excluding CR */
    if (c != '\r') *buf++ = c;
  }

  return -1;
}

/**
 * Send a printf style string with arguments to an output stream (socket).
 * @param out Output stream
 * @param format printf style format string
 * @param ... Arguments for the printf format string
 * @return Number of bytes written or -1 on error
 */
int
fluid_ostream_printf (fluid_ostream_t out, const char* format, ...)
{
  char buf[4096];
  va_list args;
  int len;

  va_start (args, format);
  len = vsnprintf (buf, 4095, format, args);
  va_end (args);

  if (len == 0)
  {
    return 0;
  }

  if (len < 0)
  {
    printf("fluid_ostream_printf: buffer overflow");
    return -1;
  }

  buf[4095] = 0;

#ifndef WIN32
  return write (out, buf, strlen (buf));
#else
  {
    int retval;

    /* Handle write differently depending on if its a socket or file descriptor */
    if (!(out & WIN32_SOCKET_FLAG))
      return write(out, buf, strlen (buf));

    /* Socket */
    retval = send (out & ~WIN32_SOCKET_FLAG, buf, strlen (buf), 0);

    return retval != SOCKET_ERROR ? retval : -1;
  }
#endif
}

int fluid_server_socket_join(fluid_server_socket_t *server_socket)
{
  return fluid_thread_join (server_socket->thread);
}


#ifndef WIN32           // Not win32?

#define SOCKET_ERROR -1

fluid_istream_t fluid_socket_get_istream (fluid_socket_t sock)
{
  return sock;
}

fluid_ostream_t fluid_socket_get_ostream (fluid_socket_t sock)
{
  return sock;
}

void fluid_socket_close(fluid_socket_t sock)
{
  if (sock != INVALID_SOCKET)
    close (sock);
}

static void
fluid_server_socket_run (void *data)
{
  fluid_server_socket_t *server_socket = (fluid_server_socket_t *)data;
  fluid_socket_t client_socket;
#ifdef IPV6_SUPPORT
  struct sockaddr_in6 addr;
  char straddr[INET6_ADDRSTRLEN];
#else
  struct sockaddr_in addr;
  char straddr[INET_ADDRSTRLEN];
#endif
  socklen_t addrlen = sizeof (addr);
  int retval;
  FLUID_MEMSET((char *)&addr, 0, sizeof(addr));

  FLUID_LOG (FLUID_DBG, "Server listening for connections");

  while (server_socket->cont)
  {
    client_socket = accept (server_socket->socket, (struct sockaddr *)&addr, &addrlen);

    FLUID_LOG (FLUID_DBG, "New client connection");

    if (client_socket == INVALID_SOCKET)
    {
      if (server_socket->cont)
        FLUID_LOG(FLUID_ERR, "Failed to accept connection");

      server_socket->cont = 0;
      return;
    } else {
#ifdef HAVE_INETNTOP
#ifdef IPV6_SUPPORT
      inet_ntop(AF_INET6, &addr.sin6_addr, straddr, sizeof(straddr));
#else
      inet_ntop(AF_INET, &addr.sin_addr, straddr, sizeof(straddr));
#endif
#endif
#ifdef HAVE_INETNTOP
      retval = server_socket->func (server_socket->data, client_socket,
                                    straddr);
#else
      retval = server_socket->func (server_socket->data, client_socket,
                                    inet_ntoa (addr.sin_addr));
#endif

      if (retval != 0)
        fluid_socket_close(client_socket);
    }
  }

  FLUID_LOG(FLUID_DBG, "Server closing");
}

fluid_server_socket_t*
new_fluid_server_socket(int port, fluid_server_func_t func, void* data)
{
  fluid_server_socket_t* server_socket;
#ifdef IPV6_SUPPORT
  struct sockaddr_in6 addr;
#else
  struct sockaddr_in addr;
#endif
  fluid_socket_t sock;

  fluid_return_val_if_fail (func != NULL, NULL);
#ifdef IPV6_SUPPORT
  sock = socket(AF_INET6, SOCK_STREAM, 0);
  if (sock == INVALID_SOCKET) {
    FLUID_LOG(FLUID_ERR, "Failed to create server socket");
    return NULL;
  }

  FLUID_MEMSET((char *)&addr, 0, sizeof(struct sockaddr_in6));
  addr.sin6_family = AF_INET6;
  addr.sin6_addr = in6addr_any;
  addr.sin6_port = htons(port);
#else

  sock = socket(AF_INET, SOCK_STREAM, 0);
  if (sock == INVALID_SOCKET) {
    FLUID_LOG(FLUID_ERR, "Failed to create server socket");
    return NULL;
  }

  FLUID_MEMSET((char *)&addr, 0, sizeof(struct sockaddr_in));
  addr.sin_family = AF_INET;
  addr.sin_addr.s_addr = htonl(INADDR_ANY);
  addr.sin_port = htons(port);
#endif
  if (bind(sock, (const struct sockaddr *) &addr, sizeof(addr)) == SOCKET_ERROR) {
    FLUID_LOG(FLUID_ERR, "Failed to bind server socket");
    fluid_socket_close(sock);
    return NULL;
  }

  if (listen(sock, 10) == SOCKET_ERROR) {
    FLUID_LOG(FLUID_ERR, "Failed listen on server socket");
    fluid_socket_close(sock);
    return NULL;
  }

  server_socket = FLUID_NEW(fluid_server_socket_t);
  if (server_socket == NULL) {
    FLUID_LOG(FLUID_ERR, "Out of memory");
    fluid_socket_close(sock);
    return NULL;
  }
  FLUID_MEMSET (server_socket, 0, sizeof (fluid_server_socket_t));

  server_socket->socket = sock;
  server_socket->func = func;
  server_socket->data = data;
  server_socket->cont = 1;

  server_socket->thread = new_fluid_thread("server", fluid_server_socket_run, server_socket,
                                           0, FALSE);
  if (server_socket->thread == NULL) {
    if( server_socket != NULL )
    {
      FLUID_FREE (server_socket);
      server_socket = NULL;
    }
    fluid_socket_close(sock);
    return NULL;
  }

  return server_socket;
}

int delete_fluid_server_socket(fluid_server_socket_t* server_socket)
{
  server_socket->cont = 0;
  if (server_socket->socket != INVALID_SOCKET) {
    fluid_socket_close(server_socket->socket);
  }
  if (server_socket->thread) {
    delete_fluid_thread(server_socket->thread);
    server_socket->thread = NULL;
  }

  if( server_socket != NULL )
  {
    FLUID_FREE(server_socket);
    server_socket = NULL;
  }
  return FLUID_OK;
}


#else           // Win32 is "special"


#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif

fluid_istream_t fluid_socket_get_istream (fluid_socket_t sock)
{
  return sock | WIN32_SOCKET_FLAG;
}

fluid_ostream_t fluid_socket_get_ostream (fluid_socket_t sock)
{
  return sock | WIN32_SOCKET_FLAG;
}

void fluid_socket_close (fluid_socket_t sock)
{
  if (sock != INVALID_SOCKET)
    closesocket (sock);
}

static void fluid_server_socket_run (void *data)
{
  fluid_server_socket_t *server_socket = (fluid_server_socket_t *)data;
  fluid_socket_t client_socket;
#ifdef IPV6_SUPPORT
  struct sockaddr_in6 addr;
  char straddr[INET6_ADDRSTRLEN];
#else
  struct sockaddr_in addr;
#ifdef HAVE_INETNTOP
  char straddr[INET_ADDRSTRLEN];
#endif
#endif
  socklen_t addrlen = sizeof (addr);
  int r;
  FLUID_MEMSET((char *)&addr, 0, sizeof(addr));

  FLUID_LOG(FLUID_DBG, "Server listening for connections");

  while (server_socket->cont)
  {
    client_socket = accept (server_socket->socket, (struct sockaddr *)&addr, &addrlen);

    FLUID_LOG (FLUID_DBG, "New client connection");

    if (client_socket == INVALID_SOCKET)
    {
      if (server_socket->cont)
        FLUID_LOG (FLUID_ERR, "Failed to accept connection: %ld", WSAGetLastError ());

      server_socket->cont = 0;
      return;
    }
    else
    {
#ifdef HAVE_INETNTOP
#ifdef IPV6_SUPPORT
      inet_ntop(AF_INET6, &addr.sin6_addr, straddr, sizeof(straddr));
#else
      inet_ntop(AF_INET, &addr.sin_addr, straddr, sizeof(straddr));
#endif
#endif
#ifdef HAVE_INETNTOP
      r = server_socket->func (server_socket->data, client_socket,
                               straddr);
#else
      r = server_socket->func (server_socket->data, client_socket,
                               inet_ntoa (addr.sin_addr));
#endif
      if (r != 0)
        fluid_socket_close (client_socket);
    }
  }

  FLUID_LOG (FLUID_DBG, "Server closing");
}

fluid_server_socket_t*
new_fluid_server_socket(int port, fluid_server_func_t func, void* data)
{
  fluid_server_socket_t* server_socket;
#ifdef IPV6_SUPPORT
  struct sockaddr_in6 addr;
#else
  struct sockaddr_in addr;
#endif

  fluid_socket_t sock;
  WSADATA wsaData;
  int retval;

  fluid_return_val_if_fail (func != NULL, NULL);

  // Win32 requires initialization of winsock
  retval = WSAStartup (MAKEWORD (2,2), &wsaData);

  if (retval != 0)
  {
    FLUID_LOG(FLUID_ERR, "Server socket creation error: WSAStartup failed: %d", retval);
    return NULL;
  }
#ifdef IPV6_SUPPORT
  sock = socket (AF_INET6, SOCK_STREAM, 0);
  if (sock == INVALID_SOCKET)
  {
    FLUID_LOG (FLUID_ERR, "Failed to create server socket: %ld", WSAGetLastError ());
    WSACleanup ();
    return NULL;
  }
  addr.sin6_family = AF_INET6;
  addr.sin6_port = htons (port);
  addr.sin6_addr = in6addr_any;
#else

  sock = socket (AF_INET, SOCK_STREAM, 0);

  if (sock == INVALID_SOCKET)
  {
    FLUID_LOG (FLUID_ERR, "Failed to create server socket: %ld", WSAGetLastError ());
    WSACleanup ();
    return NULL;
  }

  addr.sin_family = AF_INET;
  addr.sin_port = htons (port);
  addr.sin_addr.s_addr = htonl (INADDR_ANY);
#endif
  retval = bind (sock, (struct sockaddr *)&addr, sizeof (addr));

  if (retval == SOCKET_ERROR)
  {
    FLUID_LOG (FLUID_ERR, "Failed to bind server socket: %ld", WSAGetLastError ());
    fluid_socket_close (sock);
    WSACleanup ();
    return NULL;
  }

  if (listen (sock, SOMAXCONN) == SOCKET_ERROR)
  {
    FLUID_LOG (FLUID_ERR, "Failed to listen on server socket: %ld", WSAGetLastError ());
    fluid_socket_close (sock);
    WSACleanup ();
    return NULL;
  }

  server_socket = FLUID_NEW (fluid_server_socket_t);
  if (server_socket == NULL)
  {
    FLUID_LOG (FLUID_ERR, "Out of memory");
    fluid_socket_close (sock);
    WSACleanup ();
    return NULL;
  }
  FLUID_MEMSET (server_socket, 0, sizeof (fluid_server_socket_t));

  server_socket->socket = sock;
  server_socket->func = func;
  server_socket->data = data;
  server_socket->cont = 1;

  server_socket->thread = new_fluid_thread("server", fluid_server_socket_run, server_socket,
                                           0, FALSE);
  if (server_socket->thread == NULL)
  {
    FLUID_FREE (server_socket);
    fluid_socket_close (sock);
    WSACleanup ();
    return NULL;
  }

  return server_socket;
}

int delete_fluid_server_socket(fluid_server_socket_t *server_socket)
{
  server_socket->cont = 0;

  if (server_socket->socket != INVALID_SOCKET)
    fluid_socket_close (server_socket->socket);

  if (server_socket->thread) {
    fluid_thread_join(server_socket->thread);
    delete_fluid_thread (server_socket->thread);
  }

  FLUID_FREE (server_socket);

  WSACleanup ();        // Should be called the same number of times as WSAStartup

  return FLUID_OK;
}

#endif
