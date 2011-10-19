#ifndef _INC_MESSAGE_H
#define _INC_MESSAGE_H
#if !defined(PHREEQC_CLASS)
#define CLASS_STATIC
#else
#if !defined(_INC_PHREEQC_H)
#define CLASS_STATIC
#endif
#endif
#include <stdarg.h>

#ifdef USE_OLD_IO
typedef int (*PFN_OUTPUT_CALLBACK) (const int action, const int type,
									const char *err_str, const int stop,
									void *cookie, const char *, va_list args);

struct output_callback
{
	PFN_OUTPUT_CALLBACK callback;
	void *cookie;
};
int add_output_callback(PFN_OUTPUT_CALLBACK pfn, void *cookie);
int clean_up_output_callbacks(void);
int output_msg(int type, const char *format, ...);
int warning_msg(const char *err_str, ...);
int error_msg(const char *err_str, int stop, ...);
CLASS_STATIC int phreeqc_handler(int action, int type, const char *err_str,
					int stop, void *cookie, const char *, va_list args);
/*
 *  Functions for output callbacks
 */
int output_open(const int type, const char *file_name, ...);
int output_fflush(const int type, ...);
int output_rewind(const int type, ...);
int output_close(const int type, ...);
#else
int output_msg(int type, const char *format, ...);
int warning_msg(const char *err_str);
//int error_msg(const char *err_str, int stop);
/*
 *  Functions for output callbacks
 */
int output_open(int type, const char *file_name);
int output_fflush(int type);
int output_rewind(int type);
int output_close(int type);
#endif
int fpunchf(const char *name, const char *format, ...);
int fpunchf_user(int user_index, const char *format, ...);
int fpunchf_end_row(const char *format, ...);
void set_forward_output_to_log(int value);
int get_forward_output_to_log(void);
typedef enum
{
	//OUTPUT_ERROR,
	OUTPUT_WARNING,
	OUTPUT_MESSAGE,
	OUTPUT_PUNCH,
	OUTPUT_SCREEN,
	OUTPUT_LOG,
	OUTPUT_CHECKLINE,
	OUTPUT_GUI_ERROR,
	//OUTPUT_DUMP,
	//OUTPUT_SEND_MESSAGE,
	OUTPUT_ECHO,
	OUTPUT_PUNCH_END_ROW
} output_type;

typedef enum
{
	ACTION_OPEN,
	ACTION_OUTPUT,
	ACTION_CLOSE,
	ACTION_REWIND,
	ACTION_FLUSH
} action_type;

#endif /* _INC_MESSAGE_H */
