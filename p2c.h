#ifndef P2C_H
#define P2C_H

#ifdef PHREEQC_IDENT
static char const svnidp2c[] = "$Id$";
#endif

/* Header file for code generated by "p2c", the Pascal-to-C translator */

/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version 1.20.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */


#include <stdio.h>



/* If the following heuristic fails, compile -DBSD=0 for non-BSD systems,
   or -DBSD=1 for BSD systems. */

#ifdef M_XENIX
# define BSD 0
#endif

#ifdef vms
# define BSD 0
# ifndef __STDC__
#  define __STDC__ 1
# endif
#endif

#ifdef __TURBOC__
# define MSDOS 1
#endif

#ifdef MSDOS
# define BSD 0
#endif

#ifdef FILE       /* a #define in BSD, a typedef in SYSV (hp-ux, at least) */
# ifndef BSD	  /*  (a convenient, but horrible kludge!) */
#  define BSD 1
# endif
#endif

#ifdef BSD
# if !BSD
#  undef BSD
# endif
#endif


#if (defined(__STDC__) && !defined(M_XENIX)) || defined(__TURBOC__)
# include <stddef.h>
# include <stdlib.h>
# define HAS_STDLIB
# if defined(vms) || defined(__TURBOC__)
#  define __ID__(a)a
# endif
#else
# ifndef BSD
#  ifndef __TURBOC__
#   include <memory.h>
#  endif
# endif
# ifdef hpux
#  ifdef _INCLUDE__STDC__
#   include <stddef.h>
#   include <stdlib.h>
#  endif
# endif
# include <sys/types.h>
# if !defined(MSDOS) || defined(__TURBOC__)
#  define __ID__(a)a
# endif
#endif

#ifdef __ID__
# define __CAT__(a,b)__ID__(a)b
#else
# define __CAT__(a,b)a##b
#endif


#ifdef BSD
# include <strings.h>
# define memcpy(a,b,n) (bcopy(b,a,n),a)
# define memcmp(a,b,n) bcmp(a,b,n)
# define strchr(s,c) index(s,c)
# define strrchr(s,c) rindex(s,c)
#else
# include <string.h>
#endif

#include <ctype.h>
#include <math.h>
#include <setjmp.h>
#include <assert.h>


#ifndef NO_LACK
#ifdef vms

#define LACK_LABS
#define LACK_MEMMOVE
#define LACK_MEMCPY

#else

#define LACK_LABS       /* Undefine these if your library has these */
#ifdef SKIP
#define LACK_MEMMOVE
#endif

#endif
#endif


typedef struct __p2c_jmp_buf {
    struct __p2c_jmp_buf *next;
    jmp_buf jbuf;
} __p2c_jmp_buf;


/* Warning: The following will not work if setjmp is used simultaneously.
   This also violates the ANSI restriction about using vars after longjmp,
   but a typical implementation of longjmp will get it right anyway. */

#ifndef FAKE_TRY
# define TRY(x)         do { __p2c_jmp_buf __try_jb;  \
			     __try_jb.next = __top_jb;  \
			     if (!setjmp((__top_jb = &__try_jb)->jbuf)) {
# define RECOVER(x)	__top_jb = __try_jb.next; } else {
#ifdef SKIP
# define RECOVER2(x,L)  __top_jb = __try_jb.next; } else {  \
			     if (0) { L: __top_jb = __try_jb.next; }
#endif
# define RECOVER2(x,L)  __top_jb = __try_jb.next; } else {  \
			     { L: __top_jb = __try_jb.next; }
# define ENDTRY(x)      } } while (0) 
#else
# define TRY(x)         if (1) {
# define RECOVER(x)     } else do {
# define RECOVER2(x,L)  } else do { L: ;
# define ENDTRY(x)      } while (0)
#endif



#ifdef M_XENIX  /* avoid compiler bug */
# define INT_MAX  (32767)
# define INT_MIN  (-32768)
#endif


/* The following definitions work only on twos-complement machines */
#ifndef INT_MAX
# define INT_MAX  ((int)(((unsigned int) -1) >> 1))
# define INT_MIN  (~INT_MAX)
#endif

#ifndef INT_MAX
# define INT_MAX    ((int)(((unsigned int) -1) >> 1))
# define INT_MIN    (~INT_MAX)
#endif

#ifndef LONG_MAX
# define LONG_MAX   ((long)(((unsigned long) -1) >> 1))
# define LONG_MIN   (~LONG_MAX)
#endif

#ifndef SEEK_SET
# define SEEK_SET   0
# define SEEK_CUR   1
# define SEEK_END   2
#endif

#ifndef EXIT_SUCCESS
# ifdef vms
#  define EXIT_SUCCESS  1
#  define EXIT_FAILURE  (02000000000L)
# else
#  define EXIT_SUCCESS  0
#  define EXIT_FAILURE  1
# endif
#endif


#define SETBITS  32


#if defined(__STDC__) || defined(__TURBOC__)
# if !defined(vms) && !defined(M_LINT)
#  define Signed    signed
# else
#  define Signed
# endif
# define Void       void      /* Void f() = procedure */
# ifndef Const
#  define Const     const
# endif
# ifndef Volatile
# define Volatile   volatile
# endif
# ifdef M_LINT
#  define P2PP(x)     ()
#  define PV()	    ()
typedef char *Anyptr;
# else
#  define P2PP(x)     x         /* function prototype */
#  define PV()      (void)    /* null function prototype */
typedef void *Anyptr;
# endif
#else
# define Signed
# define Void       void
# ifndef Const
#  define Const
# endif
# ifndef Volatile
#  define Volatile
# endif
# define P2PP(x)      ()
# define PV()       ()
typedef char *Anyptr;
#endif

#ifdef __GNUC__
# define Inline     inline
#else
# define Inline
#endif

#define Register    register  /* Register variables */
#define Char        char      /* Characters (not bytes) */

#ifndef Static
# define Static     static    /* Private global funcs and vars */
#endif

#ifndef Local
# define Local      static    /* Nested functions */
#endif

typedef Signed   char schar;
/*typedef unsigned char uchar;*/
typedef unsigned char boolean;

#ifndef true
# define true    1
# define false   0
#endif

#ifndef TRUE
# define TRUE    1
# define FALSE   0
#endif


typedef struct {
    Anyptr proc, link;
} _PROCEDURE;

#ifndef _FNSIZE
# define _FNSIZE  120
#endif


extern Void    PASCAL_MAIN  P2PP( (int, Char **) );
/*
extern Char    **P_argv;
extern int     P_argc;
*/
extern int   P_escapecode;
extern int     P_ioresult;
extern __p2c_jmp_buf *__top_jb;


#ifdef P2C_H_PROTO   /* if you have Ansi C but non-prototyped header files */
extern Char    *strcat      P2PP( (Char *, Const Char *) );
extern Char    *strchr      P2PP( (Const Char *, int) );
extern int      strcmp      P2PP( (Const Char *, Const Char *) );
extern Char    *strcpy      P2PP( (Char *, Const Char *) );
extern size_t   strlen      P2PP( (Const Char *) );
extern Char    *strncat     P2PP( (Char *, Const Char *, size_t) );
extern int      strncmp     P2PP( (Const Char *, Const Char *, size_t) );
extern Char    *strncpy     P2PP( (Char *, Const Char *, size_t) );
extern Char    *strrchr     P2PP( (Const Char *, int) );

extern Anyptr   memchr      P2PP( (Const Anyptr, int, size_t) );
extern Anyptr   memmove     P2PP( (Anyptr, Const Anyptr, size_t) );
extern Anyptr   memset      P2PP( (Anyptr, int, size_t) );
#ifndef memcpy
extern Anyptr   memcpy      P2PP( (Anyptr, Const Anyptr, size_t) );
extern int      memcmp      P2PP( (Const Anyptr, Const Anyptr, size_t) );
#endif

extern int      atoi        P2PP( (Const Char *) );
extern LDBLE   atof        P2PP( (Const Char *) );
extern long     atol        P2PP( (Const Char *) );
extern LDBLE   strtod      P2PP( (Const Char *, Char **) );
extern long     strtol      P2PP( (Const Char *, Char **, int) );
#endif /*P2C_H_PROTO*/

/* force stdlib for DEC */
#define HAS_STDLIB
#ifndef HAS_STDLIB
extern Anyptr   malloc      P2PP( (size_t) );
extern Void     free        P2PP( (Anyptr) );
#endif

extern int      _OutMem     PV();
extern int      _CaseCheck  PV();
extern int      _NilCheck   PV();
extern int	_Escape     P2PP( (int) );
extern int	_EscIO      P2PP( (int) );

extern long     ipow        P2PP( (long, long) );
extern Char    *strsub      P2PP( (Char *, Char *, int, int) );
extern Char    *strltrim    P2PP( (Char *) );
extern Char    *strrtrim    P2PP( (Char *) );
extern Char    *strrpt      P2PP( (Char *, Char *, int) );
extern Char    *strpad      P2PP( (Char *, Char *, int, int) );
extern int      strpos2     P2PP( (Char *, Char *, int) );
extern long     memavail    PV();
extern int      P_peek      P2PP( (FILE *) );
extern int      P_eof       P2PP( (FILE *) );
extern int      P_eoln      P2PP( (FILE *) );
extern Void     P_readpaoc  P2PP( (FILE *, Char *, int) );
extern Void     P_readlnpaoc P2PP( (FILE *, Char *, int) );
extern long     P_maxpos    P2PP( (FILE *) );
extern Char    *P_trimname  P2PP( (Char *, int) );
extern long    *P_setunion  P2PP( (long *, long *, long *) );
extern long    *P_setint    P2PP( (long *, long *, long *) );
extern long    *P_setdiff   P2PP( (long *, long *, long *) );
extern long    *P_setxor    P2PP( (long *, long *, long *) );
extern int      P_inset     P2PP( (unsigned, long *) );
extern int      P_setequal  P2PP( (long *, long *) );
extern int      P_subset    P2PP( (long *, long *) );
extern long    *P_addset    P2PP( (long *, unsigned) );
extern long    *P_addsetr   P2PP( (long *, unsigned, unsigned) );
extern long    *P_remset    P2PP( (long *, unsigned) );
extern long    *P_setcpy    P2PP( (long *, long *) );
extern long    *P_expset    P2PP( (long *, long) );
extern long     P_packset   P2PP( (long *) );
extern int      P_getcmdline P2PP( (int, int, Char *) );
extern Void     TimeStamp   P2PP( (int *, int *, int *,
				 int *, int *, int *) );
extern Void	P_sun_argv  P2PP( (char *, int, int) );


/* I/O error handling */
#define _CHKIO(cond,ior,val,def)  ((cond) ? P_ioresult=0,(val)  \
					  : P_ioresult=(ior),(def))
#define _SETIO(cond,ior)          (P_ioresult = (cond) ? 0 : (ior))

/* Following defines are suitable for the HP Pascal operating system */
#define FileNotFound     10
#define FileNotOpen      13
#define FileWriteError   38
#define BadInputFormat   14
#define EndOfFile        30

#define FILENOTFOUND     10
#define FILENOTOPEN      13
#define FILEWRITEERROR   38
#define BADINPUTFORMAT   14
#define ENDOFFILE        30

/* Creating temporary files */
#if (defined(BSD) || defined(NO_TMPFILE)) && !defined(HAVE_TMPFILE)
# define tmpfile()  (fopen(tmpnam(NULL), "w+"))
#endif

/* File buffers */
#define FILEBUF(f,sc,type) sc int __CAT__(f,_BFLAGS);   \
			   sc type __CAT__(f,_BUFFER)
#define FILEBUFNC(f,type)  int __CAT__(f,_BFLAGS);   \
			   type __CAT__(f,_BUFFER)

#define RESETBUF(f,type)   (__CAT__(f,_BFLAGS) = 1)
#define SETUPBUF(f,type)   (__CAT__(f,_BFLAGS) = 0)

#define GETFBUF(f,type)    (*((__CAT__(f,_BFLAGS) == 1 &&   \
			       ((__CAT__(f,_BFLAGS) = 2),   \
				fread(&__CAT__(f,_BUFFER),  \
				      sizeof(type),1,(f)))),\
			      &__CAT__(f,_BUFFER)))
#define AGETFBUF(f,type)   ((__CAT__(f,_BFLAGS) == 1 &&   \
			     ((__CAT__(f,_BFLAGS) = 2),   \
			      fread(__CAT__(f,_BUFFER),  \
				    sizeof(type),1,(f)))),\
			    __CAT__(f,_BUFFER))

#define PUTFBUF(f,type,v)  (GETFBUF(f,type) = (v))
#define CPUTFBUF(f,v)      (PUTFBUF(f,char,v))
#define APUTFBUF(f,type,v) (memcpy(AGETFBUF(f,type), (v),  \
				   sizeof(__CAT__(f,_BUFFER))))

#define GET(f,type)        (__CAT__(f,_BFLAGS) == 1 ?   \
			    fread(&__CAT__(f,_BUFFER),sizeof(type),1,(f)) :  \
			    (__CAT__(f,_BFLAGS) = 1))

#define PUT(f,type)        (fwrite(&__CAT__(f,_BUFFER),sizeof(type),1,(f)),  \
			    (__CAT__(f,_BFLAGS) = 0))
#define CPUT(f)            (PUT(f,char))

#define BUFEOF(f)	   (__CAT__(f,_BFLAGS) != 2 && P_eof(f))
#define BUFFPOS(f)	   (ftell(f) - (__CAT__(f,_BFLAGS) == 2))

#ifdef SKIP
typedef struct {
    FILE *f;
    FILEBUFNC(f,Char);
    Char name[_FNSIZE];
} _TEXT;
#endif

/* Memory allocation */
#ifdef __GCC__
# define Malloc(n)  (PHRQ_malloc(n) ?: (Anyptr)_OutMem())
#else
extern Anyptr __MallocTemp__;
# define Malloc(n)  ((__MallocTemp__ = PHRQ_malloc(n)) ? __MallocTemp__ : (Anyptr)_OutMem())
#endif
#define FreeR(p)    (free((Anyptr)(p)))    /* used if arg is an rvalue */
#define Free(p)     (free((Anyptr)(p)), (p)=NULL)

/* sign extension */
#define SEXT(x,n)   ((x) | -(((x) & (1L<<((n)-1))) << 1))

/* packed arrays */   /* BEWARE: these are untested! */
#define P_getbits_UB(a,i,n,L)   ((int)((a)[(i)>>(L)-(n)] >>   \
				       (((~(i))&((1<<(L)-(n))-1)) << (n)) &  \
				       (1<<(1<<(n)))-1))

#define P_getbits_SB(a,i,n,L)   ((int)((a)[(i)>>(L)-(n)] <<   \
				       (16 - ((((~(i))&((1<<(L)-(n))-1))+1) <<\
					      (n)) >> (16-(1<<(n))))))

#define P_putbits_UB(a,i,x,n,L) ((a)[(i)>>(L)-(n)] |=   \
				 (x) << (((~(i))&((1<<(L)-(n))-1)) << (n)))

#define P_putbits_SB(a,i,x,n,L) ((a)[(i)>>(L)-(n)] |=   \
				 ((x) & (1<<(1<<(n)))-1) <<   \
				 (((~(i))&((1<<(L)-(n))-1)) << (n)))

#define P_clrbits_B(a,i,n,L)    ((a)[(i)>>(L)-(n)] &=   \
				 ~( ((1<<(1<<(n)))-1) <<   \
				   (((~(i))&((1<<(L)-(n))-1)) << (n))) )

/* small packed arrays */
#define P_getbits_US(v,i,n)     ((int)((v) >> ((i)<<(n)) & (1<<(1<<(n)))-1))
#define P_getbits_SS(v,i,n)     ((int)((long)(v) << (SETBITS - (((i)+1) << (n))) >> (SETBITS-(1<<(n)))))
#define P_putbits_US(v,i,x,n)   ((v) |= (x) << ((i) << (n)))
#define P_putbits_SS(v,i,x,n)   ((v) |= ((x) & (1<<(1<<(n)))-1) << ((i)<<(n)))
#define P_clrbits_S(v,i,n)      ((v) &= ~( ((1<<(1<<(n)))-1) << ((i)<<(n)) ))

#define P_max(a,b)   ((a) > (b) ? (a) : (b))
#define P_min(a,b)   ((a) < (b) ? (a) : (b))


/* Fix ANSI-isms */

#ifdef LACK_LABS
# ifndef labs
#  define labs  my_labs
   extern long my_labs P2PP( (long) );
# endif
#endif

#ifdef LACK_MEMMOVE
# ifndef memmove
#  define memmove  my_memmove
   extern Anyptr my_memmove P2PP( (Anyptr, Const Anyptr, size_t) );
# endif
#endif

#ifdef LACK_MEMCPY
# ifndef memcpy
#  define memcpy  my_memcpy
   extern Anyptr my_memcpy P2PP( (Anyptr, Const Anyptr, size_t) );
# endif
# ifndef memcmp
#  define memcmp  my_memcmp
   extern int my_memcmp P2PP( (Const Anyptr, Const Anyptr, size_t) );
# endif
# ifndef memset
#  define memset  my_memset
   extern Anyptr my_memset P2PP( (Anyptr, int, size_t) );
# endif
#endif

/* Fix toupper/tolower on Suns and other stupid BSD systems */
#ifdef toupper
# undef toupper
# undef tolower
# define toupper(c)   my_toupper(c)
# define tolower(c)   my_tolower(c)
#endif

#ifndef _toupper
# if 'A' == 65 && 'a' == 97
#  define _toupper(c)  ((c)-'a'+'A')
#  define _tolower(c)  ((c)-'A'+'a')
# else
#  ifdef toupper
#   undef toupper   /* hope these are shadowing real functions, */
#   undef tolower   /* because my_toupper calls _toupper! */
#  endif
#  define _toupper(c)  toupper(c)
#  define _tolower(c)  tolower(c)
# endif
#endif


#endif    /* P2C_H */



/* End. */


