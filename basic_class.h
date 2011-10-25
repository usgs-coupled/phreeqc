// basic.c -------------------------------

int basic_main(char *commands);
void cmd_initialize(void);
void cmd_free(void);
int basic_compile(char *commands, void **lnbase, void **vbase, void **lpbase);
int basic_run(char *commands, void *lnbase, void *vbase, void *lpbase);
int basic_init(void);
#ifdef PHREEQ98
void GridChar(char *s, char *a);
#endif
int sget_logical_line(char **ptr, int *l, char *return_line);
int free_dim_stringvar(varrec *varbase);
void exec(void);
int basic_renumber(char *commands, void **lnbase, void **vbase, void **lpbase);
void restoredata(void);
void clearloops(void);
void clearvar(varrec * v);
void clearvars(void);
Char * numtostr(Char * Result, LDBLE n);
void parse(Char * inbuf, tokenrec ** buf);
void listtokens(FILE * f, tokenrec * buf);
void disposetokens(tokenrec ** tok);
void parseinput(tokenrec ** buf);
void errormsg(const Char * s);
void snerr(const Char * s);
void tmerr(const Char * s);
void badsubscr(void);
LDBLE realfactor(struct LOC_exec *LINK);
Char * strfactor(struct LOC_exec * LINK);
Char *stringfactor(Char * Result, struct LOC_exec *LINK);
long intfactor(struct LOC_exec *LINK);
LDBLE realexpr(struct LOC_exec *LINK);
Char * strexpr(struct LOC_exec * LINK);
Char * stringexpr(Char * Result, struct LOC_exec * LINK);
long intexpr(struct LOC_exec *LINK);
void require(int k, struct LOC_exec *LINK);
void skipparen(struct LOC_exec *LINK);
varrec * findvar(struct LOC_exec *LINK);
valrec factor(struct LOC_exec *LINK);
valrec upexpr(struct LOC_exec * LINK);
valrec term(struct LOC_exec * LINK);
valrec sexpr(struct LOC_exec * LINK);
valrec relexpr(struct LOC_exec * LINK);
valrec andexpr(struct LOC_exec * LINK);
valrec expr(struct LOC_exec *LINK);
void checkextra(struct LOC_exec *LINK);
boolean iseos(struct LOC_exec *LINK);
void skiptoeos(struct LOC_exec *LINK);
linerec * findline(long n);
linerec * mustfindline(long n);
void cmdend(struct LOC_exec *LINK);
void cmdnew(struct LOC_exec *LINK);
void cmdlist(struct LOC_exec *LINK);
void cmdload(boolean merging, Char * name, struct LOC_exec *LINK);
void cmdrun(struct LOC_exec *LINK);
void cmdsave(struct LOC_exec *LINK);
void cmdput(struct LOC_exec *LINK);
void cmdchange_por(struct LOC_exec *LINK);
void cmdchange_surf(struct LOC_exec *LINK);
void cmdbye(void);
void cmddel(struct LOC_exec *LINK);
void cmdrenum(struct LOC_exec *LINK);
void cmdprint(struct LOC_exec *LINK);
void cmdpunch(struct LOC_exec *LINK);
#if defined PHREEQ98 || defined MULTICHART
void cmdgraph_x(struct LOC_exec *LINK);
void cmdgraph_y(struct LOC_exec *LINK);
void cmdgraph_sy(struct LOC_exec *LINK);
#endif
#if defined MULTICHART
void cmdplot_xy(struct LOC_exec *LINK);
#endif
void cmdlet(boolean implied, struct LOC_exec *LINK);
void cmdgoto(struct LOC_exec *LINK);
void cmdif(struct LOC_exec *LINK);
void cmdelse(struct LOC_exec *LINK);
boolean skiploop(int up, int dn, struct LOC_exec *LINK);
void cmdfor(struct LOC_exec *LINK);
void cmdnext(struct LOC_exec *LINK);
void cmdwhile(struct LOC_exec *LINK);
void cmdwend(struct LOC_exec *LINK);
void cmdgosub(struct LOC_exec *LINK);
void cmdreturn(struct LOC_exec *LINK);
void cmdread(struct LOC_exec *LINK);
void cmddata(struct LOC_exec *LINK);
void cmdrestore(struct LOC_exec *LINK);
void cmdgotoxy(struct LOC_exec *LINK);
void cmdon(struct LOC_exec *LINK);
void cmddim(struct LOC_exec *LINK);
void cmdpoke(struct LOC_exec *LINK);

void PASCAL_MAIN(int argc, Char **argv);
long my_labs(long x);
Anyptr my_memmove(Anyptr d, Const Anyptr s, size_t n);
Anyptr my_memcpy(Anyptr d, Const Anyptr s, size_t n);
int my_memcmp(Const Anyptr s1, Const Anyptr s2, size_t n);
Anyptr my_memset(Anyptr d, int c, size_t n);
int my_toupper(int c);
int my_tolower(int c);
long ipow(long a, long b);
char * strsub(register char *ret, register char *s, register int pos,
	   register int len);
int strpos2(char *s, register char *pat, register int pos);
int strcicmp(register char *s1, register char *s2);
char * strltrim(register char *s);
char * strrtrim(register char *s);
void strmove(register int len, register char *s, register int spos,
		register char *d, register int dpos);
void strinsert(register char *src, register char *dst, register int pos);
int P_peek(FILE * f);
int P_eof(void);
int P_eoln(FILE * f);
void P_readpaoc(FILE * f, char *s, int len);
void P_readlnpaoc(FILE * f, char *s, int len);
long P_maxpos(FILE * f);
Char * P_trimname(register Char * fn, register int len);
long memavail(void);
long maxavail(void);
long * P_setunion(register long *d, register long *s1, register long *s2);
long * P_setint(register long *d, register long *s1, register long *s2);
long * P_setdiff(register long *d, register long *s1, register long *s2);
long * P_setxor(register long *d, register long *s1, register long *s2);
long * P_addset(register long *s, register unsigned val);
long * P_addsetr(register long *s, register unsigned v1, register unsigned v2);
long * P_remset(register long *s, register unsigned val);
int P_setequal(register long *s1, register long *s2);
int P_subset(register long *s1, register long *s2);
long * P_setcpy(register long *d, register long *s);
long * P_expset(register long *d, register long s);
long P_packset(register long *s);	
int _OutMem(void);
int _CaseCheck(void);
int _NilCheck(void);
static char * _ShowEscape(char *buf, int code, int ior, char *prefix);
int _Escape(int code);
int _EscIO(int code);

// data members
/* basic.c ------------------------------- */

#ifdef PHREEQ98
int colnr, rownr;
#endif

//int n_user_punch_index;
Char *inbuf;
linerec *linebase;
varrec *varbase;
looprec *loopbase;
long curline;
linerec *stmtline, *dataline;
tokenrec *stmttok, *datatok, *buf;
boolean exitflag;
long EXCP_LINE;
HashTable *command_hash_table;
struct const_key *command;
int NCMDS;

Anyptr __MallocTemp__;
int P_argc;
char **P_argv;
int P_escapecode;
int P_ioresult;
__p2c_jmp_buf *__top_jb;

//jmp_buf mark;