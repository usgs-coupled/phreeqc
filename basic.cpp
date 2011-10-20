/* Output from p2c, the Pascal-to-C translator */
/* From input file "basic.p" */

typedef unsigned char boolean;
#include "Phreeqc.h"
#include "phqalloc.h"
#include "phrqproto.h"
#include "p2c.h"
#include "NameDouble.h"

/* ---------------------------------------------------------------------- */
void Phreeqc::
cmd_initialize(void)
/* ---------------------------------------------------------------------- */
{
	ENTRY item, *found_item;
	int i;
	char *token;
/*
 *   create hash table
 */
	hcreate_multi((unsigned) 2 * NCMDS, &command_hash_table);
/*
 *   fill with commands
 */
	for (i = 0; i < NCMDS; i++)
	{
		token = string_hsave(command[i].name);
		item.key = token;
		item.data = (void *) &command[i];
		found_item = hsearch_multi(command_hash_table, item, ENTER);
		if (found_item == NULL)
		{
			sprintf(error_string,
					"Hash table error in basic commands initialization.");
			error_msg(error_string, STOP);
		}
	}
	return;
}

/* ---------------------------------------------------------------------- */
void Phreeqc::
cmd_free(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   destroy hash table
 */

	hdestroy_multi(command_hash_table);
	command_hash_table = NULL;
	return;
}

int Phreeqc::
basic_compile(char *commands, void **lnbase, void **vbase, void **lpbase)
{								/*main */
	int l;
	char *ptr;

	PASCAL_MAIN(0, NULL);
	inbuf = (char *) PHRQ_calloc(max_line, sizeof(char));
	if (inbuf == NULL)
		malloc_error();
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
	exitflag = false;
	ptr = commands;
	do
	{
		TRY(try2);
		ptr = commands;
		do
		{
			if (sget_logical_line(&ptr, &l, inbuf) == EOF)
			{
				strcpy(inbuf, "bye");
			}
			parseinput(&buf);
			if (curline == 0)
			{
				stmtline = NULL;
				stmttok = buf;
				if (stmttok != NULL)
					exec();
				disposetokens(&buf);
			}
		}
		while (!(exitflag || P_eof()));
		RECOVER(try2);
		if (P_escapecode != -20)
		{
			sprintf(error_string, "%d/%d", (int) P_escapecode,
					(int) P_ioresult);
			warning_msg(error_string);
		}
		else
		{
			putchar('\n');
		}
		ENDTRY(try2);
	}
	while (!(exitflag || P_eof()));
	/*  exit(EXIT_SUCCESS); */
	PHRQ_free(inbuf);
	*lnbase = (void *) linebase;
	*vbase = (void *) varbase;
	*lpbase = (void *) loopbase;
	return (P_escapecode);
}

int Phreeqc::
basic_renumber(char *commands, void **lnbase, void **vbase, void **lpbase)
{								/*main */
	int l, i;
	char *ptr;
	PASCAL_MAIN(0, NULL);
	inbuf = (char *) PHRQ_calloc(max_line, sizeof(char));
	if (inbuf == NULL)
		malloc_error();
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
	exitflag = false;
	ptr = commands;
	do
	{
		TRY(try2);
		i = 0;
		ptr = commands;
		do
		{
			if (sget_logical_line(&ptr, &l, inbuf) == EOF)
			{
				i++;
				if (i == 1)
				{
					strcpy(inbuf, "renum");
				}
				else if (i == 2)
				{
					strcpy(inbuf, "list");
				}
				else if (i == 3)
				{
					strcpy(inbuf, "new");
				}
				else if (i == 4)
				{
					strcpy(inbuf, "bye");
				}
			}
			parseinput(&buf);
			if (curline == 0)
			{
				stmtline = NULL;
				stmttok = buf;
				if (stmttok != NULL)
					exec();
				disposetokens(&buf);
			}
		}
		while (!(exitflag || P_eof()));
		RECOVER(try2);
		if (P_escapecode != -20)
		{
			sprintf(error_string, "%d/%d", (int) P_escapecode,
					(int) P_ioresult);
			warning_msg(error_string);
		}
		else
		{
			putchar('\n');
		}
		ENDTRY(try2);
	}
	while (!(exitflag || P_eof()));
	/*  exit(EXIT_SUCCESS); */
	PHRQ_free(inbuf);
	*lnbase = (void *) linebase;
	*vbase = (void *) varbase;
	*lpbase = (void *) loopbase;

	return (P_escapecode);
}

int Phreeqc::
basic_run(char *commands, void *lnbase, void *vbase, void *lpbase)
{								/*main */
	int l;
	char *ptr;
	PASCAL_MAIN(0, NULL);
	inbuf = (char *) PHRQ_calloc(max_line, sizeof(char));
	if (inbuf == NULL)
		malloc_error();
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
	exitflag = false;
	ptr = commands;
	linebase = (linerec *) lnbase;
	varbase = (varrec *) vbase;
	loopbase = (looprec *) lpbase;
	do
	{
		TRY(try2);
		do
		{
			if (sget_logical_line(&ptr, &l, inbuf) == EOF)
			{
				strcpy(inbuf, "bye");
			}
			parseinput(&buf);
			if (curline == 0)
			{
				stmtline = NULL;
				stmttok = buf;
				if (stmttok != NULL)
					exec();
				disposetokens(&buf);
			}
		}
		while (!(exitflag || P_eof()));
		RECOVER(try2);
		if (P_escapecode != -20)
		{
			sprintf(error_string, "%d/%d", (int) P_escapecode,
					(int) P_ioresult);
			warning_msg(error_string);
		}
		else
		{
			putchar('\n');
		}
		ENDTRY(try2);
	}
	while (!(exitflag || P_eof()));

	/*  exit(EXIT_SUCCESS); */
	PHRQ_free(inbuf);
	return (P_escapecode);
}

int Phreeqc::
basic_main(char *commands)
{								/*main */
	int l;
	char *ptr;
	PASCAL_MAIN(0, NULL);
	inbuf = (char *) PHRQ_calloc(max_line, sizeof(char));
	if (inbuf == NULL)
		malloc_error();
	linebase = NULL;
	varbase = NULL;
	loopbase = NULL;
#ifdef SKIP
	printf("Chipmunk BASIC 1.0\n\n");
#endif
	exitflag = false;
	ptr = commands;
	do
	{
		TRY(try2);
		do
		{
#ifdef SKIP
			putchar('>');
#endif
			if (sget_logical_line(&ptr, &l, inbuf) == EOF)
			{
				strcpy(inbuf, "bye");
			}
			parseinput(&buf);
			if (curline == 0)
			{
				stmtline = NULL;
				stmttok = buf;
				if (stmttok != NULL)
					exec();
				disposetokens(&buf);
			}
		}
		while (!(exitflag || P_eof()));
		RECOVER(try2);
		if (P_escapecode != -20)
		{
			sprintf(error_string, "%d/%d", (int) P_escapecode,
					(int) P_ioresult);
			warning_msg(error_string);
		}
		else
		{
			putchar('\n');
		}
		ENDTRY(try2);
	}
	while (!(exitflag || P_eof()));
	return 1;
/*  exit(EXIT_SUCCESS); */
}

/* End. */
/* ---------------------------------------------------------------------- */
int Phreeqc::
sget_logical_line(char **ptr, int *l, char *return_line)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads file fp until end of line, ";", or eof
 *   stores characters in line_save
 *   reallocs line_save and line if more space is needed
 *
 *   returns:
 *	   EOF on empty line on end of file or
 *	   OK otherwise
 *	   *l returns length of line
 */
	int i;
	char c;
	i = 0;
	if (**ptr == '\0')
		return (EOF);
	for (;;)
	{
		c = **ptr;
		if (c == '\0')
			break;
		(*ptr)++;
		if (c == ';' || c == '\n')
			break;
		return_line[i++] = c;
	}
	return_line[i] = '\0';
	*l = i;
	return (1);
}
Static void Phreeqc::
restoredata(void)
{
	dataline = NULL;
	datatok = NULL;
}



Static void Phreeqc::
clearloops(void)
{
	looprec *l;

	while (loopbase != NULL)
	{
		l = loopbase->next;
		PHRQ_free(loopbase);
		loopbase = l;
	}
}



Static void Phreeqc::
clearvar(varrec * v)
{
	if (v->numdims != 0)
	{
		if (v->stringvar == 0)
		{
			PHRQ_free(v->UU.U0.arr);
			v->UU.U0.arr = NULL;
		}
		else
		{
			free_dim_stringvar(v);
		}
	}
	else if (v->stringvar && v->UU.U1.sv != NULL)
	{
		PHRQ_free(v->UU.U1.sv);
	}
	v->numdims = 0;
	if (v->stringvar)
	{
		v->UU.U1.sv = NULL;
		v->UU.U1.sval = &v->UU.U1.sv;
	}
	else
	{
		v->UU.U0.rv = 0.0;
		v->UU.U0.val = &v->UU.U0.rv;
	}
}


Static void Phreeqc::
clearvars(void)
{
	varrec *v;

	v = varbase;
	while (v != NULL)
	{
		clearvar(v);
		v = v->next;
	}
}

Static Char * Phreeqc::
numtostr(Char * Result, LDBLE n)
{
	/*string255 s; */
	char *l_s;
	long i;

	l_s = (char *) PHRQ_calloc(max_line, sizeof(char));
	if (l_s == NULL)
		malloc_error();
	l_s[max_line - 1] = '\0';
/*  if ((n != 0 && fabs(n) < 1e-2) || fabs(n) >= 1e12) { */
	if (ceil(n) == floor(n))
	{
		if (punch.high_precision == FALSE)
		{
			sprintf(l_s, "%12g", (double) n);
		}
		else
		{
			sprintf(l_s, "%19g", (double) n);
		}
	}
	else
	{
		if (punch.high_precision == FALSE)
		{
			sprintf(l_s, "%12.4e", (double) n);
		}
		else
		{
			sprintf(l_s, "%20.12e", (double) n);
		}
	}
	i = (int) strlen(l_s) + 1;
	l_s[i - 1] = '\0';
/* p2c: basic.p, line 237:
 * Note: Modification of string length may translate incorrectly [146] */
	strcpy(Result, l_s);
	free_check_null(l_s);
	return (Result);
/*  } else {
    if (punch.high_precision == FALSE) sprintf(l_s, "%30.10f", n);
      else sprintf(l_s, "%30.12f", n);
    i = strlen(l_s) + 1;
    do {
      i--;
    } while (l_s[i - 1] == '0');
    if (l_s[i - 1] == '.')
      i--;
    l_s[i] = '\0';
 * p2c: basic.p, line 248:
 * Note: Modification of string length may translate incorrectly [146] *
     return strcpy(Result, strltrim(l_s));
  } */
}

#define toklength       20


typedef long chset[9];





Static void Phreeqc::
parse(Char * l_inbuf, tokenrec ** l_buf)
{
	long i, j, begin, len, m, lp, q;
	Char token[toklength + 1] = {0};
	tokenrec *t, *tptr;
	varrec *v;
	Char ch;
	ENTRY item, *found_item;
	char *ptr;

	tptr = NULL;
	*l_buf = NULL;
	i = 1;
	lp = q = 0;
	do
	{
		ch = ' ';
		while (i <= (int) strlen(l_inbuf) && (ch == ' ' || ch == '\t'))
		{
			ch = l_inbuf[i - 1];
			i++;
		}
		if (ch != ' ')
		{
			t = (tokenrec *) PHRQ_calloc(1, sizeof(tokenrec));
			if (t == NULL)
				malloc_error();
			if (tptr == NULL)
				*l_buf = t;
			else
				tptr->next = t;
			tptr = t;
			t->next = NULL;
			switch (ch)
			{

			case '"':
			case '\'':
				q += 1;
				t->kind = tokstr;
				j = 0;
				len = (int) strlen(l_inbuf);
				begin = i;
				while (i <= len && l_inbuf[i - 1] != ch)
				{
					++j;
					++i;
				}
				if (l_inbuf[i - 1] == ch) q -= 1;
				m = 256;
				if (j + 1 > m)
					m = j + 1;
				t->UU.sp = (char *) PHRQ_calloc(m, sizeof(char));
				if (t->UU.sp == NULL)
					malloc_error();
				strncpy(t->UU.sp, l_inbuf + begin - 1, j);
				t->UU.sp[j] = '\0';
/* p2c: basic.p, line 415:
 * Note: Modification of string length may translate incorrectly [146] */
				i++;
				break;

			case '+':
				t->kind = tokplus;
				break;

			case '-':
				t->kind = tokminus;
				break;

			case '*':
				t->kind = toktimes;
				break;

			case '/':
				t->kind = tokdiv;
				break;

			case '^':
				t->kind = tokup;
				break;

			case '(':
			case '[':
				t->kind = toklp;
				lp += 1;
				break;

			case ')':
			case ']':
				t->kind = tokrp;
				lp -= 1;
				break;

			case ',':
				t->kind = tokcomma;
				break;

			case ';':
				t->kind = toksemi;
				break;

			case ':':
				t->kind = tokcolon;
				break;

			case '?':
				t->kind = tokprint;
				break;

			case '=':
				t->kind = tokeq;
				break;

			case '<':
				if (i <= (int) strlen(l_inbuf) && l_inbuf[i - 1] == '=')
				{
					t->kind = tokle;
					i++;
				}
				else if (i <= (int) strlen(l_inbuf) && l_inbuf[i - 1] == '>')
				{
					t->kind = tokne;
					i++;
				}
				else
					t->kind = toklt;
				break;

			case '>':
				if (i <= (int) strlen(l_inbuf) && l_inbuf[i - 1] == '=')
				{
					t->kind = tokge;
					i++;
				}
				else
					t->kind = tokgt;
				break;

			default:
				if (isalpha((int) ch))
				{
					i--;
					j = 0;
					token[toklength] = '\0';
					while (i <= (int) strlen(l_inbuf) &&
						   (l_inbuf[i - 1] == '$' || l_inbuf[i - 1] == '_' ||
							isalnum((int) l_inbuf[i - 1])))
					{
						if (j < toklength)
						{
							j++;
							token[j - 1] = l_inbuf[i - 1];
						}
						i++;
					}
					token[j] = '\0';
/* p2c: basic.p, line 309:
 * Note: Modification of string length may translate incorrectly [146] */
#define INT
#ifdef INT
/*
 *   Search hash list
 */
					str_tolower(token);
					item.key = token;
					item.data = NULL;
					found_item =
						hsearch_multi(command_hash_table, item, FIND);
					if (found_item != NULL)
					{
						t->kind =
							((struct key *) (found_item->data))->keycount;
						if (t->kind == tokrem)
						{
							m = (int) strlen(l_inbuf) + 1;
							if (m < 256)
								m = 256;
							t->UU.sp = (char *) PHRQ_calloc(m, sizeof(char));
							if (t->UU.sp == NULL)
								malloc_error();
							sprintf(t->UU.sp, "%.*s",
									(int) (strlen(l_inbuf) - i + 1),
									l_inbuf + i - 1);
							i = (int) strlen(l_inbuf) + 1;
						}
#endif
#ifdef LONG
						if (!strcmp(token, "and"))
							t->kind = tokand;
						else if (!strcmp(token, "or"))
							t->kind = tokor;
						else if (!strcmp(token, "xor"))
							t->kind = tokxor;
						else if (!strcmp(token, "not"))
							t->kind = toknot;
						else if (!strcmp(token, "mod"))
							t->kind = tokmod;
						else if (!strcmp(token, "sqr"))
							t->kind = toksqr;
						else if (!strcmp(token, "sqrt"))
							t->kind = toksqrt;
						else if (!strcmp(token, "ceil"))
							t->kind = tokceil;
						else if (!strcmp(token, "floor"))
							t->kind = tokfloor;
						else if (!strcmp(token, "sin"))
							t->kind = toksin;
						else if (!strcmp(token, "cos"))
							t->kind = tokcos;
						else if (!strcmp(token, "tan"))
							t->kind = toktan;
						else if (!strcmp(token, "arctan"))
							t->kind = tokarctan;
						else if (!strcmp(token, "log"))
							t->kind = toklog;
						else if (!strcmp(token, "exp"))
							t->kind = tokexp;
						else if (!strcmp(token, "abs"))
							t->kind = tokabs;
						else if (!strcmp(token, "sgn"))
							t->kind = toksgn;
						else if (!strcmp(token, "str$"))
							t->kind = tokstr_;
						else if (!strcmp(token, "val"))
							t->kind = tokval;
						else if (!strcmp(token, "chr$"))
							t->kind = tokchr_;
						else if (!strcmp(token, "eol$"))
							t->kind = tokeol_;
						else if (!strcmp(token, "asc"))
							t->kind = tokasc;
						else if (!strcmp(token, "len"))
							t->kind = toklen;
						else if (!strcmp(token, "mid$"))
							t->kind = tokmid_;
						else if (!strcmp(token, "peek"))
							t->kind = tokpeek;
						else if (!strcmp(token, "let"))
							t->kind = toklet;
						else if (!strcmp(token, "print"))
							t->kind = tokprint;
						else if (!strcmp(token, "punch"))
							t->kind = tokpunch;
#if defined PHREEQ98 
						else if (!strcmp(token, "graph_x"))
							t->kind = tokgraph_x;
						else if (!strcmp(token, "graph_y"))
							t->kind = tokgraph_y;
						else if (!strcmp(token, "graph_sy"))
							t->kind = tokgraph_sy;
#endif
						else if (!strcmp(token, "input"))
							t->kind = tokinput;
						else if (!strcmp(token, "goto"))
							t->kind = tokgoto;
						else if (!strcmp(token, "go to"))
							t->kind = tokgoto;
						else if (!strcmp(token, "if"))
							t->kind = tokif;
						else if (!strcmp(token, "end"))
							t->kind = tokend;
						else if (!strcmp(token, "stop"))
							t->kind = tokstop;
						else if (!strcmp(token, "for"))
							t->kind = tokfor;
						else if (!strcmp(token, "next"))
							t->kind = toknext;
						else if (!strcmp(token, "while"))
							t->kind = tokwhile;
						else if (!strcmp(token, "wend"))
							t->kind = tokwend;
						else if (!strcmp(token, "gosub"))
							t->kind = tokgosub;
						else if (!strcmp(token, "return"))
							t->kind = tokreturn;
						else if (!strcmp(token, "read"))
							t->kind = tokread;
						else if (!strcmp(token, "data"))
							t->kind = tokdata;
						else if (!strcmp(token, "restore"))
							t->kind = tokrestore;
						else if (!strcmp(token, "gotoxy"))
							t->kind = tokgotoxy;
						else if (!strcmp(token, "on"))
							t->kind = tokon;
						else if (!strcmp(token, "dim"))
							t->kind = tokdim;
						else if (!strcmp(token, "poke"))
							t->kind = tokpoke;
						else if (!strcmp(token, "list"))
							t->kind = toklist;
						else if (!strcmp(token, "run"))
							t->kind = tokrun;
						else if (!strcmp(token, "new"))
							t->kind = toknew;
						else if (!strcmp(token, "load"))
							t->kind = tokload;
						else if (!strcmp(token, "merge"))
							t->kind = tokmerge;
						else if (!strcmp(token, "save"))
							t->kind = toksave;
						else if (!strcmp(token, "bye"))
							t->kind = tokbye;
						else if (!strcmp(token, "quit"))
							t->kind = tokbye;
						else if (!strcmp(token, "del"))
							t->kind = tokdel;
						else if (!strcmp(token, "renum"))
							t->kind = tokrenum;
						else if (!strcmp(token, "then"))
							t->kind = tokthen;
						else if (!strcmp(token, "else"))
							t->kind = tokelse;
						else if (!strcmp(token, "to"))
							t->kind = tokto;
						else if (!strcmp(token, "step"))
							t->kind = tokstep;
						/*
						 *   dlp: added functions
						 */
						else if (!strcmp(token, "tc"))
							t->kind = toktc;
						else if (!strcmp(token, "tk"))
							t->kind = toktk;
						else if (!strcmp(token, "time"))
							t->kind = toktime;
						else if (!strcmp(token, "sim_time"))
							t->kind = toksim_time;
						else if (!strcmp(token, "total_time"))
							t->kind = toktotal_time;
						else if (!strcmp(token, "m0"))
							t->kind = tokm0;
						else if (!strcmp(token, "m"))
							t->kind = tokm;
						else if (!strcmp(token, "parm"))
							t->kind = tokparm;
						else if (!strcmp(token, "act"))
							t->kind = tokact;
						else if (!strcmp(token, "change_por"))
							t->kind = tokchange_por;
						else if (!strcmp(token, "get_por"))
							t->kind = tokget_por;
						else if (!strcmp(token, "change_surf"))
							t->kind = tokchange_surf;
						else if (!strcmp(token, "porevolume"))
							t->kind = tokporevolume;
						else if (!strcmp(token, "edl"))
							t->kind = tokedl;
						else if (!strcmp(token, "surf"))
							t->kind = toksurf;
						else if (!strcmp(token, "equi"))
							t->kind = tokequi;
						else if (!strcmp(token, "kin"))
							t->kind = tokkin;
						else if (!strcmp(token, "gas"))
							t->kind = tokgas;
						else if (!strcmp(token, "s_s"))
							t->kind = toks_s;
						else if (!strcmp(token, "misc1"))
							t->kind = tokmisc1;
						else if (!strcmp(token, "misc2"))
							t->kind = tokmisc2;
						else if (!strcmp(token, "mu"))
							t->kind = tokmu;
						else if (!strcmp(token, "osmotic"))
							t->kind = tokosmotic;
						else if (!strcmp(token, "alk"))
							t->kind = tokalk;
						else if (!strcmp(token, "lk_species"))
							t->kind = toklk_species;
						else if (!strcmp(token, "lk_named"))
							t->kind = toklk_named;
						else if (!strcmp(token, "lk_phase"))
							t->kind = toklk_phase;
						else if (!strcmp(token, "sum_species"))
							t->kind = toksum_species;
						else if (!strcmp(token, "sum_gas"))
							t->kind = toksum_gas;
						else if (!strcmp(token, "sum_s_s"))
							t->kind = toksum_s_s;
						else if (!strcmp(token, "calc_value"))
							t->kind = tokcalc_value;
						else if (!strcmp(token, "description"))
							t->kind = tokdescription;
						else if (!strcmp(token, "sys"))
							t->kind = toksys;
						else if (!strcmp(token, "instr"))
							t->kind = tokinstr;
						else if (!strcmp(token, "ltrim"))
							t->kind = tokltrim;
						else if (!strcmp(token, "rtrim"))
							t->kind = tokrtrim;
						else if (!strcmp(token, "trim"))
							t->kind = toktrim;
						else if (!strcmp(token, "pad"))
							t->kind = tokpad;
						else if (!strcmp(token, "rxn"))
							t->kind = tokrxn;
						else if (!strcmp(token, "dist"))
							t->kind = tokdist;
						else if (!strcmp(token, "mol"))
							t->kind = tokmol;
						else if (!strcmp(token, "la"))
							t->kind = tokla;
						else if (!strcmp(token, "lm"))
							t->kind = toklm;
						else if (!strcmp(token, "sr"))
							t->kind = toksr;
						else if (!strcmp(token, "step_no"))
							t->kind = tokstep_no;
						else if (!strcmp(token, "cell_no"))
							t->kind = tokcell_no;
						else if (!strcmp(token, "sim_no"))
							t->kind = toksim_no;
						else if (!strcmp(token, "si"))
							t->kind = toksi;
						else if (!strcmp(token, "tot"))
							t->kind = toktot;
						else if (!strcmp(token, "totmole"))
							t->kind = toktotmole;
						else if (!strcmp(token, "totmol"))
							t->kind = toktotmole;
						else if (!strcmp(token, "totmoles"))
							t->kind = toktotmole;
						else if (!strcmp(token, "log10"))
							t->kind = toklog10;
						else if (!strcmp(token, "put"))
							t->kind = tokput;
						else if (!strcmp(token, "get"))
							t->kind = tokget;
						else if (!strcmp(token, "exists"))
							t->kind = tokexists;
						else if (!strcmp(token, "charge_balance"))
							t->kind = tokcharge_balance;
						else if (!strcmp(token, "percent_error"))
							t->kind = tokpercent_error;
						else if (!strcmp(token, "SC"))
							t->kind = tokspcond;
						else if (!strcmp(token, "rem"))
						{
							t->kind = tokrem;
							m = strlen(l_inbuf) + 1;
							if (m < 256)
								m = 256;
							t->UU.sp = (char *) PHRQ_malloc(m);
							if (t->UU.sp == NULL)
								malloc_error();
							sprintf(t->UU.sp, "%.*s",
									(int) (strlen(l_inbuf) - i + 1),
									l_inbuf + i - 1);
							i = strlen(l_inbuf) + 1;
						}
#endif
					}
					else
					{
						t->kind = tokvar;
						v = varbase;
						while (v != NULL && strcmp(v->name, token))
							v = v->next;
						if (v == NULL)
						{
							v = (varrec *) PHRQ_calloc(1, sizeof(varrec));
							if (v == NULL)
								malloc_error();
							v->UU.U0.arr = NULL;
							v->next = varbase;
							varbase = v;
							strcpy(v->name, token);
							v->numdims = 0;
							if (token[strlen(token) - 1] == '$')
							{
								v->stringvar = true;
								v->UU.U1.sv = NULL;
								v->UU.U1.sval = &v->UU.U1.sv;
							}
							else
							{
								v->stringvar = false;
								v->UU.U0.rv = 0.0;
								v->UU.U0.val = &v->UU.U0.rv;
							}
						}
						t->UU.vp = v;
					}
				}
				else if (isdigit((int) ch) || ch == '.')
				{
					t->kind = toknum;
					i--;
					t->UU.num = strtod(&l_inbuf[i - 1], &ptr);
					if (&l_inbuf[i - 1] == ptr)
					{
						/*
						   Note: the following causes an infinite loop:
						   X = ..9
						 */
						t->kind = toksnerr;
						t->UU.snch = ch;
						i++;
						break;
					}
					i += (int) (ptr - &l_inbuf[i - 1]);
				}
				else
				{
					t->kind = toksnerr;
					t->UU.snch = ch;
				}
				break;
			}
		}
	}
	while (i <= (int) strlen(l_inbuf));
	if (q) {
		sprintf(error_string, " missing \" or \' in BASIC line\n %ld %s", curline, l_inbuf);
		error_msg(error_string, STOP);
	}
	if (lp > 0) {
		sprintf(error_string, " missing ) or ] in BASIC line\n %ld %s", curline, l_inbuf);
		error_msg(error_string, STOP);
	}
	else if (lp < 0) {
		sprintf(error_string, " missing ( or [ in BASIC line\n %ld %s", curline, l_inbuf);
		error_msg(error_string, STOP);
	}
}

#undef toklength



Static void Phreeqc::
listtokens(FILE * f, tokenrec * l_buf)
{
	boolean ltr;
	Char STR1[256] = {0};
	char *string;
	ltr = false;
	while (l_buf != NULL)
	{
		if ((l_buf->kind >= (long) toknot && l_buf->kind <= (long) tokrenum) ||
			l_buf->kind == (long) toknum || l_buf->kind == (long) tokvar ||
			l_buf->kind >= (long) toktc)
		{
			if (ltr)
				/*putc(' ', f); */
				output_temp_msg(" ");
			ltr = (boolean) (l_buf->kind != toknot);
		}
		else
			ltr = false;
		switch (l_buf->kind)
		{

		case tokvar:
			/*fputs(l_buf->UU.vp->name, f); */
			output_temp_msg(sformatf("%s", l_buf->UU.vp->name));
			break;

		case toknum:
			/*fputs(numtostr(STR1, l_buf->UU.num), f); */
			string = numtostr(STR1, l_buf->UU.num);
			string_trim(string);
			output_temp_msg(sformatf("%s", string));
			break;

		case tokstr:
			output_temp_msg(sformatf("\"%s\"", l_buf->UU.sp));
			break;

		case toksnerr:
			output_temp_msg(sformatf("{%c}", l_buf->UU.snch));
			break;

		case tokplus:
			/*putc('+', f); */
			output_temp_msg("+");
			break;

		case tokminus:
			/*putc('-', f); */
			output_temp_msg("-");
			break;

		case toktimes:
			/*putc('*', f); */
			output_temp_msg("*");
			break;

		case tokdiv:
			/*putc('/', f); */
			output_temp_msg("/");
			break;

		case tokup:
			/*putc('^', f); */
			output_temp_msg("^");
			break;

		case toklp:
			/*putc('(', f); */
			output_temp_msg("(");
			break;

		case tokrp:
			/*putc(')', f); */
			output_temp_msg(")");
			break;

		case tokcomma:
			/*putc(',', f); */
			output_temp_msg(",");
			break;

		case toksemi:
			/*putc(';', f); */
			output_temp_msg(";");
			break;

		case tokcolon:
			output_temp_msg(" : ");
			break;

		case tokeq:
			output_temp_msg(" = ");
			break;

		case toklt:
			output_temp_msg(" < ");
			break;

		case tokgt:
			output_temp_msg(" > ");
			break;

		case tokle:
			output_temp_msg(" <= ");
			break;

		case tokge:
			output_temp_msg(" >= ");
			break;

		case tokne:
			output_temp_msg(" <> ");
			break;

		case tokand:
			output_temp_msg(" AND ");
			break;

		case tokor:
			output_temp_msg(" OR ");
			break;

		case tokxor:
			output_temp_msg(" XOR ");
			break;

		case tokmod:
			output_temp_msg(" MOD ");
			break;

		case toknot:
			output_temp_msg("NOT ");
			break;

		case toksqr:
			output_temp_msg("SQR");
			break;

		case toksqrt:
			output_temp_msg("SQRT");
			break;

		case tokceil:
			output_temp_msg("CEIL");
			break;

		case tokfloor:
			output_temp_msg("FLOOR");
			break;

		case toksin:
			output_temp_msg("SIN");
			break;

		case tokcos:
			output_temp_msg("COS");
			break;

		case toktan:
			output_temp_msg("TAN");
			break;

		case tokarctan:
			output_temp_msg("ARCTAN");
			break;

		case toklog:
			output_temp_msg("LOG");
			break;

		case tokexp:
			output_temp_msg("EXP");
			break;

		case tokabs:
			output_temp_msg("ABS");
			break;

		case toksgn:
			output_temp_msg("SGN");
			break;

		case tokstr_:
			output_temp_msg("STR$");
			break;

		case tokval:
			output_temp_msg("VAL");
			break;

		case tokchr_:
			output_temp_msg("CHR$");
			break;

		case tokeol_:
			output_temp_msg("EOL$");
			break;

		case tokasc:
			output_temp_msg("ASC");
			break;

		case toklen:
			output_temp_msg("LEN");
			break;

		case tokmid_:
			output_temp_msg("MID$");
			break;

		case tokpeek:
			output_temp_msg("PEEK");
			break;

		case tokrem:
			output_temp_msg(sformatf("REM%s", l_buf->UU.sp));
			break;

		case toklet:
			output_temp_msg("LET");
			break;

		case tokprint:
			output_temp_msg("PRINT");
			break;

		case tokinput:
			output_temp_msg("INPUT");
			break;

		case tokgoto:
			output_temp_msg("GOTO");
			break;

		case tokif:
			output_temp_msg("IF");
			break;

		case tokend:
			output_temp_msg("END");
			break;

		case tokstop:
			output_temp_msg("STOP");
			break;

		case tokfor:
			output_temp_msg("FOR");
			break;

		case toknext:
			output_temp_msg("NEXT");
			break;

		case tokwhile:
			output_temp_msg("WHILE");
			break;

		case tokwend:
			output_temp_msg("WEND");
			break;

		case tokgosub:
			output_temp_msg("GOSUB");
			break;

		case tokreturn:
			output_temp_msg("RETURN");
			break;

		case tokread:
			output_temp_msg("READ");
			break;

		case tokdata:
			output_temp_msg("DATA");
			break;

		case tokrestore:
			output_temp_msg("RESTORE");
			break;

		case tokgotoxy:
			output_temp_msg("GOTOXY");
			break;

		case tokon:
			output_temp_msg("ON");
			break;

		case tokdim:
			output_temp_msg("DIM");
			break;

		case tokpoke:
			output_temp_msg("POKE");
			break;

		case toklist:
			output_temp_msg("LIST");
			break;

		case tokrun:
			output_temp_msg("RUN");
			break;

		case toknew:
			output_temp_msg("NEW");
			break;

		case tokload:
			output_temp_msg("LOAD");
			break;

		case tokmerge:
			output_temp_msg("MERGE");
			break;

		case toksave:
			output_temp_msg("SAVE");
			break;

		case tokbye:
			output_temp_msg("BYE");
			break;

		case tokdel:
			output_temp_msg("DEL");
			break;

		case tokrenum:
			output_temp_msg("RENUM");
			break;

		case tokthen:
			output_temp_msg(" THEN ");
			break;

		case tokelse:
			output_temp_msg(" ELSE ");
			break;

		case tokto:
			output_temp_msg(" TO ");
			break;

		case tokstep:
			output_temp_msg(" STEP ");
			break;

		case toktc:
			output_temp_msg("TC");
			break;

		case tokm0:
			output_temp_msg("M0");
			break;

		case tokm:
			output_temp_msg("M");
			break;

		case tokparm:
			output_temp_msg("PARM");
			break;

		case tokact:
			output_temp_msg("ACT");
			break;

		case tokchange_por:
			output_temp_msg("CHANGE_POR");
			break;

		case tokget_por:
			output_temp_msg("GET_POR");
			break;

		case tokchange_surf:
			output_temp_msg("CHANGE_SURF");
			break;

		case tokporevolume:
			output_temp_msg("POREVOLUME");
			break;

		case tokmol:
			output_temp_msg("MOL");
			break;

		case tokla:
			output_temp_msg("LA");
			break;

		case toklm:
			output_temp_msg("LM");
			break;

		case toksr:
			output_temp_msg("SR");
			break;

		case toksi:
			output_temp_msg("SI");
			break;

		case toktot:
			output_temp_msg("TOT");
			break;

		case toktotmole:
		case toktotmol:
		case toktotmoles:
			output_temp_msg("TOTMOLE");
			break;

		case toktk:
			output_temp_msg("TK");
			break;

		case toktime:
			output_temp_msg("TIME");
			break;

		case toklog10:
			output_temp_msg("LOG10");
			break;

		case toksim_time:
			output_temp_msg("SIM_TIME");
			break;

		case tokequi:
			output_temp_msg("EQUI");
			break;

		case tokgas:
			output_temp_msg("GAS");
			break;

		case tokpunch:
			output_temp_msg("PUNCH");
			break;

		case tokkin:
			output_temp_msg("KIN");
			break;

		case toks_s:
			output_temp_msg("S_S");
			break;

		case tokmu:
			output_temp_msg("MU");
			break;

		case tokosmotic:
			output_temp_msg("OSMOTIC");
			break;

		case tokalk:
			output_temp_msg("ALK");
			break;

		case toklk_species:
			output_temp_msg("LK_SPECIES");
			break;

		case toklk_named:
			output_temp_msg("LK_NAMED");
			break;

		case toklk_phase:
			output_temp_msg("LK_PHASE");
			break;

		case toksum_species:
			output_temp_msg("SUM_SPECIES");
			break;

		case toksum_gas:
			output_temp_msg("SUM_GAS");
			break;

		case toksum_s_s:
			output_temp_msg("SUM_s_s");
			break;

		case tokcalc_value:
			output_temp_msg("CALC_VALUE");
			break;

		case tokdescription:
			output_temp_msg("DESCRIPTION");
			break;

		case toksys:
			output_temp_msg("SYS");
			break;

		case tokinstr:
			output_temp_msg("INSTR");
			break;

		case tokltrim:
			output_temp_msg("LTRIM");
			break;

		case tokrtrim:
			output_temp_msg("RTRIM");
			break;

		case toktrim:
			output_temp_msg("TRIM");
			break;

		case tokpad:
			output_temp_msg("PAD");
			break;

		case tokrxn:
			output_temp_msg("RXN");
			break;

		case tokdist:
			output_temp_msg("DIST");
			break;

		case tokmisc1:
			output_temp_msg("MISC1");
			break;

		case tokmisc2:
			output_temp_msg("MISC2");
			break;

		case tokedl:
			output_temp_msg("EDL");
			break;

		case toksurf:
			output_temp_msg("SURF");
			break;

		case tokstep_no:
			output_temp_msg("STEP_NO");
			break;

		case toksim_no:
			output_temp_msg("SIM_NO");
			break;

		case toktotal_time:
			output_temp_msg("TOTAL_TIME");
			break;

		case tokput:
			output_temp_msg("PUT");
			break;

		case tokget:
			output_temp_msg("GET");
			break;

		case tokcharge_balance:
			output_temp_msg("CHARGE_BALANCE");
			break;

		case tokpercent_error:
			output_temp_msg("PERCENT_ERROR");
			break;

#if defined PHREEQ98 || defined MULTICHART
		case tokgraph_x:
			output_temp_msg("GRAPH_X");
			break;

		case tokgraph_y:
			output_temp_msg("GRAPH_Y");
			break;

		case tokgraph_sy:
			output_temp_msg("GRAPH_SY");
			break;
#endif

#if defined MULTICHART
		case tokplot_xy:
			output_temp_msg("PLOT_XY");
			break;
#endif

		case tokcell_no:
			output_temp_msg("CELL_NO");
			break;

		case tokexists:
			output_temp_msg("EXISTS");
			break;

		case toksc:
			output_temp_msg("SC");
			break;

		case tokgamma:
			output_temp_msg("GAMMA");
			break;

		case toklg:
			output_temp_msg("LG");
			break;

/* VP: Density Start */
		case tokrho:
			output_temp_msg("RHO");
			break;
/* VP: Density End */
		case tokcell_volume:
			output_temp_msg("CELL_VOLUME");
			break;
		case tokcell_pore_volume:
			output_temp_msg("CELL_PORE_VOLUME");
			break;
		case tokcell_porosity:
			output_temp_msg("CELL_POROSITY");
			break;
		case tokcell_saturation:
			output_temp_msg("CELL_SATURATION");
			break;
		case tokiso:
			output_temp_msg("ISO");
			break;
		case tokiso_unit:
			output_temp_msg("ISO_UNIT");
			break;
		case tokphase_formula:
			output_temp_msg("PHASE_FORMULA");
			break;			
		case toklist_s_s:
			output_temp_msg("LIST_S_S");
			break;
		}
		l_buf = l_buf->next;
	}
}



Static void Phreeqc::
disposetokens(tokenrec ** tok)
{
	tokenrec *tok1;

	while (*tok != NULL)
	{
		tok1 = (*tok)->next;
		if ((*tok)->kind == (long) tokrem || (*tok)->kind == (long) tokstr)
		{
			(*tok)->UU.sp = (char *) free_check_null((*tok)->UU.sp);
		}
		*tok = (tokenrec *) free_check_null(*tok);
		*tok = tok1;
	}
}



Static void Phreeqc::
parseinput(tokenrec ** l_buf)
{
	linerec *l, *l0, *l1;

	while (replace("\t", " ", inbuf));
	while (replace("\r", " ", inbuf));
	string_trim(inbuf);
	curline = 0;
	while (*inbuf != '\0' && isdigit((int) inbuf[0]))
	{
		curline = curline * 10 + inbuf[0] - 48;
		memmove(inbuf, inbuf + 1, strlen(inbuf));
	}
	parse(inbuf, l_buf);
	if (curline == 0)
		return;
	l = linebase;
	l0 = NULL;
	while (l != NULL && l->num < curline)
	{
		l0 = l;
		l = l->next;
	}
	if (l != NULL && l->num == curline)
	{
		l1 = l;
		l = l->next;
		if (l0 == NULL)
			linebase = l;
		else
			l0->next = l;
		disposetokens(&l1->txt);
		PHRQ_free(l1);
	}
	if (*l_buf != NULL)
	{
		l1 = (linerec *) PHRQ_calloc(1, sizeof(linerec));
		if (l1 == NULL)
			malloc_error();
		l1->next = l;
		if (l0 == NULL)
			linebase = l1;
		else
			l0->next = l1;
		l1->num = curline;
		l1->txt = *l_buf;
		strncpy(l1->inbuf, inbuf, MAX_LINE);
		l1->inbuf[MAX_LINE-1] = '\0';
	}
	clearloops();
	restoredata();
}

Static void Phreeqc::
errormsg(const Char * l_s)
{
	error_msg(l_s, CONTINUE);
	_Escape(42);
}


Static void Phreeqc::
snerr(const Char * l_s)
{
  char str[MAX_LENGTH] = {0};
  strcpy(str, "Syntax_error ");
  errormsg(strcat(str, l_s));
}


Static void Phreeqc::
tmerr(const Char * l_s)
{
  char str[MAX_LENGTH] = {0};
  strcpy(str, "Type mismatch error");
  errormsg(strcat(str, l_s));
}


Static void Phreeqc::
badsubscr(void)
{
	errormsg("Bad subscript");
}


Local LDBLE Phreeqc::
realfactor(struct LOC_exec *LINK)
{
	valrec n;

	n = factor(LINK);
	if (n.stringval)
		tmerr(": found characters, not a number");
	return (n.UU.val);
}

Local Char * Phreeqc::
strfactor(struct LOC_exec * LINK)
{
	valrec n;

	n = factor(LINK);
	if (!n.stringval)
		tmerr(": chemical name is not enclosed in \"  \"" );
	return (n.UU.sval);
}

Local Char * Phreeqc::
stringfactor(Char * Result, struct LOC_exec * LINK)
{
	valrec n;

	n = factor(LINK);
	if (!n.stringval)
		tmerr(": chemical name is not enclosed in \"  \"" );
	strcpy(Result, n.UU.sval);
	PHRQ_free(n.UU.sval);
	return Result;
}

Local long Phreeqc::
intfactor(struct LOC_exec *LINK)
{
	return ((long) floor(realfactor(LINK) + 0.5));
}

Local LDBLE Phreeqc::
realexpr(struct LOC_exec *LINK)
{
	valrec n;

	n = expr(LINK);
	if (n.stringval)
		tmerr(": found characters, not a number");
	return (n.UU.val);
}

Local Char * Phreeqc::
strexpr(struct LOC_exec * LINK)
{
	valrec n;

	n = expr(LINK);
	if (!n.stringval)
		tmerr(": chemical name is not enclosed in \"  \"" );
	return (n.UU.sval);
}

Local Char * Phreeqc::
stringexpr(Char * Result, struct LOC_exec * LINK)
{
	valrec n;

	n = expr(LINK);
	if (!n.stringval)
		tmerr(": chemical name is not enclosed in \"  \"" );
	strcpy(Result, n.UU.sval);
	PHRQ_free(n.UU.sval);
	return Result;
}

Local long Phreeqc::
intexpr(struct LOC_exec *LINK)
{
	return ((long) floor(realexpr(LINK) + 0.5));
}


Local void Phreeqc::
require(int k, struct LOC_exec *LINK)
{
  char str[MAX_LENGTH] = {0};
	int i;
	if (LINK->t == NULL || LINK->t->kind != k)
	{
		for (i = 0; i < NCMDS; i++)
		{
			if (command[i].keycount == k)
				break;
		}
		if (i == NCMDS)
			snerr(": missing unknown command");
		else {
			strcpy(str, ": missing ");
			snerr(strcat(str, command[i].name));
		}
	}
	LINK->t = LINK->t->next;
}


Local void Phreeqc::
skipparen(struct LOC_exec *LINK)
{
	do
	{
		if (LINK->t == NULL)
			snerr(": parenthesis missing");
		if (LINK->t->kind == tokrp || LINK->t->kind == tokcomma)
			goto _L1;
		if (LINK->t->kind == toklp)
		{
			LINK->t = LINK->t->next;
			skipparen(LINK);
		}
		LINK->t = LINK->t->next;
	}
	while (true);
  _L1:;
}


Local varrec * Phreeqc::
findvar(struct LOC_exec *LINK)
{
	varrec *v;
	long i, j, k;
	tokenrec *tok;
	long FORLIM;

	if (LINK->t == NULL || LINK->t->kind != tokvar)
		snerr(": can't find variable");
	v = LINK->t->UU.vp;
	LINK->t = LINK->t->next;
	if (LINK->t == NULL || LINK->t->kind != toklp)
	{
		if (v->numdims != 0)
			badsubscr();
		return v;
	}
	if (v->numdims == 0)
	{
		tok = LINK->t;
		i = 0;
		j = 1;
		do
		{
			if (i >= maxdims)
				badsubscr();
			LINK->t = LINK->t->next;
			skipparen(LINK);
			j *= 11;
			i++;
			v->dims[i - 1] = 11;
		}
		while (LINK->t->kind != tokrp);
		v->numdims = (char) i;
		if (v->stringvar)
		{
			v->UU.U1.sarr = (Char **) PHRQ_malloc(j * sizeof(char *));
			if (v->UU.U1.sarr == NULL)
				malloc_error();
			for (k = 0; k < j; k++)
				v->UU.U1.sarr[k] = NULL;
		}
		else
		{
			v->UU.U0.arr = (LDBLE *) PHRQ_malloc(j * sizeof(LDBLE));
			if (v->UU.U0.arr == NULL)
				malloc_error();
			for (k = 0; k < j; k++)
				v->UU.U0.arr[k] = 0.0;
		}
		LINK->t = tok;
	}
	k = 0;
	LINK->t = LINK->t->next;
	FORLIM = v->numdims;
	for (i = 1; i <= FORLIM; i++)
	{
		j = intexpr(LINK);
		if ((unsigned long) j >= (unsigned long) v->dims[i - 1])
			badsubscr();
		k = k * v->dims[i - 1] + j;
		if (i < v->numdims)
			require(tokcomma, LINK);
	}
	require(tokrp, LINK);
	if (v->stringvar)
		v->UU.U1.sval = &v->UU.U1.sarr[k];
	else
		v->UU.U0.val = &v->UU.U0.arr[k];
	return v;
}


Local valrec Phreeqc::
factor(struct LOC_exec * LINK)
{
	char string[MAX_LENGTH] = {0};
	struct solution *soln_ptr;
	int nn;
	varrec *v;
	tokenrec *facttok;
	valrec n;
	long i, j, m;
	tokenrec *tok, *tok1;
	Char *l_s;
	LDBLE l_dummy;
	int i_rate;
	union
	{
		long i;
		Char *c;
	} trick;
	struct save_values s_v, *s_v_ptr;
	int k;
	LDBLE TEMP;
	Char STR1[256] = {0}, STR2[256] = {0};
	char *elt_name, *surface_name, *mytemplate, *name;
	varrec *count_varrec = NULL, *names_varrec = NULL, *types_varrec =
		NULL, *moles_varrec = NULL;
	char **names_arg, **types_arg;
	LDBLE *moles_arg;
	int arg_num;
	LDBLE count_species;
	char *ptr, *string1, *string2;

	if (LINK->t == NULL)
		snerr(": missing variable or command");
	facttok = LINK->t;
	LINK->t = LINK->t->next;
	n.stringval = false;
	s_v.count_subscripts = 0;
	/*s_v.subscripts = (int *) PHRQ_malloc (sizeof (int)); */
	s_v.subscripts = NULL;
	switch (facttok->kind)
	{

	case toknum:
		n.UU.val = facttok->UU.num;
		break;

	case tokstr:
		n.stringval = true;
		m = (int) strlen(facttok->UU.sp) + 1;
		if (m < 256)
			m = 256;
		n.UU.sval = (char *) PHRQ_calloc(m, sizeof(char));
		if (n.UU.sval == NULL)
			malloc_error();
		strcpy(n.UU.sval, facttok->UU.sp);
		break;

	case tokvar:
		LINK->t = facttok;
		v = findvar(LINK);
		n.stringval = v->stringvar;
		if (n.stringval)
		{
			if (*v->UU.U1.sval != NULL)
			{
				m = (int) strlen(*v->UU.U1.sval) + 1;
				if (m < 256)
					m = 256;
			}
			else
			{
				m = 256;
			}
			n.UU.sval = (char *) PHRQ_calloc(m, sizeof(char));
			if (n.UU.sval == NULL)
				malloc_error();
			if (*v->UU.U1.sval != NULL)
			{
				strcpy(n.UU.sval, *v->UU.U1.sval);
			}

		}
		else
			n.UU.val = *v->UU.U0.val;
		break;

	case toklp:
		n = expr(LINK);
		require(tokrp, LINK);
		break;

	case tokminus:
		n.UU.val = -realfactor(LINK);
		break;

	case tokplus:
		n.UU.val = realfactor(LINK);
		break;

	case toknot:
		n.UU.val = ~intfactor(LINK);
		break;

	case toksqr:
		TEMP = realfactor(LINK);
		n.UU.val = TEMP * TEMP;
		break;

	case toksqrt:
		n.UU.val = sqrt(realfactor(LINK));
		break;

	case tokceil:
		n.UU.val = ceil(realfactor(LINK));
		break;

	case tokfloor:
		n.UU.val = floor(realfactor(LINK));
		break;

	case toktc:
		n.UU.val = tc_x;
		break;

	case toktk:
		n.UU.val = tc_x + 273.15;
		break;

	case toktime:
		n.UU.val = rate_time;
		break;

	case toksim_time:
		if (use.kinetics_in == FALSE)
		{
			if (state == PHAST)
			{
				n.UU.val = rate_sim_time;
			}
			else if (state == TRANSPORT)
			{
				n.UU.val = transport_step * timest;
			}
			else if (state == ADVECTION)
			{
				if (advection_kin_time_defined == TRUE)
				{
					n.UU.val = advection_step * advection_kin_time;
				}
				else
				{
					n.UU.val = advection_step;
				}
			}
			else
			{
				n.UU.val = 0;
			}
		}
		else
		{
			n.UU.val = rate_sim_time;
		}
		break;

	case toktotal_time:
		if (use.kinetics_in == FALSE)
		{
			if (state == PHAST)
			{
				n.UU.val = rate_sim_time_end;
			}
			else if (state == TRANSPORT)
			{
				n.UU.val = initial_total_time + transport_step * timest;
			}
			else if (state == ADVECTION)
			{
				n.UU.val =
					initial_total_time + advection_step * advection_kin_time;
			}
			else
			{
				n.UU.val = 0;
			}
		}
		else
		{
			n.UU.val = initial_total_time + rate_sim_time;
		}
		break;

	case tokm0:
		n.UU.val = rate_m0;
		break;

	case tokm:
		n.UU.val = rate_m;
		break;

	case tokparm:
		i_rate = intfactor(LINK);
		if (i_rate > count_rate_p || i_rate == 0)
		{
			errormsg("Parameter subscript out of range.");
		}
		n.UU.val = rate_p[i_rate - 1];
		break;

	case tokact:
		n.UU.val = activity(stringfactor(STR1, LINK));
		break;

	case tokgamma:
		n.UU.val = activity_coefficient(stringfactor(STR1, LINK));
		break;

	case toklg:
		n.UU.val = log_activity_coefficient(stringfactor(STR1, LINK));
		break;

	case tokget_por:
		i = intfactor(LINK);
		if (phast != TRUE)
		{
			if (i <= 0 || i > count_cells * (1 + stag_data->count_stag) + 1
				|| i == count_cells + 1)
			{
				/*		warning_msg("Note... no porosity for boundary solutions."); */
				n.UU.val = 0;
				break;
			}
			else
				n.UU.val = cell_data[i - 1].por;
			break;
		}
		else
		{
			n.UU.val = cell_porosity;
			break;
		}

	case tokedl:
		require(toklp, LINK);
		elt_name = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			surface_name = stringfactor(STR2, LINK);
		}
		else
		{
			surface_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = diff_layer_total(elt_name, surface_name);
		break;

	case toksurf:
		require(toklp, LINK);
		elt_name = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			surface_name = stringfactor(STR2, LINK);
		}
		else
		{
			surface_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = surf_total(elt_name, surface_name);
		break;

	case tokequi:
		n.UU.val = equi_phase(stringfactor(STR1, LINK));
		break;

	case tokkin:
		n.UU.val = kinetics_moles(stringfactor(STR1, LINK));
		break;

	case tokgas:
		n.UU.val = find_gas_comp(stringfactor(STR1, LINK));
		break;

	case toks_s:
		n.UU.val = find_s_s_comp(stringfactor(STR1, LINK));
		break;

	case tokmisc1:
		n.UU.val = find_misc1(stringfactor(STR1, LINK));
		break;

	case tokmisc2:
		n.UU.val = find_misc2(stringfactor(STR1, LINK));
		break;

	case tokmu:
		n.UU.val = mu_x;
		break;

	case tokosmotic:
		if (pitzer_model == TRUE || sit_model == TRUE)
		{
			n.UU.val = COSMOT;
		}
		else
		{
			n.UU.val = 0.0;
		}
		break;

	case tokalk:
		n.UU.val = total_alkalinity / mass_water_aq_x;
		break;

	case toklk_species:
		n.UU.val = calc_logk_s(stringfactor(STR1, LINK));
		break;

	case toklk_named:
		n.UU.val = calc_logk_n(stringfactor(STR1, LINK));
		break;

	case toklk_phase:
		n.UU.val = calc_logk_p(stringfactor(STR1, LINK));
		break;

	case toksum_species:
		require(toklp, LINK);
		mytemplate = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			elt_name = stringfactor(STR2, LINK);
		}
		else
		{
			elt_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = sum_match_species(mytemplate, elt_name);
		break;

	case toksum_gas:
		require(toklp, LINK);
		mytemplate = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			elt_name = stringfactor(STR2, LINK);
		}
		else
		{
			elt_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = sum_match_gases(mytemplate, elt_name);
		break;

	case toksum_s_s:
		require(toklp, LINK);
		mytemplate = stringfactor(STR1, LINK);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			elt_name = stringfactor(STR2, LINK);
		}
		else
		{
			elt_name = NULL;
		}
		require(tokrp, LINK);
		n.UU.val = sum_match_s_s(mytemplate, elt_name);
		break;

	case tokcalc_value:
		require(toklp, LINK);
		name = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		n.UU.val = get_calculate_value(name);
		break;

	case tokdescription:
		n.stringval = true;
		if (state == REACTION)
		{
			if (use.mix_in == TRUE)
			{
				sprintf(string, "Mix %d", use.n_mix_user);
				n.UU.sval = string_duplicate(string);
			}
			else
			{
				soln_ptr = solution_bsearch(use.n_solution_user, &nn, TRUE);
				if (soln_ptr != NULL)
				{
					n.UU.sval = string_duplicate(soln_ptr->description);
				}
				else
				{
					n.UU.sval = string_duplicate("Unknown");
				}
			}
		}
		else if (state == ADVECTION || state == TRANSPORT || state == PHAST)
		{
			sprintf(string, "Cell %d", cell_no);
			n.UU.sval = string_duplicate(string);
		}
		else
		{
			if (use.solution_ptr != NULL)
			{
				n.UU.sval = string_duplicate(use.solution_ptr->description);
			}
			else
			{
				n.UU.sval = string_duplicate("Unknown");
			}
		}
		while (replace("\t", " ", n.UU.sval));
		break;

	case tokinstr:
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokcomma, LINK);
		string2 = stringfactor(STR2, LINK);
		require(tokrp, LINK);
		ptr = strstr(string1, string2);
		if (ptr == NULL)
		{
			n.UU.val = 0;
		}
		else
		{
			n.UU.val = ((LDBLE) (ptr - string1)) + 1;
		}
		break;

	case tokltrim:
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		string_trim_left(string1);
		n.UU.sval = string_duplicate(string1);
		break;

	case tokrtrim:
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		string_trim_right(string1);
		n.UU.sval = string_duplicate(string1);
		break;

	case toktrim:
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		string_trim(string1);
		n.UU.sval = string_duplicate(string1);
		break;

	case tokiso:
		n.UU.val = iso_value(stringfactor(STR1, LINK));
		break;

	case tokiso_unit:
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokrp, LINK);
		string_trim(string1);
		n.UU.sval = iso_unit(string1);
		break;

	case tokpad:
		n.stringval = true;
		require(toklp, LINK);
		string1 = stringfactor(STR1, LINK);
		require(tokcomma, LINK);
		i = intexpr(LINK);
		require(tokrp, LINK);
		n.UU.sval = string_pad(string1, i);
		break;

	case toksys:
		require(toklp, LINK);
		elt_name = stringfactor(STR1, LINK);
		/*
		 *  Parse arguments
		 */
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			/*int c; */
			/*  struct varrec *count_varrec, *names_varrec, *types_varrec, *moles_varrec; */
			/* return number of species */
			LINK->t = LINK->t->next;
			count_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
				snerr(": can't find variable");

			/* return number of names of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			names_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || names_varrec->stringvar != 1)
				snerr(": can't find name of species");

			/* return number of types of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			types_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || types_varrec->stringvar != 1)
				snerr(": can't find type of species");

			/* return number of moles  of species */
			LINK->t = LINK->t->next;
			require(tokcomma, LINK);
			moles_varrec = LINK->t->UU.vp;
			if (LINK->t->kind != tokvar || moles_varrec->stringvar != 0)
				snerr(": can't find moles of species");
			LINK->t = LINK->t->next;
			arg_num = 4;
		}
		else
		{
			arg_num = 1;
		}
		require(tokrp, LINK);

		if (arg_num > 1)
		{
			free_dim_stringvar(names_varrec);
			free_dim_stringvar(types_varrec);
			free_check_null(moles_varrec->UU.U0.arr);
			moles_varrec->UU.U0.arr = NULL;
		}
		/*
		 *  Call subroutine
		 */
		/*
		   n.UU.val = system_total(elt_name, count_varrec->UU.U0.val, &(names_varrec->UU.U1.sarr), &(types_varrec->UU.U1.sarr), &(moles_varrec->UU.U0.arr));
		 */
		n.UU.val =
			system_total(elt_name, &count_species, &(names_arg),
						 &(types_arg), &(moles_arg));

		/*
		 *  fill in varrec structure
		 */
		if (arg_num > 1)
		{
			*count_varrec->UU.U0.val = count_species;
			names_varrec->UU.U1.sarr = names_arg;
			types_varrec->UU.U1.sarr = types_arg;
			moles_varrec->UU.U0.arr = moles_arg;

			for (i = 0; i < maxdims; i++)
			{
				names_varrec->dims[i] = 0;
				types_varrec->dims[i] = 0;
				moles_varrec->dims[i] = 0;
			}
			names_varrec->dims[0] = (long) (*count_varrec->UU.U0.val) + 1;
			types_varrec->dims[0] = (long) (*count_varrec->UU.U0.val) + 1;
			moles_varrec->dims[0] = (long) (*count_varrec->UU.U0.val) + 1;
			names_varrec->numdims = 1;
			types_varrec->numdims = 1;
			moles_varrec->numdims = 1;
		}
		else
		{
			for (i = 0; i < count_species + 1; i++)
			{
				free_check_null(names_arg[i]);
				free_check_null(types_arg[i]);
			}
			free_check_null(names_arg);
			free_check_null(types_arg);
			free_check_null(moles_arg);
		}
		break;

	case toklist_s_s:
		{
			/* list_s_s("calcite", count, name$, moles) */
			/* return total moles */
			require(toklp, LINK);
			std::string s_s_name(stringfactor(STR1, LINK));
			cxxNameDouble composition;
			/*
			*  Parse arguments
			*/
			if (LINK->t != NULL && LINK->t->kind == tokcomma)
			{
				LINK->t = LINK->t->next;
				count_varrec = LINK->t->UU.vp;
				if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
					snerr(": Cannot find count variable");

				/* return number of names of components */
				LINK->t = LINK->t->next;
				require(tokcomma, LINK);
				names_varrec = LINK->t->UU.vp;
				if (LINK->t->kind != tokvar || names_varrec->stringvar != 1)
					snerr(": Cannot find component string variable");

				/* return number of moles  of components */
				LINK->t = LINK->t->next;
				require(tokcomma, LINK);
				moles_varrec = LINK->t->UU.vp;
				if (LINK->t->kind != tokvar || moles_varrec->stringvar != 0)
					snerr(": Cannot find moles of component variable");
				LINK->t = LINK->t->next;
				arg_num = 4;
			}
			else
			{
				snerr(": Expected 4 arguments for list_s_s");
			}
			require(tokrp, LINK);

			if (arg_num > 1)
			{
				free_dim_stringvar(names_varrec);
				free_check_null(moles_varrec->UU.U0.arr);
				moles_varrec->UU.U0.arr = NULL;
			}
			/*
			*  Call subroutine
			*/
			// return total moles
			n.UU.val = list_s_s(s_s_name, composition);

			/*
			*  fill in varrec structure
			*/

			if (arg_num > 1)
			{
				size_t count = composition.size();
				*count_varrec->UU.U0.val = (LDBLE) count;
				/*
				* malloc space
				*/
				names_varrec->UU.U1.sarr = (char **) PHRQ_malloc((count + 1) * sizeof(char *));
				if (names_varrec->UU.U1.sarr == NULL)
					malloc_error();
				moles_varrec->UU.U0.arr = (LDBLE *) PHRQ_malloc((count + 1) * sizeof(LDBLE));
				if (moles_varrec->UU.U0.arr == NULL)
					malloc_error();

				// first position not used
				names_varrec->UU.U1.sarr[0] = NULL;
				moles_varrec->UU.U0.arr[0] = 0;

				// set dims for Basic array
				for (i = 0; i < maxdims; i++)
				{
					names_varrec->dims[i] = 0;
					moles_varrec->dims[i] = 0;
				}
				// set dims for first dimension and number of dims
				names_varrec->dims[0] = (long) (count + 1);
				moles_varrec->dims[0] = (long) (count + 1);
				names_varrec->numdims = 1;
				moles_varrec->numdims = 1;

				// fill in arrays
				i = 1;

				//for (cxxNameDouble::iterator it = composition.begin(); it != composition.end(); it++)
				//{
				//	names_varrec->UU.U1.sarr[i] = string_duplicate((it->first).c_str());
				//	moles_varrec->UU.U0.arr[i] = it->second;
				//	i++;
				//}

				std::vector< std::pair<std::string, LDBLE> > sort_comp = composition.sort_second();
				size_t j;
				for (j = 0; j != sort_comp.size(); j++)
				{
					names_varrec->UU.U1.sarr[i] = string_duplicate(sort_comp[j].first.c_str());
					moles_varrec->UU.U0.arr[i] = sort_comp[j].second;
					i++;
				}

			}
			break;
		}

	case tokphase_formula:
		{
			require(toklp, LINK);
			std::string phase_name(stringfactor(STR1, LINK));
			varrec *elts_varrec = NULL, *coef_varrec = NULL;
			cxxNameDouble stoichiometry;
			/*
			*  Parse arguments
			*/
			if (LINK->t != NULL && LINK->t->kind == tokcomma)
			{
				/* phase_formula("calcite", count, elt, coef) */
				/* return formula */
				/*int c; */
				/*  struct varrec *count_varrec, *names_varrec, *types_varrec, *moles_varrec; */
				/*  struct varrec *count_varrec, *elt_varrec, *coef_varrec; */
				/* return number of species */
				LINK->t = LINK->t->next;
				count_varrec = LINK->t->UU.vp;
				if (LINK->t->kind != tokvar || count_varrec->stringvar != 0)
					snerr(": Cannot find count variable");

				/* return number of names of species */
				LINK->t = LINK->t->next;
				require(tokcomma, LINK);
				elts_varrec = LINK->t->UU.vp;
				if (LINK->t->kind != tokvar || elts_varrec->stringvar != 1)
					snerr(": Cannot find element string variable");

				/* return coefficients of species */
				LINK->t = LINK->t->next;
				require(tokcomma, LINK);
				coef_varrec = LINK->t->UU.vp;
				if (LINK->t->kind != tokvar || coef_varrec->stringvar != 0)
					snerr(": Cannot find coefficient variable");
				LINK->t = LINK->t->next;
				arg_num = 4;
			}
			else
			{
				arg_num = 1;
			}
			require(tokrp, LINK);

			if (arg_num > 1)
			{
				free_dim_stringvar(elts_varrec);
				free_check_null(coef_varrec->UU.U0.arr);
				coef_varrec->UU.U0.arr = NULL;
			}
			/*
			*  Call subroutine
			*/
			std::string form = phase_formula(phase_name, stoichiometry);

			// put formula as return value
			n.stringval = true;
			n.UU.sval = string_duplicate(form.c_str());

			/*
			*  fill in varrec structure
			*/

			if (arg_num > 1)
			{
				size_t count = stoichiometry.size();
				*count_varrec->UU.U0.val = (LDBLE) count;
				/*
				* malloc space
				*/
				elts_varrec->UU.U1.sarr = (char **) PHRQ_malloc((count + 1) * sizeof(char *));
				if (elts_varrec->UU.U1.sarr == NULL)
					malloc_error();
				coef_varrec->UU.U0.arr = (LDBLE *) PHRQ_malloc((count + 1) * sizeof(LDBLE));
				if (coef_varrec->UU.U0.arr == NULL)
					malloc_error();

				// first position not used
				elts_varrec->UU.U1.sarr[0] = NULL;
				coef_varrec->UU.U0.arr[0] = 0;

				// set dims for Basic array
				for (i = 0; i < maxdims; i++)
				{
					elts_varrec->dims[i] = 0;
					coef_varrec->dims[i] = 0;
				}
				// set dims for first dimension and number of dims
				elts_varrec->dims[0] = (long) (count + 1);
				coef_varrec->dims[0] = (long) (count + 1);
				elts_varrec->numdims = 1;
				coef_varrec->numdims = 1;

				// fill in arrays
				i = 1;
				for (cxxNameDouble::iterator it = stoichiometry.begin(); it != stoichiometry.end(); it++)
				{
					elts_varrec->UU.U1.sarr[i] = string_duplicate((it->first).c_str());
					coef_varrec->UU.U0.arr[i] = it->second;
					i++;
				}

			}
			break;
		}

	case tokrxn:
		if (state == REACTION || 
			state == ADVECTION ||
			state == TRANSPORT)
		{
			n.UU.val = step_x;
		}
		else
		{
			n.UU.val = 0.0;
		}
		break;

	case tokdist:
		if (state == PHAST)
		{
			n.UU.val = 0;
		}
		else if (state == TRANSPORT)
		{
			n.UU.val = cell_data[cell - 1].mid_cell_x;
		}
		else if (state == ADVECTION)
		{
			n.UU.val = (LDBLE) use.n_solution_user;
		}
		else
		{
			n.UU.val = 0;
		}
		break;

	case tokmol:
		n.UU.val = molality(stringfactor(STR1, LINK));
		break;

	case tokla:
		n.UU.val = log_activity(stringfactor(STR1, LINK));
		break;

	case toklm:
		n.UU.val = log_molality(stringfactor(STR1, LINK));
		break;

	case toksr:
		n.UU.val = saturation_ratio(stringfactor(STR1, LINK));
		break;

	case tokstep_no:
		if (state == PHAST)
		{
			n.UU.val = 0;
		}
		else if (state == TRANSPORT)
		{
			n.UU.val = transport_step;
		}
		else if (state == ADVECTION)
		{
			n.UU.val = advection_step;
		}
		else if (state == REACTION)
		{
			n.UU.val = reaction_step;
		}
		else
		{
			n.UU.val = 0;
		}
		break;

	case tokcell_no:
		if (state == TRANSPORT)
		{
			n.UU.val = cell_no;
		}
		else if (state == PHAST)
		{
			n.UU.val = cell_no;
		}
		else if (state == ADVECTION)
		{
			n.UU.val = cell_no;
		}
		else if (state < REACTION)
		{
			n.UU.val = use.solution_ptr->n_user;
		}
		else
		{
			if (use.mix_in == TRUE)
			{
				n.UU.val = use.n_mix_user;
			}
			else
			{
				n.UU.val = use.n_solution_user;
			}
		}
		break;

	case toksim_no:
		n.UU.val = simulation;
		break;

	case tokget:
		require(toklp, LINK);

		s_v.count_subscripts = 0;
		/* get first subscript */
		if (LINK->t != NULL && LINK->t->kind != tokrp)
		{
			i = intexpr(LINK);
			if (s_v.subscripts == NULL)
			{
				s_v.subscripts = (int *) PHRQ_malloc(sizeof(int));
				if (s_v.subscripts == NULL)
					malloc_error();
			}
			s_v.subscripts =
				(int *) PHRQ_realloc(s_v.subscripts,
									 (size_t) (s_v.count_subscripts +
											   1) * sizeof(int));
			if (s_v.subscripts == NULL)
				malloc_error();
			s_v.subscripts[s_v.count_subscripts] = i;
			s_v.count_subscripts++;
		}

		/* get other subscripts */
		for (;;)
		{
			if (LINK->t != NULL && LINK->t->kind == tokcomma)
			{
				LINK->t = LINK->t->next;
				j = intexpr(LINK);
				if (s_v.subscripts == NULL)
				{
					s_v.subscripts = (int *) PHRQ_malloc(sizeof(int));
					if (s_v.subscripts == NULL)
						malloc_error();
				}
				s_v.subscripts =
					(int *) PHRQ_realloc(s_v.subscripts,
										 (size_t) (s_v.count_subscripts +
												   1) * sizeof(int));
				if (s_v.subscripts == NULL)
					malloc_error();
				s_v.subscripts[s_v.count_subscripts] = j;
				s_v.count_subscripts++;
			}
			else
			{
				/* get right parentheses */
				require(tokrp, LINK);
				break;
			}
		}
		s_v_ptr = save_values_bsearch(&s_v, &k);
		if (s_v_ptr == NULL)
		{
			n.UU.val = 0;
		}
		else
		{
			n.UU.val = s_v_ptr->value;
		}
		break;

	case tokexists:
		require(toklp, LINK);

		s_v.count_subscripts = 0;
		/* get first subscript */
		if (LINK->t != NULL && LINK->t->kind != tokrp)
		{
			i = intexpr(LINK);
			if (s_v.subscripts == NULL)
			{
				s_v.subscripts = (int *) PHRQ_malloc(sizeof(int));
				if (s_v.subscripts == NULL)
					malloc_error();
			}
			s_v.subscripts =
				(int *) PHRQ_realloc(s_v.subscripts,
									 (size_t) (s_v.count_subscripts +
											   1) * sizeof(int));
			if (s_v.subscripts == NULL)
				malloc_error();
			s_v.subscripts[s_v.count_subscripts] = i;
			s_v.count_subscripts++;
		}

		/* get other subscripts */
		for (;;)
		{
			if (LINK->t != NULL && LINK->t->kind == tokcomma)
			{
				LINK->t = LINK->t->next;
				j = intexpr(LINK);
				if (s_v.subscripts == NULL)
				{
					s_v.subscripts = (int *) PHRQ_malloc(sizeof(int));
					if (s_v.subscripts == NULL)
						malloc_error();
				}
				s_v.subscripts =
					(int *) PHRQ_realloc(s_v.subscripts,
										 (size_t) (s_v.count_subscripts +
												   1) * sizeof(int));
				if (s_v.subscripts == NULL)
					malloc_error();
				s_v.subscripts[s_v.count_subscripts] = j;
				s_v.count_subscripts++;
			}
			else
			{
				/* get right parentheses */
				require(tokrp, LINK);
				break;
			}
		}
		s_v_ptr = save_values_bsearch(&s_v, &k);
		if (s_v_ptr == NULL)
		{
			n.UU.val = 0;
		}
		else
		{
			n.UU.val = 1;
		}
		break;

	case tokcharge_balance:
		n.UU.val = cb_x;
		break;

	case tokpercent_error:
		n.UU.val = 100 * cb_x / total_ions_x;
		break;

	case toksi:
		saturation_index(stringfactor(STR1, LINK), &l_dummy, &n.UU.val);
		break;

	case toktot:
		n.UU.val = total(stringfactor(STR1, LINK));
		break;

	case toktotmole:
	case toktotmol:
	case toktotmoles:
		n.UU.val = total_mole(stringfactor(STR1, LINK));
		break;

	case tokcell_pore_volume:
	case tokporevolume:
		n.UU.val = cell_pore_volume;
		break;

/* VP : Density Start */
	case tokrho:
		n.UU.val = calc_dens();
		break;
/* VP: Density End */
	case tokcell_volume:
		n.UU.val = cell_volume;
		break;
	case tokcell_porosity:
		n.UU.val = cell_porosity;
		break;
	case tokcell_saturation:
		n.UU.val = cell_saturation;
		break;
	case toksc:
		n.UU.val = calc_SC();
		break;

	case toklog10:
		n.UU.val = log10(realfactor(LINK));
		break;

	case toksin:
		n.UU.val = sin(realfactor(LINK));
		break;

	case tokcos:
		n.UU.val = cos(realfactor(LINK));
		break;

	case toktan:
		n.UU.val = realfactor(LINK);
		n.UU.val = sin(n.UU.val) / cos(n.UU.val);
		break;

	case tokarctan:
		n.UU.val = atan(realfactor(LINK));
		break;

	case toklog:
		n.UU.val = log(realfactor(LINK));
		break;

	case tokexp:
		n.UU.val = exp(realfactor(LINK));
		break;

	case tokabs:
		n.UU.val = fabs(realfactor(LINK));
		break;

	case toksgn:
		n.UU.val = realfactor(LINK);
		n.UU.val = (n.UU.val > 0) - (n.UU.val < 0);
		break;

	case tokstr_:
		n.stringval = true;
		n.UU.sval = (char *) PHRQ_calloc(256, sizeof(char));
		if (n.UU.sval == NULL)
			malloc_error();
		numtostr(n.UU.sval, realfactor(LINK));
		break;

	case tokval:
		l_s = strfactor(LINK);
		tok1 = LINK->t;
		parse(l_s, &LINK->t);
		tok = LINK->t;
		if (tok == NULL)
			n.UU.val = 0.0;
		else
			n = expr(LINK);
		disposetokens(&tok);
		LINK->t = tok1;
		PHRQ_free(l_s);
		break;

	case tokchr_:
		n.stringval = true;
		n.UU.sval = (char *) PHRQ_calloc(256, sizeof(char));
		if (n.UU.sval == NULL)
			malloc_error();
		strcpy(n.UU.sval, " ");
		n.UU.sval[0] = (Char) intfactor(LINK);
		break;

	case tokeol_:
		n.stringval = true;
		n.UU.sval = (char *) PHRQ_calloc(256, sizeof(char));
		if (n.UU.sval == NULL)
			malloc_error();
		strcpy(n.UU.sval, "\n");
		break;

	case tokasc:
		l_s = strfactor(LINK);
		if (*l_s == '\0')
			n.UU.val = 0.0;
		else
			n.UU.val = l_s[0];
		PHRQ_free(l_s);
		break;

	case tokmid_:
		n.stringval = true;
		require(toklp, LINK);
		n.UU.sval = strexpr(LINK);
		require(tokcomma, LINK);
		i = intexpr(LINK);
		if (i < 1)
			i = 1;
		/*j = 255; */
		j = (int) strlen(n.UU.sval);
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			j = intexpr(LINK);
		}
		if (j > (int) strlen(n.UU.sval) - i + 1)
			j = (int) strlen(n.UU.sval) - i + 1;
		if (i > (int) strlen(n.UU.sval))
			*n.UU.sval = '\0';
		else
		{
			if (j + 1 > 256)
			{
				warning_msg("String too long in factor\n");
/*
      STR1 = (char *) PHRQ_realloc (STR1, j + 1);
      if (STR1 == NULL)
	malloc_error ();
*/
			}
			sprintf(STR1, "%.*s", (int) j, n.UU.sval + i - 1);
			strcpy(n.UU.sval, STR1);
		}
		require(tokrp, LINK);
		break;

	case toklen:
		l_s = strfactor(LINK);
		n.UU.val = (double) strlen(l_s);
		PHRQ_free(l_s);
		break;

	case tokpeek:
/* p2c: basic.p, line 1029: Note: Range checking is OFF [216] */
		trick.i = intfactor(LINK);
		n.UU.val = *trick.c;
/* p2c: basic.p, line 1032: Note: Range checking is ON [216] */
		break;

	default:
		snerr(": missing \" or (");
		break;
	}
	s_v.subscripts = (int *) free_check_null(s_v.subscripts);
	return n;
}

Local valrec Phreeqc::
upexpr(struct LOC_exec * LINK)
{
	valrec n, n2;

	n = factor(LINK);
	while (LINK->t != NULL && LINK->t->kind == tokup)
	{
		if (n.stringval)
			tmerr(": not a number before ^");
		LINK->t = LINK->t->next;
		n2 = upexpr(LINK);
		if (n2.stringval)
			tmerr(": not a number after ^");
		if (n.UU.val >= 0)
		{
			if (n.UU.val > 0)
			{
				n.UU.val = exp(n2.UU.val * log(n.UU.val));
			}
			continue;
		}
		if (n2.UU.val != (long) n2.UU.val)
		{
			tmerr(": negative number cannot be raised to a fractional power.");
		} 
		else
		{
			n.UU.val = exp(n2.UU.val * log(-n.UU.val));
			if (((long) n2.UU.val) & 1)
				n.UU.val = -n.UU.val;
		}
	}
	return n;
}

Local valrec Phreeqc::
term(struct LOC_exec * LINK)
{
	valrec n, n2;
	int k;

	n = upexpr(LINK);
	while (LINK->t != NULL && (unsigned long) LINK->t->kind < 32 &&
		   ((1L << ((long) LINK->t->kind)) & ((1L << ((long) toktimes)) |
											  (1L << ((long) tokdiv)) | (1L <<
																		 ((long) tokmod)))) != 0)
	{
		k = LINK->t->kind;
		LINK->t = LINK->t->next;
		n2 = upexpr(LINK);
		if (n.stringval || n2.stringval)
			tmerr(": found char, but need a number for * or /");
		if (k == tokmod)
		{
			/*      n.UU.val = (long)floor(n.UU.val + 0.5) % (long)floor(n2.UU.val + 0.5); */
			if (n.UU.val != 0)
			{
				n.UU.val =
					fabs(n.UU.val) / n.UU.val * fmod(fabs(n.UU.val) +
													 1e-14, n2.UU.val);
			}
			else
			{
				n.UU.val = 0;
			}
/* p2c: basic.p, line 1078:
 * Note: Using % for possibly-negative arguments [317] */
		}
		else if (k == toktimes)
			n.UU.val *= n2.UU.val;
		else if (n2.UU.val != 0)
		{
			n.UU.val /= n2.UU.val;
		}
		else
		{
			sprintf(error_string, "Zero divide in BASIC line\n %ld %s.\nValue set to zero.", stmtline->num, stmtline->inbuf);
			warning_msg(error_string);
			n.UU.val = 0;
		}
	}
	return n;
}

Local valrec Phreeqc::
sexpr(struct LOC_exec * LINK)
{
	valrec n, n2;
	int k, m;

	n = term(LINK);
	while (LINK->t != NULL && (unsigned long) LINK->t->kind < 32 &&
		   ((1L << ((long) LINK->t->kind)) &
			((1L << ((long) tokplus)) | (1L << ((long) tokminus)))) != 0)
	{
		k = LINK->t->kind;
		LINK->t = LINK->t->next;
		n2 = term(LINK);
		if (n.stringval != n2.stringval)
			tmerr(": found char, but need a number for + or - ");
		if (k == tokplus)
		{
			if (n.stringval)
			{
				m = (int) strlen(n.UU.sval) + (int) strlen(n2.UU.sval) + 1;
				if (m < 256)
					m = 256;

				n.UU.sval = (char *) PHRQ_realloc(n.UU.sval, (size_t) m * sizeof(char));
				if (n.UU.sval == NULL)
					malloc_error();
				strcat(n.UU.sval, n2.UU.sval);
				PHRQ_free(n2.UU.sval);
			}
			else
				n.UU.val += n2.UU.val;
		}
		else
		{
			if (n.stringval)
				tmerr(": found char, but need a number for - ");
			else
				n.UU.val -= n2.UU.val;
		}
	}
	return n;
}

Local valrec Phreeqc::
relexpr(struct LOC_exec * LINK)
{
	valrec n, n2;
	boolean f;
	int k;

	n = sexpr(LINK);
	while (LINK->t != NULL && (unsigned long) LINK->t->kind < 32 &&
		   ((1L << ((long) LINK->t->kind)) &
			((1L << ((long) tokne + 1)) - (1L << ((long) tokeq)))) != 0)
	{
		k = LINK->t->kind;
		LINK->t = LINK->t->next;
		n2 = sexpr(LINK);
		if (n.stringval != n2.stringval)
			tmerr("");
		if (n.stringval)
		{
			f = (boolean) ((!strcmp(n.UU.sval, n2.UU.sval)
							&& (unsigned long) k < 32
							&& ((1L << ((long) k)) &
								((1L << ((long) tokeq)) |
								 (1L << ((long) tokge)) | (1L <<
														   ((long) tokle))))
							!= 0) || (strcmp(n.UU.sval, n2.UU.sval) < 0
									  && (unsigned long) k < 32
									  && ((1L << ((long) k)) &
										  ((1L << ((long) toklt)) |
										   (1L << ((long) tokle)) | (1L <<
																	 ((long)
																	  tokne))))
									  != 0)
						   || (strcmp(n.UU.sval, n2.UU.sval) > 0
							   && (unsigned long) k < 32
							   && ((1L << ((long) k)) &
								   ((1L << ((long) tokgt)) |
									(1L << ((long) tokge)) | (1L <<
															  ((long)
															   tokne)))) !=
							   0));
			/* p2c: basic.p, line 2175: Note:
			 * Line breaker spent 0.0+1.00 seconds, 5000 tries on line 1518 [251] */
			PHRQ_free(n.UU.sval);
			PHRQ_free(n2.UU.sval);
		}
		else
			f = (boolean) ((n.UU.val == n2.UU.val && (unsigned long) k < 32 &&
							((1L << ((long) k)) & ((1L << ((long) tokeq)) |
												   (1L << ((long) tokge)) |
												   (1L << ((long) tokle)))) !=
							0) || (n.UU.val < n2.UU.val
								   && (unsigned long) k < 32
								   && ((1L << ((long) k)) &
									   ((1L << ((long) toklt)) |
										(1L << ((long) tokle)) | (1L <<
																  ((long)
																   tokne))))
								   != 0) || (n.UU.val > n2.UU.val
											 && (unsigned long) k < 32
											 && ((1L << ((long) k)) &
												 ((1L << ((long) tokgt)) |
												  (1L << ((long) tokge)) | (1L
																			<<
																			((long) tokne)))) != 0));
		/* p2c: basic.p, line 2175: Note:
		 * Line breaker spent 0.0+2.00 seconds, 5000 tries on line 1532 [251] */
		n.stringval = false;
		n.UU.val = f;
	}
	return n;
}

Local valrec Phreeqc::
andexpr(struct LOC_exec * LINK)
{
	valrec n, n2;

	n = relexpr(LINK);
	while (LINK->t != NULL && LINK->t->kind == tokand)
	{
		LINK->t = LINK->t->next;
		n2 = relexpr(LINK);
		if (n.stringval || n2.stringval)
			tmerr("");
		n.UU.val = ((long) n.UU.val) & ((long) n2.UU.val);
	}
	return n;
}

Local valrec Phreeqc::
expr(struct LOC_exec * LINK)
{
	valrec n, n2;
	int k;

	n = andexpr(LINK);
	while (LINK->t != NULL && (unsigned long) LINK->t->kind < 32 &&
		   ((1L << ((long) LINK->t->kind)) &
			((1L << ((long) tokor)) | (1L << ((long) tokxor)))) != 0)
	{
		k = LINK->t->kind;
		LINK->t = LINK->t->next;
		n2 = andexpr(LINK);
		if (n.stringval || n2.stringval)
			tmerr("");
		if (k == tokor)
			n.UU.val = ((long) n.UU.val) | ((long) n2.UU.val);
		else
			n.UU.val = ((long) n.UU.val) ^ ((long) n2.UU.val);
	}
	return n;
}


Local void Phreeqc::
checkextra(struct LOC_exec *LINK)
{
	if (LINK->t != NULL)
		errormsg("Extra information on line");
}


Local boolean Phreeqc::
iseos(struct LOC_exec *LINK)
{
	return ((boolean) (LINK->t == NULL || LINK->t->kind == (long) tokelse ||
					   LINK->t->kind == (long) tokcolon));
}


Local void Phreeqc::
skiptoeos(struct LOC_exec *LINK)
{
	while (!iseos(LINK))
		LINK->t = LINK->t->next;
}


Local linerec * Phreeqc::
findline(long n)
{
	linerec *l;

	l = linebase;
	while (l != NULL && l->num != n)
		l = l->next;
	return l;
}

Local linerec * Phreeqc::
mustfindline(long n)
{
	linerec *l;

	l = findline(n);
	if (l == NULL) {
		sprintf(error_string, "Undefined line %ld", n);
		errormsg(error_string);
	}
	return l;
}


Local void Phreeqc::
cmdend(struct LOC_exec *LINK)
{
	stmtline = NULL;
	LINK->t = NULL;
}


Local void Phreeqc::
cmdnew(struct LOC_exec *LINK)
{
	void *p;
	int i, k;

	cmdend(LINK);
	clearloops();
	restoredata();
	while (linebase != NULL)
	{
		p = linebase->next;
		disposetokens(&linebase->txt);
		PHRQ_free(linebase);
		linebase = (linerec *) p;
	}
	while (varbase != NULL)
	{
		p = varbase->next;
		if (varbase->stringvar)
		{
			if (varbase->numdims > 0)
			{
				k = 1;
				for (i = 0; i < varbase->numdims; i++)
				{
					k = k * (varbase->dims[i]);
				}
				for (i = 0; i < k; i++)
				{
					free_check_null(varbase->UU.U1.sarr[i]);
				}
				free_check_null(varbase->UU.U1.sarr);
			}
			else if (*varbase->UU.U1.sval != NULL)
			{
				*varbase->UU.U1.sval =
					(char *) free_check_null(*varbase->UU.U1.sval);
			}

		}
		else
		{
			free_check_null(varbase->UU.U0.arr);
			varbase->UU.U0.arr = NULL;
		}
		PHRQ_free(varbase);
		varbase = (varrec *) p;
	}
}


Local void Phreeqc::
cmdlist(struct LOC_exec *LINK)
{
	linerec *l;
	long n1, n2;

	do
	{
		n1 = 0;
		n2 = LONG_MAX;
		if (LINK->t != NULL && LINK->t->kind == toknum)
		{
			n1 = (long) LINK->t->UU.num;
			LINK->t = LINK->t->next;
			if (LINK->t == NULL || LINK->t->kind != tokminus)
				n2 = n1;
		}
		if (LINK->t != NULL && LINK->t->kind == tokminus)
		{
			LINK->t = LINK->t->next;
			if (LINK->t != NULL && LINK->t->kind == toknum)
			{
				n2 = (long) LINK->t->UU.num;
				LINK->t = LINK->t->next;
			}
			else
				n2 = LONG_MAX;
		}
		l = linebase;
		while (l != NULL && l->num <= n2)
		{
			if (l->num >= n1)
			{
				/* printf("%ld ", l->num); */
				/*	listtokens(stdout, l->txt); */
				/* putchar('\n'); */
				output_temp_msg(sformatf("%ld ", l->num));
				listtokens(NULL, l->txt);
				output_temp_msg("\n");
			}
			l = l->next;
		}
		if (!iseos(LINK))
			require(tokcomma, LINK);
	}
	while (!iseos(LINK));
}


Local void Phreeqc::
cmdload(boolean merging, Char * name, struct LOC_exec *LINK)
{
	FILE *f;
	tokenrec *l_buf;
	Char STR1[256] = {0};
	Char *TEMP;

	f = NULL;
	if (!merging)
		cmdnew(LINK);
	if (f != NULL)
	{
		sprintf(STR1, "%s.TEXT", name);
		f = freopen(STR1, "r", f);
	}
	else
	{
		sprintf(STR1, "%s.TEXT", name);
		f = fopen(STR1, "r");
	}
	if (f == NULL)
		_EscIO(FileNotFound);
	while (fgets(inbuf, 256, f) != NULL)
	{
		TEMP = strchr(inbuf, '\n');
		if (TEMP != NULL)
			*TEMP = 0;
		parseinput(&l_buf);
		if (curline == 0)
		{
			printf("Bad line in file\n");
			disposetokens(&l_buf);
		}
	}
	if (f != NULL)
		fclose(f);
	f = NULL;
	if (f != NULL)
		fclose(f);
}


Local void Phreeqc::
cmdrun(struct LOC_exec *LINK)
{
	linerec *l;
	long i;
	/*string255 s; */
	char *l_s;

	l_s = (char *) PHRQ_calloc(max_line, sizeof(char));
	if (l_s == NULL)
		malloc_error();

	l = linebase;
	if (!iseos(LINK))
	{
		if (LINK->t->kind == toknum)
			/*l = mustfindline(intexpr(LINK), LINK); */
			l = mustfindline(intexpr(LINK));
		else
		{
			stringexpr(l_s, LINK);
			i = 0;
			if (!iseos(LINK))
			{
				require(tokcomma, LINK);
				i = intexpr(LINK);
			}
			checkextra(LINK);
			cmdload(false, l_s, LINK);
			if (i == 0)
				l = linebase;
			else
			l = mustfindline(i);
		}
	}
	stmtline = l;
	LINK->gotoflag = true;
	clearvars();
	clearloops();
	restoredata();
	free_check_null(l_s);
	return;
}

/* replace basic save command with transport of rate back to calc_kinetic_rate */
Local void Phreeqc::
cmdsave(struct LOC_exec *LINK)
{
	valrec n;
	while (!iseos(LINK))
	{
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (n.stringval)
		{
			snerr(": in SAVE command");
		}
		else
		{
			rate_moles = n.UU.val;
		}
	}
}
Local void Phreeqc::
cmdput(struct LOC_exec *LINK)
{
	int j;
	struct save_values s_v;

	s_v.count_subscripts = 0;
	s_v.subscripts = (int *) PHRQ_malloc(sizeof(int));

	/* get parentheses */
	require(toklp, LINK);

	/* get first argumen */
	s_v.value = realexpr(LINK);

	for (;;)
	{
		if (LINK->t != NULL && LINK->t->kind == tokcomma)
		{
			LINK->t = LINK->t->next;
			j = intexpr(LINK);
			s_v.count_subscripts++;
			s_v.subscripts =
				(int *) PHRQ_realloc(s_v.subscripts,
									 (size_t) s_v.count_subscripts *
									 sizeof(int));
			if (s_v.subscripts == NULL)
				malloc_error();
			s_v.subscripts[s_v.count_subscripts - 1] = j;
		}
		else
		{
			/* get right parentheses */
			require(tokrp, LINK);
			break;
		}
	}

	save_values_store(&s_v);
	s_v.subscripts = (int *) free_check_null(s_v.subscripts);
}

Local void Phreeqc::
cmdchange_por(struct LOC_exec *LINK)
{
	int j;
	LDBLE TEMP;
	require(toklp, LINK);
	/* get new porosity */
	TEMP = realexpr(LINK);
	require(tokcomma, LINK);
	/* get cell_no */
	j = intexpr(LINK);
	require(tokrp, LINK);
	if (j > 0 && j <= count_cells * (1 + stag_data->count_stag) + 1
		&& j != count_cells + 1)
		cell_data[j - 1].por = TEMP;
}

Local void Phreeqc::
cmdchange_surf(struct LOC_exec *LINK)
{
/*
    change_surf("Hfo",    0.3,      "Sfo",      0,      5      )
	       (old_name, fraction, new_name, new_Dw, cell_no)
 */
	char *c1;
	int count;

	change_surf_count += 1;
	count = change_surf_count;
	if (change_surf[count - 1].next == FALSE)
		change_surf = change_surf_alloc(count + 1);

	require(toklp, LINK);
	/* get surface component name (change affects all comps of the same charge structure) */
	c1 = strexpr(LINK);
	change_surf[count - 1].comp_name = string_hsave(c1);
	PHRQ_free(c1);
	require(tokcomma, LINK);
	/* get fraction of comp to change */
	change_surf[count - 1].fraction = realexpr(LINK);
	require(tokcomma, LINK);
	/* get new surface component name */
	c1 = strexpr(LINK);
	change_surf[count - 1].new_comp_name = string_hsave(c1);
	PHRQ_free(c1);
	require(tokcomma, LINK);
	/* get new Dw (no transport if 0) */
	change_surf[count - 1].new_Dw = realexpr(LINK);
	require(tokcomma, LINK);
	/* get cell_no */
	change_surf[count - 1].cell_no = intexpr(LINK);
	require(tokrp, LINK);

	if (change_surf->cell_no == 0 || change_surf->cell_no == count_cells + 1)
		change_surf[count - 1].cell_no = -99;
}

Local void Phreeqc::
cmdbye(void)
{
	exitflag = true;
}

Local void Phreeqc::
cmddel(struct LOC_exec *LINK)
{
	linerec *l, *l0, *l1;
	long n1, n2;

	do
	{
		if (iseos(LINK))
			snerr(": no variable name after del");
		n1 = 0;
		n2 = LONG_MAX;
		if (LINK->t != NULL && LINK->t->kind == toknum)
		{
			n1 = (long) LINK->t->UU.num;
			LINK->t = LINK->t->next;
			if (LINK->t == NULL || LINK->t->kind != tokminus)
				n2 = n1;
		}
		if (LINK->t != NULL && LINK->t->kind == tokminus)
		{
			LINK->t = LINK->t->next;
			if (LINK->t != NULL && LINK->t->kind == toknum)
			{
				n2 = (long) LINK->t->UU.num;
				LINK->t = LINK->t->next;
			}
			else
				n2 = LONG_MAX;
		}
		l = linebase;
		l0 = NULL;
		while (l != NULL && l->num <= n2)
		{
			l1 = l->next;
			if (l->num >= n1)
			{
				if (l == stmtline)
				{
					cmdend(LINK);
					clearloops();
					restoredata();
				}
				if (l0 == NULL)
					linebase = l->next;
				else
					l0->next = l->next;
				disposetokens(&l->txt);
				PHRQ_free(l);
			}
			else
				l0 = l;
			l = l1;
		}
		if (!iseos(LINK))
			require(tokcomma, LINK);
	}
	while (!iseos(LINK));
}


Local void Phreeqc::
cmdrenum(struct LOC_exec *LINK)
{
	linerec *l, *l1;
	tokenrec *tok;
	long lnum, step;

	lnum = 10;
	step = 10;
	if (!iseos(LINK))
	{
		lnum = intexpr(LINK);
		if (!iseos(LINK))
		{
			require(tokcomma, LINK);
			step = intexpr(LINK);
		}
	}
	l = linebase;
	if (l == NULL)
		return;
	while (l != NULL)
	{
		l->num2 = lnum;
		lnum += step;
		l = l->next;
	}
	l = linebase;
	do
	{
		tok = l->txt;
		do
		{
			if (tok->kind == (long) tokdel || tok->kind == (long) tokrestore
				|| tok->kind == (long) toklist || tok->kind == (long) tokrun
				|| tok->kind == (long) tokelse || tok->kind == (long) tokthen
				|| tok->kind == (long) tokgosub
				|| tok->kind == (long) tokgoto)
			{
				while (tok->next != NULL && tok->next->kind == toknum)
				{
					tok = tok->next;
					lnum = (long) floor(tok->UU.num + 0.5);
					l1 = linebase;
					while (l1 != NULL && l1->num != lnum)
						l1 = l1->next;
					if (l1 == NULL)
						printf("Undefined line %ld in line %ld\n", lnum,
							   l->num2);
					else
						tok->UU.num = l1->num2;
					if (tok->next != NULL && tok->next->kind == tokcomma)
						tok = tok->next;
				}
			}
			tok = tok->next;
		}
		while (tok != NULL);
		l = l->next;
	}
	while (l != NULL);
	l = linebase;
	while (l != NULL)
	{
		l->num = l->num2;
		l = l->next;
	}
}


Local void Phreeqc::
cmdprint(struct LOC_exec *LINK)
{
	boolean semiflag;
	valrec n;
	Char STR1[256] = {0};

	semiflag = false;
	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (n.stringval)
		{
/*      fputs(n.UU.sval, stdout); */
			output_temp_msg(sformatf("%s ", n.UU.sval));
			PHRQ_free(n.UU.sval);
		}
		else
/*      printf("%s ", numtostr(STR1, n.UU.val)); */
			output_temp_msg(sformatf("%s ", numtostr(STR1, n.UU.val)));
	}
	if (!semiflag)
/*    putchar('\n');*/
		output_temp_msg("\n");
}

Local void Phreeqc::
cmdpunch(struct LOC_exec *LINK)
{
	valrec n;
	/*  Char STR1[256]; */

	while (!iseos(LINK))
	{
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (n.stringval)
		{
/*      fputs(n.UU.sval, stdout); */
			if (punch.high_precision == FALSE)
			{
				if (strlen(n.UU.sval) <= 12)
				{
					fpunchf_user(n_user_punch_index, "%12.12s\t", n.UU.sval);
				}
				else
				{
					fpunchf_user(n_user_punch_index, "%s\t", n.UU.sval);
				}
			}
			else
			{
				if (strlen(n.UU.sval) <= 20)
				{
					fpunchf_user(n_user_punch_index, "%20.20s\t", n.UU.sval);
				}
				else
				{
					fpunchf_user(n_user_punch_index, "%s\t", n.UU.sval);
				}
			}
			PHRQ_free(n.UU.sval);
		}
		else if (punch.high_precision == FALSE)
		{
			fpunchf_user(n_user_punch_index, "%12.4e\t", (double) n.UU.val);
		}
		else
		{
			fpunchf_user(n_user_punch_index, "%20.12e\t", (double) n.UU.val);
		}
		++n_user_punch_index;
	}
}

#if defined PHREEQ98 
Local void Phreeqc::
cmdgraph_x(struct LOC_exec *LINK)
{
	boolean semiflag;
	valrec n;
	Char STR1[256];
	semiflag = false;
	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (colnr == 0)
		{
			rownr++;
		}
		if (n.stringval)
		{
/*      fputs(n.UU.sval, stdout); */
			GridChar(n.UU.sval, "x");
			PHRQ_free(n.UU.sval);
		}
		else
			GridChar(numtostr(STR1, n.UU.val), "x");
		colnr++;
	}
}

Local void Phreeqc::
cmdgraph_y(struct LOC_exec *LINK)
{
	boolean semiflag;
	valrec n;
	Char STR1[256];
	semiflag = false;
	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (colnr == 0)
		{
			rownr++;
		}
		if (n.stringval)
		{
/*      fputs(n.UU.sval, stdout); */
			GridChar(n.UU.sval, "y");
			PHRQ_free(n.UU.sval);
		}
		else
			GridChar(numtostr(STR1, n.UU.val), "y");
		colnr++;
	}
}

Local void Phreeqc::
cmdgraph_sy(struct LOC_exec *LINK)
{
	boolean semiflag;
	valrec n;
	Char STR1[256];
	semiflag = false;
	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);
		if (colnr == 0)
		{
			rownr++;
		}
		if (n.stringval)
		{
/*      fputs(n.UU.sval, stdout); */
			GridChar(n.UU.sval, "s");
			PHRQ_free(n.UU.sval);
		}
		else
			GridChar(numtostr(STR1, n.UU.val), "s");
		colnr++;
	}
}
#endif

Local void Phreeqc::
cmdlet(boolean implied, struct LOC_exec *LINK)
{
	varrec *v;
	Char *old, *mynew;
	LDBLE d_value;
	LDBLE *target;
	char **starget;
	target = NULL;
	starget = NULL;
	if (implied)
		LINK->t = stmttok;
	v = findvar(LINK);
	if (v->stringvar)
	{
		starget = v->UU.U1.sval;
	}
	else
	{
		target = v->UU.U0.val;
	}
	require(tokeq, LINK);
	if (!v->stringvar)
	{
		/* in case array is used on right hand side */
		d_value = realexpr(LINK);
		v->UU.U0.val = target;
		*v->UU.U0.val = d_value;
		/*  *v->UU.U0.val = realexpr(LINK); */
		return;
	}
	mynew = strexpr(LINK);
	v->UU.U1.sval = starget;
	old = *v->UU.U1.sval;
	*v->UU.U1.sval = mynew;
	if (old != NULL)
		PHRQ_free(old);
}


Local void Phreeqc::
cmdgoto(struct LOC_exec *LINK)
{
	stmtline = mustfindline(intexpr(LINK));
	LINK->t = NULL;
	LINK->gotoflag = true;
}


Local void Phreeqc::
cmdif(struct LOC_exec *LINK)
{
	LDBLE n;
	long i;

	n = realexpr(LINK);
	require(tokthen, LINK);
	if (n == 0)
	{
		i = 0;
		do
		{
			if (LINK->t != NULL)
			{
				if (LINK->t->kind == tokif)
					i++;
				if (LINK->t->kind == tokelse)
					i--;
				LINK->t = LINK->t->next;
			}
		}
		while (LINK->t != NULL && i >= 0);
	}
	if (LINK->t != NULL && LINK->t->kind == toknum)
		cmdgoto(LINK);
	else
		LINK->elseflag = true;
}


Local void Phreeqc::
cmdelse(struct LOC_exec *LINK)
{
	LINK->t = NULL;
}


Local boolean Phreeqc::
skiploop(int up, int dn, struct LOC_exec *LINK)
{
	boolean Result;
	long i;
	linerec *saveline;

	saveline = stmtline;
	i = 0;
	do
	{
		while (LINK->t == NULL)
		{
			if (stmtline == NULL || stmtline->next == NULL)
			{
				Result = false;
				stmtline = saveline;
				goto _L1;
			}
			stmtline = stmtline->next;
			LINK->t = stmtline->txt;
		}
		if (LINK->t->kind == up)
			i++;
		if (LINK->t->kind == dn)
			i--;
		LINK->t = LINK->t->next;
	}
	while (i >= 0);
	Result = true;
  _L1:
	return Result;
}


Local void Phreeqc::
cmdfor(struct LOC_exec *LINK)
{
	looprec *l, lr;
	linerec *saveline;
	long i, j;

	lr.UU.U0.vp = findvar(LINK);
	if (lr.UU.U0.vp->stringvar)
		snerr(": error in FOR command");
	require(tokeq, LINK);
	*lr.UU.U0.vp->UU.U0.val = realexpr(LINK);
	require(tokto, LINK);
	lr.UU.U0.max = realexpr(LINK);
	if (LINK->t != NULL && LINK->t->kind == tokstep)
	{
		LINK->t = LINK->t->next;
		lr.UU.U0.step = realexpr(LINK);
	}
	else
		lr.UU.U0.step = 1.0;
	lr.homeline = stmtline;
	lr.hometok = LINK->t;
	lr.kind = forloop;
	lr.next = loopbase;
	if ((lr.UU.U0.step >= 0 && *lr.UU.U0.vp->UU.U0.val > lr.UU.U0.max) ||
		(lr.UU.U0.step <= 0 && *lr.UU.U0.vp->UU.U0.val < lr.UU.U0.max))
	{
		saveline = stmtline;
		i = 0;
		j = 0;
		do
		{
			while (LINK->t == NULL)
			{
				if (stmtline == NULL || stmtline->next == NULL)
				{
					stmtline = saveline;
					errormsg("FOR without NEXT");
				}
				stmtline = stmtline->next;
				LINK->t = stmtline->txt;
			}
			if (LINK->t->kind == tokfor)
			{
				if (LINK->t->next != NULL && LINK->t->next->kind == tokvar &&
					LINK->t->next->UU.vp == lr.UU.U0.vp)
					j++;
				else
					i++;
			}
			if (LINK->t->kind == toknext)
			{
				if (LINK->t->next != NULL && LINK->t->next->kind == tokvar &&
					LINK->t->next->UU.vp == lr.UU.U0.vp)
					j--;
				else
					i--;
			}
			LINK->t = LINK->t->next;
		}
		while (i >= 0 && j >= 0);
		skiptoeos(LINK);
		return;
	}
	l = (looprec *) PHRQ_calloc(1, sizeof(looprec));
	if (l == NULL)
		malloc_error();
	*l = lr;
	loopbase = l;
}


Local void Phreeqc::
cmdnext(struct LOC_exec *LINK)
{
	varrec *v;
	boolean found;
	looprec *l, *WITH;

	if (!iseos(LINK))
		v = findvar(LINK);
	else
		v = NULL;
	do
	{
		if (loopbase == NULL || loopbase->kind == gosubloop)
			errormsg("NEXT without FOR");
		found = (boolean) (loopbase->kind == forloop &&
						   (v == NULL || loopbase->UU.U0.vp == v));
		if (!found)
		{
			l = loopbase->next;
			PHRQ_free(loopbase);
			loopbase = l;
		}
	}
	while (!found);
	WITH = loopbase;
	*WITH->UU.U0.vp->UU.U0.val += WITH->UU.U0.step;
	if ((WITH->UU.U0.step < 0
		 || *WITH->UU.U0.vp->UU.U0.val <= WITH->UU.U0.max)
		&& (WITH->UU.U0.step > 0
			|| *WITH->UU.U0.vp->UU.U0.val >= WITH->UU.U0.max))
	{
		stmtline = WITH->homeline;
		LINK->t = WITH->hometok;
		return;
	}
	l = loopbase->next;
	PHRQ_free(loopbase);
	loopbase = l;
}


Local void Phreeqc::
cmdwhile(struct LOC_exec *LINK)
{
	looprec *l;

	l = (looprec *) PHRQ_calloc(1, sizeof(looprec));
	if (l == NULL)
		malloc_error();
	l->next = loopbase;
	loopbase = l;
	l->kind = whileloop;
	l->homeline = stmtline;
	l->hometok = LINK->t;
	if (iseos(LINK))
		return;
	if (realexpr(LINK) != 0)
		return;
	if (!skiploop(tokwhile, tokwend, LINK))
		errormsg("WHILE without WEND");
	l = loopbase->next;
	PHRQ_free(loopbase);
	loopbase = l;
	skiptoeos(LINK);
}


Local void Phreeqc::
cmdwend(struct LOC_exec *LINK)
{
	tokenrec *tok;
	linerec *tokline;
	looprec *l;
	boolean found;

	do
	{
		if (loopbase == NULL || loopbase->kind == gosubloop)
			errormsg("WEND without WHILE");
		found = (boolean) (loopbase->kind == whileloop);
		if (!found)
		{
			l = loopbase->next;
			PHRQ_free(loopbase);
			loopbase = l;
		}
	}
	while (!found);
	if (!iseos(LINK))
	{
		if (realexpr(LINK) != 0)
			found = false;
	}
	tok = LINK->t;
	tokline = stmtline;
	if (found)
	{
		stmtline = loopbase->homeline;
		LINK->t = loopbase->hometok;
		if (!iseos(LINK))
		{
			if (realexpr(LINK) == 0)
				found = false;
		}
	}
	if (found)
		return;
	LINK->t = tok;
	stmtline = tokline;
	l = loopbase->next;
	PHRQ_free(loopbase);
	loopbase = l;
}


Local void Phreeqc::
cmdgosub(struct LOC_exec *LINK)
{
	looprec *l;

	l = (looprec *) PHRQ_calloc(1, sizeof(looprec));
	if (l == NULL)
		malloc_error();
	l->next = loopbase;
	loopbase = l;
	l->kind = gosubloop;
	l->homeline = stmtline;
	l->hometok = LINK->t;
	cmdgoto(LINK);
}


Local void Phreeqc::
cmdreturn(struct LOC_exec *LINK)
{
	looprec *l;
	boolean found;

	do
	{
		if (loopbase == NULL)
			errormsg("RETURN without GOSUB");
		found = (boolean) (loopbase->kind == gosubloop);
		if (!found)
		{
			l = loopbase->next;
			PHRQ_free(loopbase);
			loopbase = l;
		}
	}
	while (!found);
	stmtline = loopbase->homeline;
	LINK->t = loopbase->hometok;
	l = loopbase->next;
	PHRQ_free(loopbase);
	loopbase = l;
	skiptoeos(LINK);
}


Local void Phreeqc::
cmdread(struct LOC_exec *LINK)
{
	varrec *v;
	tokenrec *tok;
	boolean found;

	do
	{
		v = findvar(LINK);
		tok = LINK->t;
		LINK->t = datatok;
		if (dataline == NULL)
		{
			dataline = linebase;
			LINK->t = dataline->txt;
		}
		if (LINK->t == NULL || LINK->t->kind != tokcomma)
		{
			do
			{
				while (LINK->t == NULL)
				{
					if (dataline == NULL || dataline->next == NULL)
						errormsg("Out of Data");
					dataline = dataline->next;
					LINK->t = dataline->txt;
				}
				found = (boolean) (LINK->t->kind == tokdata);
				LINK->t = LINK->t->next;
			}
			while (!found || iseos(LINK));
		}
		else
			LINK->t = LINK->t->next;
		if (v->stringvar)
		{
			if (*v->UU.U1.sval != NULL)
				*v->UU.U1.sval = (char *) free_check_null(*v->UU.U1.sval);
			*v->UU.U1.sval = strexpr(LINK);
		}
		else
			*v->UU.U0.val = realexpr(LINK);
		datatok = LINK->t;
		LINK->t = tok;
		if (!iseos(LINK))
			require(tokcomma, LINK);
	}
	while (!iseos(LINK));
}


Local void Phreeqc::
cmddata(struct LOC_exec *LINK)
{
	skiptoeos(LINK);
}


Local void Phreeqc::
cmdrestore(struct LOC_exec *LINK)
{
	if (iseos(LINK))
		restoredata();
	else
	{
		dataline = mustfindline(intexpr(LINK));
		datatok = dataline->txt;
	}
}


Local void Phreeqc::
cmdgotoxy(struct LOC_exec *LINK)
{
	intexpr(LINK);
	require(tokcomma, LINK);
}


Local void Phreeqc::
cmdon(struct LOC_exec *LINK)
{
	long i;
	looprec *l;

	i = intexpr(LINK);
	if (LINK->t != NULL && LINK->t->kind == tokgosub)
	{
		l = (looprec *) PHRQ_calloc(1, sizeof(looprec));
		if (l == NULL)
			malloc_error();
		l->next = loopbase;
		loopbase = l;
		l->kind = gosubloop;
		l->homeline = stmtline;
		l->hometok = LINK->t;
		LINK->t = LINK->t->next;
	}
	else
		require(tokgoto, LINK);
	if (i < 1)
	{
		skiptoeos(LINK);
		return;
	}
	while (i > 1 && !iseos(LINK))
	{
		require(toknum, LINK);
		if (!iseos(LINK))
			require(tokcomma, LINK);
		i--;
	}
	if (!iseos(LINK))
		cmdgoto(LINK);
}


Local void Phreeqc::
cmddim(struct LOC_exec *LINK)
{
	long i, j, k;
	varrec *v;
	boolean done;

	do
	{
		if (LINK->t == NULL || LINK->t->kind != tokvar)
			snerr(": error in DIM command");
		v = LINK->t->UU.vp;
		LINK->t = LINK->t->next;
		if (v->numdims != 0)
			errormsg("Array already dimensioned before");
		j = 1;
		i = 0;
		require(toklp, LINK);
		do
		{
			k = intexpr(LINK) + 1;
			if (k < 1)
				badsubscr();
			if (i >= maxdims)
				badsubscr();
			i++;
			v->dims[i - 1] = k;
			j *= k;
			done = (boolean) (LINK->t != NULL && LINK->t->kind == tokrp);
			if (!done)
				require(tokcomma, LINK);
		}
		while (!done);
		LINK->t = LINK->t->next;
		v->numdims = (char) i;
		if (v->stringvar)
		{
			v->UU.U1.sarr = (Char **) PHRQ_malloc(j * sizeof(char *));
			if (v->UU.U1.sarr == NULL)
				malloc_error();
			for (i = 0; i < j; i++)
				v->UU.U1.sarr[i] = NULL;
		}
		else
		{
			v->UU.U0.arr = (LDBLE *) PHRQ_malloc(j * sizeof(LDBLE));
			if (v->UU.U0.arr == NULL)
				malloc_error();
			for (i = 0; i < j; i++)
				v->UU.U0.arr[i] = 0.0;
		}
		if (!iseos(LINK))
			require(tokcomma, LINK);
	}
	while (!iseos(LINK));
}

Local void Phreeqc::
cmdpoke(struct LOC_exec *LINK)
{
	union
	{
		long i;
		Char *c;
	} trick;

/* p2c: basic.p, line 2073: Note: Range checking is OFF [216] */
	trick.i = intexpr(LINK);
	require(tokcomma, LINK);
	*trick.c = (Char) intexpr(LINK);
/* p2c: basic.p, line 2077: Note: Range checking is ON [216] */
}

Static void Phreeqc::
exec(void)
{
	struct LOC_exec V;
	Char *ioerrmsg;
	Char STR1[256] = {0};


	TRY(try1);
	do
	{
		do
		{
			V.gotoflag = false;
			V.elseflag = false;
			while (stmttok != NULL && stmttok->kind == tokcolon)
				stmttok = stmttok->next;
			V.t = stmttok;
			if (V.t != NULL)
			{
				V.t = V.t->next;
				switch (stmttok->kind)
				{

				case tokrem:
					/* blank case */
					break;

				case toklist:
					cmdlist(&V);
					break;

				case tokrun:
					cmdrun(&V);
					break;

				case toknew:
					cmdnew(&V);
					break;

				case tokload:
					cmdload(false, stringexpr(STR1, &V), &V);
					break;

				case tokmerge:
					cmdload(true, stringexpr(STR1, &V), &V);
					break;

				case toksave:
					cmdsave(&V);
					break;

				case tokbye:
					cmdbye();
					break;

				case tokdel:
					cmddel(&V);
					break;

				case tokrenum:
					cmdrenum(&V);
					break;

				case toklet:
					cmdlet(false, &V);
					break;

				case tokvar:
					cmdlet(true, &V);
					break;

				case tokprint:
					cmdprint(&V);
					break;

				case tokpunch:
					cmdpunch(&V);
					break;

				case tokput:
					cmdput(&V);
					break;

				case tokchange_por:
					cmdchange_por(&V);
					break;

				case tokchange_surf:
					cmdchange_surf(&V);
					break;

#if defined PHREEQ98 || defined MULTICHART
				case tokgraph_x:
					cmdgraph_x(&V);
					break;

				case tokgraph_y:
					cmdgraph_y(&V);
					break;

				case tokgraph_sy:
					cmdgraph_sy(&V);
					break;
#endif
#if defined MULTICHART
				case tokplot_xy:
					cmdplot_xy(&V);
					break;
#endif

				case tokinput:
					error_msg
						("Basic command INPUT is not a legal command in PHREEQC.",
						 STOP);
					break;

				case tokgoto:
					cmdgoto(&V);
					break;

				case tokif:
					cmdif(&V);
					break;

				case tokelse:
					cmdelse(&V);
					break;

				case tokend:
					cmdend(&V);
					break;

				case tokstop:
					P_escapecode = -20;
					goto _Ltry1;

				case tokfor:
					cmdfor(&V);
					break;

				case toknext:
					cmdnext(&V);
					break;

				case tokwhile:
					cmdwhile(&V);
					break;

				case tokwend:
					cmdwend(&V);
					break;

				case tokgosub:
					cmdgosub(&V);
					break;

				case tokreturn:
					cmdreturn(&V);
					break;

				case tokread:
					cmdread(&V);
					break;

				case tokdata:
					cmddata(&V);
					break;

				case tokrestore:
					cmdrestore(&V);
					break;

				case tokgotoxy:
					cmdgotoxy(&V);
					break;

				case tokon:
					cmdon(&V);
					break;

				case tokdim:
					cmddim(&V);
					break;

				case tokpoke:
					cmdpoke(&V);
					break;

				default:
					errormsg("Illegal command");
					break;
				}
			}
			if (!V.elseflag && !iseos(&V))
				checkextra(&V);
			stmttok = V.t;
		}
		while (V.t != NULL);
		if (stmtline != NULL)
		{
			if (!V.gotoflag)
				stmtline = stmtline->next;
			if (stmtline != NULL)
				stmttok = stmtline->txt;
		}
	}
	while (stmtline != NULL);
	RECOVER2(try1, _Ltry1);
	if (P_escapecode == -20)
		warning_msg("Break");
	/* printf("Break"); */
	else if (P_escapecode != 42)
	{
		switch (P_escapecode)
		{

		case -4:
			sprintf(error_string, "Integer overflow in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
			warning_msg(error_string);
			break;

		case -5:
			sprintf(error_string, "Divide by zero in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
			warning_msg(error_string);
			break;

		case -6:
			sprintf(error_string, "Real math overflow in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
			warning_msg(error_string);
			break;

		case -7:
			sprintf(error_string, "Real math underflow in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
			warning_msg(error_string);
			break;

		case -8:
		case -19:
		case -18:
		case -17:
		case -16:
		case -15:
			sprintf(error_string, "Value range error in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
			warning_msg(error_string);
			break;

		case -10:
			ioerrmsg = (char *) PHRQ_calloc(256, sizeof(char));
			if (ioerrmsg == NULL)
				malloc_error();
			sprintf(ioerrmsg, "I/O Error %d", (int) P_ioresult);
			warning_msg(ioerrmsg);
			PHRQ_free(ioerrmsg);
			break;

		default:
			if (EXCP_LINE != -1)
			{
				sprintf(error_string, "%12ld\n", EXCP_LINE);
				warning_msg(error_string);
			}
			_Escape(P_escapecode);
			break;
		}
	}
	if (stmtline != NULL)
	{
		sprintf(error_string, " in BASIC line\n %ld %s", stmtline->num, stmtline->inbuf);
		error_msg(error_string, CONTINUE);
	}
	ENDTRY(try1);
}								/*exec */
int Phreeqc::
free_dim_stringvar(varrec *l_varbase)
{
	int i, k;
	if (l_varbase->numdims > 0)
	{
		k = 1;
		for (i = 0; i < l_varbase->numdims; i++)
		{
			k = k * (l_varbase->dims[i]);
		}
		for (i = 0; i < k; i++)
		{
			free_check_null(l_varbase->UU.U1.sarr[i]);
		}
		l_varbase->UU.U1.sarr = (char **) free_check_null(l_varbase->UU.U1.sarr);
	}
	return (OK);
}
#if defined MULTICHART
Local void Phreeqc::
cmdplot_xy(struct LOC_exec *LINK)
{
	boolean semiflag;
	valrec n[2];
	Char STR[2][256];
	int i = 0;
	semiflag = false;

	while (!iseos(LINK) && i < 2)
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			i++;
			continue;
		}
		n[i] = expr(LINK);
		if (n[i].stringval)
		{
			strcpy(STR[i], n[i].UU.sval);
			PHRQ_free(n[i].UU.sval);
		}
		else
			numtostr(STR[i], n[i].UU.val);
	}

	ChartObject *chart = chart_handler.Get_current_chart();
	if (chart == NULL) return;

	std::string x_str(STR[0]), y_str(STR[1]);

// Code formerly in PlotXY, included here
	{

		bool new_sim = false, new_trans = false;

		if (chart->Get_FirstCallToUSER_GRAPH() && chart->Get_colnr() == 0)
			chart->Set_prev_sim_no(simulation);
		else
			chart->Set_prev_sim_no(simulation);

		// Add a curve if necessary
		if ((int) chart->Get_Curves().size() == chart->Get_colnr())
		{
			std::string head("");
			if (chart->Get_new_headings().size() > 0)
			{
				head = chart->Get_new_headings()[0];
				chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			}
			if (chart->Get_new_plotxy_curves().size() > 0)
			{
				// find plotxy curve definitions
				chart->Add_curve(true, head,
					chart->Get_new_plotxy_curves()[0].Get_line_w(),
					chart->Get_new_plotxy_curves()[0].Get_symbol(),
					chart->Get_new_plotxy_curves()[0].Get_symbol_size(),
					chart->Get_new_plotxy_curves()[0].Get_y_axis(),
					chart->Get_new_plotxy_curves()[0].Get_color()					
					);
				// pop plotxy curve definition
				chart->Get_new_plotxy_curves().erase(chart->Get_new_plotxy_curves().begin());
			}
			else
			{
				chart->Add_curve(true, head);
			}
			chart->Set_curve_added(true);
		}

		if (x_str.size() > 0 && y_str.size() > 0)
		{
			chart->Get_Curves()[chart->Get_colnr()]->Get_x().push_back(atof(x_str.c_str()));
			chart->Get_Curves()[chart->Get_colnr()]->Get_y().push_back(atof(y_str.c_str()));
			chart->Set_point_added(true);

			// Mark added curve for first point, might have been invisible in DefineCurves
			if (chart->Get_Curves()[chart->Get_colnr()]->Get_x().size() == 1)
				chart->Set_curve_added(true);
		}
	}
	chart->Set_colnr(chart->Get_colnr() + 1);
}
Local void Phreeqc::
cmdgraph_x(struct LOC_exec *LINK)
{
	boolean semiflag;
	valrec n;
	semiflag = false;

	ChartObject *chart = chart_handler.Get_current_chart();
	if (chart == NULL) return;

	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);

		if (n.stringval)
		{
			if (strlen(n.UU.sval) > 0)
			{
				chart->Set_graph_x(atof(n.UU.sval));
			}
			PHRQ_free(n.UU.sval);
		}
		else
			chart->Set_graph_x(n.UU.val);
	}
	if ((int) chart->Get_Curves().size() == chart->Get_colnr())
	{
		if (chart->Get_new_headings().size() > 0)
		{
			// remove x heading
			//if (chart->Get_colnr() == chart->Get_ColumnOffset())
			{
				chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			}
		}
	}
}

Local void Phreeqc::
cmdgraph_y(struct LOC_exec *LINK)
{
	boolean semiflag;
	valrec n;
	semiflag = false;

	ChartObject *chart = chart_handler.Get_current_chart();
	if (chart == NULL) return;

	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);

		// wait until x value is known before storing in curve
		if (n.stringval)
		{
			if (strlen(n.UU.sval) > 0)
			{
				chart->Get_graph_y()[chart->Get_colnr()] = atof(n.UU.sval);
			}	
			PHRQ_free(n.UU.sval);
		}
		else
		{
			chart->Get_graph_y()[chart->Get_colnr()] = n.UU.val;
		}
		chart->Set_point_added(true);

		// Add a new curve if necessary
		if ((int) chart->Get_Curves().size() == chart->Get_colnr())
		{
			//if (chart->Get_new_headings().size() > 0)
			//{
			//	// remove x heading
			//	if (chart->Get_colnr() == chart->Get_ColumnOffset())
			//	{
			//		chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			//	}
			//}
			if (chart->Get_new_headings().size() > 0)
			{
				chart->Add_curve(false, chart->Get_new_headings()[0]);
				chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			}
			else
			{
				chart->Add_curve(false);
			}
			chart->Set_curve_added(true);
		}
		chart->Set_colnr(chart->Get_colnr() + 1);
	}
}

Local void Phreeqc::
cmdgraph_sy(struct LOC_exec *LINK)
{
	boolean semiflag;
	valrec n;
	semiflag = false;

	ChartObject *chart = chart_handler.Get_current_chart();
	if (chart == NULL) return;

	while (!iseos(LINK))
	{
		semiflag = false;
		if ((unsigned long) LINK->t->kind < 32 &&
			((1L << ((long) LINK->t->kind)) &
			 ((1L << ((long) toksemi)) | (1L << ((long) tokcomma)))) != 0)
		{
			semiflag = true;
			LINK->t = LINK->t->next;
			continue;
		}
		n = expr(LINK);

		chart->Get_secondary_y()[chart->Get_colnr()] = true;

		// wait until x value is known before storing in curve
		if (n.stringval)
		{
			if (strlen(n.UU.sval) > 0)
			{
				chart->Get_graph_y()[chart->Get_colnr()] = atof(n.UU.sval);
			}
			PHRQ_free(n.UU.sval);
		}
		else
			chart->Get_graph_y()[chart->Get_colnr()] = n.UU.val;
		chart->Set_point_added(true);
		// Add a new curve if necessary
		if ((int) chart->Get_Curves().size() == chart->Get_colnr())
		{
			//if (chart->Get_new_headings().size() > 0)
			//{
			//	// remove x heading
			//	if (chart->Get_colnr() == chart->Get_ColumnOffset())
			//	{
			//		chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			//	}
			//}
			if (chart->Get_new_headings().size() > 0)
			{
				chart->Add_curve(false, chart->Get_new_headings()[0]);
				chart->Get_new_headings().erase(chart->Get_new_headings().begin());
			}
			else
			{
				chart->Add_curve(false);
			}
			chart->Set_curve_added(true);
		}
		chart->Set_colnr(chart->Get_colnr() + 1);
	}
}
#endif // MULTICHART
