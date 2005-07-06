#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] = "$Id: model.c 198 2005-03-31 18:11:06Z dlpark $";

static int initial_guesses(void);
static int revise_guesses(void);
static int remove_unstable_phases;
