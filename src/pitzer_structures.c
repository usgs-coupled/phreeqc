#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
#include "pitzer.h"

static char const svnid[] = "$Id: structures.c 269 2005-04-27 19:54:25Z dlpark $";

struct pitz_param *pitz_param_duplicate(struct pitz_param *old_ptr);
int pitz_param_copy(struct pitz_param *old_ptr, struct pitz_param *new_ptr);



/* ---------------------------------------------------------------------- */
int pitzer_clean_up(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Free all allocated memory, except strings
 */
	int i;

	if (svnid == NULL) fprintf(stderr," ");
	for (i = 0; i < count_pitz_param; i++) {
		pitz_params[i] = free_check_null(pitz_params[i]);
	}
	pitz_params = free_check_null(pitz_params);
	return OK;
}
/* **********************************************************************
 *
 *   Routines related to structure "pitz_param"
 *
 * ********************************************************************** */
/* ---------------------------------------------------------------------- */
struct pitz_param *pitz_param_alloc (void)
/* ---------------------------------------------------------------------- */
{
	struct pitz_param *pitz_param_ptr;
	pitz_param_ptr = (struct pitz_param *) PHRQ_malloc(sizeof (struct pitz_param));
	if (pitz_param_ptr == NULL) malloc_error();
	return ( pitz_param_ptr );
}
/* ---------------------------------------------------------------------- */
int pitz_param_free (struct pitz_param *pitz_param_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees all data associated with pitz_param structure.
 */
	if (pitz_param_ptr == NULL) return(ERROR);
/*
 *   Free space allocated for pitz_param structure
 */
	return(OK);
}
/* ---------------------------------------------------------------------- */
int pitz_param_init (struct pitz_param *pitz_param_ptr)
/* ---------------------------------------------------------------------- */
{
	int i;
/*
 *   Frees all data associated with pitz_param structure.
 */

	if (pitz_param_ptr == NULL) return(ERROR);
	pitz_param_ptr->species[1] = NULL;
	pitz_param_ptr->species[2] = NULL;
	pitz_param_ptr->species[3] = NULL;
	pitz_param_ptr->type = TYPE_Other;
	pitz_param_ptr->U.b0 = 0.0;
	for (i = 0; i < 5; i++) {
		pitz_param_ptr->a[i] = 0.0;
	}
	return(OK);
}
/* ---------------------------------------------------------------------- */
struct pitz_param *pitz_param_read (char *string, int n)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read pitzer parameter info from string
 *   n is number of species (character values)
 *          
 */
	int l, i, j;
	char *ptr;
	char token[2*MAX_LENGTH];
	struct pitz_param pzp, *pzp_ptr;

	if (n != 2 && n != 3) return (NULL);
	if (string == NULL) return (NULL);

	ptr = string;
	if (copy_token(token, &ptr, &l) == EMPTY) return(NULL);
	for (i = 0; i < n; i++) {
		if (copy_token(token, &ptr, &l) == EMPTY) return(NULL);
		pzp.species[i] = string_hsave(token);
	}
	l = 0;
	for (i = 0; i < 5; i++) {
		if (copy_token(token, &ptr, &l) == EMPTY) return(NULL);
		j=sscanf(token,"%le", &pzp.a[i]);
		if (j <= 0) break;
		l++;
	}
	if (l <= 0) return(NULL);
	pzp_ptr = pitz_param_duplicate(&pzp);
	return ( pzp_ptr);
}
/* ---------------------------------------------------------------------- */
struct pitz_param *pitz_param_duplicate(struct pitz_param *old_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Allocates space and makes duplicate copy of pitz_param structure
 */
	struct pitz_param *new_ptr;

	new_ptr = pitz_param_alloc();
	pitz_param_init(new_ptr);
/*
 *   Copy data
 */
	pitz_param_copy(old_ptr, new_ptr);
	return(new_ptr);
}
/* ---------------------------------------------------------------------- */
int pitz_param_copy(struct pitz_param *old_ptr, 
		    struct pitz_param *new_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Copies pitz_param data from old_ptr to new location, new_ptr.
 *   Space for the new_ptr structure must already be malloced.
 */
/*
 *   Store data for structure pitz_param
 */
	memcpy(new_ptr, old_ptr, sizeof(struct pitz_param));
	return(OK);
}
/* ---------------------------------------------------------------------- */
int pitz_param_search(struct pitz_param *pzp_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *  Does linear search of pitz_params for same type and species
 *  Returns -1 if not found, index number in pitz_params if found
 */
	int i, j;
	if (pzp_ptr == NULL) return -1;
	if (pzp_ptr->type == TYPE_Other) return -1;
	for (i = 0; i < count_pitz_param; i++) {
		if (pitz_params[i]->type != pzp_ptr->type) continue;
		for (j = 0; j < 3; j++) {
			if (pitz_params[i]->species[j] != pzp_ptr->species[j]) continue;
		}
		break;
	}
	if (i >= count_pitz_param) {
		return -1;
	} 
	return i;
}
