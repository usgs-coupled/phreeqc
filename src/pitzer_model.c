#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] = "$Id: model.c 198 2005-03-31 18:11:06Z dlpark $";

static int initial_guesses(void);
static int revise_guesses(void);

/* ---------------------------------------------------------------------- */
int set_pz(int initial)
/* ---------------------------------------------------------------------- */
{
/*
 *   Sets initial guesses for unknowns if initial == TRUE
 *   Revises guesses whether initial is true or not
 */
	int i;
	struct solution *solution_ptr;
/*
 *   Set initial log concentrations to zero
 */
	iterations = -1;
	solution_ptr = use.solution_ptr;
	for (i=0; i < count_s_x; i++) {
		s_x[i]->lm = LOG_ZERO_MOLALITY;
		s_x[i]->lg = 0.0;
	}
/*
 *   Set master species activities
 */

	tc_x=solution_ptr->tc;
	tk_x=tc_x+273.15;
/*
 *   H+, e-, H2O
 */
	mass_water_aq_x = solution_ptr->mass_water;
	mu_x = solution_ptr->mu;
	s_h2o->moles = mass_water_aq_x/gfw_water;
	s_h2o->la = log10(solution_ptr->ah2o);
	s_hplus->la = - solution_ptr->ph;
	s_hplus->lm = s_hplus->la;
	s_hplus->moles = exp(s_hplus->lm * LOG_10)*mass_water_aq_x;
	s_eminus->la= - solution_ptr->solution_pe;
	if (initial == TRUE) initial_guesses(); 
	if (diffuse_layer_x == TRUE) initial_surface_water();
	revise_guesses();
	return(OK);
}
/* ---------------------------------------------------------------------- */
int initial_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Make initial guesses for activities of master species and
 *   ionic strength
 */
	int i;
	struct solution *solution_ptr;

	solution_ptr = use.solution_ptr;
	mu_x = s_hplus->moles + exp((solution_ptr->ph - 14.) * LOG_10) * mass_water_aq_x;
	mu_x /= mass_water_aq_x;
	s_h2o->la=0.0;
	for ( i=0; i < count_unknowns; i++ ) {
		if (x[i] == ph_unknown || x[i] == pe_unknown ) continue;
		if (x[i]->type < CB) {
			mu_x += x[i]->moles / mass_water_aq_x * 0.5 * x[i]->master[0]->s->z *
				x[i]->master[0]->s->z;
			x[i]->master[0]->s->la = log10(x[i]->moles/mass_water_aq_x);
		} else if (x[i]->type == CB) {
			x[i]->master[0]->s->la = log10(0.001 * x[i]->moles/mass_water_aq_x);
		} else if (x[i]->type == SOLUTION_PHASE_BOUNDARY) {
			x[i]->master[0]->s->la = log10(0.001 * x[i]->moles/mass_water_aq_x);
		} else if (x[i]->type == EXCH) {
			if (x[i]->moles <= 0) {
				x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
			} else {
				x[i]->master[0]->s->la = log10(x[i]->moles);
			}
		} else if (x[i]->type == SURFACE) {
			if (x[i]->moles <= 0) {
				x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
			} else {
				x[i]->master[0]->s->la = log10(0.1 * x[i]->moles);
			}
		} else if (x[i]->type == SURFACE_CB) {
			x[i]->master[0]->s->la = 0.0;
		}
	}
	return(OK);
}
/* ---------------------------------------------------------------------- */
int revise_guesses(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Revise molalities species
 */
	int i;
	int iter, max_iter, repeat, fail;
	LDBLE weight, f;

	max_iter = 10;
	/* gammas(mu_x);*/
	iter = 0;
	repeat = TRUE;
 	fail = FALSE;;
	while ( repeat == TRUE ) {
		iter++;
		if (debug_set == TRUE) {
			output_msg(OUTPUT_MESSAGE,"\nBeginning set iteration %d.\n", iter);
		}			
 		if (iter == max_iter + 1) {
			output_msg(OUTPUT_LOG, "Did not converge in set, iteration %d.\n", iterations);
 			fail = TRUE;
 		}
 		if (iter > 2*max_iter) {
			output_msg(OUTPUT_LOG, "Did not converge with relaxed criteria in set.\n");
 			return(OK);
  		}
		molalities(TRUE);
		pitzer();
		/*molalities(TRUE);*/
		mb_sums();
		if (state < REACTION) {
			sum_species();
		} else {
			for (i = 0; i < count_unknowns; i++)  {
				x[i]->sum = x[i]->f;
			}
		}
		/*
		if (debug_set == TRUE) {
			pr.species = TRUE;
			pr.all = TRUE;
			print_species();
		}
		*/
		repeat=FALSE;
		for ( i=0; i < count_unknowns; i++ ) {
			if (x[i] == ph_unknown || x[i] == pe_unknown) continue;
			if (x[i]->type == MB || 
/*			    x[i]->type == ALK || */
			    x[i]->type == CB || 
			    x[i]->type == SOLUTION_PHASE_BOUNDARY || 
			    x[i]->type == EXCH || 
			    x[i]->type == SURFACE ) {
				
				if ( debug_set == TRUE ) {
					output_msg(OUTPUT_MESSAGE,"\n\t%5s  at beginning of set %d: %e\t%e\t%e\n", x[i]->description, iter, (double) x[i]->sum, (double) x[i]->moles, (double) x[i]->master[0]->s->la);
				}
				if (fabs(x[i]->moles) < 1e-30) x[i]->moles = 0;
				f = fabs(x[i]->sum);
				if (f == 0 && x[i]->moles == 0) {
					x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					continue;
				} else if (f == 0) {
					repeat = TRUE;
					x[i]->master[0]->s->la += 5;
/*!!!!*/				if (x[i]->master[0]->s->la < -999.) x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
 				} else if (fail == TRUE && f < 1.5 * fabs(x[i]->moles)) {
 					continue;
				} else if (f > 1.5 * fabs(x[i]->moles) || f < 1e-5 * fabs(x[i]->moles) ) {
					weight = (f < 1e-5 * fabs(x[i]->moles)) ? 0.3 : 1.0;
					if (x[i]->moles <= 0) {
						x[i]->master[0]->s->la = MIN_RELATED_LOG_ACTIVITY;
					} else {
						repeat = TRUE;
						x[i]->master[0]->s->la += weight * log10(fabs(x[i]->moles / x[i]->sum));
					}
					if ( debug_set == TRUE ) {
						output_msg(OUTPUT_MESSAGE,"\t%5s not converged in set %d: %e\t%e\t%e\n", x[i]->description, iter, (double) x[i]->sum, (double) x[i]->moles, (double) x[i]->master[0]->s->la);
					}
				}
			} else if (x[i]->type == ALK) {
				f = total_co2;
 				if (fail == TRUE && f < 1.5 * fabs(x[i]->moles)) {
 					continue;
 				}
				if (f > 1.5 * fabs(x[i]->moles) || f < 1e-5 * fabs(x[i]->moles) ) {
					repeat = TRUE;
					weight = (f < 1e-5 * fabs(x[i]->moles)) ? 0.3 : 1.0;
					x[i]->master[0]->s->la += weight * 
						log10(fabs(x[i]->moles / x[i]->sum));
					if ( debug_set == TRUE ) {
						output_msg(OUTPUT_MESSAGE,"%s not converged in set. %e\t%e\t%e\n", x[i]->description, (double) x[i]->sum, (double) x[i]->moles, (double) x[i]->master[0]->s->la);
					}
				}
			}
		}
	}
	output_msg(OUTPUT_LOG,"Iterations in revise_guesses: %d\n", iter);
	mu_x = mu_unknown->f * 0.5 / mass_water_aq_x;
	if (mu_x <= 1e-8) {
		mu_x = 1e-8;
	}
	/*gammas(mu_x);*/
	return(OK);
}
/* ---------------------------------------------------------------------- */
int build_pitzer_complexes(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Includes calculation of inverse saturation index in sum_mb.
 *   Puts coefficients in iap and mass balance equations for each phase.
 */
	int i;
	int stop, j, k, l;
	char token[MAX_LENGTH];
	char *ptr;
	struct master *master_ptr;
	struct rxn_token *rxn_ptr;
/*
 *   Build into sums the logic to calculate function for pitzer complexes
 */
	if (!pitzer_model) return(OK);
/*
 *   Calculate function for mass action equation
 */
	for (i = 0; i < count_unknowns; i++) {
		if (x[i]->type != COMPLEX || x[i]->s->rxn_x == NULL) continue;
		if (pitzer_complex_unknown == NULL) pitzer_complex_unknown = x[i];
		fprintf(stderr, "Species: %s\n", x[i]->s->name);

		store_mb(&(x[i]->s->lk), &(x[i]->f), 1.0);

		for (rxn_ptr = x[i]->s->rxn_x->token; rxn_ptr->s != NULL; rxn_ptr++) {
			store_mb(&(rxn_ptr->s->la), &(x[i]->f), rxn_ptr->coef);
			fprintf(stderr, "\t%s %e\n", rxn_ptr->s->name, rxn_ptr->coef);
		}
	}
	return(OK);
}
