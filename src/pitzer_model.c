#define EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"

static char const svnid[] = "$Id: model.c 198 2005-03-31 18:11:06Z dlpark $";

static int initial_guesses(void);
static int revise_guesses(void);
static int remove_unstable_phases;

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
		/*pitzer();*/
		/*s_h2o->la = 0.0;*/
		/*molalities(TRUE);*/
		mb_sums();
		if (state < REACTION) {
			sum_species();
		} else {
			for (i = 0; i < count_unknowns; i++)  {
				x[i]->sum = x[i]->f;
			}
		}
		/*n
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
	/*mu_x = mu_unknown->f * 0.5 / mass_water_aq_x;*/
	if (mu_x <= 1e-8) {
		mu_x = 1e-8;
	}
	/*gammas(mu_x);*/
	return(OK);
}
#ifdef SKIP
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
#endif
/* ---------------------------------------------------------------------- */
int jacobian_pz(void)
/* ---------------------------------------------------------------------- */
{
	double *base;
	double d, d1, d2;
	int i, j;

	base = (LDBLE *) PHRQ_malloc((size_t) count_unknowns * sizeof(LDBLE));
	if (base == NULL) malloc_error();
	for (i = 0; i < count_unknowns; i++) {
		base[i] = residual[i];
	}
	d = 0.001;
	d1 = d*log(10.0);
	for (i = 0; i < count_unknowns; i++) {
		switch (x[i]->type) {
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
			x[i]->master[0]->s->la += d;
			d2 = d1;
			break;
		case MH2O:
			mass_water_aq_x *= (1.0 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x/gfw_water;
			d2 = log(1.0 + d);
			break;
		case MU:
		case AH2O:
		case MH:
		case PP:
		case S_S_MOLES:
			continue;
			break;
		}
		molalities();
		mb_sums();
		residuals();
		for (j = 0; j < count_unknowns; j++) {
			array[j*(count_unknowns + 1) + i] = -(residual[j] - base[j])/d2;
		}
		switch (x[i]->type) {
		case MB:
		case ALK:
		case CB:
		case SOLUTION_PHASE_BOUNDARY:
		case EXCH:
		case SURFACE:
		case SURFACE_CB:
			x[i]->master[0]->s->la -= d;
			break;
		case MH2O:
			mass_water_aq_x /= (1 + d);
			x[i]->master[0]->s->moles = mass_water_aq_x/gfw_water;
			break;
		}
	}
	molalities();
	mb_sums();
	residuals();
	free_check_null(base);
	return OK;
}
/* ---------------------------------------------------------------------- */
int model_pz(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   model is called after the equations have been set up by prep
 *   and initial guesses have been made in set.
 * 
 *   Here is the outline of the calculation sequence:
 *      residuals--residuals are calculated, if small we are done
 *      sum_jacobian--jacobian is calculated 
 *      ineq--inequality solver is called
 *      reset--estimates of unknowns revised, if changes are small solution
 *         has been found, usually convergence is found in residuals.
 *      gammas--new activity coefficients
 *      molalities--calculate molalities
 *      mb_sums--calculate mass-balance sums
 *      mb_gases--decide if gas_phase exists
 *      mb_s_s--decide if solid_solutions exists
 *      switch_bases--check to see if new basis species is needed
 *         reprep--rewrite equations with new basis species if needed
 *         revise_guesses--revise unknowns to get initial mole balance
 *      check_residuals--check convergence one last time
 *         sum_species--calculate sums of elements from species concentrations
 *
 *      An additional pass through may be needed if unstable phases still exist
 *         in the phase assemblage. 
 */
	int i;
	int kode, return_kode;
	int r;
	int count_infeasible, count_basis_change;
	int debug_model_save;
	int mass_water_switch_save;
	if (svnid == NULL) fprintf(stderr," ");

/*	debug_model = TRUE; */
/*	debug_prep = TRUE; */
/*	debug_set = TRUE; */
	/* mass_water_switch == TRUE, mass of water is constant */
	mass_water_switch_save = mass_water_switch;
	if (mass_water_switch_save == FALSE && delay_mass_water == TRUE) {
		mass_water_switch = TRUE;
	}
	debug_model_save = debug_model;
	pe_step_size_now = pe_step_size;
	step_size_now = step_size;
	status(0, NULL);
	iterations=0;
	gamma_iterations = 0;
	count_basis_change = count_infeasible = 0;
	stop_program = FALSE;
	remove_unstable_phases = FALSE;
	for (; ; ) {
		mb_gases();
		mb_s_s();
		kode = 1;
		while ( ( r = residuals() ) != CONVERGED || remove_unstable_phases == TRUE) {
#if defined(PHREEQCI_GUI)
			if (WaitForSingleObject(g_hKill /*g_eventKill*/, 0) == WAIT_OBJECT_0)
				{
					error_msg("Execution canceled by user.", CONTINUE);
					RaiseException(USER_CANCELED_RUN, 0, 0, NULL);
				}
#endif
			iterations++;
			if (iterations > itmax - 1 && debug_model == FALSE && pr.logfile == TRUE) {
				set_forward_output_to_log(TRUE);
				debug_model = TRUE;
			}
			if (debug_model == TRUE) {
				output_msg(OUTPUT_MESSAGE,"\nIteration %d\tStep_size = %f\n", 
					   iterations, (double) step_size_now);
				output_msg(OUTPUT_MESSAGE,"\t\tPe_step_size = %f\n\n", (double) pe_step_size_now);
			}
			/*
			 *   Iterations exceeded
			 */
			if (iterations > itmax ) {
				sprintf(error_string,"Maximum iterations exceeded, %d\n",itmax);
				warning_msg(error_string);
				stop_program = TRUE;
				break;
			}
			/*
			 *   Calculate jacobian
			 */
			jacobian_sums();
			jacobian_pz();
			/*
			 *   Full matrix with pure phases
			 */
			if ( r == OK || remove_unstable_phases == TRUE) {
				return_kode = ineq(kode);
				if ( return_kode != OK ) {
					if (debug_model == TRUE) {
						output_msg(OUTPUT_MESSAGE, "Ineq had infeasible solution, "
							   "kode %d, iteration %d\n", 
							   return_kode, iterations);
					}
					output_msg(OUTPUT_LOG, "Ineq had infeasible solution, "
						   "kode %d, iteration %d\n", return_kode, iterations);
					count_infeasible++;
				}
				if ( return_kode == 2 ) { 
					ineq(0);
				}
				reset();
			}
			gammas_pz();
			molalities(TRUE);
			if(use.surface_ptr != NULL && 
			   use.surface_ptr->diffuse_layer == TRUE &&
			   use.surface_ptr->related_phases == TRUE)
				initial_surface_water();
			mb_sums();
			mb_gases();
			mb_s_s();
			/* debug
			   species_list_sort();
			   sum_species();
			   print_species();
			   print_exchange();
			   print_surface();
			*/
			if (stop_program == TRUE) {
				break;
			}
		}
/*
 *   Check for stop_program
 */

		if (stop_program == TRUE) {
			break;
		}
		if (check_residuals() == ERROR) {
			stop_program = TRUE;
			break;
		}
		if (remove_unstable_phases == FALSE && mass_water_switch_save == FALSE &&
		    mass_water_switch == TRUE) {
			output_msg(OUTPUT_LOG,"\nChanging water switch to FALSE. Iteration %d.\n", iterations);
			mass_water_switch = FALSE;
			continue;
		}
		gamma_iterations++;
		if (check_gammas() != TRUE) continue;
		if (remove_unstable_phases == FALSE) break;
		if (debug_model == TRUE) {
			output_msg(OUTPUT_MESSAGE,"\nRemoving unstable phases. Iteration %d.\n", iterations);
		}
		output_msg(OUTPUT_LOG,"\nRemoving unstable phases. Iteration %d.\n", iterations);
        }
	output_msg(OUTPUT_LOG,"\nNumber of infeasible solutions: %d\n",count_infeasible);
	output_msg(OUTPUT_LOG,"Number of basis changes: %d\n\n",count_basis_change);
	output_msg(OUTPUT_LOG,"Number of iterations: %d\n\n", iterations);
	debug_model = debug_model_save;
	set_forward_output_to_log(FALSE);
	if (stop_program == TRUE) {
		return(ERROR);
	}
	return(OK);
}
/* ---------------------------------------------------------------------- */
int check_gammas(void)
/* ---------------------------------------------------------------------- */
{
	double *base, old_aw, old_mu, tol;
	int converge, i;

	base = (LDBLE *) PHRQ_malloc((size_t) count_s_x * sizeof(LDBLE));
	if (base == NULL) malloc_error();

	for (i = 0; i < count_s_x; i++) {
		base[i] = s_x[i]->lg;
	}
	old_mu = mu_x;
	old_aw = s_h2o->la;
	pitzer();
	molalities();
	mb_sums();
	converge = TRUE;
	tol = convergence_tolerance*10.;
	for (i = 0; i < count_s_x; i++) {
		if (fabs(base[i] - s_x[i]->lg) > tol) converge = FALSE;
	}
	if (fabs(old_mu - mu_x) > tol) converge = FALSE;
	if (fabs(old_aw - s_h2o->la) > tol) converge = FALSE;
	
	/* underrelaxation for gammas and la water */
	if (converge == FALSE && mu_x > 1) {
		for (i = 0; i < count_s_x; i++) {
			s_x[i]->lg = base[i] + (s_x[i]->lg - base[i])*.8;
		}
		if (fabs(old_mu - mu_x) > tol) converge = FALSE;
		s_h2o->la = old_aw + (s_h2o->la - old_aw)*.8;
	}

	base = free_check_null(base);
	return converge;
}
/* ---------------------------------------------------------------------- */
int gammas_pz ()
/* ---------------------------------------------------------------------- */
{
/*
 *   Need exchange gammas for pitzer
 */
	int i, j;
	int ifirst, ilast;
	LDBLE f, a_llnl, b_llnl, bdot_llnl, log_g_co2, dln_g_co2, c2_llnl;
	LDBLE s1, s2, s3;
	LDBLE c1, c2, a, b;
	LDBLE muhalf;
	/* Initialize */
	a_llnl = b_llnl = bdot_llnl = log_g_co2 = dln_g_co2 = c2_llnl = 0;
/*
 *   Calculate activity coefficients
 */
	for (i=0; i < count_s_x; i++) {
		switch (s_x[i]->gflag) {
		    case 0:                   /* uncharged */
		    case 1:                   /* Davies */
		    case 2:                   /* Extended D-H, WATEQ D-H */
		    case 3:                   /* Always 1.0 */
			    break;
		    case 4:		   /* Exchange */
/*
 *   Find CEC
 *   z contains valence of cation for exchange species, alk contains cec
 */
/* !!!!! */
			for (j=1; s_x[i]->rxn_x->token[j].s != NULL; j++) {
				if (s_x[i]->rxn_x->token[j].s->type == EX) {
					s_x[i]->alk = s_x[i]->rxn_x->token[j].s->primary->unknown->moles;
					break;
				}
			}
			/*
			 *   Master species is a dummy variable with meaningless activity and mass
			 */
			if (s_x[i]->primary != NULL) {
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
			} else {
				if (s_x[i]->alk <= 0) {
					s_x[i]->lg = 0.0;
				} else {
					s_x[i]->lg = log10(fabs(s_x[i]->equiv)/s_x[i]->alk);
				}
				s_x[i]->dg = 0.0;
			}
			break;
		    case 5:                   /* Always 1.0 */
			    break;
		    case 6:		   /* Surface */
/*
 *   Find moles of sites. 
 *   s_x[i]->equiv is stoichiometric coefficient of sites in species
 */
			for (j=1; s_x[i]->rxn_x->token[j].s != NULL; j++) {
				if (s_x[i]->rxn_x->token[j].s->type == SURF) {
					s_x[i]->alk = s_x[i]->rxn_x->token[j].s->primary->unknown->moles;
					break;
				}
			}
			if (s_x[i]->alk > 0) {
				s_x[i]->lg = log10(s_x[i]->equiv / s_x[i]->alk);
				s_x[i]->dg = 0.0;
			} else {
				s_x[i]->lg = 0.0;
				s_x[i]->dg = 0.0;
			}
			break;
		    case 7:		   /* LLNL */
			    break;
		    case 8:		   /* LLNL CO2*/
			    break;
		    case 9:		   /* activity water */
			s_x[i]->lg = log10(exp( s_h2o->la * LOG_10) * gfw_water);
			s_x[i]->dg = 0.0;
			break;
		}
/*
		if (mu_unknown != NULL) {
			if (fabs(residual[mu_unknown->number]) > 0.1 &&
			    fabs(residual[mu_unknown->number])/mu_x > 0.5) {
				s_x[i]->dg = 0.0;
			}
		}
 */
	}
	return(OK);
}
