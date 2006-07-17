#define  EXTERNAL extern
#include "global.h"
#include "phqalloc.h"
#include "output.h"
#include "phrqproto.h"
static char const svnid[] = "$Id$";

struct spec {
	char *name;		/* name of species */
	LDBLE a;		/* activity */
	LDBLE lm;		/* log(concentration) */
	LDBLE lg;		/* log(gamma) */
	LDBLE c;		/* concentration */
	LDBLE z;		/* charge number */
	LDBLE Dp;		/* pore water diffusion coefficient, m2/s */
};
struct sol_D {
	int count_spec;	/* number of species */
	struct spec *spec;
} *sol_D;

struct J_ij {
	char *name;
	LDBLE tot;
} *J_ij;

static int multi_D(LDBLE kin_time);
static int find_J(int cell_no);
static int fill_spec(int cell_no);
static int sort_species_name(const void *ptr1, const void *ptr2);
static int multi_Dstag(int mobile_cell);
static int find_Jstag(int icell, int jcell, LDBLE mixf);

LDBLE diffc_max, J_ij_sum;
int J_ij_count_spec;

static int init_mix(void);
static int init_heat_mix(int nmix);
static int heat_mix(int heat_nmix);
static int mix_stag(int i, LDBLE stagkin_time, int punch, LDBLE step_fraction_kin);


LDBLE *heat_mix_array;
LDBLE *temp1, *temp2;
int heat_nmix;
LDBLE heat_mix_f_imm, heat_mix_f_m;

/* ---------------------------------------------------------------------- */
int transport(void)
/* ---------------------------------------------------------------------- */
{
	int i, j, k, n;
	int j_imm, n_m, n_imm;
	LDBLE b, f, mix_f_m, mix_f_imm;
	LDBLE water_m, water_imm;
	int first_c, last_c, b_c;
	int nmix;
	int max_iter;
	char token[MAX_LENGTH];
	LDBLE kin_time, stagkin_time, kin_time_save;
	struct mix *mix_ptr;
	int punch_boolean;
	LDBLE step_fraction;
	LDBLE cb_tol;
	if (svnid == NULL) fprintf(stderr," ");

	state = TRANSPORT;
	diffc_max = 0.0;
	cb_tol = 1e-9;
/*	  mass_water_switch = TRUE; */
/*
 *   Check existence of solutions
 */
	j = -1;
	/* check column solutions */
	for (i = 1; i <= count_cells; i++) {
		use.solution_ptr = solution_bsearch(i, &n, TRUE);
		if (use.solution_ptr == NULL) {
			input_error++;
			sprintf(error_string, "Solution %d is needed for transport, but is not defined.", i);
			error_msg(error_string, CONTINUE);
		}
		else
			cell_data[i - 1].temp = use.solution_ptr->tc;
	}

	if (multi_Dflag == TRUE)
		sol_D = (struct sol_D *) PHRQ_malloc((size_t) (count_cells + 2 + stag_data->count_stag * count_cells) * sizeof(struct sol_D));

	/* check solution 0 */
	if (solution_bsearch(0, &n, FALSE) == NULL) {
		if (ishift == 1) {
			input_error++;
			sprintf(error_string, "Solution 0 is needed for transport, but is not defined.");
			error_msg(error_string, CONTINUE);
		}
		else
			solution_duplicate(1, 0);
	}

	/* check solution count_cells */
	if (solution_bsearch(count_cells + 1, &n, FALSE) == NULL) {
		if (ishift == -1) {
			input_error++;
			sprintf(error_string, "Solution %d is needed for transport, but is not defined.", count_cells+1);
			error_msg(error_string, CONTINUE);
		}
		else
			solution_duplicate(count_cells, count_cells + 1);
	}
/*
 *   Initialize temperature in stagnant cells ...
 */
	for (n = 1; n <= stag_data->count_stag; n++) {
		for (i = 1; i <= count_cells; i++) {
			k = i + 1 + n * count_cells;
			use.solution_ptr = solution_bsearch(k, &use.n_solution, FALSE);
			if (use.solution_ptr != NULL)
				cell_data[k - 1].temp = use.solution_ptr->tc;
		}
	}
/*
 * First equilibrate solutions
 */
	dup_print("Equilibrating initial solutions", TRUE);
	transport_step = 0;
	/*for (i = 0; i <= count_cells + 1; i++) {*/
	for (i = 0; i <= count_cells ; i++) {
		set_initial_moles(i);
		cell_no = i;
		set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
		if (multi_Dflag == TRUE) {
			if (fabs(cb_x) > (cb_tol * mu_x)) {
/*				input_error++;
				sprintf(error_string, "Solution %d must be charge-balanced for multicomponent diffusion.", i);
				error_msg(error_string, CONTINUE);
 */
/*				sprintf(token, "Solution %d has %.2e charge imbalance in multicomponent diffusion.", i, cb_x);
				warning_msg(token);
 */			}
			fill_spec(cell_no);
		}
		if (cell_data[i - 1].punch == TRUE)
			punch_all();
		if (cell_data[i - 1].print == TRUE)
			print_all();
		saver();
	}
/*
 * Also stagnant cells
 */
	for (n = 1; n <= stag_data->count_stag; n++) {
		for (i = 1; i <= count_cells; i++) {
			k = i + 1 + n * count_cells;
			cell_no = k;
			if (solution_bsearch(k, &use.n_solution, FALSE) != 0 ) {
				set_initial_moles(k);
				set_and_run_wrapper(k, NOMIX, FALSE, k, 0.0);
				if (multi_Dflag == TRUE) {
					if (fabs(cb_x) > (cb_tol * mu_x)) {
/*						input_error++;
						sprintf(error_string, "Solution %d must be charge-balanced for multicomponent diffusion.", i);
						error_msg(error_string, CONTINUE);
 */
/*				sprintf(token, "Solution %d has %.2e charge imbalance in multicomponent diffusion.", i, cb_x);
				warning_msg(token);
 */			}
					fill_spec(cell_no);
				}
				if ((cell_data[k - 1].punch == TRUE))
					punch_all();
				if ((cell_data[k - 1].print == TRUE) && (transport_step % print_modulus == 0))
					print_all();
				saver();
			}
		}
	}
#ifdef SKIP
	if (multi_Dflag == TRUE) {
		set_and_run_wrapper(0, NOMIX, FALSE, 0, 0.0);
		if (fabs(cb_x) > (cb_tol * mu_x)) {
/*			input_error++;
			sprintf(error_string, "Solution %d must be charge-balanced for multicomponent diffusion.", i);
			error_msg(error_string, CONTINUE);
 */
/*				sprintf(token, "Solution %d has %.2e charge imbalance in multicomponent diffusion.", i, cb_x);
				warning_msg(token);
 */			}
		fill_spec(0);

		saver();
		set_and_run_wrapper(count_cells + 1, NOMIX, FALSE, count_cells + 1, 0.0);
		if (fabs(cb_x) > (cb_tol * mu_x)) {
/*			input_error++;
			sprintf(error_string, "Solution %d must be charge-balanced for multicomponent diffusion.", i);
			error_msg(error_string, CONTINUE);
 */
/*				sprintf(token, "Solution %d has %.2e charge imbalance in multicomponent diffusion.", i, cb_x);
				warning_msg(token);
 */			}
		fill_spec(count_cells + 1);
		saver();
   }
#endif
/*
 *  Initialize mixing factors, define kinetics times
 *  for multicomponent diffusion, limit mixing by diffc_max (usually from H+)
 */
	if (multi_Dflag == TRUE)
		diffc = diffc_max;
	if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1)) {
		if (multi_Dflag == TRUE)
			error_msg("In multicomponent diffusion, stagnant mixing must be defined with keyword MIX.", CONTINUE);
		mix_ptr = &mix[0];
		for (i = 0; i < count_mix; i++)
			mix_free(mix_ptr++);
		count_mix = 2 * count_cells;
/*
 * stagnant mix factors go in mix[0 .. count_cells]
 */
		mix = (struct mix *) PHRQ_realloc(mix, (size_t) count_mix * sizeof(struct mix));
		if (mix == NULL) malloc_error();
		memset(mix, 0, sizeof(struct mix) * count_mix);
	}
/*
 * mix[] is extended in init_mix(), to accommodate column mix factors
 */
	nmix = init_mix();
	heat_nmix = init_heat_mix(nmix);
	if (nmix < 2)
		stagkin_time = timest;
	else
		stagkin_time = timest / nmix;
	if (ishift != 0 )
		kin_time = timest / (1 + nmix);
	else
		kin_time = stagkin_time;
	kin_time_save = kin_time;

/* Reaction defined for a shift... */
	step_fraction = 1.0 / (1.0 + nmix);
/*
 *   Set boundary conditions, transport direction
 */
	last_model.force_prep = TRUE;
	if ((ishift == 0) || (bcon_first == 1) || (bcon_last == 1))
		b_c = 1;
	else
		b_c = 0;
	if (ishift >= 0) {
		last_c = count_cells;
		first_c = 1;
	}
	else {
		last_c = 1;
		first_c = count_cells;
	}
/*
 * Define stagnant/mobile mix structure, if not read explicitly.
 *
 * With count_stag = 1, mix factors are calculated from exchange factor à
 * (= exch_f), mobile é_m (= th_m) and immobile é_im (= th_im) porosity.
 * These variables are read under keyword TRANSPORT, after stagnant, in
 * structure stag_data.
 * MIX 'cell_no' in input file can be an alternative for the calculation here.
 */

	if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1)) {

		b = stag_data->th_m / (stag_data->th_m + stag_data->th_im);
		f = exp( -stag_data->exch_f * stagkin_time / (b * stag_data->th_im) );
		mix_f_imm = b - b * f;
		mix_f_m = mix_f_imm * stag_data->th_im / stag_data->th_m;

		n = 0;
		for (i = 1; i <= count_cells; i++) {
			j = i;
			j_imm = j + (1 + count_cells);
			if (solution_bsearch(j, &n_m, TRUE) == NULL)
				error_msg("Could not find mobile cell solution in TRANSPORT.", STOP);
			if (solution_bsearch(j_imm, &n_imm, TRUE) == NULL)
				error_msg("Could not find immobile cell solution in TRANSPORT.", STOP);
			water_m = solution[n_m]->mass_water;
			water_imm = solution[n_imm]->mass_water;
/*
 * Define C_m = (1 - mix_f_m) * C_m0  +  mix_f_m) * C_im0
 */
			mix[n].comps = (struct mix_comp *) PHRQ_malloc((size_t) 2 * sizeof(struct mix_comp));
			if (mix[n].comps == NULL) malloc_error();
			mix[n].count_comps = 2;
			mix[n].description = string_duplicate(" ");
			mix[n].n_user = j;
			mix[n].n_user_end = j;
			mix[n].comps[0].n_solution = j;
			mix[n].comps[0].fraction = 1 - mix_f_m;
			mix[n].comps[1].n_solution = j_imm;
			mix[n].comps[1].fraction = mix_f_m * water_m / water_imm;
			n++;
/*
 * Define C_im = mix_f_imm * C_m0  +  (1 - mix_f_imm) * C_im0,  or...
 */
			mix[n].comps = (struct mix_comp *) PHRQ_malloc((size_t) 2 * sizeof(struct mix_comp));
			if (mix[n].comps == NULL) malloc_error();
			mix[n].count_comps = 2;
			mix[n].description = string_duplicate(" ");
			mix[n].n_user = j_imm;
			mix[n].n_user_end = j_imm;
			mix[n].comps[0].n_solution = j_imm;
			mix[n].comps[0].fraction = 1 - mix_f_imm;
			mix[n].comps[1].n_solution = j;
			mix[n].comps[1].fraction = mix_f_imm * water_imm / water_m;
			n++;
		}

		if (heat_nmix > 0) {
/*
 * Assumption: D_e used for calculating exch_f in input file equals diffc
 */
			f = stag_data->exch_f * (heat_diffc - diffc) / diffc / tempr;
			f = exp( - f * stagkin_time / (b * stag_data->th_im) );
			heat_mix_f_imm = b - b * f;
			heat_mix_f_m = heat_mix_f_imm * stag_data->th_im / stag_data->th_m;
		}
	}
/*
 *   Stop if error
 */
	if (input_error > 0) {
		error_msg("Program terminating due to input errors.", STOP);
	}
/*
 * Now transport
 */
	max_iter = 0;
	for (transport_step = transport_start; transport_step <= count_shifts; transport_step++) {
		/*
		 *  Set initial moles of phases
		 */
		for (i = 1; i <= count_cells; i++)
			set_initial_moles(i);
/*
 * Start diffusing if boundary cond = 1, (fixed c, or closed)
 */
		if (b_c == 1) {

			/* For half of mixing steps */
			for (j = 1; j <= floor((double) nmix / 2); j++) {
				rate_sim_time_start = (transport_step - 1) * timest + (j - 1) * kin_time;
				rate_sim_time = rate_sim_time_start + kin_time;

				if (multi_Dflag)
					sprintf(token, "Transport step %3d. Multicomponent diffusion run %3d.", transport_step, j);
				else
					sprintf(token, "Transport step %3d. Mixrun %3d.", transport_step, j);
				dup_print(token, FALSE);

				if (heat_nmix > 0) {
					heat_mix(heat_nmix);
					/* equilibrate again ... */
					for (i = 1; i <= count_cells; i++) {
						cell_no = i;
						set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
						if (multi_Dflag)
							fill_spec(i);
						saver();
					}
				}
				/* Go through cells */

				if (multi_Dflag)
					multi_D(stagkin_time);

				for (i = 1; i <= count_cells; i++) {

					if (iterations > max_iter) max_iter = iterations;
					if (multi_Dflag)
						sprintf(token, "Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)", transport_step, j, i, max_iter);
					else
					sprintf(token, "Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)", transport_step, j, i, max_iter);
					status(0, token);

					cell_no = i;
					run_reactions(i, kin_time, DISP, step_fraction);
					if (multi_Dflag)
						fill_spec(i);

					/* punch and output file */
					if ((ishift == 0) && (j == nmix) && ((stag_data->count_stag == 0) ||
						solution_bsearch(i+1+count_cells, &use.n_solution, FALSE) == 0)) {
						if ((cell_data[i-1].punch == TRUE) && (transport_step % punch_modulus == 0))
							punch_all();
						if ((cell_data[i-1].print == TRUE) && (transport_step % print_modulus == 0))
							print_all();
					}
					if (i > 1) solution_duplicate(-2, i - 1);
					saver();
				}
				solution_duplicate(-2, count_cells);

				/* Stagnant zone mixing after completion of each
				   diffusive/dispersive step ...  */

				rate_sim_time_start = (transport_step - 1) * timest + (j - 1) * stagkin_time;
				rate_sim_time = rate_sim_time_start + stagkin_time;

				if (stag_data->count_stag > 0) {
					if ((ishift == 0) && (j == nmix))
						punch_boolean = TRUE;
					else
						punch_boolean = FALSE;
					for (i = 1; i <= count_cells; i++)
						mix_stag(i, stagkin_time, punch_boolean, step_fraction);
				}
			}
		}
/*
 * Advective transport
 */
		if (ishift != 0) {
			sprintf(token, "Transport step %3d.", transport_step);
			dup_print(token, FALSE);
			if (b_c == 1)
				rate_sim_time_start = (transport_step - 1) * timest + (j - 1) * kin_time;
			else
				rate_sim_time_start = (transport_step - 1) * timest;
			rate_sim_time = rate_sim_time_start + kin_time;

/* halftime kinetics for resident water in first cell ... */
			if (kinetics_bsearch(first_c, &i) != NULL && count_cells > 1) {
				cell_no = first_c;
				kin_time = kin_time_save / 2;
				run_reactions(first_c, kin_time, NOMIX, 0.0);
				saver();
				kin_time = kin_time_save;
			}

			/* for each cell in column */
/* Begin revision Dec 7, 1999 */
			for (i = last_c; i != (first_c - ishift); i -= ishift)
				solution_duplicate(i - ishift, i);
/*
 * thermal diffusion when nmix = 0...
 */
			if ((nmix == 0) && (heat_nmix > 0)) {
				heat_mix(heat_nmix);
			/* equilibrate again ... */
				for (i = 1; i <= count_cells; i++) {
					cell_no = i;
					set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
					if (multi_Dflag)
						fill_spec(i);
					saver();
				}
			}

			for (i = 1; i <= count_cells; i++) {
				if (i == first_c && count_cells > 1)
					kin_time /= 2;
				if (multi_Dflag)
					sprintf(token, "Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)", transport_step, 0, i, max_iter);
				else
					sprintf(token, "Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)", transport_step, 0, i, max_iter);
				status(0, token);
				cell_no = i;
				run_reactions(i, kin_time, NOMIX, step_fraction);
				if (multi_Dflag == TRUE)
					fill_spec(i);
				if (iterations > max_iter)
					max_iter = iterations;
/* end revision Dec 7, 1999 */

				if ((nmix == 0) && ((stag_data->count_stag == 0) ||
					(solution_bsearch(i + 1 + count_cells, &use.n_solution, FALSE) == 0))) {
					if ((cell_data[i - 1].punch == TRUE) && (transport_step % punch_modulus == 0))
						punch_all();
					if ((cell_data[i - 1].print == TRUE) && (transport_step % print_modulus == 0))
						print_all();
				}
				if (i == first_c && count_cells > 1)
					kin_time = kin_time_save;
				saver();

				/* If nmix is zero, stagnant zone mixing after
				   advective step ... */

				if ((nmix == 0) && (stag_data->count_stag > 0))
					mix_stag(i, stagkin_time, TRUE, step_fraction);
			}
		}
/*
 * Further dispersive and diffusive transport
 */
		if (b_c != 1)
			j = 1;
		for (j = j; j <= nmix; j++) {
			if (multi_Dflag)
				sprintf(token, "Transport step %3d. Multicomponent diffusion run %3d.", transport_step, j);
			else
				sprintf(token, "Transport step %3d. Mixrun %3d.", transport_step, j);
			dup_print(token, FALSE);
			rate_sim_time_start = (transport_step - 1) * timest + (j - 1) * kin_time;
			if (ishift != 0)
				rate_sim_time_start += kin_time;
			rate_sim_time = rate_sim_time_start + kin_time;

			if (heat_nmix > 0) {
				heat_mix(heat_nmix);
				/* equilibrate again ... */
				for (i = 1; i <= count_cells; i++) {
					cell_no = i;
					set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
					if (multi_Dflag)
						fill_spec(i);
					saver();
				}
			}
			/* for each cell in column */

			if (multi_Dflag == TRUE)
				multi_D(stagkin_time);

			for (i = 1; i <= count_cells; i++) {
				if (iterations > max_iter) max_iter = iterations;
				if (multi_Dflag)
					sprintf(token, "Transport step %3d. MCDrun %3d. Cell %3d. (Max. iter %3d)", transport_step, j, i, max_iter);
				else
					sprintf(token, "Transport step %3d. Mixrun %3d. Cell %3d. (Max. iter %3d)", transport_step, j, i, max_iter);
				status(0, token);
				cell_no = i;
				run_reactions(i, kin_time, DISP, step_fraction);
				if (multi_Dflag == TRUE)
					fill_spec(i);

				if ((j == nmix) && ((stag_data->count_stag == 0) ||
					 (solution_bsearch(i + 1 + count_cells, &use.n_solution, FALSE) == 0))) {
					if ((cell_data[i - 1].punch == TRUE) && (transport_step % punch_modulus == 0))
						punch_all();
					if ((cell_data[i - 1].print == TRUE) && (transport_step % print_modulus == 0))
						print_all();
				}
				if (i > 1)
					solution_duplicate(-2, i - 1);
				saver();
			}
			solution_duplicate(-2, count_cells);
			/* Stagnant zone mixing after completion of each
			   diffusive/dispersive step ... */

			rate_sim_time_start = (transport_step - 1) * timest + (j - 1) * stagkin_time;
			rate_sim_time = rate_sim_time_start + stagkin_time;

			if (stag_data->count_stag > 0) {
				if (j == nmix)
					punch_boolean = TRUE;
				else
					punch_boolean = FALSE;
				for (i = 1; i <= count_cells; i++)
					mix_stag(i, stagkin_time, punch_boolean, step_fraction);
			}
		}
		if (dump_modulus != 0 && (transport_step % dump_modulus) == 0) dump();
	}
#ifdef DOS
	output_msg(OUTPUT_SCREEN, "\n");
#else
	output_msg(OUTPUT_SCREEN, "%s%-80s", "\n", " ");
#endif
	/* free_model_allocs(); */
/*
 * free mix structures
 */
	if ((stag_data->exch_f > 0) && (stag_data->count_stag == 1)) {
		mix_ptr = &mix[0];
		for (i = 0; i < count_mix; i++) mix_free(mix_ptr++);
		count_mix =0;
	}
	else {
		if (nmix > 0) {
			mix_ptr = &mix[count_mix - count_cells];
			for (i = count_mix - count_cells; i < count_mix; i++)
				mix_free(mix_ptr++);

			count_mix -= count_cells;
			mix = (struct mix *) PHRQ_realloc(mix, (size_t) (count_mix + 1)*sizeof(struct mix));
			if (mix == NULL) malloc_error();
		}
	}
	if (heat_nmix > 0) {
		free_check_null(heat_mix_array);
		free_check_null(temp1);
				free_check_null(temp2);
	}
	if (multi_Dflag == TRUE) {
		free_check_null(sol_D);
		for (i = 0; i < count_cells + 2 + stag_data->count_stag * count_cells; i++)
			free_check_null(sol_D[i].spec);
	}
	initial_total_time += rate_sim_time;
	rate_sim_time = 0;
	mass_water_switch = FALSE;
	return(OK);
}
/* ---------------------------------------------------------------------- */
int init_mix(void)
/* ---------------------------------------------------------------------- */
{
	LDBLE dav, lav, mixf, maxmix, corr_disp, diffc_here, mD;
	int i, n, nmix, count_comps, max_mix;
	LDBLE *m;

	m = (LDBLE *) PHRQ_malloc((count_cells+1) * sizeof(LDBLE));
	if (m == NULL) malloc_error();
	if (multi_Dflag == TRUE)
		diffc_here = 0.0;
	else
		diffc_here = diffc;
/*
 * Define mixing factors among inner cells
 */
	corr_disp = 1.;
	if (correct_disp == TRUE && ishift != 0) {
	if (bcon_first == 3)
		corr_disp += 1. / count_cells;
	if (bcon_last == 3)
		corr_disp += 1. / count_cells;
	}
	maxmix = 0.0;
	for (i = 1; i < count_cells; i++) {
		lav = (cell_data[i - 1].length + cell_data[i].length) / 2;
		if (ishift != 0)
			dav = (cell_data[i - 1].disp + cell_data[i].disp) / 2;
		else
			dav = 0;

		mixf = (diffc_here * timest / lav + dav) * corr_disp / cell_data[i].length;
		if (mixf > maxmix)
			maxmix = mixf;
		m[i] = mixf;		   /* m[i] has mixf with lower cell */
		if (multi_Dflag == TRUE) {
			mD = diffc_max * timest / (lav * lav);
			if (mD > maxmix)
				maxmix = mD;
		}
	}
/*
 * Also for boundary cells
 */
	if (bcon_first == 1) {
		lav = cell_data[0].length;
		if (ishift != 0)
			dav = cell_data[0].disp;
		else
			dav = 0;

		mixf = (diffc_here * timest / lav + dav) / lav;
		if (mixf > maxmix)
			maxmix = mixf;
		m[0] = 2 * mixf;
		if (multi_Dflag == TRUE) {
			mD = diffc_max * timest / (lav * lav);
			if (mD > maxmix)
				maxmix = mD;
		}
	}
	else
		m[0] = 0;

	if (bcon_last == 1) {
		lav = cell_data[count_cells - 1].length;
		if (ishift != 0)
			dav = cell_data[count_cells-1].disp;
		else
			dav = 0;

		mixf = (diffc_here * timest / lav + dav) / lav;
		if (mixf > maxmix)
			maxmix = mixf;
		m[count_cells] = 2 * mixf;
		if (multi_Dflag == TRUE) {
			mD = diffc_max * timest / (lav * lav);
			if (mD > maxmix)
				maxmix = mD;
		}
	}
	else
		m[count_cells] = 0;

/*
 * Find number of mixes
 */
	if (maxmix == 0)
		nmix = 0;
	else {
		if ((bcon_first == 1) || (bcon_last == 1))
			nmix = 1 + (int) floor(4.5 * maxmix);
		else
			nmix = 1 + (int) floor(3.0 * maxmix);

		if ((ishift != 0) && ((bcon_first == 1) || (bcon_last == 1))) {
			if (nmix < 2)
				nmix = 2;
		}
		for (i = 0; i <= count_cells; i++)
			m[i] /= nmix;
	}
 /*
  * Fill mix structure
  */
	if (nmix != 0) {
		mix = (struct mix *) PHRQ_realloc(mix, (size_t) (count_mix + count_cells)*sizeof(struct mix));
		if (mix == NULL) malloc_error();
		count_mix += count_cells;
		for (n = count_mix - count_cells; n < count_mix; n++) {
			mix[n].description = NULL;
			mix[n].count_comps = 3;
			mix[n].comps = (struct mix_comp *) PHRQ_malloc((size_t) 3 * sizeof(struct mix_comp));
			if (mix[n].comps == NULL) malloc_error();
		}

		n = count_mix - count_cells;
/*
 * max_mix brings n_user outside range of active cells
 * mix[n].n_user = mix[n].n_user_end = -999 has same effect
 * but max_mix keeps mix in sort order in case mix_bsearch
 * is used
 */
		if (n - 1 <= 0)
			max_mix = 1;
		else
			max_mix = mix[n - 1].n_user + 1;

		if (max_mix < count_cells * (stag_data->count_stag + 1) + 1)
			max_mix = count_cells * (stag_data->count_stag + 1) + 1;

		for (i = 1; i <= count_cells; i++) {
			dav = 0;
			count_comps = 0;
			mix[n].description = (char *) free_check_null(mix[n].description);
			mix[n].description = string_duplicate(" ");
/*
 * again, max_mix brings n_user outside range of active cells, etc...
 */
			mix[n].n_user = max_mix + i ;
			mix[n].n_user_end = max_mix + i ;

			mix[n].comps[count_comps].n_solution = i - 1;
			mix[n].comps[count_comps].fraction = m[i - 1];
			dav += m[i - 1];
			count_comps++;
			mix[n].comps[count_comps].n_solution = i + 1;
			mix[n].comps[count_comps].fraction = m[i];
			dav += m[i];
			count_comps++;
			mix[n].comps[count_comps].n_solution = i;
			mix[n].comps[count_comps].fraction = 1.0 - dav;

			n++;
		}
	}
	m = (LDBLE *) free_check_null(m);
	return(nmix);
}
/* ---------------------------------------------------------------------- */
int mix_stag(int i, LDBLE kin_time, int punch, LDBLE step_fraction)
/* ---------------------------------------------------------------------- */
{
	int n, k, l;
	LDBLE t_imm;
	struct solution *ptr_imm, *ptr_m;
#ifdef SKIP
	char str[MAX_LENGTH];
#endif
/*
 * Kinetics in transport cell is done while transporting
 */
	for (n = 1; n <= stag_data->count_stag; n++) {
		k = i + 1 + n * count_cells;
		if ((ptr_imm = solution_bsearch(k, &use.n_solution, FALSE)) != NULL) {
			if (n == 1) {
				if (heat_nmix > 0) {
					ptr_m = solution_bsearch(i, &use.n_solution, FALSE);
					t_imm = heat_mix_f_imm * ptr_m->tc + (1 - heat_mix_f_imm) * ptr_imm->tc;
					ptr_m->tc = heat_mix_f_m * ptr_imm->tc + (1 - heat_mix_f_m) * ptr_m->tc;
					cell_data[i - 1].temp = ptr_m->tc;
					cell_data[k - 1].temp = ptr_imm->tc = t_imm;
					/* equilibrate again ... */
					cell_no = i;
					set_and_run_wrapper(i, NOMIX, FALSE, i, 0.0);
					if (multi_Dflag == TRUE)
						fill_spec(cell_no);
					saver();
					cell_no = k;
					set_and_run_wrapper(k, NOMIX, FALSE, k, 0.0);
					if (multi_Dflag == TRUE)
						fill_spec(cell_no);
					saver();
				}
/*
 * Mobile cell, kinetics already done ...
 */
				cell_no = i;
				if (multi_Dflag == TRUE)
					multi_Dstag(i);
				set_and_run_wrapper(i, STAG, FALSE, -2, 0.0);
				if (multi_Dflag == TRUE)
					fill_spec(cell_no);
				if ((use.kinetics_ptr = kinetics_bsearch(i, &l)) != NULL) {
					use.n_kinetics_user = i;
					use.kinetics_in = TRUE;
				}
				if ((punch == TRUE) && (cell_data[i - 1].punch == TRUE) &&
					(transport_step % punch_modulus == 0))
					punch_all();
				if ((punch == TRUE) && (cell_data[i - 1].print == TRUE) &&
					(transport_step % print_modulus == 0))
					print_all();
				saver();
			}
			cell_no = k;
			run_reactions(k, kin_time, STAG, step_fraction);
			if (multi_Dflag == TRUE)
				fill_spec(cell_no);
			if ((cell_data[k - 1].punch == TRUE) && (punch == TRUE) &&
				(transport_step % punch_modulus == 0))
				punch_all();
			if ((cell_data[k - 1].print == TRUE) && (punch == TRUE) &&
				(transport_step % print_modulus == 0))
				print_all();
			saver();
		}
	}
	for (n = 1; n <= stag_data->count_stag; n++) {
		k = i + 1 + n * count_cells;
		if (solution_bsearch(k, &use.n_solution, FALSE) != 0) {
			solution_duplicate(-2 - k, k);
			if (n == 1)
				solution_duplicate(-2, i);
		}
	}
	return(OK);
}
/* ---------------------------------------------------------------------- */
int init_heat_mix(int nmix)
/* ---------------------------------------------------------------------- */
{
	LDBLE lav, mixf, maxmix, corr_disp;
	int i, j, k, n;
	int heat_nmix;
	LDBLE t0;
/*
 * Check for need to model thermal diffusion...
 */
	if (heat_diffc <= diffc)
		return(0);
	if (count_cells < 2)
		return(0);

	heat_nmix = 0;
	t0 = solution_bsearch(0, &n, FALSE)->tc;
	for (i = 0; i < count_cells; i++) {
		if (fabs(cell_data[i].temp - t0) > 1.0) {
			heat_nmix = 1;
			break;
		}
	}
	if (heat_nmix == 0) {
		if (fabs(solution_bsearch(count_cells + 1, &n, FALSE)->tc - t0) > 1.0)
			heat_nmix = 1;
		for (n = 1; n <= stag_data->count_stag; n++) {
			for (i = 1; i < count_cells; i++) {
				k = i + 1 + n * count_cells;
				if (solution_bsearch(k, &j, FALSE) != 0) {
					if (fabs(cell_data[k - 1].temp - t0) > 1.0) {
						heat_nmix = 1;
						break;
					}
				}
			}
		}
	}
	if (heat_nmix == 0) return(0);
/*
 * Initialize arrays...
 */
	heat_mix_array = (LDBLE *) PHRQ_malloc((count_cells+2)*sizeof(LDBLE));
	if (heat_mix_array == NULL) malloc_error();

	temp1 = (LDBLE *) PHRQ_malloc((count_cells + 2)*sizeof(LDBLE));
	if (temp1 == NULL) malloc_error();

	temp2 = (LDBLE *) PHRQ_malloc((count_cells + 2)*sizeof(LDBLE));
	if (temp2 == NULL) malloc_error();
/*
 * Define mixing factors among inner cells...
 */
	corr_disp = 1.;
	if (correct_disp == TRUE && ishift != 0) {
		if (bcon_first == 3)
			corr_disp += 1. / count_cells;
		if (bcon_last == 3)
			corr_disp += 1. / count_cells;
	}
	if (nmix > 0) corr_disp /= nmix;
		maxmix = 0.0;
	for (i = 1; i < count_cells; i++) {
		lav = (cell_data[i - 1].length + cell_data[i].length) / 2;
		mixf = (heat_diffc - diffc) * timest * corr_disp / tempr / (lav *lav);
		if (mixf > maxmix)
			maxmix = mixf;
		heat_mix_array[i + 1] = mixf;	/* m[i] has mixf with lower cell */
	}
/*
 * Also for boundary cells
 */
	if (bcon_first == 1) {
		lav = cell_data[0].length;
		mixf = (heat_diffc - diffc) * timest * corr_disp / tempr / (lav * lav);
		if (2 * mixf > maxmix)
			maxmix = 2 * mixf;
		heat_mix_array[1] = 2 * mixf;
	}
	else
		heat_mix_array[1] = 0;

	if (bcon_last == 1) {
		lav = cell_data[count_cells - 1].length;
		mixf = (heat_diffc - diffc) * timest * corr_disp / tempr / (lav * lav);
		if (2 * mixf > maxmix) maxmix = 2 * mixf;
		heat_mix_array[count_cells + 1] = 2 * mixf;
	}
	else
		heat_mix_array[count_cells+1] = 0;
/*
 * Find number of mixes
 */
	if (maxmix == 0)
		heat_nmix = 0;
	else {
		heat_nmix = 1 + (int) floor(3.0 * maxmix);
		for (i = 1; i <= count_cells + 1; i++)
			heat_mix_array[i] /= heat_nmix;
	}

	return(heat_nmix);
}
/* ---------------------------------------------------------------------- */
int heat_mix(int heat_nmix)
/* ---------------------------------------------------------------------- */
{
	int i, j;

	for (i = 1; i <= count_cells; i++) temp1[i] = solution_bsearch(i, &j, FALSE)->tc;
	temp1[0] = solution_bsearch(0, &j, FALSE)->tc;
	temp1[count_cells+1] = solution_bsearch((count_cells+1), &j, FALSE)->tc;

	for (i = 1; i <= heat_nmix; i++) {
		for (j = 1; j <= count_cells; j++)
			temp2[j] = heat_mix_array[j] * temp1[j-1] + heat_mix_array[j+1] * temp1[j+1] +
					(1 - heat_mix_array[j] - heat_mix_array[j+1]) * temp1[j];
		for (j = 1; j <= count_cells; j++)
			temp1[j] = temp2[j];
	}

	for (i = 1; i <= count_cells; i++) {
		cell_data[i-1].temp = temp1[i];
		solution_bsearch(i, &j, FALSE)->tc = temp1[i];
	}

	return(OK);
}
/* ---------------------------------------------------------------------- */
int set_initial_moles(int i)
/* ---------------------------------------------------------------------- */
{
	struct pp_assemblage *pp_assemblage_ptr;
	struct gas_phase *gas_phase_ptr;
	struct kinetics *kinetics_ptr;
	struct s_s_assemblage *s_s_assemblage_ptr;
	int j, k, n;
	/*
	 *   Pure phase assemblage
	 */
	pp_assemblage_ptr = pp_assemblage_bsearch (i, &n);
	if (pp_assemblage_ptr != NULL) {
		for (j = 0; j < pp_assemblage_ptr->count_comps; j++) {
			pp_assemblage_ptr->pure_phases[j].initial_moles = pp_assemblage_ptr->pure_phases[j].moles;
			if (pp_assemblage_ptr->pure_phases[j].initial_moles < 0)
				pp_assemblage_ptr->pure_phases[j].initial_moles = 0;
		}
	}
	/*
	 *   Gas phase
	 */
	gas_phase_ptr = gas_phase_bsearch (i, &n);
	if (gas_phase_ptr != NULL) {
		for (j = 0; j < gas_phase_ptr->count_comps; j++)
			gas_phase_ptr->comps[j].initial_moles = gas_phase_ptr->comps[j].moles;
	}
	/*
	 *   Kinetics
	 */
	kinetics_ptr = kinetics_bsearch (i, &n);
	if (kinetics_ptr != NULL) {
		for (j = 0; j < kinetics_ptr->count_comps; j++)
			kinetics_ptr->comps[j].initial_moles = kinetics_ptr->comps[j].m;
	}
	/*
	 *   Solid solutions
	 */
	s_s_assemblage_ptr = s_s_assemblage_bsearch (i, &n);
	if (s_s_assemblage_ptr != NULL) {
		for (k = 0; k < s_s_assemblage_ptr->count_s_s; k++) {
			for (j = 0; j < s_s_assemblage_ptr->s_s[k].count_comps; j++)
				s_s_assemblage_ptr->s_s[k].comps[j].init_moles = s_s_assemblage_ptr->s_s[k].comps[j].moles;
		}
	}
	return(OK);
}
/* ---------------------------------------------------------------------- */
int multi_D(LDBLE DDt)
/* ---------------------------------------------------------------------- */
{
/* basic scheme:
 * 1. determine mole transfer (mol/s) of all solute species > 1e-20 mol/L
 * 2. sum up as mole transfer of master_species
 * 3. add moles of master_species to solutions for mixing timestep DDt
 */
	int i, j, k, l, length;
	char *ptr;
	char token[MAX_LENGTH];
	struct m_s {
		char *name;
		LDBLE tot;
	} *m_s;
	int count_m_s;
	LDBLE tot_h, tot_o, temp;

	m_s = (struct m_s *) PHRQ_malloc((size_t) count_elements * sizeof(struct m_s));
	if (m_s == NULL) malloc_error();

	for (i = 0; i <= count_cells; i++) {
/*
 * 1. obtain J_ij...
 */
		find_J(i);
/*
 * 2. sum up the primary or secondary master_species from all the solute species
 *	H and O go in total_h and total_o
 */
		tot_h = tot_o = 0.0;
		count_m_s = 0;
		for (j = 0; j < J_ij_count_spec; j++) {
			ptr = J_ij[j].name;
			count_elts = 0;
			get_elts_in_species(&ptr, 1);
			for (k = 0; k < count_elts; k++) {
				if (strcmp(elt_list[k].elt->name, "H") == 0)
					tot_h += elt_list[k].coef * J_ij[j].tot;
				else if (strcmp(elt_list[k].elt->name, "O") == 0)
					tot_o += elt_list[k].coef * J_ij[j].tot;
				else {
					for (l = 0; l < count_m_s; l++) {
						length = strlen(elt_list[k].elt->name);
						if (strncmp(m_s[l].name, elt_list[k].elt->name, length) == 0) {
							m_s[l].tot += elt_list[k].coef * J_ij[j].tot;
							break;
						}
					}
					if (l == count_m_s) {
						m_s[l].name = string_hsave(elt_list[k].elt->name);
						m_s[l].tot = elt_list[k].coef * J_ij[j].tot;
						count_m_s++;
					}
				}
			}
		}
/*
 * multiply with timestep...
*/
		for (l = 0; l < count_m_s; l++)
			m_s[l].tot *= DDt;
		tot_h *= DDt;
		tot_o *= DDt;
		J_ij_sum *= DDt;
/*
 * 3. find the solutions, multiply with A/V, and add or subtract mol/L...
 */
		if (i > 0) {
			if (i == count_cells && bcon_last != 1)
				continue;
			use.solution_ptr = solution_bsearch(i, &use.n_solution, FALSE);
			use.solution_ptr->total_h -= tot_h / cell_data[i - 1].length;
			use.solution_ptr->total_o -= tot_o / cell_data[i - 1].length;
			use.solution_ptr->cb -= J_ij_sum / cell_data[i - 1].length;
			for (l = 0; l < count_m_s; l++) {
				temp = 0.0;
				length = strlen(m_s[l].name);
				for (j = 0; use.solution_ptr->totals[j].description != NULL; j++) {
					if (strncmp(m_s[l].name, use.solution_ptr->totals[j].description, length) == 0) {
						if (use.solution_ptr->totals[j].moles < m_s[l].tot / cell_data[i - 1].length) {
							temp = use.solution_ptr->totals[j].moles;
							use.solution_ptr->totals[j].moles = 0;
							/* see if other redox states have more moles... */
							for (k = 1; use.solution_ptr->totals[j + k].description != NULL; k++) {
								if (strncmp(m_s[l].name, use.solution_ptr->totals[j + k].description, length) == 0) {
									temp += use.solution_ptr->totals[j + k].moles;
									if (temp < m_s[l].tot / cell_data[i - 1].length) {
										use.solution_ptr->totals[j + k].moles = 0;
									}
									else {
										use.solution_ptr->totals[j + k].moles = temp - m_s[l].tot / cell_data[i - 1].length;
										temp = 0.0;
										break;
									}
								}
							}
							if (temp != 0.0) {
								sprintf(token,"Negative concentration in MCD: added %.1e moles %s in cell %d.",
									m_s[l].tot / cell_data[i - 1].length - temp, m_s[l].name, i);
								warning_msg(token);
							}
						}
						else
							use.solution_ptr->totals[j].moles -= m_s[l].tot / cell_data[i - 1].length;
						break;
					}
#ifdef SKIP
					if (strncmp(m_s[l].name, use.solution_ptr->totals[j].description, length) == NULL) {
						use.solution_ptr->totals[j].moles -= m_s[l].tot / cell_data[i - 1].length;
						if (use.solution_ptr->totals[j].moles < 0) {
							sprintf(token,"Negative concentration in MCD: added %.1e moles %s in cell %d.",
								-use.solution_ptr->totals[j].moles, m_s[l].name, i);
							warning_msg(token);
							use.solution_ptr->totals[j].moles = 0;
						}
						break;
					}
#endif
				}
				if (use.solution_ptr->totals[j].description == NULL) {
					use.solution_ptr->totals = (struct conc *) PHRQ_realloc(use.solution_ptr->totals,
						(size_t) (j + 2) * sizeof(struct conc));
					use.solution_ptr->totals[j].description = string_hsave(m_s[l].name);
					use.solution_ptr->totals[j].moles = -m_s[l].tot / cell_data[i - 1].length;
					if (use.solution_ptr->totals[j].moles < 0) {
						sprintf(token,"Negative concentration in MCD: added %.2e moles %s in cell %d.",
							-use.solution_ptr->totals[j].moles, m_s[l].name, i);
						warning_msg(token);
						use.solution_ptr->totals[j].moles = 0;
					}
					use.solution_ptr->totals[j + 1].description = NULL;
				}
			}
		}
		if (i < count_cells) {
			if (i == 0 && bcon_first != 1)
				continue;
			use.solution_ptr = solution_bsearch(i + 1, &use.n_solution, FALSE);
			use.solution_ptr->total_h += tot_h / cell_data[i].length;
			use.solution_ptr->total_o += tot_o / cell_data[i].length;
			use.solution_ptr->cb += J_ij_sum / cell_data[i].length;
			for (l = 0; l < count_m_s; l++) {
				temp = 0.0;
				length = strlen(m_s[l].name);
				for (j = 0; use.solution_ptr->totals[j].description != NULL; j++) {
					if (strncmp(m_s[l].name, use.solution_ptr->totals[j].description, length) == 0) {
						if (use.solution_ptr->totals[j].moles < -m_s[l].tot / cell_data[i].length) {
							temp = use.solution_ptr->totals[j].moles;
							use.solution_ptr->totals[j].moles = 0;
							/* see if other redox states have more moles... */
							for (k = 1; use.solution_ptr->totals[j + k].description != NULL; k++) {
								if (strncmp(m_s[l].name, use.solution_ptr->totals[j + k].description, length) == 0) {
									temp += use.solution_ptr->totals[j + k].moles;
									if (temp < -m_s[l].tot / cell_data[i].length) {
										use.solution_ptr->totals[j + k].moles = 0;
									}
									else {
										use.solution_ptr->totals[j + k].moles = temp + m_s[l].tot / cell_data[i].length;
										temp = 0.0;
										break;
									}
								}
							}
							if (temp != 0.0) {
								sprintf(token,"Negative concentration in MCD: added %.3e moles %s in cell %d",
									-m_s[l].tot / cell_data[i].length - temp, m_s[l].name, i + 1);
								warning_msg(token);
							}
						}
						else
							use.solution_ptr->totals[j].moles += m_s[l].tot / cell_data[i].length;
						break;
					}
#ifdef SKIP
					if (strncmp(m_s[l].name, use.solution_ptr->totals[j].description, length) == NULL) {
						use.solution_ptr->totals[j].moles += m_s[l].tot / cell_data[i].length;
						if (use.solution_ptr->totals[j].moles < 0) {
							sprintf(token,"Negative concentration in MCD: added %.3e moles %s in cell %d.",
								-use.solution_ptr->totals[j].moles, m_s[l].name, i + 1);
							warning_msg(token);
							use.solution_ptr->totals[j].moles = 0;
						}
						break;
					}
#endif
				}
				if (use.solution_ptr->totals[j].description == NULL) {
					use.solution_ptr->totals = (struct conc *) PHRQ_realloc(use.solution_ptr->totals,
						(size_t) (j + 2) * sizeof(struct conc));
					use.solution_ptr->totals[j].description = string_hsave(m_s[l].name);
					use.solution_ptr->totals[j].moles = m_s[l].tot / cell_data[i].length;
					if (use.solution_ptr->totals[j].moles < 0) {
							sprintf(token,"Negative concentration in MCD: added %.4e moles %s in cell %d.",
								-use.solution_ptr->totals[j].moles, m_s[l].name, i + 1);
						warning_msg(token);
						use.solution_ptr->totals[j].moles = 0;
					}
					use.solution_ptr->totals[j + 1].description = NULL;
				}
			}
		}
	}
	m_s = (struct m_s *) free_check_null(m_s);
	J_ij = (struct J_ij *) free_check_null(J_ij);
	return(OK);
}
/* ---------------------------------------------------------------------- */
int find_J(int cell_no)
/* ---------------------------------------------------------------------- */
{
/* mole transfer (mol/s) of the individual master_species:
 * Eqn 1:
 * J_ij = A_ij * (-D_i*grad(a) + D_i*z_i*c_i * SUM(D_i*z_i*grad(a)) / SUM(D_i*(z_i)^2*c_i))
 */
	int i, i_max, j, j_max, k, l;
	LDBLE lav, ddlm;
	LDBLE *grad, *D, *z, *Dz, *Dzc, *Dzc_dl, Dz2c, Dz2c_dl, c, A_ij;
	struct surface *s_ptr1, *s_ptr2;
	LDBLE dl_s, dl_aq1, dl_aq2, visc1, visc2, c_dl;

	Dzc_dl = NULL;
	if (cell_no == 0) {
		lav = cell_data[0].length / 2;
		A_ij = cell_data[0].por;
	}
	else if (cell_no == count_cells) {
		lav = cell_data[count_cells - 1].length / 2;
		A_ij = cell_data[count_cells - 1].por;
	} else {
		lav = (cell_data[cell_no - 1].length + cell_data[cell_no].length) / 2;
		A_ij = (cell_data[cell_no - 1].por < cell_data[cell_no].por ? cell_data[cell_no - 1].por : cell_data[cell_no].por);
	}
/*
 * check if DL calculations must be made...
 */
	dl_s = dl_aq1 = dl_aq2 = 0.0;
	visc1 = visc2 = 1.0;
	s_ptr1 = surface_bsearch(cell_no, &i);
	if (s_ptr1 != NULL ) {
		if (s_ptr1->diffuse_layer == TRUE) {
			dl_aq1 = s_ptr1->charge->mass_water;
			visc1 = s_ptr1->DDL_viscosity;
		}
	}
	s_ptr2 = surface_bsearch(cell_no + 1, &i);
	if (s_ptr2 != NULL ) {
		if (s_ptr2->diffuse_layer == TRUE) {
			dl_aq2 = s_ptr2->charge->mass_water;
			visc2 = s_ptr2->DDL_viscosity;
		}
	}
	if (cell_no == 0)
		visc1 = visc2;
	else if (cell_no == count_cells)
		visc2 = visc1;

	/* in each cell: DL surface = mass_water_DL / (cell_length * tortuosity)
					 free pore surface = mass_water_free / (cell_length * tortuosity)
	   determine DL surface as a fraction of the total pore surface... */
	if (dl_aq1 > 0) {
		dl_aq1 /= (dl_aq1 + solution[cell_no]->mass_water);
		dl_s = dl_aq1;
	}
	if (dl_aq2 > 0)
		dl_aq2 /= (dl_aq2 + solution[cell_no + 1]->mass_water);
	if (dl_aq1 > 0 && dl_aq2 > 0)
		/* average the 2... */
		dl_s = (dl_aq1 + dl_aq2) / 2;
	else if (dl_aq2 > 0)
		/* there is one DL surface... */
		dl_s = dl_aq2;
/*
 * malloc sufficient space...
 */
	k = sol_D[cell_no].count_spec + sol_D[cell_no + 1].count_spec;

	J_ij = (struct J_ij *) free_check_null(J_ij);
	J_ij = (struct J_ij *) PHRQ_malloc((size_t) k * sizeof(struct J_ij));
	if (J_ij == NULL) malloc_error();

	grad = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (grad == NULL) malloc_error();

	D = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (D == NULL) malloc_error();

	z = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (z == NULL) malloc_error();

	Dz = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (Dz == NULL) malloc_error();

	Dzc = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (Dzc == NULL) malloc_error();

	if (dl_s > 0) {
		Dzc_dl = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
		if (Dzc_dl == NULL) malloc_error();
	}

	Dz2c = Dz2c_dl = 0.0;

/*
 * coefficients in Eqn (1)...
 */
	i = j = k = 0;
	i_max = sol_D[cell_no].count_spec;
	j_max = sol_D[cell_no + 1].count_spec;

	while (i < i_max || j < j_max) {
		if (j == j_max || (strcmp(sol_D[cell_no].spec[i].name, sol_D[cell_no + 1].spec[j].name) < 0 &&
				i < i_max)) {
			/* species 'name' is only in cell_no */
			J_ij[k].name = string_hsave(sol_D[cell_no].spec[i].name);
			D[k] = sol_D[cell_no].spec[i].Dp;
			z[k] = sol_D[cell_no].spec[i].z;
			Dz[k] = D[k] * z[k];
			Dzc[k] = Dz[k] * sol_D[cell_no].spec[i].c / 2;
			if (dl_s > 0 && dl_aq1 > 0) {
				for (l = 0; l < s_ptr1->charge->count_g; l++) {
					if (equal(s_ptr1->charge->g[l].charge, z[k], G_TOL) == TRUE) {
						Dzc_dl[k] = Dz[k] * sol_D[cell_no].spec[i].c / 2 * s_ptr1->charge->g[l].g;
						break;
					}
				}
				Dz2c_dl += Dzc_dl[k] * z[k];
			}
			Dz2c += Dzc[k] * z[k];
			grad[k] = -sol_D[cell_no].spec[i].a / lav;
			if (i < i_max) i++;
		}
		else if (i == i_max || (strcmp(sol_D[cell_no].spec[i].name, sol_D[cell_no + 1].spec[j].name) > 0 &&
				j < j_max)) {
			/* species 'name' is only in (cell_no + 1) */
			J_ij[k].name = string_hsave(sol_D[cell_no + 1].spec[j].name);
			D[k] = sol_D[cell_no + 1].spec[j].Dp;
			z[k] = sol_D[cell_no + 1].spec[j].z;
			Dz[k] = D[k] * z[k];
			Dzc[k] = Dz[k] * sol_D[cell_no + 1].spec[j].c / 2;
			if (dl_s > 0 && dl_aq2 > 0) {
				for (l = 0; l < s_ptr2->charge->count_g; l++) {
					if (equal(s_ptr2->charge->g[l].charge, z[k], G_TOL) == TRUE) {
						Dzc_dl[k] = Dz[k] * sol_D[cell_no + 1].spec[j].c / 2 * s_ptr2->charge->g[l].g;
						break;
					}
				}
				Dz2c_dl += Dzc_dl[k] * z[k];
			}
			Dz2c += Dzc[k] * z[k];
			grad[k] = sol_D[cell_no + 1].spec[j].a / lav;
			if (j < j_max) j++;
		}
		else if (strcmp(sol_D[cell_no].spec[i].name, sol_D[cell_no + 1].spec[j].name) == 0) {
			/* species 'name' is in both cells */
			J_ij[k].name = string_hsave(sol_D[cell_no].spec[i].name);
			if (sol_D[cell_no].spec[i].Dp == 0 || sol_D[cell_no + 1].spec[j].Dp == 0)
				D[k] = 0.0;
			else
				D[k] = (sol_D[cell_no].spec[i].Dp + sol_D[cell_no + 1].spec[j].Dp) / 2;
			z[k] = sol_D[cell_no].spec[i].z;
			Dz[k] = D[k] * z[k];
			Dzc[k] = Dz[k] * (sol_D[cell_no].spec[i].c + sol_D[cell_no + 1].spec[j].c) / 2;
/*			Dzc[k] = Dz[k] * (sol_D[cell_no].spec[i].c > sol_D[cell_no + 1].spec[j].c ? sol_D[cell_no].spec[i].c : sol_D[cell_no + 1].spec[j].c);
 */
			if (dl_s > 0) {
				c_dl = 0.0;
				if (dl_aq1 > 0) {
					for (l = 0; l < s_ptr1->charge->count_g; l++) {
						if (equal(s_ptr1->charge->g[l].charge, z[k], G_TOL) == TRUE) {
							c_dl = sol_D[cell_no].spec[i].c / 2 * s_ptr1->charge->g[l].g;
							break;
						}
					}
				}
				else c_dl = sol_D[cell_no].spec[i].c / 2;

				if (dl_aq2 > 0) {
					for (l = 0; l < s_ptr2->charge->count_g; l++) {
						if (equal(s_ptr2->charge->g[l].charge, z[k], G_TOL) == TRUE) {
							c_dl += sol_D[cell_no + 1].spec[j].c / 2 * s_ptr2->charge->g[l].g;
							break;
						}
					}
				}
				else c_dl += sol_D[cell_no + 1].spec[j].c / 2;

				Dzc_dl[k] = Dz[k] * c_dl;
				Dz2c_dl += Dzc_dl[k] * z[k];
			}
			Dz2c += Dzc[k] * z[k];
			grad[k] = (sol_D[cell_no + 1].spec[j].c - sol_D[cell_no].spec[i].c) / lav;
			ddlm = sol_D[cell_no + 1].spec[j].lm - sol_D[cell_no].spec[i].lm;
			if (fabs(ddlm) > 1e-10)
				grad[k] *= (1 + (sol_D[cell_no + 1].spec[j].lg - sol_D[cell_no].spec[i].lg) / ddlm);
			if (i < i_max) i++;
			if (j < j_max) j++;
		}
		k++;
	}
/*
 * fill in J_ij...
 */
	if (Dz2c == 0)
		k = 0;
	J_ij_count_spec = k;
	J_ij_sum = 0;
	c = 0.0;
	for (j = 0; j < k; j++)
		c += Dz[j] * grad[j];
	for (i = 0; i < k; i++) {
		J_ij[i].tot = -D[i] * grad[i] + c * Dzc[i] / Dz2c;
		if (Dz2c_dl > 0)
			J_ij[i].tot = J_ij[i].tot * (1 - dl_s) +
				(-D[i] * grad[i] + c * Dzc_dl[i] / Dz2c_dl) * dl_s * 2 / (visc1 + visc2);
		J_ij[i].tot *= A_ij;
		J_ij_sum += z[i] * J_ij[i].tot;
	}

	D = (LDBLE *) free_check_null(D);
	z = (LDBLE *) free_check_null(z);
	Dz = (LDBLE *) free_check_null(Dz);
	Dzc = (LDBLE *) free_check_null(Dzc);
	if (dl_s > 0)
		Dzc_dl = (LDBLE *) free_check_null(Dzc_dl);
	grad = (LDBLE *) free_check_null(grad);

	return(OK);
}
/* ---------------------------------------------------------------------- */
int fill_spec(int cell_no)
/* ---------------------------------------------------------------------- */
{
/* copy species activities into sol_D.spec... */

	int i, count_spec;
	LDBLE lm;
	LDBLE por, por_factor, viscos;

	sol_D[cell_no].spec = (struct spec *) free_check_null(sol_D[cell_no].spec);
	sol_D[cell_no].spec = (struct spec *) PHRQ_malloc((size_t) count_species_list * sizeof(struct spec));
	if (sol_D[cell_no].spec == NULL) malloc_error();
	if (cell_no == 0) {
		por = cell_data[0].por;
		if (por < multi_Dpor_lim)
			por = por_factor = 0.0;
		else
			por_factor = pow(por, multi_Dn);
	}
	else if (cell_no == count_cells + 1) {
		por = cell_data[count_cells - 1].por;
		if (por < multi_Dpor_lim)
			por = por_factor = 0.0;
		else
			por_factor = pow(por, multi_Dn);
	}
	else {
		por = cell_data[cell_no - 1].por;
		if (por < multi_Dpor_lim)
			por = por_factor = 0.0;
		else
			por_factor = pow(por, multi_Dn);
	}
/*
 * correct diffusion coefficient for temperature and viscosity, D_T = D_298 * Tk * viscos_298 / (298 * viscos)
 */
	viscos = pow(10, -(1.37023 * (tc_x - 20) + 0.000836 * (tc_x - 20) * (tc_x - 20)) / (109 + tc_x));
/*
 * put temperature factor in por_factor which corrects for porous medium...
 */
	por_factor *= tk_x * 0.88862 / (298.15 * viscos);

	count_spec = 0;
/*
 * sort species by name...
 */
	if (count_species_list > 0)
		qsort (&species_list[0], (size_t) count_species_list,
				(size_t) sizeof(struct species_list), sort_species_name);

	for (i = 0; i < count_species_list; i++) {
		if (species_list[i].s->type == EX) continue;
		if (species_list[i].s->type == SURF) continue;
/*
 *   copy species data
 */
		if (i > 0 && strcmp(species_list[i].s->name, species_list[i-1].s->name) == 0) continue;
		if (species_list[i].s == s_h2o) continue;
/* appt			lm = s_h2o->la - species_list[i].s->lg;
		else
 */			lm = species_list[i].s->lm;

		if (lm > -20) {
			sol_D[cell_no].spec[count_spec].name = string_hsave(species_list[i].s->name);
			sol_D[cell_no].spec[count_spec].c = species_list[i].s->moles / mass_water_aq_x;
			sol_D[cell_no].spec[count_spec].a = under(lm + species_list[i].s->lg);
			sol_D[cell_no].spec[count_spec].lm = lm;
			sol_D[cell_no].spec[count_spec].lg = species_list[i].s->lg;
			sol_D[cell_no].spec[count_spec].z = species_list[i].s->z;
			if (species_list[i].s->dw == 0)
				sol_D[cell_no].spec[count_spec].Dp = default_Dw * por_factor;
			else
				sol_D[cell_no].spec[count_spec].Dp = species_list[i].s->dw * por_factor;
			if (/*strcmp(species_list[i].s->name, "H+") != NULL &&
				strcmp(species_list[i].s->name, "OH-") != NULL && */
				por * sol_D[cell_no].spec[count_spec].Dp > diffc_max)
				diffc_max = por * sol_D[cell_no].spec[count_spec].Dp;
			count_spec++;
		}
	}
	sol_D[cell_no].spec = (struct spec *) PHRQ_realloc(sol_D[cell_no].spec, (size_t) count_spec * sizeof(struct spec));
	if (sol_D[cell_no].spec == NULL) malloc_error();

	sol_D[cell_no].count_spec = count_spec;

	return(OK);
}
/* ---------------------------------------------------------------------- */
int sort_species_name (const void *ptr1, const void *ptr2)
/* ---------------------------------------------------------------------- */
{
	const struct species_list *nptr1, *nptr2;

	nptr1 = (const struct species_list *) ptr1;
	nptr2 = (const struct species_list *) ptr2;

	return(strcmp(nptr1->s->name, nptr2->s->name));
}

/* ---------------------------------------------------------------------- */
int multi_Dstag(int mobile_cell)
/* ---------------------------------------------------------------------- */
{
/* basic scheme follows multi_D, but uses the mix factors predefined with MIX
 * 1. determine mole transfer (mol/s) of all solute species > 1e-20 mol/L
		for the interface between 2 cells.
 * 2. sum up as mole transfer of master_species
 * 3. add moles of master_species to the 2 cells
 *	NOTE. Define the water content of stagnant cells relative to the
 *	mobile cell (with, for example, 1 kg water)
 *	Define properties of each interface only 1 time with MIX.
 */
	int icell, jcell, i, j, k, l, n, length;
	char *ptr;
	char token[MAX_LENGTH];
	struct m_s {
		char *name;
		LDBLE tot;
	} *m_s;
	int count_m_s;
	LDBLE tot_h, tot_o, mixf, temp;

	m_s = (struct m_s *) PHRQ_malloc((size_t) count_elements * sizeof(struct m_s));
	if (m_s == NULL) malloc_error();

	for (n = 0; n <= stag_data->count_stag; n++) {
		icell = mobile_cell + 1 + n * count_cells;
		if (n == 0) icell -= 1;

/*
 *	find the mix ptr for icell and go along the cells that mix with it
 */
		use.mix_ptr = mix_search (icell, &use.n_mix, FALSE);
		if (use.mix_ptr == NULL) continue;
		for (i = 0; i < use.mix_ptr->count_comps; i++) {
			if (use.mix_ptr->comps[i].n_solution == icell) continue;
			jcell = use.mix_ptr->comps[i].n_solution;
			mixf = use.mix_ptr->comps[i].fraction;
/*
 * 1. obtain J_ij...
 */
			find_Jstag(icell, jcell, mixf);
/*
 * 2. sum up the primary or secondary master_species from all the solute species
 *	H and O go in total_h and total_o
 */
			tot_h = tot_o = 0.0;
			count_m_s = 0;
			for (j = 0; j < J_ij_count_spec; j++) {
				ptr = J_ij[j].name;
				count_elts = 0;
				get_elts_in_species(&ptr, 1);
				for (k = 0; k < count_elts; k++) {
					if (strcmp(elt_list[k].elt->name, "H") == 0)
						tot_h += elt_list[k].coef * J_ij[j].tot;
					else if (strcmp(elt_list[k].elt->name, "O") == 0)
						tot_o += elt_list[k].coef * J_ij[j].tot;
					else {
						for (l = 0; l < count_m_s; l++) {
							length = strlen(elt_list[k].elt->name);
							if (strncmp(m_s[l].name, elt_list[k].elt->name, length) == 0) {
								m_s[l].tot += elt_list[k].coef * J_ij[j].tot;
								break;
							}
						}
						if (l == count_m_s) {
							m_s[l].name = string_hsave(elt_list[k].elt->name);
							m_s[l].tot = elt_list[k].coef * J_ij[j].tot;
							count_m_s++;
						}
					}
				}
			}
/*
 * timestep is in mixf.
 * NOTE. The timestep calculated in init_mix for MCD (by PHREEQC) must be equal or smaller than
 *	  the timestep taken (by the user) for calculating mixf in MIX.
 *	  Make this timestep small enough, consider the largest Dw in phreeqd.dat (usually H+).
 *	  Dw used for calculating mixf must be given as default_Dw in the input file.
 */
/*
 * 3. find the solutions, add or subtract the moles...
 */
			use.solution_ptr = solution_bsearch(icell, &use.n_solution, FALSE);
			use.solution_ptr->total_h -= tot_h;
			use.solution_ptr->total_o -= tot_o;
			use.solution_ptr->cb -= J_ij_sum;
			for (l = 0; l < count_m_s; l++) {
				temp = 0.0;
				length = strlen(m_s[l].name);
				for (j = 0; use.solution_ptr->totals[j].description != NULL; j++) {
					if (strncmp(m_s[l].name, use.solution_ptr->totals[j].description, length) == 0) {
						if (use.solution_ptr->totals[j].moles < m_s[l].tot) {
							temp = use.solution_ptr->totals[j].moles;
							use.solution_ptr->totals[j].moles = 0;
							/* see if other redox states have more moles... */
							for (k = 1; use.solution_ptr->totals[j + k].description != NULL; k++) {
								if (strncmp(m_s[l].name, use.solution_ptr->totals[j + k].description, length) == 0) {
									temp += use.solution_ptr->totals[j + k].moles;
									if (temp < m_s[l].tot) {
										use.solution_ptr->totals[j + k].moles = 0;
									}
									else {
										use.solution_ptr->totals[j + k].moles = temp - m_s[l].tot;
										temp = 0.0;
										break;
									}
								}
							}
							if (temp != 0.0) {
								sprintf(token,"Negative concentration in MCD: added %.1e moles %s in cell %d.",
									m_s[l].tot - temp, m_s[l].name, icell);
								warning_msg(token);
							}
						}
						else
							use.solution_ptr->totals[j].moles -= m_s[l].tot;
						break;
					}
#ifdef SKIP
						use.solution_ptr->totals[j].moles -= m_s[l].tot;
						if (use.solution_ptr->totals[j].moles < 0) {
							sprintf(token,"Negative concentration in MCD: added %e moles %s",
								-use.solution_ptr->totals[j].moles, m_s[l].name);
							warning_msg(token);
							use.solution_ptr->totals[j].moles = 0;
						}
						break;
					}
#endif
				}
				if (use.solution_ptr->totals[j].description == NULL) {
					use.solution_ptr->totals = (struct conc *) PHRQ_realloc(use.solution_ptr->totals,
						(size_t) (j + 2) * sizeof(struct conc));
					use.solution_ptr->totals[j].description = string_hsave(m_s[l].name);
					use.solution_ptr->totals[j].moles = -m_s[l].tot;
					if (use.solution_ptr->totals[j].moles < 0) {
						sprintf(token,"Negative concentration in MCD: added %.2e moles %s in cell %d",
							-use.solution_ptr->totals[j].moles, m_s[l].name, icell);
						warning_msg(token);
						use.solution_ptr->totals[j].moles = 0;
					}
					use.solution_ptr->totals[j + 1].description = NULL;
				}
			}
			use.solution_ptr = solution_bsearch(jcell, &use.n_solution, FALSE);
			use.solution_ptr->total_h += tot_h;
			use.solution_ptr->total_o += tot_o;
			use.solution_ptr->cb += J_ij_sum;
			for (l = 0; l < count_m_s; l++) {
				temp = 0.0;
				length = strlen(m_s[l].name);
				for (j = 0; use.solution_ptr->totals[j].description != NULL; j++) {
					if (strncmp(m_s[l].name, use.solution_ptr->totals[j].description, length) == 0) {
						if (use.solution_ptr->totals[j].moles < -m_s[l].tot) {
							temp = use.solution_ptr->totals[j].moles;
							use.solution_ptr->totals[j].moles = 0;
							/* see if other redox states have more moles... */
							for (k = 1; use.solution_ptr->totals[j + k].description != NULL; k++) {
								if (strncmp(m_s[l].name, use.solution_ptr->totals[j + k].description, length) == 0) {
									temp += use.solution_ptr->totals[j + k].moles;
									if (temp < -m_s[l].tot) {
										use.solution_ptr->totals[j + k].moles = 0;
									}
									else {
										use.solution_ptr->totals[j + k].moles = temp + m_s[l].tot;
										temp = 0.0;
										break;
									}
								}
							}
							if (temp != 0.0) {
								sprintf(token,"Negative concentration in MCD: added %.3e moles %s in cell %d",
									-m_s[l].tot - temp, m_s[l].name, jcell);
								warning_msg(token);
							}
						}
						else
							use.solution_ptr->totals[j].moles += m_s[l].tot / cell_data[i].length;
						break;
					}
#ifdef SKIP
						use.solution_ptr->totals[j].moles += m_s[l].tot;
						if (use.solution_ptr->totals[j].moles < 0) {
							sprintf(token,"Negative concentration in MCD: added %e moles %s",
								-use.solution_ptr->totals[j].moles, m_s[l].name);
							warning_msg(token);
							use.solution_ptr->totals[j].moles = 0;
						}
						break;
					}
#endif
				}
				if (use.solution_ptr->totals[j].description == NULL) {
					use.solution_ptr->totals = (struct conc *) PHRQ_realloc(use.solution_ptr->totals,
						(size_t) (j + 2) * sizeof(struct conc));
					use.solution_ptr->totals[j].description = string_hsave(m_s[l].name);
					use.solution_ptr->totals[j].moles = m_s[l].tot;
					if (use.solution_ptr->totals[j].moles < 0) {
						sprintf(token,"Negative concentration in MCD: added %.4e moles %s in cell %d",
							-use.solution_ptr->totals[j].moles, m_s[l].name, jcell);
						warning_msg(token);
						use.solution_ptr->totals[j].moles = 0;
					}
					use.solution_ptr->totals[j + 1].description = NULL;
				}
			}
		}
	}
	m_s = (struct m_s *) free_check_null(m_s);
	J_ij = (struct J_ij *) free_check_null(J_ij);
	return(OK);
}

/* ---------------------------------------------------------------------- */
int find_Jstag(int icell, int jcell, LDBLE mixf)
/* ---------------------------------------------------------------------- */
{
/* mole transfer (mol/s) of the individual master_species:
 * Eqn 1:
 * J_ij = mixf_ij * (-D_i*grad(a) + D_i*z_i*c_i * SUM(D_i*z_i*grad(a)) / SUM(D_i*(z_i)^2*c_i))
 *		mixf_ij = mixf / (Dw * init_pf) * new_por / init_por
 *		mixf is defined in MIX; Dw is default multicomponent diffusion coefficient;
 *		init_pf equals multi_Dpor^multi_Dn;
 */
	int i, i_max, j, j_max, k, l;
	LDBLE ddlm;
	LDBLE *grad, *D, *z, *Dz, *Dzc, *Dzc_dl, Dz2c, Dz2c_dl, c;
	struct surface *s_ptr1, *s_ptr2;
	LDBLE dl_s, dl_aq1, dl_aq2, c_dl, visc1, visc2;

        Dzc_dl = NULL;
	if (cell_data[icell - 1].por < multi_Dpor_lim || cell_data[jcell - 1].por < multi_Dpor_lim)
		mixf = 0.0;
	else
		mixf *= cell_data[icell - 1].por / (default_Dw * pow(multi_Dpor, multi_Dn) * multi_Dpor);
/*
 * check if DL calculations must be made...
 */
	dl_s = dl_aq1 = dl_aq2 = 0.0;
	visc1 = visc2 = 1.0;
	s_ptr1 = surface_bsearch(icell, &i);
	if (s_ptr1 != NULL ) {
		if (s_ptr1->diffuse_layer == TRUE) {
			dl_aq1 = s_ptr1->charge->mass_water;
			visc1 = s_ptr1->DDL_viscosity;
		}
	}
	s_ptr2 = surface_bsearch(jcell, &i);
	if (s_ptr2 != NULL ) {
		if (s_ptr2->diffuse_layer == TRUE) {
			dl_aq2 = s_ptr2->charge->mass_water;
			visc2 = s_ptr2->DDL_viscosity;
		}
	}
	/* in each cell: DL surface = mass_water_DL / (cell_length * tortuosity)
					 free pore surface = mass_water_free / (cell_length * tortuosity)
	   determine DL surface as a fraction of the total pore surface... */
	if (dl_aq1 > 0) {
		dl_aq1 /= (dl_aq1 + solution[icell]->mass_water);
		dl_s = dl_aq1;
	}
	if (dl_aq2 > 0)
		dl_aq2 /= (dl_aq2 + solution[jcell]->mass_water);
	if (dl_aq1 > 0 && dl_aq2 > 0)
		/* average the 2... */
		dl_s = (dl_aq1 + dl_aq2) / 2;
	else if (dl_aq2 > 0)
		/* there is one DL surface... */
		dl_s = dl_aq2;
/*
 * malloc sufficient space...
 */
	k = sol_D[icell].count_spec + sol_D[jcell].count_spec;

	J_ij = (struct J_ij *) free_check_null(J_ij);
	J_ij = (struct J_ij *) PHRQ_malloc((size_t) k * sizeof(struct J_ij));
	if (J_ij == NULL) malloc_error();

	grad = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (grad == NULL) malloc_error();

	D = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (D == NULL) malloc_error();

	z = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (z == NULL) malloc_error();

	Dz = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (Dz == NULL) malloc_error();

	Dzc = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
	if (Dzc == NULL) malloc_error();

	if (dl_s > 0) {
		Dzc_dl = (LDBLE *) PHRQ_malloc((size_t) k * sizeof(LDBLE));
		if (Dzc_dl == NULL) malloc_error();
	}

	Dz2c = Dz2c_dl = 0.0;

/*
 * coefficients in Eqn (1)...
 */
	i = j = k = 0;
	i_max = sol_D[icell].count_spec;
	j_max = sol_D[jcell].count_spec;

	while (i < i_max || j < j_max) {
		if (j == j_max || (strcmp(sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name) < 0 &&
				i < i_max)) {
			J_ij[k].name = string_hsave(sol_D[icell].spec[i].name);
			D[k] = sol_D[icell].spec[i].Dp;
			z[k] = sol_D[icell].spec[i].z;
			Dz[k] = D[k] * z[k];
			Dzc[k] = Dz[k] * sol_D[icell].spec[i].c / 2;
			if (dl_s > 0 && dl_aq1 > 0) {
				for (l = 0; l < s_ptr1->charge->count_g; l++) {
					if (equal(s_ptr1->charge->g[l].charge, z[k], G_TOL) == TRUE) {
						Dzc_dl[k] = Dz[k] * sol_D[icell].spec[i].c / 2 * s_ptr1->charge->g[l].g;
						break;
					}
				}
				Dz2c_dl += Dzc_dl[k] * z[k];
			}
			Dz2c += Dzc[k] * z[k];
			grad[k] = -sol_D[icell].spec[i].a;
			if (i < i_max) i++;
		}
		else if (i == i_max || (strcmp(sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name) > 0 &&
				j < j_max)) {
			J_ij[k].name = string_hsave(sol_D[jcell].spec[j].name);
			D[k] = sol_D[jcell].spec[j].Dp;
			z[k] = sol_D[jcell].spec[j].z;
			Dz[k] = D[k] * z[k];
			Dzc[k] = Dz[k] * sol_D[jcell].spec[j].c / 2;
			if (dl_s > 0 && dl_aq2 > 0) {
				for (l = 0; l < s_ptr2->charge->count_g; l++) {
					if (equal(s_ptr2->charge->g[l].charge, z[k], G_TOL) == TRUE) {
						Dzc_dl[k] = Dz[k] * sol_D[jcell].spec[j].c / 2 * s_ptr2->charge->g[l].g;
						break;
					}
				}
				Dz2c_dl += Dzc_dl[k] * z[k];
			}
			Dz2c += Dzc[k] * z[k];
			grad[k] = sol_D[jcell].spec[j].a;
			if (j < j_max) j++;
		}
		else if (strcmp(sol_D[icell].spec[i].name, sol_D[jcell].spec[j].name) == 0) {
			J_ij[k].name = string_hsave(sol_D[icell].spec[i].name);
			if (sol_D[icell].spec[i].Dp == 0 || sol_D[jcell].spec[j].Dp == 0)
				D[k] = 0.0;
			else
				D[k] = (sol_D[icell].spec[i].Dp + sol_D[jcell].spec[j].Dp) / 2;
			z[k] = sol_D[icell].spec[i].z;
			Dz[k] = D[k] * z[k];
			Dzc[k] = Dz[k] * (sol_D[icell].spec[i].c + sol_D[jcell].spec[j].c) / 2;
/*			Dzc[k] = Dz[k] * (sol_D[icell].spec[i].c > sol_D[jcell].spec[j].c ? sol_D[icell].spec[i].c : sol_D[jcell].spec[j].c);
 */
			if (dl_s > 0) {
				c_dl = 0.0;
				if (dl_aq1 > 0) {
					for (l = 0; l < s_ptr1->charge->count_g; l++) {
						if (equal(s_ptr1->charge->g[l].charge, z[k], G_TOL) == TRUE) {
							c_dl = sol_D[icell].spec[i].c / 2 * s_ptr1->charge->g[l].g;
							break;
						}
					}
				}
				else c_dl = sol_D[icell].spec[i].c / 2;

				if (dl_aq2 > 0) {
					for (l = 0; l < s_ptr2->charge->count_g; l++) {
						if (equal(s_ptr2->charge->g[l].charge, z[k], G_TOL) == TRUE) {
							c_dl += sol_D[jcell].spec[j].c / 2 * s_ptr2->charge->g[l].g;
							break;
						}
					}
				}
				else c_dl += sol_D[jcell].spec[j].c / 2;

				Dzc_dl[k] = Dz[k] * c_dl;
				Dz2c_dl += Dzc_dl[k] * z[k];
			}
			Dz2c += Dzc[k] * z[k];
			grad[k] = (sol_D[jcell].spec[j].c - sol_D[icell].spec[i].c);
			ddlm = sol_D[jcell].spec[j].lm - sol_D[icell].spec[i].lm;
			if (fabs(ddlm) > 1e-10)
				grad[k] *= (1 + (sol_D[jcell].spec[j].lg - sol_D[icell].spec[i].lg) / ddlm);
			if (i < i_max) i++;
			if (j < j_max) j++;
		}
		k++;
	}
/*
 * fill in J_ij...
 */
	if (Dz2c == 0)
		k = 0;
	J_ij_count_spec = k;
	J_ij_sum = 0;
	c = 0.0;
	for (j = 0; j < k; j++)
		c += Dz[j] * grad[j];
	for (i = 0; i < k; i++) {
		J_ij[i].tot = -D[i] * grad[i] + c * Dzc[i] / Dz2c;
		if (Dz2c_dl > 0)
			J_ij[i].tot = J_ij[i].tot * (1 - dl_s) +
				(-D[i] * grad[i] + c * Dzc_dl[i] / Dz2c_dl) * dl_s * 2 / (visc1 + visc2);
		J_ij[i].tot *= mixf;
		J_ij_sum += z[i] * J_ij[i].tot;
	}

	D = (LDBLE *) free_check_null(D);
	z = (LDBLE *) free_check_null(z);
	Dz = (LDBLE *) free_check_null(Dz);
	Dzc = (LDBLE *) free_check_null(Dzc);
	if (dl_s > 0)
		Dzc_dl = (LDBLE *) free_check_null(Dzc_dl);
	grad = (LDBLE *) free_check_null(grad);

	return(OK);
}

