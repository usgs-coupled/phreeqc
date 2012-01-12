#include <complex>
#include <assert.h>
#include <time.h>

#include "Utils.h"	
#include "Phreeqc.h"
#include "phqalloc.h"
#include "Pressure.h"
#include "Temperature.h"
#include "Parser.h"
#include "cxxMix.h"
#include "Exchange.h"
#include "GasPhase.h"
#include "Reaction.h"
#include "PPassemblage.h"

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_input(void)
/* ---------------------------------------------------------------------- */
{
	int i, j, l;
	char *ptr;
	char token[2 * MAX_LENGTH];
#define LAST_C_KEYWORD 61

	parse_error = 0;
	input_error = 0;
	next_keyword = Keywords::KEY_NONE;
	count_warnings = 0;
/*
 *  Initialize keyword counters
 */
	for (i = 0; i < Keywords::KEY_COUNT_KEYWORDS; i++)
	{
		keycount[i] = 0;
	}
/*
 *  Initialize use and save pointers
 */
	use.init();

	save.solution = FALSE;
	save.mix = FALSE;
	save.reaction = FALSE;
	save.kinetics = FALSE;
	save.pp_assemblage = FALSE;
	save.exchange = FALSE;
	save.surface = FALSE;
	save.gas_phase = FALSE;
	save.ss_assemblage = FALSE;
	title_x = (char *) free_check_null(title_x);

	while ((i =	check_line("Subroutine Read", FALSE, TRUE, TRUE, TRUE)) != KEYWORD)
	{
		/* empty, eof, keyword, print */

		if (i == EOF)
			return (EOF);
		error_string = sformatf(
				"Unknown input, no keyword has been specified.");
		warning_msg(error_string);
	}

	for (;;)
	{
		if (next_keyword > 0)
		{
			if (next_keyword != Keywords::KEY_DATABASE && !reading_database())
			{
				first_read_input = FALSE;
			}
		}
		// mark keyword read
		if (next_keyword > 0 && next_keyword < Keywords::KEY_COUNT_KEYWORDS)
		{
			keycount[next_keyword]++;
		}
		switch (next_keyword)
		{
		case Keywords::KEY_NONE:				/* Have not read line with keyword */
			error_string = sformatf( "Unknown input, no keyword has been specified.");
			warning_msg(error_string);
			while ((j =	check_line("No keyword", FALSE, TRUE, TRUE, TRUE)) != KEYWORD && j != EOF)
			{
				warning_msg(error_string);
			}
			break;
		case Keywords::KEY_END:
			goto END_OF_SIMULATION_INPUT;
		case Keywords::KEY_SOLUTION_SPECIES:				/* Read aqueous model */
			read_species();
			break;
		case Keywords::KEY_SOLUTION_MASTER_SPECIES:		/* Read master species */
			read_master_species();
			break;
		case Keywords::KEY_SOLUTION:				/* Read solution data */
			read_solution();
			solution_sort();
			break;
		case Keywords::KEY_PHASES:
			read_phases();
			break;
		case Keywords::KEY_EQUILIBRIUM_PHASES:
			//read_pure_phases();
			read_pp_assemblage();
			break;
		case Keywords::KEY_REACTION:
			read_reaction();
			break;
		case Keywords::KEY_MIX:
			read_mix();
			break;
		case Keywords::KEY_USE:
			read_use();
			break;
		case Keywords::KEY_SAVE:
			read_save();
			break;
		case Keywords::KEY_EXCHANGE_SPECIES:
			read_exchange_species();
			break;
		case Keywords::KEY_EXCHANGE_MASTER_SPECIES:
			read_exchange_master_species();
			break;
		case Keywords::KEY_EXCHANGE:
			read_exchange();
			break;
		case Keywords::KEY_SURFACE_SPECIES:
			read_surface_species();
			break;
		case Keywords::KEY_SURFACE_MASTER_SPECIES:
			read_surface_master_species();
			break;
		case Keywords::KEY_SURFACE:
			read_surf();
			break;
		case Keywords::KEY_REACTION_TEMPERATURE:
			read_temperature();
			break;
		case Keywords::KEY_INVERSE_MODELING:
			read_inverse();
			break;
		case Keywords::KEY_GAS_PHASE:
			read_gas_phase();
			break;
		case Keywords::KEY_TRANSPORT:
			read_transport();
			break;
		case Keywords::KEY_KNOBS:
			read_debug();
			break;
		case Keywords::KEY_SELECTED_OUTPUT:
			read_selected_output();
			break;
		case Keywords::KEY_PRINT:
			read_print();
			break;
		case Keywords::KEY_TITLE:
			read_title();
			break;
		case Keywords::KEY_ADVECTION:
			read_advection();
			break;
		case Keywords::KEY_KINETICS:
			read_kinetics();
			break;
		case Keywords::KEY_INCREMENTAL_REACTIONS:
			read_incremental_reactions();
			break;
		case Keywords::KEY_RATES:
			read_rates();
			break;
		case Keywords::KEY_SOLUTION_SPREAD:
			read_solution_spread();
			break;
		case Keywords::KEY_USER_PRINT:
			read_user_print();
			break;
		case Keywords::KEY_USER_PUNCH:
			read_user_punch();
			break;
		case Keywords::KEY_SOLID_SOLUTIONS:
			read_solid_solutions();
			break;
		case Keywords::KEY_USER_GRAPH:
#if defined PHREEQ98 
			read_user_graph();
#elif defined MULTICHART
			read_user_graph_handler();
# else
			for (;;)
			{
				j = check_line("Reading user_graph", FALSE, TRUE, TRUE, TRUE);
				if (j == EOF || j == KEYWORD)
				{
					break;
				}
			}
#endif
			break;
		case Keywords::KEY_LLNL_AQUEOUS_MODEL_PARAMETERS:
			read_llnl_aqueous_model_parameters();
			break;
		case Keywords::KEY_DATABASE:
			if (reading_database())
			{
				/* warning_msg("DATABASE is ignored in the database file."); */
			}
			else if (first_read_input == FALSE)
			{
				error_msg("DATABASE must be the first keyword in the input file.", CONTINUE);
				input_error++;
			}
			else
			{
				ptr = line;
				copy_token(token, &ptr, &l);
#if defined(SWIG_SHARED_OBJ)
				warning_msg("DATABASE keyword is ignored by IPhreeqc.");
#else
				user_database = string_duplicate(ptr);
				if (string_trim(user_database) == EMPTY)
				{
					error_msg("DATABASE file name is missing.", CONTINUE);
					input_error++;
					user_database = (char *) free_check_null(user_database);
				}
				first_read_input = FALSE;
#endif
			}
			j = check_line("Reading after DATABASE", FALSE, TRUE, TRUE, TRUE);
			break;
		case Keywords::KEY_NAMED_EXPRESSIONS:
			read_named_logk();
			break;
		case Keywords::KEY_ISOTOPES:
			read_isotopes();
			break;
		case Keywords::KEY_CALCULATE_VALUES:
			read_calculate_values();
			break;
		case Keywords::KEY_ISOTOPE_RATIOS:
			read_isotope_ratios();
			break;
		case Keywords::KEY_ISOTOPE_ALPHAS:
			read_isotope_alphas();
			break;
		case Keywords::KEY_COPY:
			read_copy();
			break;
		case Keywords::KEY_PITZER:
			read_pitzer();
			break;
		case Keywords::KEY_SIT:
			read_sit();
			break;
		case Keywords::KEY_SOLUTION_RAW:		
			read_solution_raw();
			break;
		case Keywords::KEY_EXCHANGE_RAW:		
			//read_exchange_raw();
			Utilities::Rxn_read_raw(Rxn_exchange_map, this);
			break;
		case Keywords::KEY_SURFACE_RAW:		
			read_surface_raw();
			break;
		case Keywords::KEY_EQUILIBRIUM_PHASES_RAW:		
			//read_equilibrium_phases_raw();
			Utilities::Rxn_read_raw(Rxn_pp_assemblage_map, this);
			break;
		case Keywords::KEY_KINETICS_RAW:
			read_kinetics_raw();
			break;
		case Keywords::KEY_SOLID_SOLUTIONS_RAW:
			read_solid_solutions_raw();
			break;
		case Keywords::KEY_GAS_PHASE_RAW:		
			//read_gas_phase_raw();
			Utilities::Rxn_read_raw(Rxn_gas_phase_map, this);
			break;
		case Keywords::KEY_REACTION_RAW:		
			//read_reaction_raw();
			Utilities::Rxn_read_raw(Rxn_reaction_map, this);
			break;
		case Keywords::KEY_MIX_RAW:		
			//read_mix_raw();
			Utilities::Rxn_read_raw(Rxn_mix_map, this);
			break;
		case Keywords::KEY_REACTION_TEMPERATURE_RAW:
			//read_temperature_raw();
			Utilities::Rxn_read_raw(Rxn_temperature_map, this);
			break;
		case Keywords::KEY_DUMP:		
			read_dump();
			break;
		case Keywords::KEY_SOLUTION_MODIFY:
			read_solution_modify();
			break;
		case Keywords::KEY_EQUILIBRIUM_PHASES_MODIFY:
			//read_equilibrium_phases_modify();
			Utilities::Rxn_read_modify(Rxn_pp_assemblage_map, this);
			break;
		case Keywords::KEY_EXCHANGE_MODIFY:
			//read_exchange_modify();
			Utilities::Rxn_read_modify(Rxn_exchange_map, this);
			break;
		case Keywords::KEY_SURFACE_MODIFY:
			read_surface_modify();
			break;
		case Keywords::KEY_SOLID_SOLUTIONS_MODIFY:
			read_solid_solutions_modify();
			break;
		case Keywords::KEY_GAS_PHASE_MODIFY:
			//read_gas_phase_modify();
			Utilities::Rxn_read_modify(Rxn_gas_phase_map, this);
			break;
		case Keywords::KEY_KINETICS_MODIFY:
			read_kinetics_modify();
			break;
		case Keywords::KEY_DELETE:
			read_delete();
			break;
		case Keywords::KEY_RUN_CELLS:
			read_run_cells();
			break;
		case Keywords::KEY_REACTION_MODIFY:
			//read_reaction_modify();
			Utilities::Rxn_read_modify(Rxn_reaction_map, this);
			break;
		//case LAST_C_KEYWORD + 22:		//reaction_temperature_modify
		//	keyword[LAST_C_KEYWORD + 22].keycount++;
		//	read_reaction_temperature_modify();
		//	break;
		case Keywords::KEY_SOLID_SOLUTION_MODIFY:
			read_solid_solutions_modify();
			break;
		case Keywords::KEY_REACTION_PRESSURE:
			read_reaction_pressure();
			break;
		case Keywords::KEY_REACTION_PRESSURE_RAW:
			//read_reaction_pressure_raw();
			Utilities::Rxn_read_raw(Rxn_pressure_map, this);
			break;
#ifdef SKIP
		//case Keywords::KEY_REACTION_PRESSURE_MODIFY:
		//	read_reaction_pressure_modify();
		//	break;
#endif
		default:
			error_msg("Error in keyword switch", STOP);
			break;		  
		}
	}
  END_OF_SIMULATION_INPUT:
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_conc(int n, int count_mass_balance, char *str)
/* ---------------------------------------------------------------------- */
{
	int j, l;
	int alk;
	int count_redox_states;

	char *ptr;
	char token[MAX_LENGTH], token1[MAX_LENGTH];
/*
 *   Set defaults
 */
	/*
	   solution[n]->totals[count_mass_balance].equation_name = NULL;
	   solution[n]->totals[count_mass_balance].phase = NULL;
	   solution[n]->totals[count_mass_balance].phase_si = 0.0;
	   solution[n]->totals[count_mass_balance].units=NULL;
	   solution[n]->totals[count_mass_balance].n_pe=-1;
	   solution[n]->totals[count_mass_balance].as=NULL;
	   solution[n]->totals[count_mass_balance].gfw= 0.0;
	 */
	conc_init(&(solution[n]->totals[count_mass_balance]));

/*
 *   Remove space between "kg" and "solution" or "water" in units
 */
	replace("Kg", "kg", str);
	replace("KG", "kg", str);
	while (replace("kg ", "kg", str) == TRUE);
	ptr = str;
/*
 *   Read master species list for mass balance equation
 */
	token1[0] = '\0';
	count_redox_states = 0;
	while (((j = copy_token(token, &ptr, &l)) == UPPER) ||
		   (token[0] == '[') ||
		   (strcmp_nocase_arg1(token, "ph") == 0) ||
		   (strcmp_nocase_arg1(token, "pe") == 0))
	{
		count_redox_states++;
		replace("(+", "(", token);
		if (count_redox_states > 1)
			strcat(token1, " ");
		strcat(token1, token);
	}
	if (count_redox_states == 0)
	{
		input_error++;
		error_msg
			("No element or master species given for concentration input.",
			 CONTINUE);
		return (ERROR);
	}
	solution[n]->totals[count_mass_balance].description =
		string_hsave(token1);
/*
 *   Determine if reading alkalinity, allow equivalents for units
 */
	str_tolower(token1);
	if (strstr(token1, "alk") == token1)
	{
		alk = TRUE;
	}
	else
	{
		alk = FALSE;
	}
/*
 *   Read concentration
 */

	j = sscanf(token, SCANFORMAT,
			   &(solution[n]->totals[count_mass_balance].input_conc));
	if (j == 0)
	{
		error_string = sformatf(
				"Concentration data error for %s in solution input.", token1);
		error_msg(error_string, CONTINUE);
		return (ERROR);
	}
	if ((j = copy_token(token, &ptr, &l)) == EMPTY)
		return (OK);
/*
 *   Read optional data
 */
	strcpy(token1, token);
/*
 *   Check for units info
 */
	if (check_units(token1, alk, FALSE, solution[n]->units, FALSE) == OK)
	{
		if (check_units(token1, alk, FALSE, solution[n]->units, TRUE) == OK)
		{
			solution[n]->totals[count_mass_balance].units =
				string_hsave(token1);
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
				return (OK);
		}
		else
		{
			return (ERROR);
		}
	}
/*
 *   Check for "as" followed by formula to be used for gfw
 */
	strcpy(token1, token);
	str_tolower(token1);
	if (strcmp(token1, "as") == 0)
	{
		copy_token(token, &ptr, &l);
		solution[n]->totals[count_mass_balance].as = string_hsave(token);
		if ((j = copy_token(token, &ptr, &l)) == EMPTY)
			return (OK);
/*
 *   Check for "gfw" followed by gram formula weight
 */
	}
	else if (strcmp(token1, "gfw") == 0 || strcmp(token1, "gfm") == 0)
	{
		if (copy_token(token, &ptr, &l) != DIGIT)
		{
			error_msg("Expecting gram formula weight.", CONTINUE);
			return (ERROR);
		}
		else
		{
			sscanf(token, SCANFORMAT,
				   &(solution[n]->totals[count_mass_balance].gfw));
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
				return (OK);
		}
	}
/*
 *   Check for redox couple for pe
 */
	if (strcmp_nocase_arg1(token, "pe") == 0)
	{
		solution[n]->totals[count_mass_balance].n_pe =
			pe_data_store(&(solution[n]->pe), token);
		if ((j = copy_token(token, &ptr, &l)) == EMPTY)
			return (OK);
	}
	else if (strstr(token, "/") != NULL)
	{
		if (parse_couple(token) == OK)
		{
			solution[n]->totals[count_mass_balance].n_pe =
				pe_data_store(&(solution[n]->pe), token);
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
				return (OK);
		}
		else
		{
			return (ERROR);
		}
	}
/*
 *   Must have phase
 */
	solution[n]->totals[count_mass_balance].equation_name =
		string_hsave(token);
	if ((j = copy_token(token, &ptr, &l)) == EMPTY)
		return (OK);
/*
 *   Check for saturation index
 */
	j = sscanf(token, SCANFORMAT,
			   &(solution[n]->totals[count_mass_balance].phase_si));
	if (j != 1)
	{
		error_msg("Expected saturation index.", CONTINUE);
		return (ERROR);
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_exchange_species(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read data for exchange species, parse equations
 */
	int i;
	int association;
	char token[MAX_LENGTH];
	char *ptr;
	struct phase *phase_ptr;

	struct species *s_ptr;
	struct elt_list *next_elt;
	struct rxn_token *token_ptr;
	LDBLE exchange_coef;
	LDBLE offset;

	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"no_check",				/* 0 */
		"check",				/* 1 */
		"mb",					/* 2 */
		"mass_balance",			/* 3 */
		"log_k",				/* 4 */
		"logk",					/* 5 */
		"delta_h",				/* 6 */
		"deltah",				/* 7 */
		"analytical_expression",	/* 8 */
		"a_e",					/* 9 */
		"ae",					/* 10 */
		"mole_balance",			/* 11 */
		"gamma",				/* 12 */
		"davies",				/* 13 */
		"offset",				/* 14 */
		"llnl_gamma",			/* 15 */
		"add_logk",				/* 16 */
		"add_log_k",			/* 17 */
		"add_constant",			/* 18 */
		"delta_v",				/* 19 */   
		"deltav",				/* 20 */
		"vm"	/* 21, molar volume, must replace delta_v */
	};
	int count_opt_list = 22;

	association = TRUE;
	s_ptr = NULL;
/*
 *   Read eqn from file and call parser
 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in EXCHANGE_SPECIES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* no_check */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->check_equation = FALSE;
			break;
		case 1:				/* check */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->check_equation = TRUE;
			break;
		case 2:				/* mb */
		case 3:				/* mass_balance */
		case 11:				/* mole_balance */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			count_elts = 0;
			paren_count = 0;
			copy_token(token, &next_char, &i);
			s_ptr->mole_balance = string_hsave(token);
			ptr = token;
			get_secondary_in_species(&ptr, 1.0);
			s_ptr->next_secondary =
				(struct elt_list *) free_check_null(s_ptr->next_secondary);
			s_ptr->next_secondary = elt_list_save();
/* debug
			for (i = 0; i < count_elts; i++) {
				output_msg(sformatf("%s\t%f\n", elt_list[i].elt->name,
					elt_list[i].coef));
			}
 */
			opt_save = OPTION_DEFAULT;
			break;
		case 4:				/* log_k */
		case 5:				/* logk */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_log_k_only(next_char, &s_ptr->logk[0]);
			opt_save = OPTION_DEFAULT;
			break;
		case 6:				/* delta_h */
		case 7:				/* deltah */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_delta_h_only(next_char, &s_ptr->logk[1],
							  &s_ptr->original_units);
			opt_save = OPTION_DEFAULT;
			break;
		case 8:				/* analytical_expression */
		case 9:				/* a_e */
		case 10:				/* ae */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_analytical_expression_only(next_char, &(s_ptr->logk[T_A1]));
			opt_save = OPTION_DEFAULT;
			break;
		case 12:				/* gamma */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->exch_gflag = 2;
			s_ptr->a_f = 0;
			i = sscanf(next_char, SCANFORMAT SCANFORMAT SCANFORMAT, &s_ptr->dha,
					   &s_ptr->dhb, &s_ptr->a_f);
			if (i < 2)
			{
				error_string = sformatf( "Expecting 2 activity-"
						"coefficient parameters, a and b.");
				warning_msg(error_string);
			}
			if (s_ptr->dha == 0 && s_ptr->dhb == 0)
			{
				s_ptr->dhb = 99.9;
				s_ptr->exch_gflag = 1;
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 13:				/* davies */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->exch_gflag = 1;
			s_ptr->dha = 0;
			s_ptr->dhb = 99.9;
			opt_save = OPTION_DEFAULT;
			break;
		case 14:				/* offset */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (sscanf(next_char, SCANFORMAT, &offset) != 1)
			{
				error_msg("No offset for log K given", STOP);
			}
			s_ptr->logk[0] += offset;
			opt_save = OPTION_DEFAULT;
			break;
		case 15:				/* llnl_gamma */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->exch_gflag = 7;	/* llnl D-H */
			i = sscanf(next_char, SCANFORMAT, &s_ptr->dha);
			if (i < 1)
			{
				error_string = sformatf(
						"Expecting activity-coefficient parameter, a.");
				warning_msg(error_string);
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 16:				/* add_logk */
		case 17:				/* add_log_k */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (s_ptr->count_add_logk == 0)
			{
				s_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				s_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(s_ptr->add_logk,
													  (size_t) ((s_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			/* read name */
			if (copy_token(token, &next_char, &i) == EMPTY)
			{
				input_error++;
				error_string = sformatf(
						"Expected the name of a NAMED_EXPRESSION.");
				error_msg(error_string, CONTINUE);
				break;
			}
			s_ptr->add_logk[s_ptr->count_add_logk].name = string_hsave(token);
			/* read coef */
			i = sscanf(next_char, SCANFORMAT,
					   &s_ptr->add_logk[s_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				s_ptr->add_logk[s_ptr->count_add_logk].coef = 1;
			}
			s_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 18:				/* add_constant */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (s_ptr->count_add_logk == 0)
			{
				s_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				s_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(s_ptr->add_logk,
													  (size_t) ((s_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			i = sscanf(next_char, SCANFORMAT,
					   &s_ptr->add_logk[s_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				input_error++;
				error_string = sformatf(
						"Expected the constant to add for log_K definition.");
				error_msg(error_string, CONTINUE);
				break;
			}
			/* set name */
			s_ptr->add_logk[s_ptr->count_add_logk].name =
				string_hsave("XconstantX");
			/* read coef */
			s_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 19:            /* delta_v */
		case 20:            /* deltav */
		case 21:            /* vm, molar volume */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
					"No reaction defined before option, %s.",
				opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_delta_v_only(next_char, &s_ptr->logk[vm0],
				&s_ptr->original_deltav_units);
			opt_save = OPTION_DEFAULT;
			break;

		case OPTION_DEFAULT:
/*
 *   Get exchange species information and parse equation
 */
			s_ptr = NULL;
			if (parse_eq(line, &next_elt, association) == ERROR)
			{
				parse_error++;
				error_msg("Parsing equation.", CONTINUE);
				error_msg(line_save, CONTINUE);
				break;
			}
/*
 *   Get pointer to each species in the reaction, store new species if necessary
 */
			trxn.token[0].s =
				s_store(trxn.token[0].name, trxn.token[0].z, TRUE);
			for (i = 1; i < count_trxn; i++)
			{
				trxn.token[i].s =
					s_store(trxn.token[i].name, trxn.token[i].z, FALSE);
			}
/*
 *   Save element list and carbon, hydrogen, and oxygen in species
 */
			trxn.token[0].s->next_elt = next_elt;
			for (; next_elt->elt != NULL; next_elt++)
			{
				if (strcmp(next_elt->elt->name, "C") == 0)
				{
					trxn.token[0].s->carbon = next_elt->coef;
				}
				if (strcmp(next_elt->elt->name, "H") == 0)
				{
					trxn.token[0].s->h = next_elt->coef;
				}
				if (strcmp(next_elt->elt->name, "O") == 0)
				{
					trxn.token[0].s->o = next_elt->coef;
				}
			}
/*
 *   Find valence of cation from coefficients of reaction components
 *   Changed to be coefficient of exchanger
 */
			exchange_coef = 0.0;
			for (i = 1; i < count_trxn; i++)
			{
				if (trxn.token[i].s->type == EX)
				{
					exchange_coef = trxn.token[i].coef;
				}
			}
			trxn.token[0].s->equiv = exchange_coef;
/*
 *   Malloc space for species reaction
 */
			trxn.token[0].s->rxn = rxn_alloc(count_trxn + 1);
/*
 *   Copy reaction to reaction for species
 */
			token_ptr = trxn.token[0].s->rxn->token;
			for (i = 0; i < count_trxn; i++)
			{
				token_ptr[i].s = trxn.token[i].s;
				token_ptr[i].coef = trxn.token[i].coef;
			}
			token_ptr[i].s = NULL;
/*
 *   Set type for species
 */
			trxn.token[0].s->type = EX;
			s_ptr = trxn.token[0].s;
/*
 *   Set gamma data
 */
			s_ptr->gflag = 4;
			s_ptr->exch_gflag = 3;
			s_ptr->dha = 0.0;
			s_ptr->dhb = 0.0;
			opt_save = OPTION_DEFAULT;
/*
 *  Save as a phase for inverse modeling only
 */
			phase_ptr = phase_store(s_ptr->name);
			if (phase_ptr == NULL)
			{
				input_error++;
				error_string = sformatf( "Copying exchange to phases.");
				error_msg(error_string, CONTINUE);
			}
			phase_ptr->formula = s_ptr->name;
			phase_ptr->check_equation = FALSE;
			phase_ptr->type = EX;
			phase_ptr->next_elt = elt_list_dup(s_ptr->next_elt);
			phase_ptr->rxn = rxn_dup(s_ptr->rxn);
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_exchange(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads exchange data
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int i, j, l;
	int n_user, n_user_end;
	LDBLE conc;
	char *ptr;
	char *description;
	char token[MAX_LENGTH], token1[MAX_LENGTH];
	//struct exchange *exchange_ptr;

	int return_value, opt;
	char *next_char;
	const char *opt_list[] = {
		"equilibrate",			/* 0 */
		"equil",				/* 1 */
		"pitzer_exchange_gammas",	/* 2 */
		"exchange_gammas",		/* 3 */
		"gammas",				/* 4 */
		"equilibrium"           /* 5 */
	};
	int count_opt_list = 6;
/*
 * kin_exch is for exchangers, related to kinetically reacting minerals
 *    they are defined if "sites" is followed by mineral name:
 *    Z	 Manganite		('equi' or 'kine')      0.25
 *    ^Name     ^equi or kinetic mineral ^switch		  ^prop.factor
 */
/*
 *   Read exchange number and description
 */

	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);
/*
 *   Default values + n_user, description
 */
	cxxExchange temp_exchange;
	cxxExchComp *comp = NULL;
	temp_exchange.Set_new_def(true);
	temp_exchange.Set_n_user(n_user);
	temp_exchange.Set_n_user_end(n_user_end);
	temp_exchange.Set_description(description);
	free_check_null(description);
/*
 *   Set use data
 */
	if (use.Get_exchange_in() == FALSE)
	{
		use.Set_exchange_in(true);
		use.Set_n_exchange_user(n_user);
	}
/*
 *   Read exchange data
 */
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in EXCHANGE keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* equilibrate */
		case 1:             /* equil */
		case 5:             /* equilibrium */
			/*
			 *   Read solution to equilibrate with
			 */
			for (;;)
			{
				i = copy_token(token, &next_char, &l);
				if (i == DIGIT)
				{
					int n_solution;
					sscanf(token, "%d", &n_solution);
					temp_exchange.Set_n_solution(n_solution);
					temp_exchange.Set_new_def(true);
					temp_exchange.Set_solution_equilibria(true);
					break;
				}
				if (i == EMPTY)
				{
					error_msg
						("Expected a solution number with which to equilibrate exchanger.",
						 CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
					break;
				}
			}
			break;
		case 2:				/* pitzer_exchange_gammas */
		case 3:				/* exchange_gammas */
		case 4:				/* gammas */
			temp_exchange.Set_pitzer_exchange_gammas(
				get_true_false(next_char, TRUE) == TRUE);
			break;
		case OPTION_DEFAULT:
			ptr = line;
			i = copy_token(token, &ptr, &l);
			/*
			 *   Species formula is stored in token
			 */
			if (i != UPPER && token[0] != '[')
			{ /* maybe a bit more clear? */
				error_string = sformatf(
					"Expected exchanger name to begin with a capital letter, but found:\n %s",
					line_save);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			delete comp;
			comp = new cxxExchComp;
			comp->Set_formula(token);
			//exchange[n].comps[count_comps].formula = string_hsave(token);
			prev_next_char = ptr;
			i = copy_token(token1, &ptr, &l);
			if (i == DIGIT)
			{
				/*
				 *   Read exchange concentration
				 */

				/* exchanger conc. is read directly .. */
				if (sscanf(token1, SCANFORMAT, &conc) < 1)
				{
				error_string = sformatf(
					"Expected concentration of exchanger, but found:\n %s",
					line_save);
					error_msg(error_string, CONTINUE);
					input_error++;
					break;
				}
				prev_next_char = ptr;
				j = copy_token(token1, &ptr, &l);
				if (j == UPPER || j == LOWER)
				{
					comp->Set_rate_name(token1);
					if (copy_token(token1, &ptr, &l) != DIGIT)
					{
						error_string = sformatf(
							"Expected a coefficient to relate exchange to kinetic reaction, but found:\n %s",
							prev_next_char);
						error_msg(error_string, CONTINUE);
						input_error++;
						break;
					}
					LDBLE p;
					sscanf(token1, SCANFORMAT, &p);
					comp->Set_phase_proportion(p);
				}
				/*
				 *   Read equilibrium phase name or kinetics rate name
				 */
			}
			else if (i != EMPTY)
			{

				/* exchanger conc. is related to mineral or kinetics */
				comp->Set_phase_name(token1);
				prev_next_char = ptr;
				j = copy_token(token1, &ptr, &l);
				if (j != DIGIT)
				{
					if (token1[0] == 'K' || token1[0] == 'k')
					{
						comp->Set_rate_name(comp->Get_phase_name().c_str());
						comp->Set_phase_name("");
					}
					else if (token1[0] != 'E' && token1[0] != 'e')
					{
						error_string = sformatf(
							"Character string expected to be 'equilibrium_phase' or 'kinetics'\n to relate exchange to mineral or kinetic reaction, but found:\n %s",
							prev_next_char);
						error_msg(error_string, CONTINUE);
						input_error++;
						break;
					}
					prev_next_char = ptr;
					j = copy_token(token1, &ptr, &l);
				}


				if (j != DIGIT)
				{
					error_string = sformatf(
						"Expected a coefficient to relate exchanger to mineral or kinetic reaction, but found:\n %s",
						prev_next_char);
					error_msg(error_string, CONTINUE);
					input_error++;
					break;
				}
				LDBLE p;
				sscanf(token1, SCANFORMAT, &p);
				comp->Set_phase_proportion(p);
				/* real conc must be defined in tidy_model */
				conc = 1.0;
			}
			else
			{
				error_msg
					("Expected concentration of exchanger, mineral name, or kinetic reaction name.",
					 CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				break;
			}
			/*
			 *   Accumulate elements in elt_list
			 */
			count_elts = 0;
			paren_count = 0;
			ptr = token;
			get_elts_in_species(&ptr, conc);
			/*
			 *   save formula for adjusting number of exchange sites
			 */
			ptr = token;
			LDBLE z;
			get_token(&ptr, token1, &z, &l);
			comp->Set_formula_z(z);
			comp->Set_formula_totals(elt_list_NameDouble()); 
			/*
			 *   Save elt_list
			 */
			comp->Set_moles(conc);
			comp->Set_totals(elt_list_NameDouble());
			comp->Set_charge_balance(0.0);
			temp_exchange.Get_exchComps()[comp->Get_formula()] = *comp;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	Rxn_exchange_map[n_user] = temp_exchange;

	//if (n_user_end > n_user)
	//{
	//	for (int i = n_user + 1; i <= n_user_end; i++)
	//	{
	//		Utilities::Rxn_copy(Rxn_exchange_map, n_user, i);
	//	}
	//}
	delete comp;
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_exchange_master_species(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads master species data from data file or input file
 */
	int j, l;
	char *ptr, *ptr1;
	LDBLE l_z;
	struct element *elts_ptr;
	struct species *s_ptr;
	char token[MAX_LENGTH], token1[MAX_LENGTH];
	for (;;)
	{
		j = check_line("Exchange species equation", FALSE, TRUE, TRUE, TRUE);
		if (j == EOF || j == KEYWORD)
		{
			break;
		}
/*
 *   Get element name with valence, allocate space, store
 */
		ptr = line;
/*
 *   Get element name and save pointer to character string
 */
		if (copy_token(token, &ptr, &l) != UPPER && token[0] != '[')
		{
			parse_error++;
			error_msg("Reading element for master species.", CONTINUE);
			error_msg(line_save, CONTINUE);
			continue;
		}
		/*
		   if (token[0] == '[') {
		   ptr1 = token;
		   get_elt(&ptr, element, &l);
		   strcpy(token, element);
		   }
		 */
		replace("(+", "(", token);
/*
 *   Delete master if it exists
 */
		master_delete(token);
/*
 *   Increase pointer array, if necessary,  and malloc space
 */
		if (count_master >= max_master)
		{
			space((void **) ((void *) &master), count_master + 1,
				  &max_master, sizeof(struct master *));
		}
		master[count_master] = master_alloc();
/*
 *   Set type to EX
 */
		master[count_master]->type = EX;
/*
 *   Save element name
 */
		master[count_master]->elt = element_store(token);
/*
 *   Save pointer to species data for master species
 */
		if ((copy_token(token, &ptr, &l) != UPPER) &&
			token[0] != '[' && (strcmp_nocase_arg1(token, "e-") != 0))
		{
			parse_error++;
			error_msg("Reading master species name.", CONTINUE);
			error_msg(line_save, CONTINUE);
			continue;
		}
		s_ptr = s_search(token);
		if (s_ptr != NULL)
		{
			master[count_master]->s = s_ptr;
		}
		else
		{
			ptr1 = token;
			get_token(&ptr1, token1, &l_z, &l);
			master[count_master]->s = s_store(token1, l_z, FALSE);
		}
/*
 *   MAKE LISTS OF PRIMARY AND SECONDARY MASTER SPECIES
 */
		master[count_master]->primary = TRUE;
		if (strcmp(master[count_master]->elt->name, "E") != 0)
		{
			elts_ptr = element_store(master[count_master]->elt->name);
			elts_ptr->gfw = 0.0;
		}

		count_master++;
		if (count_master >= max_master)
		{
			space((void **) ((void *) &master), count_master, &max_master,
				  sizeof(struct master *));
		}
	}
	return (j);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_gas_phase(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads gas phase data
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int i, j, l;
	int n_user, n_user_end;
	char *ptr;
	char *description;
	char token[MAX_LENGTH];
	//struct gas_phase *gas_phase_ptr;
	cxxGasPhase temp_gas_phase;
	int return_value, opt;
	char *next_char;
	const char *opt_list[] = {
		"pressure",				/* 0 */
		"volume",				/* 1 */
		"temp",					/* 2 */
		"temperature",			/* 3 */
		"fixed_pressure",		/* 4 */
		"fixed_volume",			/* 5 */
		"equilibrium",			/* 6 */
		"equilibrate",			/* 7 */
		"equil"					/* 8 */
	};
	int count_opt_list = 9;
/*
 *   Read gas_phase number
 */
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);

	temp_gas_phase.Set_n_user(n_user);
	temp_gas_phase.Set_n_user_end(n_user_end);
	temp_gas_phase.Set_description(description);
	temp_gas_phase.Set_new_def(true);
	free_check_null(description);
/*
 *   Set use data to first read
 */
	if (use.Get_gas_phase_in() == FALSE)
	{
		use.Set_gas_phase_in(true);
		use.Set_n_gas_phase_user(n_user);
	}
/*
 *   Read phases
 */
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in GAS_PHASE keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* pressure */
			sscanf(next_char, SCANFORMAT, &dummy);
			temp_gas_phase.Set_total_p(dummy);
			break;
		case 1:				/* Volume */
			sscanf(next_char, SCANFORMAT, &dummy);
			temp_gas_phase.Set_volume(dummy);
			break;
		case 2:				/* Temperature */
		case 3:
			j = sscanf(next_char, SCANFORMAT, &dummy);
			if (j == 1)
			{
				temp_gas_phase.Set_temperature(dummy + 273.15);
			}
			break;
		case 4:				/* fixed_pressure */
			temp_gas_phase.Set_type(cxxGasPhase::GP_PRESSURE);
			break;
		case 5:				/* fixed_volume */
			temp_gas_phase.Set_type(cxxGasPhase::GP_VOLUME);
			break;
		case 6:				/* equilibrate */
		case 7:				/* equilibrium */
		case 8:				/* equil */
/*
 *   Read solution to equilibrate with
 */
			for (;;)
			{
				i = copy_token(token, &next_char, &l);
				if (i == DIGIT)
				{
					sscanf(token, "%d", &l);
					temp_gas_phase.Set_n_solution(l);
					temp_gas_phase.Set_new_def(true);
					temp_gas_phase.Set_solution_equilibria(true);
					break;
				}
				if (i == EMPTY)
				{
					error_msg
						("Expected a solution number with which to equilibrate gas phase.",
						 CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
					break;
				}
			}
			break;
		case OPTION_DEFAULT:
			{
				cxxGasComp temp_comp;
				/*
				*   Read name
				*/
				ptr = line;
				copy_token(token, &ptr, &l);
				temp_comp.Set_phase_name(token);
				if ((j = copy_token(token, &ptr, &l)) == EMPTY)
				{
					temp_comp.Set_p_read(NAN);
					temp_gas_phase.Get_gas_comps().push_back(temp_comp);
					break;
				}
				/*
				*   Read initial partial pressure of gas
				*/

				j = sscanf(token, SCANFORMAT, &dummy);
				
				if (j != 1)
				{
					error_msg("Expected partial pressure of gas in gas phase.",
						CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
				}
				else
				{
					temp_comp.Set_p_read(dummy);
					temp_gas_phase.Get_gas_comps().push_back(temp_comp);
				}
			}
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}

	// sort components
	std::map<std::string, cxxGasComp> gc;
	for (size_t i = 0; i < temp_gas_phase.Get_gas_comps().size(); i++)
	{
		cxxGasComp *gc_ptr = &(temp_gas_phase.Get_gas_comps()[i]);
		gc[gc_ptr->Get_phase_name()] = *gc_ptr;
	}
	std::vector<cxxGasComp> vgc;
	std::map<std::string, cxxGasComp>::iterator it = gc.begin();
	for ( ; it != gc.end(); it++)
	{
		vgc.push_back(it->second);
	}
	temp_gas_phase.Set_gas_comps(vgc);


	Rxn_gas_phase_map[n_user] = temp_gas_phase;

	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_inverse(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads data for mass_balance calculations
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int n, j;
	int n_user, n_user_end;
	char *ptr;
	char *description;
	LDBLE range_max, inv_tol, water_uncertainty;

	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"solutions",			/* 0 */
		"uncertainty",			/* 1 */
		"uncertainties",		/* 2 */
		"balances",				/* 3 */
		"phase_data",			/* 4 */
		"range",				/* 5 */
		"minimal",				/* 6 */
		"minimum",				/* 7 */
		"balance",				/* 8 */
		"bal",					/* 9 */
		"sol",					/* 10 */
		"phases",				/* 11 */
		"ranges",				/* 12 */
		"tolerance",			/* 13 */
		"u_water",				/* 14 */
		"uncertainty_water",	/* 15 */
		"force",				/* 16 */
		"force_solution",		/* 17 */
		"force_solutions",		/* 18 */
		"isotopes",				/* 19 */
		"mineral_water",		/* 20 */
		"phase",				/* 21 */
		"multiple_precision",	/* 22 */
		"mp_tolerance",			/* 23 */
		"censor_mp",			/* 24 */
		"lon_netpath",			/* 25 */
		"pat_netpath"			/* 26 */
	};
	int count_opt_list = 27;

	ptr = line;
/*
 *   Read solution number and description
 */
	read_number_description(ptr, &n_user, &n_user_end, &description);
/*
 *   Malloc space for solution data
 */
	if (inverse_search(n_user, &n) != NULL)
	{
		inverse_delete(n);
	}
	inverse_alloc();
	n = count_inverse - 1;
/*
 *   Initialize structure and use
 */
	inverse[n].new_def = TRUE;
	inverse[n].n_user = n_user;
	inverse[n].range = FALSE;
	inverse[n].range_max = 1000.;
	inverse[n].tolerance = 1e-10;
	inverse[n].minimal = FALSE;
	inverse[n].description = description;
	inverse[n].count_uncertainties = 1;
	inverse[n].uncertainties[0] = 0.05;
	inverse[n].count_ph_uncertainties = 1;
	inverse[n].ph_uncertainties[0] = 0.05;
	inverse[n].water_uncertainty = 0.0;
	inverse[n].mineral_water = TRUE;
	inverse[n].mp = FALSE;
	inverse[n].mp_tolerance = 1e-12;
	inverse[n].mp_censor = 1e-20;
	inverse[n].netpath = NULL;
	inverse[n].pat = NULL;
/*
 *   Read data for inverse modeling
 */
	opt_save = OPTION_ERROR;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		opt_save = opt;
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in INVERSE_MODELING keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* solutions */
		case 10:				/* solution */
			inverse[n].solns =
				read_list_ints(&next_char, &inverse[n].count_solns, TRUE);
			opt_save = OPTION_ERROR;
			break;
		case 1:				/* uncertainty */
		case 2:				/* uncertainties */
			inverse[n].uncertainties =
				(LDBLE *) free_check_null(inverse[n].uncertainties);
			inverse[n].uncertainties =
				read_list_doubles(&next_char,
								  &inverse[n].count_uncertainties);
			opt_save = OPTION_ERROR;
			break;
		case 3:				/* balances */
		case 8:				/* balance */
		case 9:				/* bal */
			read_inv_balances(&(inverse[n]), next_char);
			break;
		case 4:				/* phase_data */
		case 11:				/* phases */
		case 21:				/* phase */
			read_inv_phases(&(inverse[n]), next_char);
			break;
		case 5:				/* range */
		case 12:				/* ranges */
			inverse[n].range = TRUE;
			j = sscanf(next_char, SCANFORMAT, &range_max);
			if (j == 1)
			{
				inverse[n].range_max = range_max;
			}
			opt_save = OPTION_ERROR;
			break;
		case 6:				/* minimal */
		case 7:				/* minimum */
			inverse[n].minimal = TRUE;
			opt_save = OPTION_ERROR;
			break;
		case 13:				/* tolerance */
			j = sscanf(next_char, SCANFORMAT, &inv_tol);
			if (j == 1)
			{
				inverse[n].tolerance = inv_tol;
			}
			opt_save = OPTION_ERROR;
			break;
		case 14:				/* u_water */
		case 15:				/* uncertainty_water */
			j = sscanf(next_char, SCANFORMAT, &water_uncertainty);
			if (j == 1)
			{
				inverse[n].water_uncertainty = water_uncertainty;
			}
			opt_save = OPTION_ERROR;
			break;
		case 16:				/* force */
		case 17:				/* force_solution */
		case 18:				/* force_solutions */
			inverse[n].force_solns =
				(int *) free_check_null(inverse[n].force_solns);
			inverse[n].force_solns =
				read_list_t_f(&next_char, &inverse[n].count_force_solns);
			opt_save = OPTION_ERROR;
			break;
		case 19:				/* isotope values */
			read_inv_isotopes(&(inverse[n]), next_char);
			break;
		case 20:				/* mineral_water */
			inverse[n].mineral_water = get_true_false(next_char, TRUE);
			break;
		case 22:				/* multiple_precision */
			inverse[n].mp = get_true_false(next_char, TRUE);
			break;
		case 23:				/* mp_tolerance */
			j = sscanf(next_char, SCANFORMAT, &inv_tol);
			if (j == 1)
			{
				inverse[n].mp_tolerance = fabs(inv_tol);
			}
			break;
		case 24:				/* censor_mp */
			j = sscanf(next_char, SCANFORMAT, &inv_tol);
			if (j == 1)
			{
				inverse[n].mp_censor = fabs(inv_tol);
			}
			break;
		case 25:				/* lon_netpath */
			/*copy_token(file_name, &next_char, &l); */
			if (string_trim(next_char) != EMPTY)
			{
				inverse[n].netpath = string_hsave(next_char);
			}
			else
			{
				inverse[n].netpath = string_hsave("netpath");
			}
			opt_save = OPTION_ERROR;
			break;
		case 26:				/* pat_netpath */
			/*copy_token(file_name, &next_char, &l); */
			if (string_trim(next_char) != EMPTY)
			{
				inverse[n].pat = string_hsave(next_char);
			}
			else
			{
				inverse[n].pat = string_hsave("netpath");
			}
			opt_save = OPTION_ERROR;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
/*
 *  Default: soln 1 -> soln 2
 */
	if (inverse[n].count_solns == 0)
	{
		inverse[n].solns = (int *) PHRQ_malloc(2 * sizeof(int));
		if (inverse[n].solns == NULL)
			malloc_error();
		inverse[n].solns[0] = 1;
		inverse[n].solns[1] = 2;
		inverse[n].count_solns = 2;
	}
/*
 *   Sort isotopes
 */
	if (inverse[n].count_isotopes > 0)
	{
		qsort(inverse[n].isotopes,
			  (size_t) inverse[n].count_isotopes,
			  (size_t) sizeof(struct inv_isotope), inverse_isotope_compare);
	}

	if (inverse[n].count_i_u > 0)
	{
		qsort(inverse[n].i_u,
			  (size_t) inverse[n].count_i_u,
			  (size_t) sizeof(struct inv_isotope), inverse_isotope_compare);
	}

	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_inv_balances(struct inverse *inverse_ptr, char *ptr)
/* ---------------------------------------------------------------------- */
{
	int j, l, count;
	char token[MAX_LENGTH];
/*
 *   Read element name
 */
	j = copy_token(token, &ptr, &l);
	if (j == EMPTY)
	{
		return (OK);
	}
	else if (j == LOWER && strcmp_nocase_arg1(token, "ph") != 0)
	{
		error_msg("Expecting element name.", CONTINUE);
		error_msg(line_save, CONTINUE);
		input_error++;
	}
	else if (strcmp_nocase_arg1(token, "ph") != 0)
	{
		inverse_ptr->elts =
			(struct inv_elts *) PHRQ_realloc(inverse_ptr->elts,
											 (size_t) (inverse_ptr->
													   count_elts +
													   1) *
											 sizeof(struct inv_elts));
		if (inverse_ptr->elts == NULL)
			malloc_error();
		replace("(+", "(", token);
		inverse_ptr->elts[inverse_ptr->count_elts].name = string_hsave(token);
/*
 *   Read element uncertainties
 */
		inverse_ptr->elts[inverse_ptr->count_elts].uncertainties =
			read_list_doubles(&ptr, &count);
		inverse_ptr->elts[inverse_ptr->count_elts].count_uncertainties =
			count;
		inverse_ptr->count_elts++;
	}
	else if (strcmp_nocase_arg1(token, "ph") == 0)
	{
		inverse_ptr->ph_uncertainties =
			(LDBLE *) free_check_null(inverse_ptr->ph_uncertainties);
		inverse_ptr->ph_uncertainties = read_list_doubles(&ptr, &count);
		inverse_ptr->count_ph_uncertainties = count;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_inv_isotopes(struct inverse *inverse_ptr, char *ptr)
/* ---------------------------------------------------------------------- */
{
	int i, j, l, l1, l2, count;
	LDBLE isotope_number;
	char token[MAX_LENGTH], token1[MAX_LENGTH];
	char *ptr1, *ptr2;
	const char * redox_name, *element_name;
/*
 *   Read element name
 */
	ptr1 = ptr;
	j = copy_token(token, &ptr1, &l);
/*
 *   ptr1 is start of uncertainties
 */
	if (j == EMPTY)
	{
		return (OK);
	}
	else if (j != DIGIT)
	{
		error_msg("Expecting isotope to begin with isotope number.",
				  CONTINUE);
		error_msg(line_save, CONTINUE);
		input_error++;
		return (ERROR);
	}
/*
 *   Read isotope name
 */
	ptr2 = token;
	get_num(&ptr2, &isotope_number);
	if (ptr2[0] == '\0' || isupper((int) ptr2[0]) == FALSE)
	{
		error_msg("Expecting element name.", CONTINUE);
		error_msg(line_save, CONTINUE);
		input_error++;
		return (ERROR);
	}

	/* redox state name with parentheses */
	redox_name = string_hsave(ptr2);

	copy_token(token, &ptr2, &l1);
	replace("(", " ", token);
	ptr2 = token;

	/* element name, without parentheses */
	copy_token(token1, &ptr2, &l2);
	element_name = string_hsave(token1);

/*
 *  add element name to inv_ptr->isotopes
 */
	for (i = 0; i < inverse_ptr->count_isotopes; i++)
	{
		if (element_name == inverse_ptr->isotopes[i].elt_name)
			break;
	}
	if (i == inverse_ptr->count_isotopes)
	{
		inverse_ptr->isotopes =
			(struct inv_isotope *) PHRQ_realloc(inverse_ptr->isotopes,
												(size_t) (inverse_ptr->
														  count_isotopes +
														  1) *
												sizeof(struct inv_isotope));
		if (inverse_ptr->isotopes == NULL)
			malloc_error();
		inverse_ptr->isotopes[inverse_ptr->count_isotopes].isotope_number =
			isotope_number;
		inverse_ptr->isotopes[inverse_ptr->count_isotopes].elt_name =
			element_name;
		inverse_ptr->isotopes[inverse_ptr->count_isotopes].uncertainties =
			(LDBLE *) PHRQ_malloc((size_t) sizeof(LDBLE));
		if (inverse_ptr->isotopes[inverse_ptr->count_isotopes].
			uncertainties == NULL)
			malloc_error();
		inverse_ptr->count_isotopes++;
	}
/*
 *  add redox state name to inv_ptr->i_u
 */
	inverse_ptr->i_u =
		(struct inv_isotope *) PHRQ_realloc(inverse_ptr->i_u,
											(size_t) (inverse_ptr->
													  count_i_u +
													  1) *
											sizeof(struct inv_isotope));
	if (inverse_ptr->i_u == NULL)
		malloc_error();
	inverse_ptr->i_u[inverse_ptr->count_i_u].elt_name = redox_name;
	inverse_ptr->i_u[inverse_ptr->count_i_u].isotope_number = isotope_number;
/*
 *   Read isotope uncertainties
 */
	inverse_ptr->i_u[inverse_ptr->count_i_u].uncertainties =
		read_list_doubles(&ptr1, &count);
	inverse_ptr->i_u[inverse_ptr->count_i_u].count_uncertainties = count;
	inverse_ptr->count_i_u++;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_inv_phases(struct inverse *inverse_ptr, char *ptr)
/* ---------------------------------------------------------------------- */
{
	int j, l;
	int count_isotopes;
	char token[MAX_LENGTH], token1[MAX_LENGTH];
	char *ptr1;
	struct isotope *isotopes;
/*
 *   Read phase name
 */
	j = copy_token(token, &ptr, &l);
	if (j == EMPTY)
		return (OK);
	inverse_ptr->phases =
		(struct inv_phases *) PHRQ_realloc(inverse_ptr->phases,
										   (size_t) (inverse_ptr->
													 count_phases +
													 1) *
										   sizeof(struct inv_phases));
	if (inverse_ptr->phases == NULL)
		malloc_error();
	inverse_ptr->phases[inverse_ptr->count_phases].name = string_hsave(token);
/*
 *   Read constraint, force, and isotopes
 */
	inverse_ptr->phases[inverse_ptr->count_phases].constraint = EITHER;
	inverse_ptr->phases[inverse_ptr->count_phases].force = FALSE;
	count_isotopes = 0;
	isotopes = (struct isotope *) PHRQ_malloc(sizeof(struct isotope));
	if (isotopes == NULL)
		malloc_error();

	for (;;)
	{
		j = copy_token(token, &ptr, &l);
		if (j == EMPTY)
			break;
		strcpy(token1, token);
		str_tolower(token1);
		if (token1[0] == 'p')
		{
			inverse_ptr->phases[inverse_ptr->count_phases].constraint =
				PRECIPITATE;
		}
		else if (token1[0] == 'd')
		{
			inverse_ptr->phases[inverse_ptr->count_phases].constraint =
				DISSOLVE;
		}
		else if (token[0] == 'f')
		{
			inverse_ptr->phases[inverse_ptr->count_phases].force = TRUE;
		}
		else if (j == DIGIT)
		{
/*
 *   read isotope data
 */
			isotopes =
				(struct isotope *) PHRQ_realloc(isotopes,
												(size_t) (count_isotopes +
														  1) *
												sizeof(struct isotope));
			if (isotopes == NULL)
				malloc_error();
			ptr1 = token;

			/* isotope number */
			get_num(&ptr1, &(isotopes[count_isotopes].isotope_number));
			if (ptr1[0] == '\0' || isupper((int) ptr1[0]) == FALSE)
			{
				error_string = sformatf( "Expecting element name: %s.", ptr1);
				error_msg(error_string, CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				break;
			}

			/* element name */
			isotopes[count_isotopes].elt_name = string_hsave(ptr1);

			/* ratio */
			j = copy_token(token, &ptr, &l);
			if (j != DIGIT)
			{
				error_msg("Expecting isotope ratio for phase.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				break;
			}
			sscanf(token, SCANFORMAT, &(isotopes[count_isotopes].ratio));

			/* read and store isotope ratio uncertainty */
			prev_next_char = ptr;
			if (copy_token(token, &ptr, &l) != DIGIT)
			{
				input_error++;
				error_string = sformatf(
					"Expected numeric value for uncertainty in isotope ratio, but found:\n %s",
					prev_next_char);
				error_msg(error_string, CONTINUE);
				continue;
			}
			sscanf(token, SCANFORMAT,
				   &(isotopes[count_isotopes].ratio_uncertainty));

			count_isotopes++;
		}
		else
		{
			error_string = sformatf(
					"Unknown option for inverse modeling phase.");
			warning_msg(error_string);
		}
	}
	if (count_isotopes > 0)
	{
		inverse_ptr->phases[inverse_ptr->count_phases].isotopes = isotopes;
		inverse_ptr->phases[inverse_ptr->count_phases].count_isotopes =
			count_isotopes;
	}
	else
	{
		inverse_ptr->phases[inverse_ptr->count_phases].isotopes = NULL;
		inverse_ptr->phases[inverse_ptr->count_phases].count_isotopes = 0;
		isotopes = (struct isotope *) free_check_null(isotopes);
	}
	inverse_ptr->count_phases++;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_kinetics(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads kinetics data
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
/*
 *   Read kinetics
 */
	int i, j, k, l, count_comps, count_steps, count_list;
	char *ptr;
	char *description;
	char token[MAX_LENGTH];
	int n;
	int n_user, n_user_end;
	struct kinetics *kinetics_ptr;
	struct kinetics_comp *kinetics_comp_ptr;
	LDBLE step, coef;

	int return_value, opt;
	char *next_char;
	const char *opt_list[] = {
		"tol",					/* 0 */
		"m",					/* 1 */
		"m0",					/* 2 */
		"parms",				/* 3 */
		"formula",				/* 4 */
		"steps",				/* 5 */
		"step_divide",			/* 6 */
		"parameters",			/* 7 */
		"runge-kutta",			/* 8 */
		"runge_kutta",			/* 9 */
		"rk",					/* 10 */
		"bad_step_max",			/* 11 */
		"cvode",				/* 12 */
		"cvode_steps",			/* 13 */
		"cvode_order",			/* 14 */
		"time_steps"			/* 15 */
	};
	int count_opt_list = 16;

/*
 *   Read kinetics number
 */
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);

/*
 *   Find space for kinetics data
 */
	kinetics_ptr = kinetics_search(n_user, &n, FALSE);
	if (kinetics_ptr != NULL)
	{
		kinetics_free(kinetics_ptr);
	}
	else
	{
		space((void **) ((void *) &kinetics), count_kinetics, &max_kinetics,
			  sizeof(struct kinetics));
		n = count_kinetics++;
	}
/*
 *   Set use data to first read
 */
	if (use.Get_kinetics_in() == FALSE)
	{
		use.Set_kinetics_in(true);
		use.Set_n_kinetics_user(n_user);
	}
/*
 *   Initialize
 */
	kinetics_init(&(kinetics[n]), n_user, n_user_end, description);
	free_check_null(description);
	kinetics_ptr = &kinetics[n];

	count_steps = 0;
/*
 *   Read kinetics data
 */
	return_value = UNKNOWN;
	kinetics_comp_ptr = NULL;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:	/* allocate space, read new name */
			count_comps = kinetics_ptr->count_comps++;
			kinetics_ptr->comps =
				(struct kinetics_comp *) PHRQ_realloc(kinetics_ptr->comps,
													  (size_t) (count_comps +
																1) *
													  sizeof(struct
															 kinetics_comp));
			if (kinetics_ptr->comps == NULL)
				malloc_error();
			kinetics_ptr->comps[count_comps].moles = 0;
			ptr = line;
			copy_token(token, &ptr, &l);
			kinetics_ptr->comps[count_comps].rate_name = string_hsave(token);
			kinetics_ptr->comps[count_comps].tol = 1e-8;
			kinetics_ptr->comps[count_comps].m0 = -1.0;
			kinetics_ptr->comps[count_comps].m = -1.0;
			kinetics_ptr->comps[count_comps].count_c_params = 0;

			kinetics_comp_ptr = &kinetics_ptr->comps[count_comps];
			kinetics_comp_ptr->d_params =
				(LDBLE *) PHRQ_malloc(sizeof(LDBLE));
			if (kinetics_comp_ptr->d_params == NULL)
				malloc_error();
			kinetics_comp_ptr->count_d_params = 0;

			kinetics_comp_ptr->c_params =
				(const char **) PHRQ_malloc(sizeof(char *));
			if (kinetics_comp_ptr->c_params == NULL)
				malloc_error();
			kinetics_comp_ptr->count_c_params = 0;

			kinetics_comp_ptr->count_list = 0;
			kinetics_comp_ptr->list =
				(struct name_coef *) PHRQ_malloc(sizeof(struct name_coef));
			if (kinetics_comp_ptr->list == NULL)
				malloc_error();
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in KINETICS keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* tolerance */
			if (kinetics_comp_ptr == NULL)
			{
				error_string = sformatf( "No rate name has been defined.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			else
			{
				prev_next_char = next_char;
				if (copy_token(token, &next_char, &l) == DIGIT)
				{
					kinetics_comp_ptr->tol = strtod(token, &ptr);
				}
				else
				{
					error_string = sformatf(
						"Expecting numerical value for tolerance, but found:\n %s",
						prev_next_char);
					error_msg(error_string, CONTINUE);
					input_error++;
				}
			}
			break;
		case 1:				/* m */
			if (kinetics_comp_ptr == NULL)
			{
				error_string = sformatf( "No rate name has been defined.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			else
			{
				prev_next_char = next_char;
				if (copy_token(token, &next_char, &l) == DIGIT)
				{
					kinetics_comp_ptr->m = strtod(token, &ptr);
				}
				else
				{
					error_string = sformatf(
						"Expecting numerical value for moles of reactant (m), but found:\n %s",
						prev_next_char);
					error_msg(error_string, CONTINUE);
					input_error++;
				}
			}
			break;
		case 2:				/* m0 */
			if (kinetics_comp_ptr == NULL)
			{
				error_string = sformatf( "No rate name has been defined.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			else
			{
				prev_next_char = next_char;
				if (copy_token(token, &next_char, &l) == DIGIT)
				{
					kinetics_comp_ptr->m0 = strtod(token, &ptr);
				}
				else
				{
					error_string = sformatf(
					"Expecting numerical value for initial moles of reactant (m0), but found:\n %s",
					prev_next_char);
				error_msg(error_string, CONTINUE);
					input_error++;
				}
			}
			break;
		case 3:				/* parms */
		case 7:				/* parameters */
			if (kinetics_comp_ptr == NULL)
			{
				error_string = sformatf( "No rate name has been defined.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			else
			{
				while ((j = copy_token(token, &next_char, &l)) != EMPTY)
				{
					/*
					 *   Store a LDBLE parameter
					 */
					if (j == DIGIT)
					{
						kinetics_comp_ptr->d_params =
							(LDBLE *) PHRQ_realloc(kinetics_comp_ptr->
												   d_params,
												   (size_t)
												   (kinetics_comp_ptr->
													count_d_params +
													1) * sizeof(LDBLE));
						if (kinetics_comp_ptr->d_params == NULL)
							malloc_error();
						kinetics_comp_ptr->d_params[kinetics_comp_ptr->
													count_d_params] =
							strtod(token, &ptr);
						kinetics_comp_ptr->count_d_params++;
					}
					else
					{
						/*
						 *   Store a character parameter
						 */
						kinetics_comp_ptr->c_params =
							(const char **) PHRQ_realloc(kinetics_comp_ptr->
												   c_params,
												   (size_t)
												   (kinetics_comp_ptr->
													count_c_params +
													1) * sizeof(char *));
						if (kinetics_comp_ptr->c_params == NULL)
							malloc_error();
						kinetics_comp_ptr->c_params[kinetics_comp_ptr->
													count_c_params] =
							string_hsave(token);
						kinetics_comp_ptr->count_c_params++;
					}
				}
			}
			break;
		case 4:				/* formula */
			if (kinetics_comp_ptr == NULL)
			{
				error_string = sformatf( "No rate name has been defined.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			else
			{
				/*
				 *   Store reactant name, default coefficient
				 */
				ptr = next_char;
				while (copy_token(token, &ptr, &l) != EMPTY)
				{
					if (isalpha((int) token[0]) || (token[0] == '(')
						|| (token[0] == '['))
					{
						count_list = kinetics_comp_ptr->count_list++;
						kinetics_comp_ptr->list =
							(struct name_coef *)
							PHRQ_realloc(kinetics_comp_ptr->list,
										 (size_t) (count_list +
												   1) *
										 sizeof(struct name_coef));
						if (kinetics_comp_ptr->list == NULL)
							malloc_error();
						kinetics_comp_ptr->list[count_list].name =
							string_hsave(token);
						kinetics_comp_ptr->list[count_list].coef = 1.0;
					}
					else
					{
						/*
						 *   Store relative coefficient
						 */
						j = sscanf(token, SCANFORMAT, &coef);
						if (j == 1)
						{
							count_list = kinetics_comp_ptr->count_list - 1;
							kinetics_comp_ptr->list[count_list].coef = coef;
						}
						else
						{
							error_msg
								("Reading relative coefficient of reactant.",
								 CONTINUE);
							error_msg(line_save, CONTINUE);
							input_error++;
						}
					}
				}
			}
			break;
		case 5:				/* steps */
		case 15:			/* time_steps */
			/*
			 *   Read one or more kinetics time increments
			 */
			while ((j = copy_token(token, &next_char, &l)) == DIGIT)
			{
				/*  Read next step increment(s) */
/* multiple, equal timesteps 15 aug. 2005 */
				if (replace("*", " ", token) == TRUE)
				{
					if (sscanf(token, "%d" SCANFORMAT, &k, &step) == 2)
					{
						for (i = 0; i < k; i++)
						{
							count_steps++;
							kinetics_ptr->steps =
								(LDBLE *) PHRQ_realloc(kinetics_ptr->steps,
													   (size_t) count_steps *
													   sizeof(LDBLE));
							if (kinetics_ptr->steps == NULL)
								malloc_error();
							kinetics_ptr->steps[kinetics_ptr->count_steps] =
								step;
							kinetics_ptr->count_steps = count_steps;
						}
					}
					else
					{
						input_error++;
						error_msg
							("Format error in multiple, equal KINETICS timesteps.\nCorrect is (for example): 20 4*10 2*5 3\n",
							 CONTINUE);
					}
				}
				else
				{
					step = strtod(token, &ptr);
					count_steps++;
					kinetics_ptr->steps =
						(LDBLE *) PHRQ_realloc(kinetics_ptr->steps,
											   (size_t) count_steps *
											   sizeof(LDBLE));
					if (kinetics_ptr->steps == NULL)
						malloc_error();
					kinetics_ptr->steps[kinetics_ptr->count_steps] = step;
					kinetics_ptr->count_steps = count_steps;
				}
			}
			if (j == EMPTY)
				break;
			/*
			 *   Read number of increments
			 */
			if (kinetics_ptr->count_steps != 1)
			{
				error_msg
					("To define equal time increments, only one total time should be defined.",
					 CONTINUE);
				input_error++;
				break;
			}
			do
			{
				j = sscanf(token, "%d", &i);
				if (j == 1 && i > 0)
				{
					kinetics_ptr->count_steps = -i;
					break;
				}
				else if (j == 1 && i <= 0)
				{
					error_msg
						("Expecting positive number for number of equal "
						 "time increments for kinetics.", CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
					break;
				}
			}
			while (copy_token(token, &next_char, &l) != EMPTY);
			break;
		case 6:				/* step_divide */
			if (copy_token(token, &next_char, &l) == DIGIT)
			{
				sscanf(token, SCANFORMAT, &kinetics_ptr->step_divide);
			}
			else
			{
				error_string = sformatf(
						"Expecting numerical value for step_divide.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			break;
		case 8:				/* runge-kutta */
		case 9:				/* runge_kutta */
		case 10:				/* rk */
			j = copy_token(token, &next_char, &l);
			if (j == DIGIT)
			{
				kinetics_ptr->rk = (int) strtod(token, &ptr);
			}
			else if (j == EMPTY)
			{
			}
			else
			{
				error_string = sformatf(
						"Expecting order for Runge-Kutta method.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			break;
		case 11:				/* bad_step_max */
			j = copy_token(token, &next_char, &l);
			if (j == DIGIT)
			{
				kinetics_ptr->bad_step_max = (int) strtod(token, &ptr);
			}
			else if (j == EMPTY)
			{
			}
			else
			{
				error_string = sformatf( "Expecting maximal bad steps number.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			break;
		case 12:				/* cvode */
			kinetics[n].use_cvode = get_true_false(next_char, TRUE);
			break;
		case 13:				/* cvode_steps */
			j = copy_token(token, &next_char, &l);
			if (j == DIGIT)
			{
				kinetics_ptr->cvode_steps = (int) strtod(token, &ptr);
			}
			else if (j == EMPTY)
			{
			}
			else
			{
				error_string = sformatf(
						"Expecting maximum number of steps for one call to cvode.");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			break;
		case 14:				/* cvode_order */
			j = copy_token(token, &next_char, &l);
			if (j == DIGIT)
			{
				kinetics_ptr->cvode_order = (int) strtod(token, &ptr);
			}
			else if (j == EMPTY)
			{
			}
			else
			{
				error_string = sformatf(
						"Expecting number of terms (order) used in cvode (1-5).");
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}

/*
 *   Default reactant
 */
	for (i = 0; i < kinetics[n].count_comps; i++)
	{
		if (kinetics[n].comps[i].count_list == 0)
		{
			kinetics[n].comps[i].list[0].name =
				kinetics_ptr->comps[i].rate_name;
			kinetics[n].comps[i].list[0].coef = 1;
			kinetics[n].comps[i].count_list = 1;
		}
	}
/*
 *   Default 1 sec
 */
	if (kinetics[n].count_steps == 0)
	{
		kinetics[n].count_steps = 1;
		kinetics[n].steps[0] = 1.0;
	}
/*
 *   set defaults for moles
 */
	for (i = 0; i < kinetics[n].count_comps; i++)
	{
		if (kinetics[n].comps[i].m0 < 0)
		{
			if (kinetics[n].comps[i].m < 0)
			{
				kinetics[n].comps[i].m0 = 1;
			}
			else
			{
				kinetics[n].comps[i].m0 = kinetics[n].comps[i].m;
			}
		}
		if (kinetics[n].comps[i].m < 0)
		{
			kinetics[n].comps[i].m = kinetics[n].comps[i].m0;
		}
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
LDBLE * Phreeqc::
read_list_doubles(char **ptr, int *count_doubles)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads a list of LDBLE numbers until end of line is reached or
 *   a LDBLE can not be read from a token.
 *
 *      Arguments:
 *	 ptr    entry: points to line to read from
 *		exit:  points to next non-LDBLE token or end of line
 *
 *	 count_doubles exit: number of LDBLEs read
 *
 *      Returns:
 *	 pointer to a list of count_doubles LDBLEs.
 */

	LDBLE *LDBLE_list;
	char token[MAX_LENGTH];
	LDBLE value;
	char *ptr_save;
	int l;

	LDBLE_list = (LDBLE *) PHRQ_malloc(sizeof(LDBLE));
	if (LDBLE_list == NULL)
		malloc_error();
	*count_doubles = 0;

	ptr_save = *ptr;
	while (copy_token(token, ptr, &l) != EMPTY)
	{
		if (sscanf(token, SCANFORMAT, &value) == 1)
		{
			*count_doubles = *count_doubles + 1;
			LDBLE_list =
				(LDBLE *) PHRQ_realloc(LDBLE_list,
									   (size_t) (*count_doubles) *
									   sizeof(LDBLE));
			if (LDBLE_list == NULL)
				malloc_error();
			LDBLE_list[(*count_doubles) - 1] = value;
			ptr_save = *ptr;
		}
		else
		{
			*ptr = ptr_save;
			break;
		}
	}
	return (LDBLE_list);
}

/* ---------------------------------------------------------------------- */
int * Phreeqc::
read_list_ints(char **ptr, int *count_ints, int positive)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads a list of int numbers until end of line is reached or
 *   an int can not be read from a token.
 *
 *      Arguments:
 *	 ptr    entry: points to line to read from
 *		exit:  points to next non-int token or end of line
 *
 *	 count_ints exit: number of LDBLEs read
 *
 *	 positive  entry: if TRUE, expects to read only positive integers
 *
 *      Returns:
 *	 pointer to a list of count_ints ints.
 */
	int *int_list;
	char token[MAX_LENGTH];
	int value;
	int l;
	char *ptr_save;

	int_list = (int *) PHRQ_malloc(sizeof(int));
	if (int_list == NULL)
		malloc_error();
	*count_ints = 0;

	ptr_save = *ptr;
	while (copy_token(token, ptr, &l) != EMPTY)
	{
		if (sscanf(token, "%d", &value) == 1)
		{
			(*count_ints)++;
			int_list =
				(int *) PHRQ_realloc(int_list,
									 (size_t) (*count_ints) * sizeof(int));
			if (int_list == NULL)
				malloc_error();
			int_list[(*count_ints) - 1] = value;
			if (value <= 0 && positive == TRUE)
			{
				error_msg("Expected an integer greater than zero.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
			}
			ptr_save = *ptr;
		}
		else
		{
			*ptr = ptr_save;
			break;
		}
	}
	return (int_list);
}

/* ---------------------------------------------------------------------- */
int * Phreeqc::
read_list_ints_range(char **ptr, int *count_ints, int positive, int *int_list)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads a list of int numbers until end of line is reached or
 *   an int can not be read from a token.
 *
 *      Arguments:
 *	 ptr    entry: points to line to read from
 *		exit:  points to next non-int token or end of line
 *
 *	 count_ints entry: number of ints already in list
 *
 *	 positive  entry: if TRUE, expects to read only positive integers
 *
 *      Returns:
 *	 pointer to a list of count_ints ints
 */
	char token[MAX_LENGTH];
	int value, value1, value2;
	int i, l;
	char *ptr_save;

	if (int_list == NULL)
	{
		int_list = (int *) PHRQ_malloc(sizeof(int));
		if (int_list == NULL)
			malloc_error();
		*count_ints = 0;
	}
	ptr_save = *ptr;
	while (copy_token(token, ptr, &l) != EMPTY)
	{
		if (sscanf(token, "%d", &value) == 1)
		{
			/* Read an integer */
			(*count_ints)++;
			int_list =
				(int *) PHRQ_realloc(int_list,
									 (size_t) (*count_ints) * sizeof(int));
			if (int_list == NULL)
				malloc_error();
			int_list[(*count_ints) - 1] = value;
			if (value <= 0 && positive == TRUE)
			{
				error_msg("Expected an integer greater than zero.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
			}
			/* Read range of integers */
			if (replace("-", " ", token) == TRUE)
			{
				if (sscanf(token, "%d %d", &value1, &value2) != 2)
				{
					error_msg("Expected an integer range n-m.", CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
				}
				else if (value2 < value1)
				{
					error_msg("Expected an integer range n-m, with n <= m.",
							  CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
				}
				else if (value2 <= 0 && positive == TRUE)
				{
					error_msg("Expected an integer greater than zero.",
							  CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
				}
				else
				{
					for (i = value1 + 1; i <= value2; i++)
					{
						(*count_ints)++;
						int_list =
							(int *) PHRQ_realloc(int_list,
												 (size_t) (*count_ints) *
												 sizeof(int));
						if (int_list == NULL)
							malloc_error();
						int_list[(*count_ints) - 1] = i;
					}
				}
			}
			ptr_save = *ptr;
		}
		else
		{
			*ptr = ptr_save;
			break;
		}
	}
	return (int_list);
}

/* ---------------------------------------------------------------------- */
int * Phreeqc::
read_list_t_f(char **ptr, int *count_ints)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads a list of true and false until end of line is reached or
 *   until non- t or f is found
 *
 *      Arguments:
 *	 ptr    entry: points to line to read from
 *		exit:  points to next non-int token or end of line
 *
 *	 count_ints exit: number of LDBLEs read
 *
 *	 positive  entry: if TRUE, expects to read only positive integers
 *
 *      Returns:
 *	 pointer to a list of count_ints ints.
 */
	int *int_list;
	char token[MAX_LENGTH];
	int value;
	int l;

	int_list = (int *) PHRQ_malloc(sizeof(int));
	if (int_list == NULL)
		malloc_error();
	*count_ints = 0;

	while (copy_token(token, ptr, &l) != EMPTY)
	{
		str_tolower(token);
		if (token[0] == 't')
		{
			value = TRUE;
		}
		else if (token[0] == 'f')
		{
			value = FALSE;
		}
		else
		{
			error_msg("Expected TRUE or FALSE.", CONTINUE);
			error_msg(line_save, CONTINUE);
			input_error++;
			break;
		}
		(*count_ints)++;
		int_list =
			(int *) PHRQ_realloc(int_list,
								 (size_t) (*count_ints) * sizeof(int));
		if (int_list == NULL)
			malloc_error();
		int_list[(*count_ints) - 1] = value;
	}
	return (int_list);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_log_k_only(char *ptr, LDBLE * log_k)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read log k
 */
	*log_k = 0.0;
	replace("=", " ", ptr);
	if (sscanf(ptr, SCANFORMAT, log_k) < 1)
	{
		input_error++;
		error_msg("Expecting log k.", CONTINUE);
		return (ERROR);
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_t_c_only(char *ptr, LDBLE *t_c)
/* ---------------------------------------------------------------------- */
{
	*t_c = 0.0;
	replace("=", " ", ptr);
	if (sscanf(ptr, SCANFORMAT, t_c) < 1)
	{
		input_error++;
		error_msg("Expecting numeric value for critical temperature T_c (K)", CONTINUE);
		return (ERROR);
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_p_c_only(char *ptr, LDBLE * p_c)
/* ---------------------------------------------------------------------- */
{
	*p_c = 0.0;
	replace("=", " ", ptr);
	if (sscanf(ptr, SCANFORMAT, p_c) < 1)
	{
		input_error++;
		error_msg("Expecting numeric value for critical pressure P_c (atm)", CONTINUE);
		return (ERROR);
	}
	return OK;
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_omega_only(char *ptr, LDBLE *omega)
/* ---------------------------------------------------------------------- */
{
	*omega = 0.0;
	replace("=", " ", ptr);
	if (sscanf(ptr, SCANFORMAT, omega) < 1)
	{
		input_error++;
		error_msg("Expecting numeric value for acentric factor Omega", CONTINUE);
		return (ERROR);
	}
	return OK;
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_delta_v_only(char *ptr, LDBLE * delta_v, DELTA_V_UNIT * units)
/* ---------------------------------------------------------------------- */
{
	int j, l;
	char token[MAX_LENGTH];
	/*
	*   Read analytical expression
	*/
	for (j = 0; j < 8; j++)
	{
		delta_v[j] = 0.0;
	}
	j = sscanf(ptr, SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT /*SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT*/,
		&(delta_v[0]), &(delta_v[1]), &(delta_v[2]), &(delta_v[3])/*,
		&(delta_v[4]), &(delta_v[5]), &(delta_v[6]), &(delta_v[7])*/);
	if (j < 1)
	{
		input_error++;
		error_msg("Expecting numeric values for the species' molar volume, vm.",
			CONTINUE);
		return (ERROR);
	}
	/*
	*   Read delta V units
	*/
	*units = cm3_per_mol;
	do
	{
		j = copy_token(token, &ptr, &l);
	} while (j == DIGIT); 

	if (j == EMPTY)
	{
		return (OK);
	}

	LDBLE factor = 1.0;
	if (j == UPPER || j == LOWER)
	{
		str_tolower(token);
		if (strstr(token, "cm3") != NULL)
		{
			/* cm3/mol */
			;
		}
		else if (strstr(token, "dm3") != NULL)
		{
			/* Convert dm3/mol to cm3/mol */
			//*delta_v *= 1E3;
			factor = 1e3;
		}
		else if (strstr(token, "m3") != NULL)
		{
			/* Convert m3/mol to cm3/mol */
			//*delta_v *= 1E6;
			factor = 1e6;
		}

		for (int i = 0; i < 8; i++)
		{
			delta_v[i] *= factor;
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_delta_h_only(char *ptr, LDBLE * delta_h, DELTA_H_UNIT * units)
/* ---------------------------------------------------------------------- */
{
	int j, l, kilo, joul;
	char token[MAX_LENGTH];
/*
 *   Read delta H
 */
	*delta_h = 0.0;
	replace("=", " ", ptr);
	j = copy_token(token, &ptr, &l);
	if (j == EMPTY)
	{
		input_error++;
		error_msg("Expecting numeric value for delta H.", CONTINUE);
		return (ERROR);
	}
	if (sscanf(token, SCANFORMAT, delta_h) < 1)
	{
		input_error++;
		error_msg("Expecting numeric value for delta H.", CONTINUE);
		return (ERROR);
	}
/*
 *   Read delta H units
 */
	j = copy_token(token, &ptr, &l);
	*units = kjoules;
	kilo = TRUE;
	joul = TRUE;
	if (j == EMPTY)
	{
		return (OK);
	}
	if (j == UPPER || j == LOWER)
	{
		str_tolower(token);
		if (strstr(token, "k") != token)
		{
			/* convert to kilo */
			kilo = FALSE;
			*delta_h /= 1000.;
		}
		if (strstr(token, "c") != NULL)
		{
			/* convert to joules */
			*delta_h *= JOULES_PER_CALORIE;
			joul = FALSE;
		}
	}
	if (kilo == TRUE && joul == TRUE)
	{
		*units = kjoules;
	}
	else if (kilo == FALSE && joul == TRUE)
	{
		*units = joules;
	}
	else if (kilo == TRUE && joul == FALSE)
	{
		*units = kcal;
	}
	else if (kilo == FALSE && joul == FALSE)
	{
		*units = cal;
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_analytical_expression_only(char *ptr, LDBLE * log_k)
/* ---------------------------------------------------------------------- */
{
	int j;
	int num_terms = T_A6 - T_A1 + 1;
/*
 *   Read analytical expression
 */
	for (j = 0; j < num_terms; j++)
	{
		log_k[j] = 0.0;
	}
	j = sscanf(ptr, SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT,
			   &(log_k[0]), &(log_k[1]), &(log_k[2]), &(log_k[3]),
			   &(log_k[4]), &(log_k[5]));
	if (j < 1)
	{
		input_error++;
		error_msg("Expecting numeric values for analytical expression.",
				  CONTINUE);
		return (ERROR);
	}
	return (OK);
}

/* VP: Density Start */
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_millero_abcdef (char *ptr, LDBLE * abcdef)
/* ---------------------------------------------------------------------- */
{
  int j;
/*
 *   Read a, b, c, d, e and f parameters for Millero density model.
 */
  for (j = 0; j < 6; j++)
  {
    abcdef[j] = 0.0;
  }
  j =
    sscanf (ptr, SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT SCANFORMAT,
	    &(abcdef[0]), &(abcdef[1]), &(abcdef[2]), &(abcdef[3]), &(abcdef[4]), &(abcdef[5]));
  if (j < 1)
  {
    input_error++;
    error_msg ("Expecting numeric values for analytical expression.",
	       CONTINUE);
    return (ERROR);
  }
  return (OK);
}
/* VP: Density End */

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_incremental_reactions(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Define flow only
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int j, l;
	char *ptr;
	char token[MAX_LENGTH];

	ptr = line;
	/* read keyword */
	copy_token(token, &ptr, &l);

	/* read true or false */
	incremental_reactions = get_true_false(ptr, TRUE);
/*
 *   find next keyword
 */
	while ((j =
			check_line("Subroutine Read", FALSE, TRUE, TRUE,
					   FALSE)) != KEYWORD)
	{
		/* empty, eof, keyword, print */
		if (j == EOF)
			return (EOF);
		error_string = sformatf( "Unknown input: %s", line);
		error_msg(error_string, CONTINUE);
		input_error++;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_master_species(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads master species data from data file or input file
 */
	int j, i, l;
	char *ptr, *ptr1;
	LDBLE l_z;
	struct element *elts_ptr;
	struct species *s_ptr;
	char token[MAX_LENGTH], token1[MAX_LENGTH];

	elts_ptr = NULL;
	for (;;)
	{
		j = check_line("Master species", FALSE, TRUE, TRUE, TRUE);
		if (j == EOF || j == KEYWORD)
		{
			break;
		}
/*
 *   Get element name with valence, allocate space, store
 */
		ptr = line;
/*
 *   Get element name and save pointer to character string
 */
		if (copy_token(token, &ptr, &l) != UPPER && token[0] != '[')
		{
			parse_error++;
			error_msg("Reading element for master species.", CONTINUE);
			error_msg(line_save, CONTINUE);
			continue;
		}
		/*
		   if (token[0] == '[') {
		   ptr1 = token;
		   get_elt(&ptr, element, &l);
		   strcpy(token, element);
		   }
		 */
		replace("(+", "(", token);
/*
 *   Delete master if it exists
 */
		master_delete(token);
/*
 *   Increase pointer array, if necessary,  and malloc space
 */
		if (count_master >= max_master)
		{
			space((void **) ((void *) &master), count_master + 1,
				  &max_master, sizeof(struct master *));
		}
		master[count_master] = master_alloc();
/*
 *   Set type to AQ
 */
		master[count_master]->type = AQ;
/*
 *   Save element name
 */
		master[count_master]->elt = element_store(token);
/*
 *   Save pointer to species data for master species
 */
		if ((copy_token(token, &ptr, &l) != UPPER) &&
			token[0] != '[' && (strcmp_nocase_arg1(token, "e-") != 0))
		{
			parse_error++;
			error_msg("Reading master species name.", CONTINUE);
			error_msg(line_save, CONTINUE);
			continue;
		}

		s_ptr = s_search(token);
		if (s_ptr != NULL)
		{
			master[count_master]->s = s_ptr;
		}
		else
		{
			ptr1 = token;
			get_token(&ptr1, token1, &l_z, &l);
			master[count_master]->s = s_store(token1, l_z, FALSE);
		}

/*
 *   Read alkalinity for species
 */
		copy_token(token, &ptr, &l);
		i = sscanf(token, SCANFORMAT, &master[count_master]->alk);
		if (i != 1)
		{
			input_error++;
			if (elts_ptr != NULL)
			{
				error_string = sformatf(
						"Expected alkalinity for master species, %s, in master species input.",
						elts_ptr->name);
			}
			else
			{
				error_string = sformatf(
						"Expected alkalinity for master species in master species input.");
			}
			error_msg(error_string, CONTINUE);
			continue;
		}
/*
 *   Read default gfw for species
 */
		i = copy_token(token, &ptr, &l);
		if (i == DIGIT)
		{
			sscanf(token, SCANFORMAT, &master[count_master]->gfw);
		}
		else if (i == UPPER)
		{
			master[count_master]->gfw_formula = string_hsave(token);
		}
		else
		{
			input_error++;
			if (elts_ptr != NULL)
			{
				error_string = sformatf(
						"Expected gram formula weight for master species, %s, in master species input.",
						elts_ptr->name);
			}
			else
			{
				error_string = sformatf(
						"Expected gram formula weight for master species in master species input.");
			}
			error_msg(error_string, CONTINUE);
			continue;
		}
/*
 *   MAKE LISTS OF PRIMARY AND SECONDARY MASTER SPECIES
 */
		if (strchr(master[count_master]->elt->name, '(') == NULL)
		{
			master[count_master]->primary = TRUE;
			/* Read gram formula weight for primary */
			if (strcmp(master[count_master]->elt->name, "E") != 0)
			{
				elts_ptr = master[count_master]->elt;
				i = copy_token(token, &ptr, &l);
				if (i == DIGIT)
				{
					sscanf(token, SCANFORMAT, &elts_ptr->gfw);
				}
				else
				{
					input_error++;
					if (elts_ptr != NULL)
					{
						error_string = sformatf(
								"Expected gram formula weight for element, %s.",
								elts_ptr->name);
					}
					else
					{
						error_string = sformatf(
								"Expected gram formula weight for element.");
					}

					error_msg(error_string, CONTINUE);
					continue;
				}
			}
		}
		else
		{
			master[count_master]->primary = FALSE;
		}
		count_master++;
		if (count_master >= max_master)
		{
			space((void **) ((void *) &master), count_master, &max_master,
				  sizeof(struct master *));
		}

	}
	return (j);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_mix(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads mixing fractions
 */
	int n_user, n_user_end;
	int return_value;
	int n_solution;
	LDBLE fraction;
	int j, i, l;
	char *ptr;
	char token[MAX_LENGTH];
	char *description;
	//struct mix *mix_ptr;
	cxxMix temp_mix;

/*
 *   Read mix number
 */
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);

	temp_mix.Set_n_user(n_user);
	temp_mix.Set_n_user_end(n_user);
	temp_mix.Set_description(description);
/*
 *   Set use data to first read
 */
	if (use.Get_mix_in() == FALSE)
	{
		use.Set_mix_in(true);
		use.Set_n_mix_user(n_user);
	}
/*
 *   Read mixture data
 */
	for (;;)
	{
		return_value = check_line("Mixture data", FALSE, TRUE, TRUE, TRUE);
		/* empty, eof, keyword, print */
		if (return_value == EOF || return_value == KEYWORD)
		{
			break;
		}
		ptr = line;
/*
 *   Read n_user
 */
		i = copy_token(token, &ptr, &l);
		if (i == DIGIT)
		{
			sscanf(token, "%d ", &n_solution);
		}
		else
		{
			input_error++;
			error_msg("Expected a solution number in mix input.", CONTINUE);
			error_msg(line_save, CONTINUE);
			continue;
		}
/*
 *   Read fraction for solution
 */
		copy_token(token, &ptr, &l);
		j = sscanf(token, SCANFORMAT, &fraction);
		if (j != 1)
		{
			input_error++;
			error_msg("Expected a mixing fraction.", CONTINUE);
			error_msg(line_save, CONTINUE);
			continue;
		}
		
/*
 *   Save data
 */
		temp_mix.Add(n_solution ,fraction);
	}
	if (temp_mix.Get_mixComps().size() == 0)
	{
		input_error++;
		error_msg
			("Must define at least one solution number and mixing fraction for MIX input.",
			 CONTINUE);
	}
	Rxn_mix_map[n_user] = temp_mix;

	// copy if needed
	if (n_user_end > n_user)
	{
		int i;
		for (i = n_user + 1; i <= n_user_end; i++)
		{
			Utilities::Rxn_copy(Rxn_mix_map, n_user, i);
		}
	}

	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_number_description(char *ptr, int *n_user,
						int *n_user_end, char **description, int allow_negative)
/* ---------------------------------------------------------------------- */
{
	int l, n;
	char token[MAX_LENGTH];
	char *ptr1;
/*
 *   Read user number, allow negative numbers Oct 3, 2011
 */
	copy_token(token, &ptr, &l);  // keyword
	ptr1 = ptr;
	copy_token(token, &ptr, &l);

	if (!isdigit(token[0]) && token[0] != '-')
	{
		*n_user = 1;
		*n_user_end = 1;
	}
	else
	{
		if (replace("-", " ", &token[1]))
		{
			n = sscanf(token, "%d%d", n_user, n_user_end);
			if (n != 2)
			{
				if (next_keyword >= 0)
				{
					error_string = sformatf( "Reading number range for %s.", Keywords::Keyword_name_search(next_keyword).c_str());
				}
				else
				{
					error_string = sformatf( "Reading number range for keyword.");
				}
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ptr1 = ptr;
		}
		else
		{
			n = sscanf(token, "%d", n_user);
			if (n != 1)
			{
				if (next_keyword >= 0)
				{
					error_string = sformatf( "Reading number range for %s.", Keywords::Keyword_name_search(next_keyword).c_str());
				}
				else
				{
					error_string = sformatf( "Reading number range for keyword.");
				}
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			*n_user_end = *n_user;
			ptr1 = ptr;
		};
	}
	if (*n_user < 0 && allow_negative == FALSE)
	{
		error_string = sformatf( "Negative number in number range not allowed for keyword.");
		error_msg(error_string, CONTINUE);
		input_error++;
	}
/*
 *   Read description
 */
	for (; isspace((int) ptr1[0]); ptr1++);
	*description = string_duplicate(ptr1);
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_phases(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read data for phases, parse equations
 */
	int j, i, l;
	int association;
	char *ptr;
	char token[MAX_LENGTH];
	char token1[MAX_LENGTH];
	struct phase *phase_ptr;
	struct elt_list *next_elt;
	struct rxn_token *token_ptr;

	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"no_check",				/* 0 */
		"check",				/* 1 */
		"log_k",				/* 2 */
		"logk",					/* 3 */
		"delta_h",				/* 4 */
		"deltah",				/* 5 */
		"analytical_expression",	/* 6 */
		"a_e",					/* 7 */
		"ae",					/* 8 */
		"add_logk",				/* 9 */
		"add_log_k",			/* 10 */
		"add_constant",			/* 11 */
		"t_c",					/* 12 */
		"p_c",					/* 13 */
		"omega",				/* 14 */
		"delta_v",				/* 15 */
		"deltav",				/* 16 */
		"vm"	/* 17, molar volume, must replace delta_v */
	};
	int count_opt_list = 18;

	association = FALSE;
/*
 *   Read eqn from file and call parser
 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	phase_ptr = NULL;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in PHASES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* no_check */
			if (phase_ptr == NULL)
				break;
			phase_ptr->check_equation = FALSE;
			break;
		case 1:				/* check */
			if (phase_ptr == NULL)
				break;
			phase_ptr->check_equation = TRUE;
			break;
		case 2:				/* log_k */
		case 3:				/* logk */
			if (phase_ptr == NULL)
				break;
			read_log_k_only(next_char, &phase_ptr->logk[0]);
			opt_save = OPTION_DEFAULT;
			break;
		case 4:				/* delta_h */
		case 5:				/* deltah */
			if (phase_ptr == NULL)
				break;
			read_delta_h_only(next_char, &phase_ptr->logk[1],
							  &phase_ptr->original_units);
			opt_save = OPTION_DEFAULT;
			break;
		case 6:				/* analytical_expression */
		case 7:				/* a_e */
		case 8:				/* ae */
			if (phase_ptr == NULL)
				break;
			read_analytical_expression_only(next_char, &(phase_ptr->logk[T_A1]));
			opt_save = OPTION_DEFAULT;
			break;
		case 9:				/* add_logk */
		case 10:				/* add_log_k */
			if (phase_ptr == NULL)
				break;
			if (phase_ptr->count_add_logk == 0)
			{
				phase_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (phase_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				phase_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(phase_ptr->add_logk,
													  (size_t) ((phase_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (phase_ptr->add_logk == NULL)
					malloc_error();
			}
			/* read name */
			if (copy_token(token, &next_char, &i) == EMPTY)
			{
				input_error++;
				error_string = sformatf(
						"Expected the name of a NAMED_EXPRESSION.");
				error_msg(error_string, CONTINUE);
				break;
			}
			phase_ptr->add_logk[phase_ptr->count_add_logk].name =
				string_hsave(token);
			/* read coef */
			i = sscanf(next_char, SCANFORMAT,
					   &phase_ptr->add_logk[phase_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				phase_ptr->add_logk[phase_ptr->count_add_logk].coef = 1;
			}
			phase_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 11:				/* add_constant */
			if (phase_ptr == NULL)
				break;
			if (phase_ptr->count_add_logk == 0)
			{
				phase_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (phase_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				phase_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(phase_ptr->add_logk,
													  (size_t) ((phase_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (phase_ptr->add_logk == NULL)
					malloc_error();
			}
			i = sscanf(next_char, SCANFORMAT,
					   &phase_ptr->add_logk[phase_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				input_error++;
				error_string = sformatf(
						"Expected the constant to add for log_K definition.");
				error_msg(error_string, CONTINUE);
				break;
			}
			/* set name */
			phase_ptr->add_logk[phase_ptr->count_add_logk].name =
				string_hsave("XconstantX");
			/* read coef */
			phase_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 12:				/* T_c */
			if (phase_ptr == NULL)
				break;
			read_t_c_only(next_char, &phase_ptr->t_c);
			opt_save = OPTION_DEFAULT;
			break;
		case 13:				/* P_c */
			if (phase_ptr == NULL)
				break;
			read_p_c_only(next_char, &phase_ptr->p_c);
			opt_save = OPTION_DEFAULT;
			break;
		case 14:				/* Omega */
			if (phase_ptr == NULL)
				break;
			read_omega_only(next_char, &phase_ptr->omega);
			opt_save = OPTION_DEFAULT;
			break;
		case 15:				/* delta_v */
		case 16:				/* delta_v */
		case 17:            /* vm, molar volume */
			if (phase_ptr == NULL)
				break;
			//read_delta_v_only(next_char, &phase_ptr->delta_v[1],
			read_delta_v_only(next_char, &(phase_ptr->logk[vm0]),
				&phase_ptr->original_deltav_units);
			phase_ptr->delta_v[1] = phase_ptr->logk[vm0];
			opt_save = OPTION_DEFAULT;
			break;
		case OPTION_DEFAULT:
/*
 *   Get element name and save pointer to character string
 */
			phase_ptr = NULL;
			ptr = line;
			copy_token(token, &ptr, &l);
/*
 *   Get and parse equation
 */
			j = check_line("Phase equation", FALSE, TRUE, TRUE, TRUE);
			if (j == EOF || j == KEYWORD)
			{
				return_value = j;
				break;
			}
			if (parse_eq(line, &next_elt, association) == ERROR)
			{
				parse_error++;
				error_msg("Parsing equation.", CONTINUE);
				error_msg(line_save, CONTINUE);
				break;
			}
			phase_ptr = phase_store(token);
/*
 *   Get pointer to each species in the reaction, store new species if necessary
 */
			strcpy(token1, trxn.token[0].name);
			replace("(g)", "", token1);
			replace("(s)", "", token1);
			replace("(G)", "", token1);
			replace("(S)", "", token1);
			phase_ptr->formula = string_hsave(token1);
			for (i = 1; i < count_trxn; i++)
			{
				if ((strstr(trxn.token[i].name, "(s)") == NULL) &&
					(strstr(trxn.token[i].name, "(g)") == NULL) &&
					(strstr(trxn.token[i].name, "(S)") == NULL) &&
					(strstr(trxn.token[i].name, "(G)") == NULL))
				{
					strcpy(token1, trxn.token[i].name);
					replace("(aq)", "", token1);
					replace("(AQ)", "", token1);
					replace("H2O(l)", "H2O", token1);
					replace("(H2O(L)", "H2O", token1);
					trxn.token[i].s = s_store(token1, trxn.token[i].z, FALSE);
				}
				else
				{
					trxn.token[i].s = NULL;
				}
			}
/*
 *   Save element list
 */
			phase_ptr->next_elt = next_elt;
/*
 *   Malloc space for phase reaction
 */
			phase_ptr->rxn = rxn_alloc(count_trxn + 1);
/*
 *   Copy reaction to reaction for phase, first token (token[0]) is not used
 *   except to check that coef of phase formula = 1.0
 */
			token_ptr = phase_ptr->rxn->token;
			/* token_ptr[0].coef=0; */
			token_ptr[0].coef = trxn.token[0].coef;
			token_ptr[0].s = trxn.token[1].s;
			for (i = 1; i < count_trxn; i++)
			{
				token_ptr[i].name = NULL;
				token_ptr[i].s = trxn.token[i].s;
				token_ptr[i].coef = trxn.token[i].coef;
				if (token_ptr[i].s == NULL)
				{
					token_ptr[i].name = trxn.token[i].name;
				}
			}
			token_ptr[0].name = trxn.token[1].name;
			/*
			   token_ptr[0].name=phase_ptr->name;
			   token_ptr[0].s=NULL;
			 */
			token_ptr[i].s = NULL;
			token_ptr[i].name = NULL;
/*
 *   Set type for phase
 */
			phase_ptr->type = SOLID;
			opt_save = OPTION_DEFAULT;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_pp_assemblage(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads pp_assemblage data
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	//int j, n, l, return_value;
	int j;
	int return_value;
	//int count_pure_phases;
	int n_user, n_user_end;
	char *ptr;
	char *description;
	//char token[MAX_LENGTH];
	std::string token;
	int opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"force_equality"		/* 0 */
	};
	int count_opt_list = 1;

	ptr = line;
	/*
	 *   Read pp_assemblage number
	 */
	read_number_description(ptr, &n_user, &n_user_end, &description);
	/*
	 *   Find pp_assemblage or realloc space for pp_assemblage
	 */
	cxxPPassemblage temp_pp_assemblage;
	cxxPPassemblageComp *comp = NULL;
	std::map<std::string, cxxPPassemblageComp> comps;
	temp_pp_assemblage.Set_new_def(true);
	temp_pp_assemblage.Set_n_user(n_user);
	temp_pp_assemblage.Set_n_user_end(n_user_end);
	temp_pp_assemblage.Set_description(description);
	free_check_null(description);
	/*
	 *   Set use data to first read
	 */
	if (use.Get_pp_assemblage_in() == FALSE)
	{
		use.Set_pp_assemblage_in(true);
		use.Set_n_pp_assemblage_user(n_user);
	}
	/*
	 *  Read equilibrium phase data
	 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in EQUILIBRIUM_PHASES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* force_equality */
			if (comp == NULL)
			{
				error_msg
					("Force_equality defined before equilibrium phase has been defined.",
					 CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
			}
			else
			{
				comp->Set_force_equality(get_true_false(next_char, TRUE) == TRUE);
			}
			break;
		case OPTION_DEFAULT:
			/*
			 *   Make space, set default
			 */
			if (comp)
			{
				comps[comp->Get_name()] = *comp;
			}
			delete comp;
			comp = new cxxPPassemblageComp;
			/*
			 *   Read name
			 */
			ptr = line;
			copy_token(token, &ptr);
			comp->Set_name(token.c_str());

			if ((j = copy_token(token, &ptr)) == EMPTY)
				continue;
			/*
			 *   Read saturation index
			 */
			j = sscanf(token.c_str(), SCANFORMAT, &dummy);
			comp->Set_si(dummy);
			comp->Set_si_org(dummy);
			if (j != 1)
			{
				error_msg("Expected saturation index.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				continue;
			}
			/*
			 *   Adding a reaction to the phase boundary
			 */
			if ((j = copy_token(token, &ptr)) == EMPTY)
				continue;
			if (j == UPPER || j == LOWER)
			{
				comp->Set_add_formula(token.c_str());
				j = copy_token(token, &ptr);
			}
			/*
			 *   Read amount
			 */
			if (j == EMPTY)
				continue;
			j = sscanf(token.c_str(), SCANFORMAT, &dummy);
			comp->Set_moles(dummy);
			if (j != 1)
			{
				error_msg("Expected amount of mineral.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				continue;
			}
			if ((j = copy_token(token, &ptr)) == EMPTY)
				continue;
			Utilities::str_tolower(token);
			if (strstr(token.c_str(), "d") == token.c_str())
			{
				comp->Set_dissolve_only(true);
				comp->Set_precipitate_only(false);
			} else if (strstr(token.c_str(), "p") == token.c_str())
			{
				comp->Set_precipitate_only(true);
				comp->Set_dissolve_only(false);
			}
			else
			{
				error_msg
					("Unexpected data at end of equilibrium-phase definition.",
					 CONTINUE);
				input_error++;
				continue;
			}
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	if (comp)
	{
		comps[comp->Get_name()] = *comp;
		delete comp;
		comp = NULL;
	}
	temp_pp_assemblage.Set_pp_assemblage_comps(comps);
	Rxn_pp_assemblage_map[n_user] = temp_pp_assemblage;

	return (return_value);
}
#ifdef SKIP
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_pure_phases(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads pure phase data
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int j, n, l, return_value;
	int count_pure_phases;
	int n_user, n_user_end;
	char *ptr;
	char *description;
	char token[MAX_LENGTH];
	int opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"force_equality"		/* 0 */
	};
	int count_opt_list = 1;

	ptr = line;
	/*
	 *   Read pp_assemblage number
	 */
	read_number_description(ptr, &n_user, &n_user_end, &description);
	/*
	 *   Find pp_assemblage or realloc space for pp_assemblage
	 */

	if (pp_assemblage_search(n_user, &n) != NULL)
	{
		pp_assemblage_free(&pp_assemblage[n]);
	}
	else
	{
		n = count_pp_assemblage++;
		space((void **) ((void *) &pp_assemblage), count_pp_assemblage,
			  &max_pp_assemblage, sizeof(struct pp_assemblage));
	}
	/*
	 *   Set use data to first read
	 */
	if (use.Get_pp_assemblage_in() == FALSE)
	{
		use.Get_pp_assemblage_in() = TRUE;
		use.Get_n_pp_assemblage_user() = n_user;
	}

	pp_assemblage_init(&(pp_assemblage[n]), n_user, n_user_end, description);
	free_check_null(description);
	count_pure_phases = 0;
	/*
	 *  Read equilibrium phase data
	 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in EQUILIBRIUM_PHASES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* force_equality */
			if (count_pure_phases == 0)
			{
				error_msg
					("Force_equality defined before equilibrium phase has been defined.",
					 CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
			}
			else
			{
				pp_assemblage[n].pure_phases[count_pure_phases -
											 1].force_equality =
					get_true_false(next_char, TRUE);
			}
			break;
		case OPTION_DEFAULT:
			/*
			 *   Make space, set default
			 */
			count_pure_phases++;
			pp_assemblage[n].pure_phases =
				(struct pure_phase *) PHRQ_realloc(pp_assemblage[n].
												   pure_phases,
												   (size_t)
												   (count_pure_phases) *
												   sizeof(struct pure_phase));
			if (pp_assemblage[n].pure_phases == NULL)
				malloc_error();
			pp_assemblage[n].pure_phases[count_pure_phases - 1].si = 0.0;
			pp_assemblage[n].pure_phases[count_pure_phases - 1].si_org = 0.0;
			pp_assemblage[n].pure_phases[count_pure_phases - 1].add_formula =
				NULL;
			pp_assemblage[n].pure_phases[count_pure_phases - 1].moles = 10.0;
			pp_assemblage[n].pure_phases[count_pure_phases - 1].delta = 0.0;
			pp_assemblage[n].pure_phases[count_pure_phases -
										 1].initial_moles = 0.0;
			pp_assemblage[n].pure_phases[count_pure_phases -
										 1].force_equality = FALSE;
			pp_assemblage[n].pure_phases[count_pure_phases -
										 1].dissolve_only = FALSE;
			pp_assemblage[n].pure_phases[count_pure_phases -
										 1].precipitate_only = FALSE;
			/*
			 *   Read name
			 */
			ptr = line;
			copy_token(token, &ptr, &l);
			pp_assemblage[n].pure_phases[count_pure_phases - 1].name =
				string_hsave(token);
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
				continue;
			/*
			 *   Read saturation index
			 */
			j = sscanf(token, SCANFORMAT, &dummy);
			pp_assemblage[n].pure_phases[count_pure_phases - 1].si =
				(LDBLE) dummy;
			pp_assemblage[n].pure_phases[count_pure_phases - 1].si_org = 
				(LDBLE) dummy;
			if (j != 1)
			{
				error_msg("Expected saturation index.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				continue;
			}
			/*
			 *   Adding a reaction to the phase boundary
			 */
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
				continue;
			if (j == UPPER || j == LOWER)
			{
				pp_assemblage[n].pure_phases[count_pure_phases -
											 1].add_formula =
					string_hsave(token);
				j = copy_token(token, &ptr, &l);
			}
			/*
			 *   Read amount
			 */
			if (j == EMPTY)
				continue;
			j = sscanf(token, SCANFORMAT, &dummy);
			pp_assemblage[n].pure_phases[count_pure_phases - 1].moles =
				(LDBLE) dummy;
			if (j != 1)
			{
				error_msg("Expected amount of mineral.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				continue;
			}
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
				continue;
			str_tolower(token);
			if (strstr(token, "d") == token)
			{
				pp_assemblage[n].pure_phases[count_pure_phases -
											 1].dissolve_only = TRUE;
				pp_assemblage[n].pure_phases[count_pure_phases -
											 1].precipitate_only = FALSE;
			} else if (strstr(token, "p") == token)
			{
				pp_assemblage[n].pure_phases[count_pure_phases -
											 1].precipitate_only = TRUE;
				pp_assemblage[n].pure_phases[count_pure_phases -
											 1].dissolve_only = FALSE;
			}
			else
			{
				error_msg
					("Unexpected data at end of equilibrium-phase definition.",
					 CONTINUE);
				input_error++;
				continue;
			}
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	pp_assemblage[n].count_comps = count_pure_phases;
	/*
	 *   Sort phases by name (lowercase)
	 */
	qsort(pp_assemblage[n].pure_phases,
		  (size_t) count_pure_phases,
		  (size_t) sizeof(struct pure_phase), pure_phase_compare);

	return (return_value);
}
#endif
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_reaction(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads reaction data
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
/*
 *   Read reaction
 */
	int l;
	char *ptr;
	char *description;
	char token[MAX_LENGTH];
	int return_value;
	int n_user, n_user_end;
	
/*
 *   Read reaction number
 */
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);

/*
 *   Set use data to first read
 */
	if (use.Get_reaction_in() == FALSE)
	{
		use.Set_reaction_in(true);
		use.Set_n_reaction_user(n_user);
	}
/*
 *   Defaults
 */
	cxxReaction temp_reaction;
	temp_reaction.Set_n_user(n_user);
	temp_reaction.Set_n_user_end(n_user_end);
	temp_reaction.Set_description(description);
	free_check_null(description);
/*
 *   Read reaction data
 */
	for (;;)
	{
/*
 *   Read line
 */
		return_value = check_line("Reaction data", FALSE, TRUE, TRUE, TRUE);
		/* empty, eof, keyword, print */
		if (return_value == EOF || return_value == KEYWORD)
		{
			break;
		}
		ptr = line;
		copy_token(token, &ptr, &l);
		if (isalpha((int) token[0]) || (token[0] == '(') || (token[0] == '['))
		{
/*
 *   Read reactant information
 */
			read_reaction_reactants(&temp_reaction);
		}
		else
		{
/*
 *   Read steps information
 */
			read_reaction_steps(&temp_reaction);
		}
	}
/*
 *   Default 1 mol of reaction
 */
	if (temp_reaction.Get_steps().size() == 0)
	{
		std::vector<LDBLE> v;
		v.push_back(1.0);
		temp_reaction.Set_steps(v);
	}
	if (temp_reaction.Get_equalIncrements())
	{
		if (temp_reaction.Get_countSteps() == 0)
		{
			temp_reaction.Set_countSteps(1);
		}
	}
	Rxn_reaction_map[n_user] = temp_reaction;
	// copy if needed
	Utilities::Rxn_copies(Rxn_reaction_map, n_user, n_user_end);

	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_reaction_reactants(cxxReaction *reaction_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read reactants, may be a chemical formula or a pure_phase name
 *   followed by relative reaction coefficient, default 1.0.
 *
 */
	std::string token, last_token;
	LDBLE coef;
	char *ptr;
/*
 *   Read one or more reactants
 */
	ptr = line;
	while (copy_token(token, &ptr) != EMPTY)
	{
/*
 *   Store reactant name, default coefficient
 */
		if (isalpha((int) token[0]) || (token[0] == '(') || (token[0] == '['))
		{

			reaction_ptr->Get_reactantList()[token] = 1.0;
			last_token = token;
/*
 *   Store relative coefficient
 */
		}
		else
		{
			int j = sscanf(token.c_str(), SCANFORMAT, &coef);
			if (j == 1 && last_token.size() > 0)
			{
				reaction_ptr->Get_reactantList()[last_token] = coef;
			}
			else
			{
				error_msg("Reading relative coefficient of reactant.",
						  CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
			}
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_reaction_steps(cxxReaction *reaction_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read amount(s) of irrev reactions in one of three forms:
 *
 *   6 millimoles in 6 steps   or
 *
 *   1 2 3 4 5 6 millimoles    or
 *
 *   6*1 millimoles
 *   INCREMENTAL_REACTIONS
 */

	char *ptr;
	std::string token, token1;

	ptr = line;
/*
 *   Read one or more reaction increments
 */
	for (;;)
	{
		if (copy_token(token, &ptr) == EMPTY)
		{
			return (OK);
		}
/*
 *   Read next step increment
 */
/* begin modif 29 july 2005... */
		if (replace("*", " ", token))
		{
			int n;
			LDBLE value;
			if (sscanf(token.c_str(), "%d" SCANFORMAT, &n, &value) == 2)
			{
				for (int i = 0; i < n; i++)
				{
					reaction_ptr->Get_steps().push_back(value);
				}
			}
			else
			{
				input_error++;
				error_msg
					("Format error in multiple, equal REACTION steps.\nCorrect is (for example): 0.2 4*0.1 2*0.5 0.3\n",
					 CONTINUE);
			}
		}
		else
		{
			LDBLE step;
			int j = sscanf(token.c_str(), SCANFORMAT, &step);
			if (j == 1)
			{
				reaction_ptr->Get_steps().push_back(step);
			}
			else
			{
				break;
			}
		}
/* ...end modif 29 july 2005 */
	}
/*
 *   Read units
 */
	token1 = token;
	token1.append("/l");
	char * t1 = string_duplicate(token1.c_str());
	if (check_units(t1, FALSE, FALSE, NULL, FALSE) == OK)
	{
		replace("/l", "", t1);
		if (strstr(t1, "Mol") == NULL)
		{
			error_string = sformatf( "Units of steps not in moles, %s.", token.c_str());
			error_msg(error_string, CONTINUE);
			input_error++;
			return (ERROR);
		}
		else
		{
			reaction_ptr->Set_units(t1);
		}
		if (copy_token(token, &ptr) == EMPTY)
		{
			return (OK);
		}
	}
/*
 *  Read number of equal increments, store as negative integer
 */
	if (reaction_ptr->Get_actualSteps() != 1)
	{
		error_msg
			("To define equal increments, only one reaction increment should be defined.",
			 CONTINUE);
		input_error++;
		return (ERROR);
	}
	do
	{
		int i;
		int j = sscanf(token.c_str(), "%d", &i);
		if (j == 1 && i > 0)
		{
			reaction_ptr->Set_countSteps(i);
			reaction_ptr->Set_equalIncrements(true);
			return (OK);
		}
		else if (j == 1 && i <= 0)
		{
			break;
		}
	}
	while (copy_token(token, &ptr) != EMPTY);

	error_msg("Expecting positive number for number of equal "
			  "increments to add.", CONTINUE);
	error_msg(line_save, CONTINUE);
	input_error++;
	return (ERROR);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_save(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution, mix, irreversible reaction, and pure phases to use
 *   in reaction calculation
 */
	int i, l, n, n_user, n_user_end;
	char *ptr;
	char token[MAX_LENGTH];
/*
 *   Read "save"
 */
	ptr = line;
	copy_token(token, &ptr, &l);
/*
 *   Read keyword
 */
	copy_token(token, &ptr, &l);
	check_key(token);
/*
 *   Read number
 */
	for (;;)
	{
		i = copy_token(token, &ptr, &l);
		if (i == DIGIT)
		{
			replace("-", " ", token);
			n = sscanf(token, "%d%d", &n_user, &n_user_end);
			if (n == 1)
			{
				n_user_end = n_user;
			}
			if (n_user < 0)
			{
				error_msg("Number must be a positive integer.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
			}
			break;
		}
		else if (i == EMPTY)
		{
			error_string = sformatf( "No number given, 1 assumed.");
			warning_msg(error_string);
			n_user = 1;
			n_user_end = 1;
			break;
		}
	}

	switch (next_keyword)
	{
	case Keywords::KEY_SOLUTION:								/* Solution */
		save.solution = TRUE;
		save.n_solution_user = n_user;
		save.n_solution_user_end = n_user_end;
		break;
	case Keywords::KEY_EQUILIBRIUM_PHASES:					/* Pure phases */
		save.pp_assemblage = TRUE;
		save.n_pp_assemblage_user = n_user;
		save.n_pp_assemblage_user_end = n_user_end;
		break;
	case Keywords::KEY_EXCHANGE:								/* exchange */
		save.exchange = TRUE;
		save.n_exchange_user = n_user;
		save.n_exchange_user_end = n_user_end;
		break;
	case Keywords::KEY_SURFACE:								/* surface */
		save.surface = TRUE;
		save.n_surface_user = n_user;
		save.n_surface_user_end = n_user_end;
		break;
	case Keywords::KEY_GAS_PHASE:								/* gas_phase */
		save.gas_phase = TRUE;
		save.n_gas_phase_user = n_user;
		save.n_gas_phase_user_end = n_user_end;
		break;
	case Keywords::KEY_SOLID_SOLUTIONS:						/* solid_solutions */
		save.ss_assemblage = TRUE;
		save.n_ss_assemblage_user = n_user;
		save.n_ss_assemblage_user_end = n_user_end;
		break;
	default:
		input_error++;
		error_msg
			("Expecting keyword solution, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions.",
			CONTINUE);
		error_msg(line_save, CONTINUE);
		check_line("End of save", FALSE, TRUE, TRUE, TRUE);
		/* empty, eof, keyword, print */
		return (ERROR);
	}

	check_line("End of save", FALSE, TRUE, TRUE, TRUE);
	/* empty, eof, keyword, print */
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_selected_output(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read data for to output to flat file
 */
	int value;
	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"file",					/* 0 */
		"totals",				/* 1 */
		"molalities",			/* 2 */
		"activities",			/* 3 */
		"pure_phases",			/* 4 */
		"si",					/* 5 */
		"saturation_indices",	/* 6 */
		"gases",				/* 7 */
		"equilibrium_phases",	/* 8 */
		"equilibria",       	/* 9 */
		"equilibrium",			/* 10 */
		"pure",					/* 11 */
		"inverse",				/* 12 */
		"kinetic_reactants",	/* 13 */
		"kinetics",				/* 14 */
		"solid_solutions",		/* 15 */
		"inverse_modeling",		/* 16 */
		"reset",				/* 17 */
		"simulation",			/* 18 */
		"sim",					/* 19 */
		"state",				/* 20 */
		"solution",				/* 21 */
		"soln",					/* 22 */
		"distance",				/* 23 */
		"dist",					/* 24 */
		"time",					/* 25 */
		"step",					/* 26 */
		"reaction",				/* 27 */
		"rxn",					/* 28 */
		"temperature",			/* 29 */
		"temp",					/* 30 */
		"ph",					/* 31 */
		"pe",					/* 32 */
		"alkalinity",			/* 33 */
		"alk",					/* 34 */
		"ionic_strength",		/* 35 */
		"mu",					/* 36 */
		"water",				/* 37 */
		"high_precision",		/* 38 */
		"user_punch",			/* 39 */
		"mol",					/* 40 */
		"kin",					/* 41 */
		"charge_balance",		/* 42 */
		"percent_error",		/* 43 */
		"selected_out",			/* 44 */
		"selected_output",		/* 45 */
		"isotopes",				/* 46 */
		"calculate_values",		/* 47 */
		"equilibrium_phase"     /* 48 */
	};
	int count_opt_list = 49;

	int i, l;
	char file_name[MAX_LENGTH], token[MAX_LENGTH];

	punch.in = TRUE;
	punch.new_def = TRUE;
	punch.count_totals = 0;
	punch.count_molalities = 0;
	punch.count_activities = 0;
	punch.count_pure_phases = 0;
	punch.count_si = 0;
	punch.count_gases = 0;
	punch.count_kinetics = 0;
	punch.count_s_s = 0;
/*
 *   Read eqn from file and call parser
 */
	opt_save = OPTION_ERROR;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		opt_save = opt;
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SELECTED_OUTPUT keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* file name */
			/* copy_token(file_name, &next_char, &l); */
			if (string_trim(next_char) != EMPTY)
			{
				strcpy(file_name, next_char);
				have_punch_name = TRUE;
				punch_close();
				if (punch_open(file_name) != OK)
				{
					error_string = sformatf( "Can't open file, %s.", file_name);
					input_error++;
					error_msg(error_string, CONTINUE);
				}
			}
			opt_save = OPTION_ERROR;
			break;
		case 1:				/* totals */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				if (i != UPPER && token[0] != '[')
				{
					error_string = sformatf( "Expected element name to"
							" begin with upper case letter.");
					warning_msg(error_string);
				}
				else
				{
					punch.count_totals++;
					punch.totals =
						(struct name_master *) PHRQ_realloc(punch.totals,
															(size_t) punch.
															count_totals *
															sizeof(struct
																   name_master));
					if (punch.totals == NULL)
						malloc_error();
					punch.totals[punch.count_totals - 1].name =
						string_hsave(token);
				}
			}
			break;
		case 2:				/* molalities */
		case 40:				/* mol */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				if (i != UPPER && token[0] != '(' && (token[0] != '['))
				{
					error_string = sformatf( "Expected species name to"
							" begin with upper case letter.");
					warning_msg(error_string);
				}
				else
				{
					punch.count_molalities++;
					punch.molalities =
						(struct name_species *) PHRQ_realloc(punch.
															 molalities,
															 (size_t) punch.
															 count_molalities
															 *
															 sizeof(struct
																	name_species));
					if (punch.molalities == NULL)
						malloc_error();
					punch.molalities[punch.count_molalities - 1].name =
						string_hsave(token);
				}
			}
			break;
		case 3:				/* activities */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				if (i != UPPER && token[0] != '(' && (token[0] != '['))
				{
					error_string = sformatf( "Expected species name to"
							" begin with upper case letter.");
					warning_msg(error_string);
				}
				else
				{
					punch.count_activities++;
					punch.activities =
						(struct name_species *) PHRQ_realloc(punch.
															 activities,
															 (size_t) punch.
															 count_activities
															 *
															 sizeof(struct
																	name_species));
					if (punch.activities == NULL)
						malloc_error();
					punch.activities[punch.count_activities - 1].name =
						string_hsave(token);
				}
			}
			break;
		case 4:				/* pure_phases */
		case 8:				/* equilibrium_phases */
		case 9:				/* equilibria */
		case 10:			/* equilibrium */
		case 11:			/* pure */
		case 48:            /* equilibrium_phase */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				punch.count_pure_phases++;
				punch.pure_phases =
					(struct name_phase *) PHRQ_realloc(punch.pure_phases,
													   (size_t) punch.
													   count_pure_phases *
													   sizeof(struct
															  name_phase));
				if (punch.pure_phases == NULL)
					malloc_error();
				punch.pure_phases[punch.count_pure_phases - 1].name =
					string_hsave(token);
			}
			break;
		case 5:				/* si */
		case 6:				/* saturation_index */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				punch.count_si++;
				punch.si =
					(struct name_phase *) PHRQ_realloc(punch.si,
													   (size_t) punch.
													   count_si *
													   sizeof(struct
															  name_phase));
				if (punch.si == NULL)
					malloc_error();
				punch.si[punch.count_si - 1].name = string_hsave(token);
			}
			break;
		case 7:				/* gases */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				punch.count_gases++;
				punch.gases =
					(struct name_phase *) PHRQ_realloc(punch.gases,
													   (size_t) punch.
													   count_gases *
													   sizeof(struct
															  name_phase));
				if (punch.gases == NULL)
					malloc_error();
				punch.gases[punch.count_gases - 1].name = string_hsave(token);
			}
			break;
		case 12:				/* inverse */
		case 16:				/* inverse_modeling */
			punch.inverse = get_true_false(next_char, TRUE);
			break;
		case 13:				/* kinetic_reactants */
		case 14:				/* kinetics */
		case 41:				/* kin */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				punch.count_kinetics++;
				punch.kinetics =
					(struct name_phase *) PHRQ_realloc(punch.kinetics,
													   (size_t) punch.
													   count_kinetics *
													   sizeof(struct
															  name_phase));
				if (punch.kinetics == NULL)
					malloc_error();
				punch.kinetics[punch.count_kinetics - 1].name =
					string_hsave(token);
			}
			break;
		case 15:				/* solid_solutions */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				punch.count_s_s++;
				punch.s_s =
					(struct name_phase *) PHRQ_realloc(punch.s_s,
													   (size_t) punch.
													   count_s_s *
													   sizeof(struct
															  name_phase));
				if (punch.s_s == NULL)
					malloc_error();
				punch.s_s[punch.count_s_s - 1].name = string_hsave(token);
			}
			break;
		case 46:				/* isotopes */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				if (i != UPPER && token[0] != '[')
				{
					error_string = sformatf( "Expected element name to"
							" begin with upper case letter.");
					warning_msg(error_string);
				}
				else
				{
					punch.count_isotopes++;
					punch.isotopes =
						(struct name_master *) PHRQ_realloc(punch.isotopes,
															(size_t) punch.
															count_isotopes *
															sizeof(struct
																   name_master));
					if (punch.isotopes == NULL)
						malloc_error();
					punch.isotopes[punch.count_isotopes - 1].name =
						string_hsave(token);
				}
			}
			break;
		case 47:				/* calculate_values */
			while ((i = copy_token(token, &next_char, &l)) != EMPTY)
			{
				punch.count_calculate_values++;
				punch.calculate_values =
					(struct name_master *) PHRQ_realloc(punch.
														calculate_values,
														(size_t) punch.
														count_calculate_values
														*
														sizeof(struct
															   name_master));
				if (punch.calculate_values == NULL)
					malloc_error();
				punch.calculate_values[punch.count_calculate_values -
									   1].name = string_hsave(token);
			}
			break;
		case 17:				/* reset */
			value = get_true_false(next_char, TRUE);
			/* matches print order */
			punch.sim = value;
			punch.state = value;
			punch.soln = value;
			punch.dist = value;
			punch.time = value;
			punch.step = value;
			punch.ph = value;
			punch.pe = value;
			punch.rxn = value;
			punch.temp = value;
			punch.alk = value;
			punch.mu = value;
			punch.water = value;
			punch.charge_balance = value;
			punch.percent_error = value;
			break;
		case 18:				/* simulation */
		case 19:				/* sim */
			punch.sim = get_true_false(next_char, TRUE);
			break;
		case 20:				/* state */
			punch.state = get_true_false(next_char, TRUE);
			break;
		case 21:				/* solution */
		case 22:				/* soln */
			punch.soln = get_true_false(next_char, TRUE);
			break;
		case 23:				/* distance */
		case 24:				/* dist */
			punch.dist = get_true_false(next_char, TRUE);
			break;
		case 25:				/* time */
			punch.time = get_true_false(next_char, TRUE);
			break;
		case 26:				/* step */
			punch.step = get_true_false(next_char, TRUE);
			break;
		case 27:				/* reaction */
		case 28:				/* rxn */
			punch.rxn = get_true_false(next_char, TRUE);
			break;
		case 29:				/* temperature */
		case 30:				/* temp */
			punch.temp = get_true_false(next_char, TRUE);
			break;
		case 31:				/* ph */
			punch.ph = get_true_false(next_char, TRUE);
			break;
		case 32:				/* pe */
			punch.pe = get_true_false(next_char, TRUE);
			break;
		case 33:				/* alkalinity */
		case 34:				/* alk */
			punch.alk = get_true_false(next_char, TRUE);
			break;
		case 35:				/* ionic strength */
		case 36:				/* mu */
			punch.mu = get_true_false(next_char, TRUE);
			break;
		case 37:				/* water */
			punch.water = get_true_false(next_char, TRUE);
			break;
		case 38:				/* high_precision */
			punch.high_precision = get_true_false(next_char, TRUE);
			if (punch.high_precision == TRUE /*&& convergence_tolerance > 1e-12*/)
			{
				convergence_tolerance = 1e-12;
			}
			break;
		case 39:				/* user_punch */
			punch.user_punch = get_true_false(next_char, TRUE);
			break;
		case 42:				/* charge_balance */
			punch.charge_balance = get_true_false(next_char, TRUE);
			break;
		case 43:				/* percent_error */
			punch.percent_error = get_true_false(next_char, TRUE);
			break;
		case 44:				/* selected_out */
		case 45:				/* selected_output */
			pr.punch = get_true_false(next_char, TRUE);
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	if (!have_punch_name)
	{
		punch_close();
		if (punch_open("selected.out") != OK)
		{
			error_string = sformatf( "Can't open file, %s.", "selected.out");
			input_error++;
			error_msg(error_string, CONTINUE);
		}
	}

	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_solution(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads solution data
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int i, j, n, l;
	int n_user, n_user_end;
	int default_pe, alk;
	int count_isotopes;
	int max_mass_balance, count_mass_balance;
	char *ptr, *ptr1;
	char *description;
	char token[MAX_LENGTH], token1[MAX_LENGTH];

	int return_value, opt;
	char *next_char;
	const char *opt_list[] = {
		"temp",					/* 0 */
		"temperature",			/* 1 */
		"dens",					/* 2 */
		"density",				/* 3 */
		"units",				/* 4 */
		"redox",				/* 5 */
		"ph",					/* 6 */
		"pe",					/* 7 */
		"unit",					/* 8 */
		"isotope",				/* 9 */
		"water",				/* 10 */
		"press",				/* 11 */
		"pressure"				/* 12 */
	};
	int count_opt_list = 13;
/*
 *   Read solution number and description
 */
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);
/*
 *   Malloc space for solution data
 */
	if (solution_bsearch(n_user, &n, FALSE) != NULL)
	{
		solution_free(solution[n]);
	}
	else
	{
		n = count_solution++;
		if (count_solution >= max_solution)
		{
			space((void **) ((void *) &(solution)), count_solution,
				  &max_solution, sizeof(struct solution *));
		}
	}
	solution[n] = solution_alloc();

	solution[n]->n_user = n_user;
	solution[n]->n_user_end = n_user_end;
	if (use.Get_solution_in() == FALSE)
	{
		use.Set_solution_in(true);
		use.Set_n_solution_user(n_user);
	}
	max_mass_balance = MAX_MASS_BALANCE;
/*
 *   Set default ph, temp, density, pe, units
 */
	solution[n]->description = description;
	solution[n]->tc = 25.0;
	solution[n]->patm = 1;
	solution[n]->ph = 7.0;
	solution[n]->density = 1.0;
	solution[n]->solution_pe = 4.0;
	solution[n]->mass_water = 1.0;
	solution[n]->ah2o = 1.0;
	solution[n]->mu = 1e-7;
	solution[n]->cb = 0.0;
	solution[n]->count_master_activity = 0;
	default_pe = 0;
	solution[n]->units = string_hsave("mMol/kgw");
	solution[n]->totals[0].description = NULL;
	count_mass_balance = 0;
	count_isotopes = 0;
/*
 *   Read concentration data
 */
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			ptr = next_char;
			if (copy_token(token, &ptr, &l) == DIGIT)
			{
				opt = 9;
			}
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SOLUTION keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* temperature */
		case 1:
			if (sscanf(next_char, SCANFORMAT, &(solution[n]->tc)) != 1)
			{
				solution[n]->tc = 25;
			}
			break;
		case 2:				/* density */
		case 3:
			sscanf(next_char, SCANFORMAT, &(solution[n]->density));
			break;
		case 4:				/* units */
		case 8:				/* unit */
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			{
				if (check_units(token, FALSE, FALSE, solution[n]->units, TRUE) ==
					OK)
				{
					solution[n]->units = string_hsave(token);
				}
				else
				{
					input_error++;
				}
			}
			break;
		case 5:				/* redox */
			if (copy_token(token, &next_char, &l) == EMPTY)
				break;
			if (parse_couple(token) == OK)
			{
				default_pe = pe_data_store(&(solution[n]->pe), token);
			}
			else
			{
				input_error++;
			}
			break;
		case 6:				/* ph */
			if (read_conc(n, count_mass_balance, line) == ERROR)
			{
				input_error++;
				break;
			}
			solution[n]->ph =
				solution[n]->totals[count_mass_balance].input_conc;
			if (solution[n]->totals[count_mass_balance].equation_name == NULL)
			{
				break;
			}
			solution[n]->totals[count_mass_balance].description =
				string_hsave("H(1)");
			count_mass_balance++;
			break;
		case 7:				/* pe */
			if (read_conc(n, count_mass_balance, line) == ERROR)
			{
				input_error++;
				break;
			}
			solution[n]->solution_pe =
				solution[n]->totals[count_mass_balance].input_conc;
			if (solution[n]->totals[count_mass_balance].equation_name == NULL)
			{
				break;
			}
			solution[n]->totals[count_mass_balance].description =
				string_hsave("E");
			count_mass_balance++;
			break;
		case 9:				/* isotope */
			if (copy_token(token, &next_char, &l) != DIGIT)
			{
				input_error++;
				error_string = sformatf( "Expected isotope name to"
						" begin with an isotopic number.");
				error_msg(error_string, CONTINUE);
				continue;
			}
			solution[n]->isotopes =
				(struct isotope *) PHRQ_realloc(solution[n]->isotopes,
												(size_t) (count_isotopes +
														  1) *
												sizeof(struct isotope));
			if (solution[n]->isotopes == NULL)
				malloc_error();

			/* read and save element name */
			ptr1 = token;
			get_num(&ptr1,
					&(solution[n]->isotopes[count_isotopes].isotope_number));
			if (ptr1[0] == '\0' || isupper((int) ptr1[0]) == FALSE)
			{
				error_msg("Expecting element name.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				return (ERROR);
			}
			solution[n]->isotopes[count_isotopes].elt_name =
				string_hsave(ptr1);

			/* read and store isotope ratio */
			if (copy_token(token, &next_char, &l) != DIGIT)
			{
				input_error++;
				error_string = sformatf(
						"Expected numeric value for isotope ratio.");
				error_msg(error_string, CONTINUE);
				continue;
			}
			sscanf(token, SCANFORMAT,
				   &(solution[n]->isotopes[count_isotopes].ratio));

			solution[n]->isotopes[count_isotopes].ratio_uncertainty = NAN;

			/* read and store isotope ratio uncertainty */
			if ((j = copy_token(token, &next_char, &l)) != EMPTY)
			{
				if (j != DIGIT)
				{
					input_error++;
					error_string = sformatf(
							"Expected numeric value for uncertainty in isotope ratio.");
					error_msg(error_string, CONTINUE);
					continue;
				}
				sscanf(token, SCANFORMAT,
					   &(solution[n]->isotopes[count_isotopes].
						 ratio_uncertainty));
			}
			count_isotopes++;
			break;
		case 10:				/* water */
			j = copy_token(token, &next_char, &l);
			if (j == EMPTY)
			{
				solution[n]->mass_water = 1.0;
			}
			else if (j != DIGIT)
			{
				input_error++;
				error_string = sformatf(
						"Expected numeric value for mass of water in solution.");
				error_msg(error_string, CONTINUE);
			}
			else
			{
				sscanf(token, SCANFORMAT, &dummy);
				solution[n]->mass_water = (LDBLE) dummy;
			}
			break;
		case 11: /* pressure */
		case 12:
			if (sscanf(next_char, SCANFORMAT, &(solution[n]->patm)) != 1)
			{
				solution[n]->patm = 1;
			}
			break;
		case OPTION_DEFAULT:
/*
 *   Read concentration
 */
			if (read_conc(n, count_mass_balance, line) == ERROR)
			{
				input_error++;
				break;
			}
			count_mass_balance++;
			break;
		}
		if (count_mass_balance + 1 >= max_mass_balance)
		{
			space((void **) ((void *) &(solution[n]->totals)),
				  count_mass_balance + 1, &max_mass_balance,
				  sizeof(struct conc));
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
/*
 *   fix up default units and default pe
 */
	for (i = 0; i < count_mass_balance; i++)
	{
		strcpy(token, solution[n]->totals[i].description);
		str_tolower(token);
		if (solution[n]->totals[i].units == NULL)
		{
			solution[n]->totals[i].units = solution[n]->units;
		}
		else
		{
			alk = FALSE;
			if (strstr(token, "alk") == token)
				alk = TRUE;
			strcpy(token1, solution[n]->totals[i].units);
			if (check_units(token1, alk, TRUE, solution[n]->units, TRUE) ==
				ERROR)
			{
				input_error++;
			}
			else
			{
				solution[n]->totals[i].units = string_hsave(token1);
			}
		}
		if (solution[n]->totals[i].n_pe < 0)
		{
			solution[n]->totals[i].n_pe = default_pe;
		}
	}
	solution[n]->default_pe = default_pe;
/*
 *   Mark end of solution
 */
	solution[n]->totals[count_mass_balance].description = NULL;
	solution[n]->count_isotopes = count_isotopes;
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_species(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read data for aqueous species, parse equations
 */
	int i;
	int association;
	struct species *s_ptr;
	struct elt_list *next_elt;
	char *ptr, token[MAX_LENGTH];
	bool vm_read;
	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"no_check",				/* 0 */
		"check",				/* 1 */
		"gamma",				/* 2 */
		"mb",					/* 3 */
		"mass_balance",			/* 4 */
		"log_k",				/* 5 */
		"logk",					/* 6 */
		"delta_h",				/* 7 */
		"deltah",				/* 8 */
		"analytical_expression",	/* 9 */
		"a_e",					/* 10 */
		"ae",					/* 11 */
		"mole_balance",			/* 12 */
		"llnl_gamma",			/* 13 */
		"co2_llnl_gamma",		/* 14 */
		"activity_water",		/* 15 */
		"add_logk",				/* 16 */
		"add_log_k",			/* 17 */
		"add_constant",			/* 18 */
		"dw",					/* 19 */
		"erm_ddl",				/* 20 */
/* VP: Density Start */
		"millero",				/* 21 */
/* VP: Density End */
		"delta_v",				/* 22 */
		"deltav",				/* 23 */
		"vm"	/* 24, molar volume, must replace delta_v */
	};
	int count_opt_list = 25;

	association = TRUE;
	s_ptr = NULL;
/*
 *   Read eqn from file and call parser
 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		vm_read = false;
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SOLUTION_SPECIES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* no_check */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->check_equation = FALSE;
			break;
		case 1:				/* check */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->check_equation = TRUE;
			break;
		case 2:				/* gamma data */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->gflag = 2;	/* Wateq D-H */
			i = sscanf(next_char, SCANFORMAT SCANFORMAT, &s_ptr->dha,
					   &s_ptr->dhb);
			if (i < 2)
			{
				error_string = sformatf( "Expecting 2 activity-"
						"coefficient parameters, a and b.");
				warning_msg(error_string);
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 3:				/* mb */
		case 4:				/* mass_balance */
		case 12:				/* mole_balance */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			count_elts = 0;
			paren_count = 0;
			copy_token(token, &next_char, &i);
			s_ptr->mole_balance = string_hsave(token);
			ptr = token;
			s_ptr->next_secondary =
				(struct elt_list *) free_check_null(s_ptr->next_secondary);
			get_secondary_in_species(&ptr, 1.0);
			s_ptr->next_secondary = elt_list_save();
/* debug
			for (i = 0; i < count_elts; i++) {
				output_msg(sformatf("%s\t%f\n", elt_list[i].elt->name,
					elt_list[i].coef));
			}
 */
			opt_save = OPTION_DEFAULT;
			break;
		case 5:				/* log_k */
		case 6:				/* logk */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_log_k_only(next_char, &s_ptr->logk[0]);
			opt_save = OPTION_DEFAULT;
			break;
		case 7:				/* delta_h */
		case 8:				/* deltah */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_delta_h_only(next_char, &s_ptr->logk[1],
							  &s_ptr->original_units);
			opt_save = OPTION_DEFAULT;
			break;
		case 9:				/* analytical_expression */
		case 10:				/* a_e */
		case 11:				/* ae */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_analytical_expression_only(next_char, &(s_ptr->logk[T_A1]));
			opt_save = OPTION_DEFAULT;
			break;
		case 13:				/* llnl_gamma */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->gflag = 7;	/* llnl D-H */
			i = sscanf(next_char, SCANFORMAT, &s_ptr->dha);
			if (i < 1)
			{
				error_string = sformatf(
						"Expecting activity-coefficient parameter, a.");
				warning_msg(error_string);
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 14:				/* co2_llnl_gamma */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->gflag = 8;	/* llnl CO2 D-H */
			opt_save = OPTION_DEFAULT;
			break;
		case 15:				/* activity water */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->gflag = 9;	/* activity_water/55.5 for HDO, D2O, H2[O18], etc */
			opt_save = OPTION_DEFAULT;
			break;
		case 16:				/* add_logk */
		case 17:				/* add_log_k */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (s_ptr->count_add_logk == 0)
			{
				s_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				s_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(s_ptr->add_logk,
													  (size_t) ((s_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			/* read name */
			if (copy_token(token, &next_char, &i) == EMPTY)
			{
				input_error++;
				error_string = sformatf(
						"Expected the name of a NAMED_EXPRESSION.");
				error_msg(error_string, CONTINUE);
				break;
			}
			s_ptr->add_logk[s_ptr->count_add_logk].name = string_hsave(token);
			/* read coef */
			i = sscanf(next_char, SCANFORMAT,
					   &s_ptr->add_logk[s_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				s_ptr->add_logk[s_ptr->count_add_logk].coef = 1;
			}
			s_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 18:				/* add_constant */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (s_ptr->count_add_logk == 0)
			{
				s_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				s_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(s_ptr->add_logk,
													  (size_t) ((s_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			i = sscanf(next_char, SCANFORMAT,
					   &s_ptr->add_logk[s_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				input_error++;
				error_string = sformatf(
						"Expected the constant to add for log_K definition.");
				error_msg(error_string, CONTINUE);
				break;
			}
			/* set name */
			s_ptr->add_logk[s_ptr->count_add_logk].name =
				string_hsave("XconstantX");
			/* read coef */
			s_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 19:				/* tracer diffusion coefficient */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			i = sscanf(next_char, SCANFORMAT, &s_ptr->dw);
			opt_save = OPTION_DEFAULT;
			break;
		case 20:				/* enrichment factor in the DDL */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			i = sscanf(next_char, SCANFORMAT, &s_ptr->erm_ddl);
			if (s_ptr->erm_ddl < 0)
			{
				error_string = sformatf( "Expecting enrichment factor > 0, "
						"resetting to erm_ddl = 1.");
				warning_msg(error_string);
				s_ptr->erm_ddl = 1.0;
			}
			opt_save = OPTION_DEFAULT;
			break;
/* VP: Density Start */
		case 21:			/* Millero a, b, c, d, e and f coefficients */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			print_density = read_millero_abcdef (next_char, &(s_ptr->millero[0]));
			if (!vm_read)
			{
			/* copy millero parms into logk, for calculating pressure dependency... */
				for (i = 0; i < 3; i++)
					s_ptr->logk[vm0 + i] = s_ptr->millero[i];
			}
			opt_save = OPTION_DEFAULT;
			break;
/* VP: Density End */

		case 22:            /* delta_v */
		case 23:            /* deltav */
		case 24:            /* vm, molar volume */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
					"No reaction defined before option, %s.",
				opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_delta_v_only(next_char, &s_ptr->logk[vm0],
				&s_ptr->original_deltav_units);
			vm_read = true;
			opt_save = OPTION_DEFAULT;
			break;
		case OPTION_DEFAULT:
/*
 *   Get space for species information and parse equation
 */
			s_ptr = NULL;
			if (parse_eq(line, &next_elt, association) == ERROR)
			{
				parse_error++;
				error_msg("Parsing equation.", CONTINUE);
				error_msg(line_save, CONTINUE);
				break;
			}
/*
 *   Get pointer to each species in the reaction, store new species if necessary
 */
			trxn.token[0].s =
				s_store(trxn.token[0].name, trxn.token[0].z, TRUE);
			for (i = 1; i < count_trxn; i++)
			{
				trxn.token[i].s =
					s_store(trxn.token[i].name, trxn.token[i].z, FALSE);
			}
/*
 *   Save element list and carbon, hydrogen, and oxygen in species
 */
			trxn.token[0].s->next_elt = next_elt;
			trxn.token[0].s->next_secondary = NULL;
			for (; next_elt->elt != NULL; next_elt++)
			{
				if (strcmp(next_elt->elt->name, "C") == 0)
				{
					trxn.token[0].s->carbon = next_elt->coef;
				}
				if (strcmp(next_elt->elt->name, "H") == 0)
				{
					trxn.token[0].s->h = next_elt->coef;
				}
				if (strcmp(next_elt->elt->name, "O") == 0)
				{
					trxn.token[0].s->o = next_elt->coef;
				}
			}
/*
 *   Malloc space for species reaction
 */
			trxn.token[0].s->rxn = rxn_alloc(count_trxn + 1);
/*
 *   Copy reaction to reaction for species
 */
			trxn_copy(trxn.token[0].s->rxn);
			s_ptr = trxn.token[0].s;
/*
 *   Default gamma data
 */
			s_ptr->dha = 0.0;
			s_ptr->dhb = 0.0;
			if (equal(s_ptr->z, 0.0, TOL) == TRUE)
			{
				s_ptr->gflag = 0;	/* Uncharged */
				s_ptr->dhb = 0.1;
			}
			else
			{
				s_ptr->gflag = 1;	/* Davies */
			}
/*
 *   Set type for species
 */
			if (strcmp(trxn.token[0].s->name, "H+") == 0)
			{
				s_hplus = trxn.token[0].s;
				s_hplus->type = HPLUS;
			}
			else if (strcmp(trxn.token[0].s->name, "H3O+") == 0)
			{
				s_h3oplus = trxn.token[0].s;
				s_h3oplus->type = HPLUS;
			}
			else if (strcmp(trxn.token[0].s->name, "e-") == 0)
			{
				s_eminus = trxn.token[0].s;
				s_eminus->type = EMINUS;
				s_eminus->gflag = 3;	/* Always 1 */
			}
			else if (strcmp(trxn.token[0].s->name, "H2O") == 0)
			{
				s_h2o = trxn.token[0].s;
				s_h2o->type = H2O;
				s_h2o->gflag = 3;	/* Always 1 */
			}
			else if (strstr(trxn.token[0].s->name, "(s)") != NULL)
			{
				trxn.token[0].s->type = SOLID;
			}
			else if (strcmp(trxn.token[0].s->name, "H2") == 0)
			{
				trxn.token[0].s->type = AQ;
				s_h2 = trxn.token[0].s;
			}
			else if (strcmp(trxn.token[0].s->name, "O2") == 0)
			{
				trxn.token[0].s->type = AQ;
				s_o2 = trxn.token[0].s;
			}
			else
			{
				trxn.token[0].s->type = AQ;
			}
			opt_save = OPTION_DEFAULT;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_use(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution, mix, irreversible reaction, and pure phases to use
 *   in reaction calculation
 */
	int i, l, n_user, return_value;
	char *ptr;
	char token[MAX_LENGTH], token1[MAX_LENGTH];;
/*
 *   Read "use"
 */
	ptr = line;
	copy_token(token, &ptr, &l);
/*
 *   Read keyword
 */
	copy_token(token, &ptr, &l);
	check_key(token);
	if (next_keyword != Keywords::KEY_SOLUTION				&&
		next_keyword != Keywords::KEY_MIX					&&
		next_keyword != Keywords::KEY_KINETICS				&&
		next_keyword != Keywords::KEY_REACTION				&&
		next_keyword != Keywords::KEY_REACTION_TEMPERATURE	&&
		next_keyword != Keywords::KEY_REACTION_PRESSURE		&&
		next_keyword != Keywords::KEY_EQUILIBRIUM_PHASES	&&
		next_keyword != Keywords::KEY_EXCHANGE				&&
		next_keyword != Keywords::KEY_SURFACE				&&
		next_keyword != Keywords::KEY_GAS_PHASE				&&
		next_keyword != Keywords::KEY_SOLID_SOLUTIONS)
	{
		input_error++;
		error_msg("Unknown item in USE keyword", CONTINUE);
		error_msg(line_save, CONTINUE);
		check_line("End of use", FALSE, TRUE, TRUE, TRUE);
		/* empty, eof, keyword, print */
		return (ERROR);
	}
/*
 *   Read number
 */
	strcpy(token1, token);
	for (;;)
	{
		i = copy_token(token, &ptr, &l);
		if (i == DIGIT)
		{
			sscanf(token, "%d", &n_user);
			if (n_user < 0)
			{
				error_msg("Number must be a positive integer.", CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
			}
			if (strstr(token, "-") != NULL)
			{
				error_string = sformatf(
						"USE does not accept a range of numbers, %s.", token);
				warning_msg(error_string);
				error_string = sformatf(
						"Only %s %d will be used in the batch-reaction calculation.",
						token1, n_user);
				warning_msg(error_string);
				error_string = sformatf(
						"NOTE--USE is not needed for ADVECTION and TRANSPORT calculations.");
				warning_msg(error_string);

			}
			break;
		}
		else if (i == EMPTY)
		{
			error_string = sformatf( "No number given, 1 assumed.");
			warning_msg(error_string);
			n_user = 1;
			break;
		}
		else if (token[0] == 'N' || token[0] == 'n')
		{
			n_user = -2;
			break;
		}
	}
	switch (next_keyword)
	{
	case Keywords::KEY_SOLUTION:					/* Solution */
		use.Set_n_solution_user(n_user);
		if (n_user >= 0)
		{
			use.Set_solution_in(true);
		}
		else
		{
			use.Set_solution_in(false);
		}
		break;
	case Keywords::KEY_EQUILIBRIUM_PHASES:					/* Pure phases */
		use.Set_n_pp_assemblage_user(n_user);
		if (n_user >= 0)
		{
			use.Set_pp_assemblage_in(true);
		}
		else
		{
			use.Set_pp_assemblage_in(false);
		}
		break;
	case Keywords::KEY_REACTION:					/* Reaction */
		use.Set_n_reaction_user(n_user);
		if (n_user >= 0)
		{
			use.Set_reaction_in(true);
		}
		else
		{
			use.Set_reaction_in(false);
		}
		break;
	case Keywords::KEY_MIX:					/* Mix */
		use.Set_n_mix_user(n_user);
		if (n_user >= 0)
		{
			use.Set_mix_in(true);
		}
		else
		{
			use.Set_mix_in(false);
		}
		break;
	case Keywords::KEY_EXCHANGE:					/* Ex */
		use.Set_n_exchange_user(n_user);
		if (n_user >= 0)
		{
			use.Set_exchange_in(true);
		}
		else
		{
			use.Set_exchange_in(false);
		}
		break;
	case Keywords::KEY_SURFACE:					/* Surface */
		use.Set_n_surface_user(n_user);
		if (n_user >= 0)
		{
			use.Set_surface_in(true);
		}
		else
		{
			use.Set_surface_in(false);
		}
		break;
	case Keywords::KEY_REACTION_TEMPERATURE:					/* Temperature */
		use.Set_n_temperature_user(n_user);
		if (n_user >= 0)
		{
			use.Set_temperature_in(true);
		}
		else
		{
			use.Set_temperature_in(false);
		}
		break;
	case Keywords::KEY_REACTION_PRESSURE:					/* pressure */
		use.Set_n_pressure_user(n_user);
		if (n_user >= 0)
		{
			use.Set_pressure_in(true);
		}
		else
		{
			use.Set_pressure_in(false);
		}
		break;
	case Keywords::KEY_GAS_PHASE:					/* Gas */
		use.Set_n_gas_phase_user(n_user);
		if (n_user >= 0)
		{
			use.Set_gas_phase_in(true);
		}
		else
		{
			use.Set_gas_phase_in(false);
		}
		break;
	case Keywords::KEY_KINETICS:					/* Kinetics */
		use.Set_n_kinetics_user(n_user);
		if (n_user >= 0)
		{
			use.Set_kinetics_in(true);
		}
		else
		{
			use.Set_kinetics_in(false);
		}
		break;
	case Keywords::KEY_SOLID_SOLUTIONS:					/* solid_solutions */
		use.Set_n_ss_assemblage_user(n_user);
		if (n_user >= 0)
		{
			use.Set_ss_assemblage_in(true);
		}
		else
		{
			use.Set_ss_assemblage_in(false);
		}
		break;
	default:
		input_error++;
		error_msg(line_save, CONTINUE);
		error_msg("Error in switch for USE.", STOP);
		break;
	}
	return_value = check_line("End of use", FALSE, TRUE, TRUE, TRUE);
	/* empty, eof, keyword, print */
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_surface_species(void)
/* ---------------------------------------------------------------------- */
{
	/*
	 *   Read data for surface species, parse equations
	 */
	int i, j;
	int association;
	char token[MAX_LENGTH];
	char *ptr;
	LDBLE offset;

	struct species *s_ptr;
	struct elt_list *next_elt;
	struct rxn_token *token_ptr;

	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"no_check",				/* 0 */
		"check",				/* 1 */
		"mb",					/* 2 */
		"mass_balance",			/* 3 */
		"log_k",				/* 4 */
		"logk",					/* 5 */
		"delta_h",				/* 6 */
		"deltah",				/* 7 */
		"analytical_expression",	/* 8 */
		"a_e",					/* 9 */
		"ae",					/* 10 */
		"mole_balance",			/* 11 */
		"offset",				/* 12 */
		"add_logk",				/* 13 */
		"add_log_k",			/* 14 */
		"add_constant",			/* 15 */
		"cd_music",				/* 16 */
		"music",				/* 17 */
		"delta_v",				/* 18 */
		"deltav",				/* 19 */
		"vm"
	};
	int count_opt_list = 21;

	association = TRUE;
	/*
	 *   Read eqn from file and call parser
	 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	s_ptr = NULL;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SURFACE_SPECIES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* no_check */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->check_equation = FALSE;
			break;
		case 1:				/* check */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			s_ptr->check_equation = TRUE;
			break;
		case 2:				/* mb */
		case 3:				/* mass_balance */
		case 11:				/* mole_balance */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			count_elts = 0;
			paren_count = 0;
			copy_token(token, &next_char, &i);
			s_ptr->mole_balance = string_hsave(token);
			ptr = token;
			s_ptr->next_secondary =
				(struct elt_list *) free_check_null(s_ptr->next_secondary);
			get_secondary_in_species(&ptr, 1.0);
			s_ptr->next_secondary = elt_list_save();
			/* debug
			   for (i = 0; i < count_elts; i++) {
			   output_msg(sformatf("%s\t%f\n", elt_list[i].elt->name,
			   elt_list[i].coef));
			   }
			 */
			opt_save = OPTION_DEFAULT;
			break;
		case 4:				/* log_k */
		case 5:				/* logk */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_log_k_only(next_char, &s_ptr->logk[0]);
			opt_save = OPTION_DEFAULT;
			break;
		case 6:				/* delta_h */
		case 7:				/* deltah */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_delta_h_only(next_char, &s_ptr->logk[1],
							  &s_ptr->original_units);
			opt_save = OPTION_DEFAULT;
			break;
		case 8:				/* analytical_expression */
		case 9:				/* a_e */
		case 10:				/* ae */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_analytical_expression_only(next_char, &(s_ptr->logk[T_A1]));
			opt_save = OPTION_DEFAULT;
			break;
		case 12:				/* offset */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (sscanf(next_char, SCANFORMAT, &offset) != 1)
			{
				error_msg("No offset for log K given", STOP);
			}
			s_ptr->logk[0] += offset;
			opt_save = OPTION_DEFAULT;
			break;
		case 13:				/* add_logk */
		case 14:				/* add_log_k */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (s_ptr->count_add_logk == 0)
			{
				s_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				s_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(s_ptr->add_logk,
													  (size_t) ((s_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			/* read name */
			if (copy_token(token, &next_char, &i) == EMPTY)
			{
				input_error++;
				error_string = sformatf(
						"Expected the name of a NAMED_EXPRESSION.");
				error_msg(error_string, CONTINUE);
				break;
			}
			s_ptr->add_logk[s_ptr->count_add_logk].name = string_hsave(token);
			/* read coef */
			i = sscanf(next_char, SCANFORMAT,
					   &s_ptr->add_logk[s_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				s_ptr->add_logk[s_ptr->count_add_logk].coef = 1;
			}
			s_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 15:				/* add_constant */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (s_ptr->count_add_logk == 0)
			{
				s_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				s_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(s_ptr->add_logk,
													  (size_t) ((s_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (s_ptr->add_logk == NULL)
					malloc_error();
			}
			i = sscanf(next_char, SCANFORMAT,
					   &s_ptr->add_logk[s_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				input_error++;
				error_string = sformatf(
						"Expected the constant to add for log_K definition.");
				error_msg(error_string, CONTINUE);
				break;
			}
			/* set name */
			s_ptr->add_logk[s_ptr->count_add_logk].name =
				string_hsave("XconstantX");
			/* read coef */
			s_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 16:				/* cd_music */
		case 17:				/* music */
			if (s_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			for (j = 0; j < 5; j++)
			{
				if (copy_token(token, &next_char, &i) == EMPTY)
					break;
				if (sscanf(token, SCANFORMAT, &s_ptr->cd_music[j]) != 1)
					break;
			}
			s_ptr->dz[0] = s_ptr->cd_music[0] + s_ptr->cd_music[3] * s_ptr->cd_music[4];
			s_ptr->dz[1] = s_ptr->cd_music[1] + (1 - s_ptr->cd_music[3]) *	s_ptr->cd_music[4];
			s_ptr->dz[2] = s_ptr->cd_music[2];
			for (j = 0; j < 3; j++)
			{
				s_ptr->rxn->dz[j] = s_ptr->dz[j];
			}
			opt_save = OPTION_DEFAULT;
			break;
		case 18:            /* delta_v */
		case 19:            /* deltav */
		case 20:            /* vm, molar volume */
		if (s_ptr == NULL)
			{
				error_string = sformatf(
					"No reaction defined before option, %s.",
				opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_delta_v_only(next_char, &s_ptr->logk[vm0],
				&s_ptr->original_deltav_units);
			opt_save = OPTION_DEFAULT;
			break;
		case OPTION_DEFAULT:
			/*
			 *   Get surface species information and parse equation
			 */
			s_ptr = NULL;
			if (parse_eq(line, &next_elt, association) == ERROR)
			{
				parse_error++;
				error_msg("Parsing equation.", CONTINUE);
				error_msg(line_save, CONTINUE);
				break;
			}
			/*
			 *   Get pointer to each species in the reaction, store new species if necessary
			 */
			trxn.token[0].s =
				s_store(trxn.token[0].name, trxn.token[0].z, TRUE);
			for (i = 1; i < count_trxn; i++)
			{
				trxn.token[i].s =
					s_store(trxn.token[i].name, trxn.token[i].z, FALSE);
			}
			/*
			 *   Save element list and carbon, hydrogen, and oxygen in species
			 */
			trxn.token[0].s->next_elt = next_elt;
			for (; next_elt->elt != NULL; next_elt++)
			{
				if (strcmp(next_elt->elt->name, "C") == 0)
				{
					trxn.token[0].s->carbon = next_elt->coef;
				}
				if (strcmp(next_elt->elt->name, "H") == 0)
				{
					trxn.token[0].s->h = next_elt->coef;
				}
				if (strcmp(next_elt->elt->name, "O") == 0)
				{
					trxn.token[0].s->o = next_elt->coef;
				}
			}
			/*
			 *   Find coefficient of surface in rxn, store in equiv
			 */
			trxn.token[0].s->equiv = 0.0;
			for (i = 1; i < count_trxn; i++)
			{
				if (trxn.token[i].s->type == SURF)
				{
					trxn.token[0].s->equiv = trxn.token[i].coef;
				}
			}
			if (trxn.token[0].s->equiv == 0.0)
				trxn.token[0].s->equiv = 1.0;
			/*
			 *   Malloc space for species reaction
			 */
			trxn.token[0].s->rxn = rxn_alloc(count_trxn + 1);
			/*
			 *   Copy reaction to reaction for species
			 */
			token_ptr = trxn.token[0].s->rxn->token;
			for (i = 0; i < count_trxn; i++)
			{
				token_ptr[i].s = trxn.token[i].s;
				token_ptr[i].coef = trxn.token[i].coef;
			}
			token_ptr[i].s = NULL;
			/*
			 *   Set type for species
			 */
			trxn.token[0].s->type = SURF;
			s_ptr = trxn.token[0].s;
			/*
			 *   Read gamma data
			 */
			s_ptr->gflag = 6;
			s_ptr->dha = 0.0;
			s_ptr->dhb = 0.0;
			opt_save = OPTION_DEFAULT;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_surf(void)
/* ---------------------------------------------------------------------- */
{
	/*
	 *      Reads surface data
	 *
	 *      Arguments:
	 *   none
	 *
	 *      Returns:
	 *   KEYWORD if keyword encountered, input_error may be incremented if
	 *	  a keyword is encountered in an unexpected position
	 *   EOF     if eof encountered while reading mass balance concentrations
	 *   ERROR   if error occurred reading data
	 *
	 */
	int i, j, n, l;
	int count_comps, count_charge;
	int n_user, n_user_end;
	LDBLE conc, area, grams, thickness;
	char *ptr, *ptr1;
	char *description;
	char token[MAX_LENGTH], token1[MAX_LENGTH], name[MAX_LENGTH];
	struct surface *surface_ptr;
	struct surface_charge *charge_ptr;
	struct surface_comp *comp_ptr;

	int return_value, opt;
	char *next_char;
	const char *opt_list[] = {
		"equilibrate",          /* 0 */
		"equil",                /* 1 */
		"diff",                 /* 2 */
		"diffuse_layer",        /* 3 */
		"no_edl",               /* 4 */
		"no_electrostatic",     /* 5 */
		"only_counter_ions",    /* 6 */
		"donnan",               /* 7 */
		"cd_music",             /* 8 */
		"capacitances",         /* 9 */
		"sites",                /* 10 */
		"sites_units",			/* 11 */
		"constant_capacitance", /* 12 */
		"ccm",                  /* 13 */
        "equilibrium",          /* 14 */
		"site_units"            /* 15 */
	};
	int count_opt_list = 16;
	/*
	 * kin_surf is for Surfaces, related to kinetically reacting minerals
	 *    they are defined if "sites" is followed by mineral name:
	 *    Surf_wOH  Manganite  [equilibrium_phases or kinetics]      0.25    4000
	 *    ^Name     mineral    ^switch		 ^prop.factor ^m2/mol
	 */
	/*
	 *   Read surface number and description
	 */
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);
	/*
	 *   Find space for surface data
	 */
	surface_ptr = surface_search(n_user, &n, FALSE);
	if (surface_ptr != NULL)
	{
		surface_free(&surface[n]);
	}
	else
	{
		n = count_surface++;
		space((void **) ((void *) &surface), count_surface, &max_surface,
			  sizeof(struct surface));
	}
	/*
	 *   Initialize
	 */
	surface_init(&(surface[n]), n_user, n_user_end, description);
	free_check_null(description);

	if (use.Get_surface_in() == FALSE)
	{
		use.Set_surface_in(true);
		use.Set_n_surface_user(n_user);
	}
	/*
	 *   Read surface data
	 */
	count_charge = 0;
	count_comps = 0;
	return_value = UNKNOWN;
	comp_ptr = NULL;
	charge_ptr = NULL;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SURFACE keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* equilibrate */
		case 1:             /* equil */
		case 14:			/* equilibrium */
			for (;;)
			{
				i = copy_token(token, &next_char, &l);
				if (i == DIGIT)
				{
					sscanf(token, "%d", &surface[n].n_solution);
					surface[n].solution_equilibria = TRUE;
					break;
				}
				if (i == EMPTY)
				{
					error_msg
						("Expected a solution number with which to equilibrate surface.",
						 CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
					break;
				}
			}
			break;
		case 2:				/* diffuse_layer */
		case 3:
			surface[n].thickness = 1e-8;
			surface[n].dl_type = BORKOVEK_DL;
			sscanf(next_char, SCANFORMAT, &surface[n].thickness);
			/*		surface[n].thickness = thickness;
			   }
			 */ break;
		case 4:				/* no electrostatic */
		case 5:
			/*surface[n].edl = FALSE; */
			surface[n].type = NO_EDL;
			break;
		case 6:
			surface[n].only_counter_ions = TRUE;
			break;
		case 7:				/* donnan for DL conc's */
			surface[n].dl_type = DONNAN_DL;
			/*surface[n].diffuse_layer = TRUE; */
			thickness = 0.0;
			for (;;)
			{
				i = copy_token(token, &next_char, &l);
				if (i == DIGIT)
				{
					sscanf(token, SCANFORMAT, &surface[n].thickness);
					thickness = 1;
					continue;
				}
				else if (i != EMPTY)
				{
					if (token[0] == 'D' || token[0] == 'd')
					{
						if (thickness != 0)
						{
							error_msg
								("You must enter EITHER thickness OR Debye lengths (1/k),\n	   and relative DDL viscosity, DDL limit.\nCorrect is (for example): -donnan 1e-8 viscosity 0.5\n or (default values):     -donnan debye_lengths 1 viscosity 1 limit 0.8",
								 CONTINUE);
							error_msg(line_save, CONTINUE);
							input_error++;
							break;
						}
						j = copy_token(token1, &next_char, &l);
						if (j == DIGIT)
						{
							sscanf(token1, SCANFORMAT,
								   &surface[n].debye_lengths);
							continue;
						}
						else if (j != EMPTY)
						{
							error_msg
								("Expected number of Debye lengths (1/k).",
								 CONTINUE);
							error_msg(line_save, CONTINUE);
							input_error++;
							break;
						}
					}
					else if (token[0] == 'V' || token[0] == 'v')
					{
						j = copy_token(token1, &next_char, &l);
						if (j == DIGIT)
						{
							sscanf(token1, SCANFORMAT,
								   &surface[n].DDL_viscosity);
							continue;
						}
						else if (j != EMPTY)
						{
							error_msg
								("Expected number for relative DDL viscosity.",
								 CONTINUE);
							error_msg(line_save, CONTINUE);
							input_error++;
							break;
						}
					}
					else if (token[0] == 'L' || token[0] == 'l')
					{
						j = copy_token(token1, &next_char, &l);
						if (j == DIGIT)
						{
							sscanf(token1, SCANFORMAT, &surface[n].DDL_limit);
							continue;
						}
						else if (j != EMPTY)
						{
							error_msg
								("Expected number for maximum of DDL water as fraction of bulk water.",
								 CONTINUE);
							error_msg(line_save, CONTINUE);
							input_error++;
							break;
						}
					}
					else
					{
						error_msg
							("Expected diffuse layer thickness (m) or Debye lengths (1/k) for \n\tcalculating the thickness, and relative DDL viscosity and DDL limit.\nCorrect is (for example): -donnan 1e-8 visc 0.5\n or (default values):     -donnan debye_lengths 1 visc 1 lim 0.8",
							 CONTINUE);
						error_msg(line_save, CONTINUE);
						input_error++;
						break;
					}
				}
				else
					break;
			}
			break;
		case 8:				/* cd_music */
			/*surface[n].edl = FALSE; */
			surface[n].type = CD_MUSIC;
			break;
		case 9:				/* capacitances */
			/*
			 *   Read capacitance for CD_MUSIC
			 */
			if (charge_ptr == NULL)
			{
				error_msg
					("Surface component has not been defined for capacitances.",
					 CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				break;
			}
			copy_token(token1, &next_char, &l);
			if (sscanf(token1, SCANFORMAT, &grams) == 1)
			{
				charge_ptr->capacitance[0] = grams;
			}
			copy_token(token1, &next_char, &l);
			if (sscanf(token1, SCANFORMAT, &grams) == 1)
			{
				charge_ptr->capacitance[1] = grams;
			}
			break;
		case 10:				/* sites */
		case 11:				/* sites_units */
		case 15:				/* site_units */
			j = copy_token(token1, &next_char, &l);
			if (token1[0] == 'A' || token1[0] == 'a')
			{
				surface[n].sites_units = SITES_ABSOLUTE;
			}
			else if (token1[0] == 'D' || token1[0] == 'd')
			{
				surface[n].sites_units = SITES_DENSITY;
			}
			else
			{
				error_msg
					("Character string expected to be 'Absolute' or 'Density' to define the units of the first item in the definition of a surface component.\n",
					 CONTINUE);
				input_error++;
				break;
			}
			break;
		case 12:			    /* constant_capacitance */
		case 13:			    /* ccm */
			surface[n].type = CCM;
			copy_token(token1, &next_char, &l);
			if (sscanf(token1, SCANFORMAT, &(charge_ptr->capacitance[0])) != 1)
			{
				error_msg("Expected capacitance for constant_capacitance model.\n",
					 CONTINUE);
				input_error++;
				break;
			}

			/* constant capacitance model not implemented yet */
			error_msg("Constant capacitance model not implemented.", CONTINUE);
			input_error++;

			break;
		case OPTION_DEFAULT:
			/*
			 *   Read surface component
			 */
			ptr = line;
			i = copy_token(token, &ptr, &l);
			if (i != UPPER && token[0] != '[')
			{
				error_msg
					("Expected surface name to begin with a capital letter.",
					 CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				break;
			}
			/*
			 *   Realloc space for new surface component
			 */
			surface[n].comps =
				(struct surface_comp *) PHRQ_realloc(surface[n].comps,
													 (size_t) (count_comps +
															   1) *
													 sizeof(struct
															surface_comp));
			if (surface[n].comps == NULL)
				malloc_error();
			surface[n].comps[count_comps].formula = string_hsave(token);
			surface[n].comps[count_comps].formula_totals = NULL;
			surface[n].comps[count_comps].formula_z = 0.0;
			surface[n].comps[count_comps].moles = 0;
			surface[n].comps[count_comps].la = 0;
			surface[n].comps[count_comps].charge = 0;
			surface[n].comps[count_comps].cb = 0;
			surface[n].comps[count_comps].phase_name = NULL;
			surface[n].comps[count_comps].phase_proportion = 0;
			surface[n].comps[count_comps].rate_name = NULL;
			comp_ptr = &(surface[n].comps[count_comps]);
			surface[n].comps[count_comps].Dw = 0;
			i = copy_token(token1, &ptr, &l);
			if (i == DIGIT)
			{
				/*
				 *   Read surface concentration
				 */
				/* surface concentration is read directly */
				if (sscanf(token1, SCANFORMAT, &conc) < 1)
				{
					error_msg("Expected number of surface sites in moles.",
							  CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
					break;
				}
				surface[n].comps[count_comps].moles = conc;
				/*
				 *   Read equilibrium phase name or kinetics rate name
				 */
			}
			else if (i != EMPTY)
			{

				/* surface conc. is related to mineral or kinetics */
				surface[n].comps[count_comps].phase_name =
					string_hsave(token1);
				j = copy_token(token1, &ptr, &l);

				/* read optional 'equilibrium_phases' or 'kinetics' */
				if (j == DIGIT)
				{
					surface[n].related_phases = TRUE;
				}
				else
				{
					if (token1[0] == 'K' || token1[0] == 'k')
					{
						surface[n].comps[count_comps].rate_name =
							surface[n].comps[count_comps].phase_name;
						surface[n].comps[count_comps].phase_name = NULL;
						surface[n].related_rate = TRUE;
					}
					else if (token1[0] == 'E' || token1[0] == 'e')
					{
						surface[n].related_phases = TRUE;
						surface[n].comps[count_comps].rate_name = NULL;
					}
					else
					{
						error_msg
							("Character string expected to be 'equilibrium_phase' or 'kinetics' to relate surface to mineral or kinetic reaction.\n",
							 CONTINUE);
						input_error++;
						break;
					}
					j = copy_token(token1, &ptr, &l);
				}

				/* read proportion */
				if (j != DIGIT)
				{
					error_msg
						("Expected a coefficient to relate surface to mineral or kinetic reaction.\n",
						 CONTINUE);
					input_error++;
					break;
				}
				sscanf(token1, SCANFORMAT,
					   &surface[n].comps[count_comps].phase_proportion);
				/* real conc must be defined in tidy_model */
				conc = 1.0;
			}
			else
			{
				error_msg
					("Expected concentration of surface, mineral name, or kinetic reaction name.",
					 CONTINUE);
				error_msg(line_save, CONTINUE);
				input_error++;
				break;
			}
			/*
			 *   Accumulate elements in elt_list
			 */
			count_elts = 0;
			paren_count = 0;
			ptr1 = token;
			get_elts_in_species(&ptr1, conc);

			/*
			 *   save formula for adjusting number of exchange sites
			 */

			ptr1 = token;
			get_token(&ptr1, token1,
					  &surface[n].comps[count_comps].formula_z, &l);
			surface[n].comps[count_comps].formula_totals = elt_list_save();
			surface[n].comps[count_comps].totals = elt_list_save();
			/*
			 *   Search for charge structure
			 */
			ptr1 = token;
			get_elt(&ptr1, name, &l);
			ptr1 = strchr(name, '_');
			if (ptr1 != NULL)
				ptr1[0] = '\0';
			for (i = 0; i < count_charge; i++)
			{
				if (strcmp(surface[n].charge[i].name, name) == 0)
					break;
			}
			if (i >= count_charge)
			{
				surface[n].charge =
					(struct surface_charge *) PHRQ_realloc(surface[n].charge,
														   (size_t)
														   (count_charge +
															1) *
														   sizeof(struct
																  surface_charge));
				if (surface[n].charge == NULL)
					malloc_error();
				i = count_charge;
				surface[n].charge[i].name = string_hsave(name);
				if (surface[n].comps[count_comps].phase_name == NULL
					&& surface[n].comps[count_comps].rate_name == NULL)
				{
					surface[n].charge[i].specific_area = 600.0;
					surface[n].charge[i].grams = 0.0;
				}
				else
				{
					surface[n].charge[i].specific_area = 0.0;
					surface[n].charge[i].grams = 1.0;
				}
				surface[n].charge[i].charge_balance = 0.0;
				surface[n].charge[i].mass_water = 0.0;
				surface[n].charge[i].diffuse_layer_totals = NULL;
				surface[n].charge[i].count_g = 0;
				surface[n].charge[i].g = NULL;
				surface[n].charge[i].la_psi = 0;
				surface[n].charge[i].la_psi1 = 0;
				surface[n].charge[i].la_psi2 = 0;
				surface[n].charge[i].psi = 0;
				surface[n].charge[i].psi1 = 0;
				surface[n].charge[i].psi2 = 0;
				surface[n].charge[i].capacitance[0] = 1.0;
				surface[n].charge[i].capacitance[1] = 5.0;
				surface[n].charge[i].sigma0 = 0;
				surface[n].charge[i].sigma1 = 0;
				surface[n].charge[i].sigma2 = 0;
				surface[n].charge[i].sigmaddl = 0;
				count_charge++;
			}
			charge_ptr = &(surface[n].charge[i]);
			surface[n].comps[count_comps].charge = i;
			count_comps++;
			/*
			 *   Read surface area (m2/g)
			 */
			copy_token(token1, &ptr, &l);
			if (sscanf(token1, SCANFORMAT, &area) == 1)
			{
				surface[n].charge[i].specific_area = area;
			}
			else
			{
				break;
			}
			/*
			 *   Read grams of solid (g)
			 */
			copy_token(token1, &ptr, &l);
			if (sscanf(token1, SCANFORMAT, &grams) == 1)
			{
				surface[n].charge[i].grams = grams;
			}
			/* read Dw */
			copy_token(token1, &ptr, &l);
			str_tolower(token1);
			if (strcmp(token1, "dw") == 0)
			{
				j = copy_token(token1, &ptr, &l);
				if (j != DIGIT)
				{
					error_msg
						("Expected surface diffusion coefficient (m2/s).\n",
						 CONTINUE);
					input_error++;
					break;
				}
				else
				{
					sscanf(token1, SCANFORMAT,
						   &surface[n].comps[count_comps - 1].Dw);
					if (surface[n].comps[count_comps - 1].Dw > 0)
					{
						surface[n].transport = TRUE;
						if (surface[n].related_rate == TRUE
							|| surface[n].related_phases == TRUE)
						{
							error_msg
								("Can't transport surfaces related to phases or rates (yet).",
								 CONTINUE);
							input_error++;
						}
					}
					break;
				}
			}
			break;

		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	surface[n].count_comps = count_comps;
	surface[n].count_charge = count_charge;

	/* sort charge */
	qsort(surface[n].charge,
		  (size_t) surface[n].count_charge,
		  (size_t) sizeof(struct surface_charge), surface_charge_compare);

	/* renumber charge variable in comps */
	for (i = 0; i < surface[n].count_comps; i++)
	{
		char * temp_formula = string_duplicate(surface[n].comps[i].formula);
		ptr = temp_formula;
		copy_token(token, &ptr, &l);
		ptr1 = token;
		get_elt(&ptr1, name, &l);
		ptr1 = strchr(name, '_');
		if (ptr1 != NULL)
			ptr1[0] = '\0';
		for (j = 0; j < surface[n].count_charge; j++)
		{
			if (strcmp(surface[n].charge[j].name, name) == 0)
			{
				surface[n].comps[i].charge = j;
				break;
			}
		}
		free_check_null(temp_formula);
		assert(j < surface[n].count_charge);
	}

	/*
	 * copy Dw > 0 to all comps with the same charge structure, check edl & dif for surface transport
	 */
	if (surface[n].transport == TRUE)
	{
		/*
		   if (surface[n].edl == FALSE || surface[n].diffuse_layer == FALSE) {
		   surface[n].edl = TRUE;
		   surface[n].diffuse_layer = TRUE;
		   warning_msg("Can't transport without edl and diffuse layer, set to true");
		   }
		 */
		if (surface[n].type <= NO_EDL)
		{
			input_error++;
			error_msg
				("Must use default Dzombak and Morel or -cd_music for surface transport.",
				 CONTINUE);
		}
		if (surface[n].dl_type <= NO_DL)
		{
			input_error++;
			error_msg
				("Must use -donnan or -diffuse_layer for surface transport.",
				 CONTINUE);
		}
		for (i = 0; i < count_comps; i++)
		{
			if (surface[n].comps[i].Dw > 0)
			{
				for (j = 0; j < count_comps; j++)
				{
					if (surface[n].comps[j].charge ==
						surface[n].comps[i].charge)
						surface[n].comps[j].Dw = surface[n].comps[i].Dw;
				}
			}
		}
	}
	/*
	 *   Make sure surface area is defined
	 */
	/*if (surface[n].edl == TRUE) { */
	if (surface[n].type == DDL || surface[n].type == CD_MUSIC)
	{
		for (i = 0; i < count_charge; i++)
		{
			if (surface[n].charge[i].grams *
				surface[n].charge[i].specific_area <= 0)
			{
				error_string = sformatf( "Surface area not defined for %s.\n",
						surface[n].charge[i].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
		}
	}
	/*
	 *  Logical checks
	 */
	if (surface[n].type == UNKNOWN_DL)
	{
		error_string = sformatf( "Unknown surface type.\n");
		error_msg(error_string, CONTINUE);
		input_error++;
	}
	else if (surface[n].type == NO_EDL)
	{
		if (surface[n].dl_type != NO_DL)
		{
			error_string = sformatf(
					"No electrostatic term calculations do not allow calculation of the diffuse layer composition.\n");
			warning_msg(error_string);
			surface[n].dl_type = NO_DL;
		}
	}
	else if (surface[n].type == DDL)
	{
		/* all values of dl_type are valid */
	}
	else if (surface[n].type == CD_MUSIC)
	{
		if (surface[n].dl_type == BORKOVEK_DL)
		{
			error_string = sformatf(
					"Borkovec and Westall diffuse layer calculation is not allowed with a CD_MUSIC surface.\n\tUsing Donnan diffuse layer calculation.");
			warning_msg(error_string);
			surface[n].dl_type = DONNAN_DL;
		}
		if (surface[n].debye_lengths > 0)
		{
			error_msg
				("Variable DDL thickness is not permitted in CD_MUSIC. Fix DDL thickness\n\tfor example (default value): -donnan 1e-8",
				 CONTINUE);
			input_error++;
		}
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_surface_master_species(void)
/* ---------------------------------------------------------------------- */
{
	/*
	 *   Reads master species data from data file or input file
	 */
	int l, return_value;
	char *ptr, *ptr1;
	LDBLE l_z;
	struct master *m_ptr;
	struct species *s_ptr;
	char token[MAX_LENGTH], token1[MAX_LENGTH];
	int opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"capacitance",			/* 0 */
		"cd_music_capacitance"	/* 1 */
	};
	int count_opt_list = 0;
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	m_ptr = NULL;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SURFACE_MASTER_SPECIES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case OPTION_DEFAULT:
			/*
			 *   Get "element" name with valence, allocate space, store
			 */
			ptr = line;
			/*
			 *   Get element name and save pointer to character string
			 */
			if (copy_token(token, &ptr, &l) != UPPER && token[0] != '[')
			{
				parse_error++;
				error_msg("Reading element for master species.", CONTINUE);
				error_msg(line_save, CONTINUE);
				continue;
			}
			replace("(+", "(", token);
			/*
			 *   Delete master if it exists
			 */
			master_delete(token);
			/*
			 *   Increase pointer array, if necessary,  and malloc space
			 */
			if (count_master + 4 >= max_master)
			{
				space((void **) ((void *) &master), count_master + 4,
					  &max_master, sizeof(struct master *));
			}
			/*
			 *   Save values in master and species structure for surface sites
			 */
			master[count_master] = master_alloc();
			master[count_master]->type = SURF;
			master[count_master]->elt = element_store(token);
			m_ptr = master[count_master];
			if (copy_token(token, &ptr, &l) != UPPER && token[0] != '[')
			{
				parse_error++;
				error_msg("Reading surface master species name.", CONTINUE);
				error_msg(line_save, CONTINUE);
				continue;
			}
			s_ptr = s_search(token);
			if (s_ptr != NULL)
			{
				master[count_master]->s = s_ptr;
			}
			else
			{
				ptr1 = token;
				get_token(&ptr1, token1, &l_z, &l);
				master[count_master]->s = s_store(token1, l_z, FALSE);
			}
			master[count_master]->primary = TRUE;
			strcpy(token, master[count_master]->elt->name);
			count_master++;
			/*
			 *   Save values in master and species structure for surface psi
			 */
			strcpy(token1, token);
			replace("_", " ", token1);
			ptr1 = token1;
			copy_token(token, &ptr1, &l);
			strcat(token, "_psi");
			add_psi_master_species(token);
			opt_save = OPTION_DEFAULT;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
add_psi_master_species(char *token)
/* ---------------------------------------------------------------------- */
{
	struct species *s_ptr;
	struct master *master_ptr;
	char *ptr;
	char token1[MAX_LENGTH];
	int i, n, plane;

	strcpy(token1, token);
	for (plane = SURF_PSI; plane <= SURF_PSI2; plane++)
	{
		strcpy(token, token1);
		switch (plane)
		{
		case SURF_PSI:
			break;
		case SURF_PSI1:
			strcat(token, "b");
			break;
		case SURF_PSI2:
			strcat(token, "d");
			break;
		}
		master_ptr = master_search(token, &n);
		if (master_ptr == NULL)
		{
			master[count_master] = master_alloc();
			master[count_master]->type = plane;
			master[count_master]->elt = element_store(token);
			s_ptr = s_search(token);
			if (s_ptr != NULL)
			{
				master[count_master]->s = s_ptr;
			}
			else
			{
				master[count_master]->s = s_store(token, 0.0, FALSE);
			}
			count_elts = 0;
			paren_count = 0;
			ptr = token;
			get_elts_in_species(&ptr, 1.0);
			master[count_master]->s->next_elt = elt_list_save();
			master[count_master]->s->type = plane;
			master[count_master]->primary = TRUE;
			master[count_master]->s->rxn = rxn_alloc(3);
			/*
			 *   Define reaction for psi
			 */
			for (i = 0; i < MAX_LOG_K_INDICES; i++)
			{
				master[count_master]->s->rxn->logk[i] = 0.0;
			}
			master[count_master]->s->rxn->token[0].s =
				master[count_master]->s;
			master[count_master]->s->rxn->token[0].coef = -1.0;
			master[count_master]->s->rxn->token[1].s =
				master[count_master]->s;
			master[count_master]->s->rxn->token[1].coef = 1.0;
			master[count_master]->s->rxn->token[2].s = NULL;
			count_master++;
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_title(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads title for simulation
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	char *ptr, *ptr1;
	int l, title_x_length, line_length;
	int return_value;
	char token[MAX_LENGTH];
/*
 *   Read anything after keyword
 */
	ptr = line;
	copy_token(token, &ptr, &l);
	ptr1 = ptr;
	title_x = (char *) free_check_null(title_x);
	if (copy_token(token, &ptr, &l) != EMPTY)
	{
		title_x = string_duplicate(ptr1);
	}
	else
	{
		title_x = (char *) PHRQ_malloc(sizeof(char));
		if (title_x == NULL)
			malloc_error();
		title_x[0] = '\0';
	}

/*
 *   Read additonal lines
 */
	for (;;)
	{
		return_value = check_line("title", TRUE, TRUE, TRUE, TRUE);
		/* empty, eof, keyword, print */
		if (return_value == EOF || return_value == KEYWORD)
			break;
/*
 *   append line to title_x
 */
		title_x_length = (int) strlen(title_x);
		line_length = (int) strlen(line);
		title_x =
			(char *) PHRQ_realloc(title_x,
								  (size_t) (title_x_length + line_length +
											2) * sizeof(char));
		if (title_x == NULL)
			malloc_error();
		if (title_x_length > 0)
		{
			title_x[title_x_length] = '\n';
			title_x[title_x_length + 1] = '\0';
		}
		strcat(title_x, line);
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_advection(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads advection information
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
/*
 *   Read advection parameters:
 *	number of cells;
 *	number of shifts;
 */
	char *ptr;
	char *description;
	int n_user, n_user_end, i;

	int count_punch, count_print;
	int *punch_temp, *print_temp;
	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"cells",				/* 0 */
		"shifts",				/* 1 */
		"print",				/* 2 */
		"selected_output",		/* 3 */
		"punch",				/* 4 */
		"print_cells",			/* 5 */
		"selected_cells",		/* 6 */
		"time_step",			/* 7 */
		"timest",				/* 8 */
		"output",				/* 9 */
		"output_frequency",		/* 10 */
		"selected_output_frequency",	/* 11 */
		"punch_frequency",		/* 12 */
		"print_frequency",		/* 13 */
		"punch_cells",			/* 14 */
		"initial_time",			/* 15 */
		"warning",				/* 16 */
		"warnings"				/* 17 */
	};
	int count_opt_list = 18;
/*
 *   Read advection number (not currently used)
 */
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);
	description = (char *) free_check_null(description);
/*
 *   Set use data
 */
	use.Set_advect_in(true);
	count_ad_cells = 0;
	count_ad_shifts = 0;
	print_ad_modulus = 1;
	punch_ad_modulus = 1;
	count_punch = 0;
	count_print = 0;
	punch_temp = (int *) PHRQ_malloc(sizeof(int));
	if (punch_temp == NULL)
		malloc_error();
	print_temp = (int *) PHRQ_malloc(sizeof(int));
	if (print_temp == NULL)
		malloc_error();
/*
 *   Read lines
 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in ADVECTION keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* cells */
			sscanf(next_char, "%d", &count_ad_cells);
			opt_save = OPTION_DEFAULT;
			break;
		case 1:				/* shifts */
			sscanf(next_char, "%d", &count_ad_shifts);
			opt_save = OPTION_DEFAULT;
			break;
		case 2:				/* print */
		case 5:				/* print_cells */
			print_temp =
				read_list_ints_range(&next_char, &count_print, TRUE,
									 print_temp);
			opt_save = 2;
			break;
		case 3:				/* selected_output */
		case 11:				/* selected_output_frequency */
		case 12:				/* punch_frequency */
			sscanf(next_char, "%d", &punch_ad_modulus);
			opt_save = OPTION_DEFAULT;
			if (punch_ad_modulus <= 0)
			{
				error_string = sformatf(
						"Punch frequency must be greater than 0. Frequency set to 1000.");
				warning_msg(error_string);
				punch_ad_modulus = 1000;
			}
			break;
		case 4:				/* punch */
		case 14:				/* punch_cells */
		case 6:				/* selected_cells */
			punch_temp =
				read_list_ints_range(&next_char, &count_punch, TRUE,
									 punch_temp);
			opt_save = 4;
			break;
		case 7:				/* time_step */
		case 8:				/* timest */
			sscanf(next_char, SCANFORMAT, &advection_kin_time);
			advection_kin_time_defined = TRUE;
			opt_save = OPTION_DEFAULT;
			break;
		case 9:				/* output */
		case 10:				/* output_frequency */
		case 13:				/* print_frequency */
			sscanf(next_char, "%d", &print_ad_modulus);
			opt_save = OPTION_DEFAULT;
			if (print_ad_modulus <= 0)
			{
				error_string = sformatf(
						"Print frequency must be greater than 0. Frequency set to 1000.");
				warning_msg(error_string);
				print_ad_modulus = 1000;
			}
			break;
		case 15:				/* initial_time */
			sscanf(next_char, SCANFORMAT, &initial_total_time);
			opt_save = OPTION_DEFAULT;
			break;
		case 16:				/* warning */
		case 17:				/* warnings */
			advection_warnings = get_true_false(next_char, TRUE);
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
/*
 *   Fill in data for punch
 */
	advection_punch =
		(int *) PHRQ_realloc(advection_punch,
							 (size_t) (count_ad_cells + 1) * sizeof(int));
	if (advection_punch == NULL)
		malloc_error();
	if (count_punch != 0)
	{
		for (i = 0; i < count_ad_cells; i++)
			advection_punch[i] = FALSE;
		for (i = 0; i < count_punch; i++)
		{
			if (punch_temp[i] > count_ad_cells || punch_temp[i] < 1)
			{
				error_string = sformatf(
						"Cell number for punch is out of range, %d. Request ignored.",
						punch_temp[i]);
				warning_msg(error_string);
			}
			else
			{
				advection_punch[punch_temp[i] - 1] = TRUE;
			}
		}
	}
	else
	{
		for (i = 0; i < count_ad_cells; i++)
			advection_punch[i] = TRUE;
	}
	punch_temp = (int *) free_check_null(punch_temp);
/*
 *   Fill in data for print
 */
	advection_print =
		(int *) PHRQ_realloc(advection_print,
							 (size_t) (count_ad_cells + 1) * sizeof(int));
	if (advection_print == NULL)
		malloc_error();
	if (count_print != 0)
	{
		for (i = 0; i < count_ad_cells; i++)
			advection_print[i] = FALSE;
		for (i = 0; i < count_print; i++)
		{
			if (print_temp[i] > count_ad_cells || print_temp[i] < 1)
			{
				error_string = sformatf(
						"Cell number for print is out of range, %d. Request ignored.",
						print_temp[i]);
				warning_msg(error_string);
			}
			else
			{
				advection_print[print_temp[i] - 1] = TRUE;
			}
		}
	}
	else
	{
		for (i = 0; i < count_ad_cells; i++)
			advection_print[i] = TRUE;
	}
	print_temp = (int *) free_check_null(print_temp);
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_debug(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads knobs and debugging info
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int return_value, opt;
	char *next_char;
	const char *opt_list[] = {
		"iterations",			/* 0 */
		"tolerance",			/* 1 */
		"step_size",			/* 2 */
		"pe_step_size",			/* 3 */
		"scale_pure_phases",	/* 4 */
		"diagonal_scale",		/* 5 */
		"debug_model",			/* 6 */
		"debug_prep",			/* 7 */
		"debug_set",			/* 8 */
		"debug_inverse",		/* 9 */
		"logfile",				/* 10 */
		"log_file",				/* 11 */
		"debug_diffuse_layer",	/* 12 */
		"delay_mass_water",		/* 13 */
		"convergence_tolerance",	/* 14 */
		"numerical_derivatives",	/* 15 */
		"tries",					/* 16 */
		"try",						/* 17 */
		"numerical_fixed_volume",    /* 18 */
		"force_numerical_fixed_volume"    /* 19 */
	};
	int count_opt_list = 20;
/*
 *   Read parameters:
 *	ineq_tol;
 *	step_size;
 *	pe_step_size;
 *	pp_scale;
 *	diagonal_scale;
 */
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in KNOBS keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* iterations */
			sscanf(next_char, "%d", &itmax);
			break;
		case 1:				/* tolerance */
			sscanf(next_char, SCANFORMAT, &ineq_tol);
			break;
		case 2:				/* step_size */
			sscanf(next_char, SCANFORMAT, &step_size);
			break;
		case 3:				/* pe_step_size */
			sscanf(next_char, SCANFORMAT, &pe_step_size);
			break;
		case 4:				/* pp_scale */
			sscanf(next_char, SCANFORMAT, &pp_scale);
			break;
		case 5:				/* diagonal_scale */
			diagonal_scale = get_true_false(next_char, TRUE);
			break;
		case 6:				/* debug_model */
			debug_model = get_true_false(next_char, TRUE);
			break;
		case 7:				/* debug_prep */
			debug_prep = get_true_false(next_char, TRUE);
			break;
		case 8:				/* debug_set */
			debug_set = get_true_false(next_char, TRUE);
			break;
		case 9:				/* debug_inverse */
			debug_inverse = get_true_false(next_char, TRUE);
			break;
		case 10:				/* logfile */
		case 11:				/* log_file */
			pr.logfile = get_true_false(next_char, TRUE);
			if (phast == TRUE)
			{
				pr.logfile = FALSE;
				warning_msg("PHREEQC log file is disabled in PHAST");
			}
			phrq_io->Set_log_on(pr.logfile == TRUE);
			break;
		case 12:				/* debug_diffuse_layer */
			debug_diffuse_layer = get_true_false(next_char, TRUE);
			break;
		case 13:				/* delay_mass_water */
			delay_mass_water = get_true_false(next_char, TRUE);
			break;
		case 14:				/* convergence_tolerance */
			{
				LDBLE ct;
				sscanf(next_char, SCANFORMAT, &ct);
				convergence_tolerance = ct;
				//if (punch.high_precision == TRUE)
				//{
				//	if (ct < 1e-12)
				//	{
				//		convergence_tolerance = ct;
				//	}
				//}
				//else
				//{
				//	{
				//		convergence_tolerance = ct;
				//	}
				//}
			}
			break;
		case 15:				/* numerical_derivatives */
			numerical_deriv = get_true_false(next_char, TRUE);
			break;
		case 16:				/* tries */
		case 17:				/* try */
			sscanf(next_char, "%d", &max_tries);
			break;
		case 18:				/* debug_inverse */
			numerical_fixed_volume = (get_true_false(next_char, TRUE) == TRUE);
			break;
		case 19:				/* debug_inverse */
			force_numerical_fixed_volume = (get_true_false(next_char, TRUE) == TRUE);
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_print(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads printing info
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int return_value, opt, l;
	char *next_char;
	char token[MAX_LENGTH];
	LDBLE num;
	const char *opt_list[] = {
		"reset",				/* 0 */
		"gas_phase",			/* 1 */
		"pure_phase",			/* 2 */
		"surface",				/* 3 */
		"exchange",				/* 4 */
		"totals",				/* 5 */
		"eh",					/* 6 */
		"species",				/* 7 */
		"saturation_indices",	/* 8 */
		"si",					/* 9 same a 8 */
		"reaction",				/* 10 irrev, not used, pr.use */
		"mix",					/* 11 */
		"use",					/* 12 */
		"selected_output",		/* 13 */
		"equilibrium_phases",	/* 14 same as 2 */
		"equilibria",	        /* 15 same as 2 */
		"equilibrium",			/* 16 same as 2 */
		"pure",					/* 17 same as 2 */
		"other",				/* 18 same as 12 */
		"status",				/* 19 */
		"inverse",				/* 20 */
		"kinetics",				/* 21 */
		"dump",					/* 22 */
		"user_print",			/* 23 */
		"user_pr",				/* 24 */
		"solid_solution",		/* 25 */
		"solid_solutions",		/* 26 */
		"inverse_modeling",		/* 27 */
		"headings",				/* 28 */
		"heading",				/* 29 */
		"user_graph",			/* 30 */
		"echo_input",			/* 31 */
		"warning",				/* 32 */
		"warnings",				/* 33 */
		"initial_isotopes",		/* 34 */
		"isotope_ratios",		/* 35 */
		"isotope_alphas",		/* 36 */
		"censor_species",		/* 37 */
		"alkalinity",			/* 38 */
		"equilibrium_phase"     /* 39 */
	};

	int count_opt_list = 40;
	int value;
/*
 *   Read flags:
 */
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in PRINT keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* reset */
			value = get_true_false(next_char, TRUE);
			pr.kinetics = value;
			pr.gas_phase = value;
			pr.pp_assemblage = value;
			pr.surface = value;
			pr.exchange = value;
			pr.totals = value;
			pr.eh = value;
			pr.species = value;
			pr.saturation_indices = value;
			pr.irrev = value;
			pr.mix = value;
			pr.reaction = value;
			pr.use = value;
			pr.inverse = value;
			pr.user_print = value;
			pr.ss_assemblage = value;
			pr.headings = value;
			pr.initial_isotopes = value;
			pr.isotope_ratios = value;
			pr.isotope_alphas = value;
			pr.echo_input = value;
			break;
		case 1:				/* gas_phase */
			pr.gas_phase = get_true_false(next_char, TRUE);
			break;
		case 2:				/* pure_phase */
		case 14:				/* equilibrium_phases */
		case 15:				/* equilibria */
		case 16:				/* equilibrium */
		case 17:				/* pure */
		case 39:                /* equilibrium_phase */
			pr.pp_assemblage = get_true_false(next_char, TRUE);
			break;
		case 3:				/* surface */
			pr.surface = get_true_false(next_char, TRUE);
			break;
		case 4:				/* exchange */
			pr.exchange = get_true_false(next_char, TRUE);
			break;
		case 5:				/* totals */
			pr.totals = get_true_false(next_char, TRUE);
			break;
		case 6:				/* eh */
			pr.eh = get_true_false(next_char, TRUE);
			break;
		case 7:				/* species */
			pr.species = get_true_false(next_char, TRUE);
			break;
		case 8:
		case 9:				/* saturation_indices */
			pr.saturation_indices = get_true_false(next_char, TRUE);
			break;
		case 10:				/* reaction, not used, pr.use controls */
			pr.irrev = get_true_false(next_char, TRUE);
			break;
		case 11:				/* mix, not used, pr.use controls */
			pr.mix = get_true_false(next_char, TRUE);
			break;
		case 12:				/* use */
		case 18:				/* other */
			pr.use = get_true_false(next_char, TRUE);
			break;
		case 13:				/* selected_output */
			pr.punch = get_true_false(next_char, TRUE);
			phrq_io->Set_punch_on(pr.punch == TRUE);
			break;
		case 19:				/* status */
			{
				int j;
				pr.status = get_true_false(next_char, TRUE);
				j = copy_token(token, &next_char, &l);
				if  (j == UPPER || j == LOWER)
				{
					j = copy_token(token, &next_char, &l);
				}
				if (j == DIGIT)
				{
					char * tptr = token;
					get_num(&tptr, &num);
					status_interval = (int) floor(num);
				}
			}
			if (status_interval < 0)
				status_interval = 0;
			break;
		case 20:				/* inverse */
		case 27:				/* inverse_modeling */
			pr.inverse = get_true_false(next_char, TRUE);
			break;
		case 21:				/* kinetics */
			pr.kinetics = get_true_false(next_char, TRUE);
			break;
		case 22:				/* dump */
			pr.dump = get_true_false(next_char, TRUE);
			phrq_io->Set_dump_on(pr.dump == TRUE);
			break;
		case 23:				/* user_print */
		case 24:				/* user_pr */
			pr.user_print = get_true_false(next_char, TRUE);
			break;
		case 25:				/* solid_solution */
		case 26:				/* solid_solutions */
			pr.ss_assemblage = get_true_false(next_char, TRUE);
			break;
		case 28:				/* headings */
		case 29:				/* heading */
			pr.headings = get_true_false(next_char, TRUE);
			break;
		case 30:				/* user_graph */
			pr.user_graph = get_true_false(next_char, TRUE);
			break;
		case 31:				/* echo_input */
			pr.echo_input = get_true_false(next_char, TRUE);
			break;
		case 32:				/* warning */
		case 33:				/* warnings */
			sscanf(next_char, "%d", &pr.warnings);
			break;
		case 34:				/* initial_isotopes */
			pr.initial_isotopes = get_true_false(next_char, TRUE);
			break;
		case 35:				/* isotope_ratios */
			pr.isotope_ratios = get_true_false(next_char, TRUE);
			break;
		case 36:				/* isotope_alphas */
			pr.isotope_alphas = get_true_false(next_char, TRUE);
			break;
		case 37:				/* censor_species */
			if (copy_token(token, &next_char, &l) != EMPTY)
			{
				sscanf(token, SCANFORMAT, &censor);
			}
			else
			{
				censor = 0;
			}
			break;
		case 38:				/* alkalinity */
			pr.alkalinity = get_true_false(next_char, TRUE);
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
/*
 *   Utilitity routines for read
 */
/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_key(const char *str)
/* ---------------------------------------------------------------------- */
{
/*
 *   Check if string begins with a keyword, returns TRUE or FALSE
 *
 *   Arguments:
 *      *str is pointer to character string to be checked for keywords
 *   Returns:
 *      TRUE,
 *      FALSE.
 */
	char *ptr;
	std::string stdtoken;
	char * token1;
	token1 = string_duplicate(str);

	ptr = token1;
	int j = copy_token(stdtoken, &ptr);
	Utilities::str_tolower(stdtoken);
	std::string key(stdtoken);

	if (j == EMPTY)
	{
		next_keyword = Keywords::KEY_END;
	}
	else
	{
		next_keyword = Keywords::Keyword_search(key);
	}

	free_check_null(token1);
	if (next_keyword > 0)
	{
		return TRUE;
	}
	return (FALSE);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
check_units(char *tot_units, int alkalinity, int check_compatibility,
			const char *default_units, int print)
/* ---------------------------------------------------------------------- */
{
#define NUNITS (sizeof(units) / sizeof(char *))
/*
 *   Check if legitimate units
 *   Input:
 *	   tot_units	   character string to check,
 *	   alkalinity	  TRUE if alkalinity, FALSE if any other total,
 *	   check_compatibility TRUE check alk and default units, FALSE otherwise
 *	   default_units       character string of default units (check /L, /kg, etc)
 *	   print	       TRUE print warning messages
 *   Output:
 *	   tot_units	   standard form for unit
 */
	int i, found;
	char *end;
	char string[MAX_LENGTH];
	const char *units[] = {
		"Mol/l",				/* 0 */
		"mMol/l",				/* 1 */
		"uMol/l",				/* 2 */
		"g/l",					/* 3 */
		"mg/l",					/* 4 */
		"ug/l",					/* 5 */
		"Mol/kgs",				/* 6 */
		"mMol/kgs",				/* 7 */
		"uMol/kgs",				/* 8 */
		"g/kgs",				/* 9 = ppt */
		"mg/kgs",				/* 10 = ppm */
		"ug/kgs",				/* 11 = ppb */
		"Mol/kgw",				/* 12 = mol/kg H2O */
		"mMol/kgw",				/* 13 = mmol/kg H2O */
		"uMol/kgw",				/* 14 = umol/kg H2O */
		"g/kgw",				/* 15 = mol/kg H2O */
		"mg/kgw",				/* 16 = mmol/kg H2O */
		"ug/kgw",				/* 17 = umol/kg H2O */
		"eq/l",					/* 18 */
		"meq/l",				/* 19 */
		"ueq/l",				/* 20 */
		"eq/kgs",				/* 21 */
		"meq/kgs",				/* 22 */
		"ueq/kgs",				/* 23 */
		"eq/kgw",				/* 24 */
		"meq/kgw",				/* 25 */
		"ueq/kgw",				/* 26 */
	};

	squeeze_white(tot_units);
	str_tolower(tot_units);
	replace("milli", "m", tot_units);
	replace("micro", "u", tot_units);
	replace("grams", "g", tot_units);
	replace("gram", "g", tot_units);
	replace("moles", "Mol", tot_units);
	replace("mole", "Mol", tot_units);
	replace("mol", "Mol", tot_units);
	replace("liter", "l", tot_units);
	replace("kgh", "kgw", tot_units);
	replace("ppt", "g/kgs", tot_units);
	replace("ppm", "mg/kgs", tot_units);
	replace("ppb", "ug/kgs", tot_units);
	replace("equivalents", "eq", tot_units);
	replace("equivalent", "eq", tot_units);
	replace("equiv", "eq", tot_units);

	if ((end = strstr(tot_units, "/l")) != NULL)
	{
		*(end + 2) = '\0';
	}
	if ((end = strstr(tot_units, "/kgs")) != NULL)
	{
		*(end + 4) = '\0';
	}
	if ((end = strstr(tot_units, "/kgw")) != NULL)
	{
		*(end + 4) = '\0';
	}
/*
 *   Check if unit in list
 */
	found = FALSE;
	for (i = 0; i < (int) NUNITS; i++)
	{
		if (strcmp(tot_units, units[i]) == 0)
		{
			found = TRUE;
			break;
		}
	}
	if (found == FALSE)
	{
		if (print == TRUE)
		{
			error_string = sformatf( "Unknown unit, %s.", tot_units);
			error_msg(error_string, CONTINUE);
		}
		return (ERROR);
	}

/*
 *   Check if units are compatible with default_units
 */
	if (check_compatibility == FALSE)
		return (OK);
/*
 *   Special cases for alkalinity
 */
	if (alkalinity == TRUE && strstr(tot_units, "Mol") != NULL)
	{
		if (print == TRUE)
		{
			error_string = sformatf(
					"Alkalinity given in moles, assumed to be equivalents.");
			warning_msg(error_string);
		}
		replace("Mol", "eq", tot_units);
	}
	if (alkalinity == FALSE && strstr(tot_units, "eq") != NULL)
	{
		if (print == TRUE)
		{
			error_msg("Only alkalinity can be entered in equivalents.",
					  CONTINUE);
		}
		return (ERROR);
	}
/*
 *   See if default_units are compatible with tot_units
 */
	if (strstr(default_units, "/l") && strstr(tot_units, "/l"))
		return (OK);
	if (strstr(default_units, "/kgs") && strstr(tot_units, "/kgs"))
		return (OK);
	if (strstr(default_units, "/kgw") && strstr(tot_units, "/kgw"))
		return (OK);

	strcpy(string, default_units);
	replace("kgs", "kg solution", string);
	replace("kgs", "kg solution", tot_units);
	replace("kgw", "kg water", string);
	replace("kgw", "kg water", tot_units);
	replace("/l", "/L", string);
	replace("Mol", "mol", string);
	replace("/l", "/L", tot_units);
	replace("Mol", "mol", tot_units);
	if (print == TRUE)
	{
		error_string = sformatf(
				"Units for master species, %s, are not compatible with default units, %s.",
				tot_units, string);
		error_msg(error_string, CONTINUE);
	}
	return (ERROR);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
find_option(const char *item, int *n, const char **list, int count_list, int exact)
/* ---------------------------------------------------------------------- */
{
/*
 *   Compares a string value to match beginning letters of a list of options
 *
 *      Arguments:
 *	 item    entry: pointer to string to compare
 *	 n       exit:  item in list that was matched
 *	 list    entry: pointer to list of character values, assumed to
 *		 be lower case
 *	 count_list entry: number of character values in list
 *
 *      Returns:
 *	 OK      item matched
 *	 ERROR   item not matched
 *	 n       -1      item not matched
 *		 i       position of match in list
 */
	int i;
	std::string stdtoken(item);

	Utilities::str_tolower(stdtoken);
	for (i = 0; i < count_list; i++)
	{
		if (exact == TRUE)
		{
			if (strcmp(list[i], stdtoken.c_str()) == 0)
			{
				*n = i;
				return (OK);
			}
		}
		else
		{
			if (strstr(list[i], stdtoken.c_str()) == list[i])
			{
				*n = i;
				return (OK);
			}
		}
	}
	*n = -1;
	return (ERROR);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_true_false(char *string, int default_value)
/* ---------------------------------------------------------------------- */
{
/*
 *   Returns true unless string starts with "F" or "f"
 */
	int l;
	char token[MAX_LENGTH];
	char *ptr;

	ptr = string;

	if (copy_token(token, &ptr, &l) == EMPTY)
	{
		return (default_value);
	}
	else
	{
		if (token[0] == 'F' || token[0] == 'f')
		{
			return (FALSE);
		}
	}
	return (TRUE);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
get_option(const char **opt_list, int count_opt_list, char **next_char)
/* ---------------------------------------------------------------------- */
{
/*
 *   Read a line and check for options
 */
	int j;
	int opt;
	char *opt_ptr;
	std::string stdoption;
/*
 *   Read line
 */
	j = check_line("get_option", FALSE, TRUE, TRUE, FALSE);
	if (j == EOF)
	{
		j = OPTION_EOF;
	}
	else if (j == KEYWORD)
	{
		j = OPTION_KEYWORD;
	}
	else if (j == OPTION)
	{
		opt_ptr = line;
		copy_token(stdoption, &opt_ptr);
		if (find_option(&(stdoption.c_str()[1]), &opt, opt_list, count_opt_list, FALSE) == OK)
		{
			j = opt;
			replace(stdoption.c_str(), opt_list[j], line_save);
			replace(stdoption.c_str(), opt_list[j], line);
			opt_ptr = line;
			copy_token(stdoption, &opt_ptr);
			*next_char = opt_ptr;
			if (pr.echo_input == TRUE)
			{
				if (!reading_database())
					output_msg(sformatf( "\t%s\n", line_save));
			}
		}
		else
		{
			if (!reading_database())
				output_msg(sformatf( "\t%s\n", line_save));
			error_msg("Unknown option.", CONTINUE);
			error_msg(line_save, CONTINUE);
			input_error++;
			j = OPTION_ERROR;
			*next_char = line;
		}
	}
	else
	{
		opt_ptr = line;
		copy_token(stdoption, &opt_ptr);
		if (find_option(&(stdoption.c_str()[0]), &opt, opt_list, count_opt_list, TRUE) == OK)
		{
			j = opt;
			*next_char = opt_ptr;
		}
		else
		{
			j = OPTION_DEFAULT;
			*next_char = line;
		}
		if (pr.echo_input == TRUE)
		{
			if (!reading_database())
				output_msg(sformatf( "\t%s\n", line_save));
		}
	}
	return (j);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_rates(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate rates
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	char *ptr;
	int l, length, line_length, n;
	int return_value, opt, opt_save;
	char token[MAX_LENGTH];
	struct rate *rate_ptr;
	char *description;
	int n_user, n_user_end;
	char *next_char;
	const char *opt_list[] = {
		"start",				/* 0 */
		"end"					/* 1 */
	};
	int count_opt_list = 2;
/*
 *   Read advection number (not currently used)
 */
	n = -1;
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);
	description = (char *) free_check_null(description);
	opt_save = OPTION_DEFAULT;
/*
 *   Read lines
 */
	return_value = UNKNOWN;
	rate_ptr = NULL;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in RATES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* start */
			opt_save = OPT_1;
			break;
		case 1:				/* end */
			opt_save = OPTION_DEFAULT;
			break;
		case OPTION_DEFAULT:	/* read rate name */
			ptr = line;
			copy_token(token, &ptr, &l);
			rate_ptr = rate_search(token, &n);
			if (rate_ptr == NULL)
			{
				rates =
					(struct rate *) PHRQ_realloc(rates,
												 (size_t) (count_rates +
														   1) *
												 sizeof(struct rate));
				if (rates == NULL)
					malloc_error();
				rate_ptr = &rates[count_rates];
				count_rates++;
			}
			else
			{
				rate_free(rate_ptr);
			}
			rate_ptr->new_def = TRUE;
			rate_ptr->commands = (char *) PHRQ_malloc(sizeof(char));
			if (rate_ptr->commands == NULL)
				malloc_error();
			rate_ptr->commands[0] = '\0';
			rate_ptr->name = string_hsave(token);
			rate_ptr->linebase = NULL;
			rate_ptr->varbase = NULL;
			rate_ptr->loopbase = NULL;
			opt_save = OPT_1;
			break;
		case OPT_1:			/* read command */
			if (rate_ptr == NULL)
			{
				input_error++;
				error_string = sformatf( "No rate name has been defined.");
				error_msg(error_string, CONTINUE);
				opt_save = OPT_1;
				break;
			}
			length = (int) strlen(rate_ptr->commands);
			line_length = (int) strlen(line);
			rate_ptr->commands =
				(char *) PHRQ_realloc(rate_ptr->commands,
									  (size_t) (length + line_length +
												2) * sizeof(char));
			if (rate_ptr->commands == NULL)
				malloc_error();
			rate_ptr->commands[length] = ';';
			rate_ptr->commands[length + 1] = '\0';
			strcat((rate_ptr->commands), line);
			opt_save = OPT_1;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
/*	output_msg(sformatf( "%s", rates[0].commands));
 */ return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_user_print(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate rates
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int length, line_length;
	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"start",				/* 0 */
		"end"					/* 1 */
	};
	int count_opt_list = 2;

	opt_save = OPTION_DEFAULT;
/*
 *   Read lines
 */
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in USER_PRINT keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* start */
			opt_save = OPTION_DEFAULT;
			break;
		case 1:				/* end */
			opt_save = OPTION_DEFAULT;
			break;
		case OPTION_DEFAULT:	/* read first command */
			rate_free(user_print);
			user_print->new_def = TRUE;
			user_print->commands = (char *) PHRQ_malloc(sizeof(char));
			if (user_print->commands == NULL)
				malloc_error();
			user_print->commands[0] = '\0';
			user_print->linebase = NULL;
			user_print->varbase = NULL;
			user_print->loopbase = NULL;
			user_print->name =
				string_hsave("user defined Basic print routine");
		case OPT_1:			/* read command */
			length = (int) strlen(user_print->commands);
			line_length = (int) strlen(line);
			user_print->commands =
				(char *) PHRQ_realloc(user_print->commands,
									  (size_t) (length + line_length +
												2) * sizeof(char));
			if (user_print->commands == NULL)
				malloc_error();
			user_print->commands[length] = ';';
			user_print->commands[length + 1] = '\0';
			strcat((user_print->commands), line);
			opt_save = OPT_1;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
/*	output_msg(sformatf( "%s", rates[0].commands));
 */ return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_user_punch(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate rates
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int length, line_length;
	int return_value, opt, opt_save;
	std::string stdtoken;
	char *next_char;
	const char *opt_list[] = {
		"start",				/* 0 */
		"end",					/* 1 */
		"heading",				/* 2 */
		"headings"				/* 3 */
	};
	int count_opt_list = 4;

	opt_save = OPTION_DEFAULT;
/*
 *   Read lines
 */
	user_punch_count_headings = 0;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in USER_PUNCH keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* start */
			opt_save = OPTION_DEFAULT;
			break;
		case 1:				/* end */
			opt_save = OPTION_DEFAULT;
			break;
		case 2:				/* headings */
		case 3:				/* heading */
			while (copy_token(stdtoken, &next_char) != EMPTY)
			{
				user_punch_headings =
					(const char **) PHRQ_realloc(user_punch_headings,
										   (size_t)
										   (user_punch_count_headings +
											1) * sizeof(char *));
				if (user_punch_headings == NULL)
					malloc_error();
				user_punch_headings[user_punch_count_headings] =
					string_hsave(stdtoken.c_str());
				user_punch_count_headings++;
			}
			break;
		case OPTION_DEFAULT:	/* read first command */
			rate_free(user_punch);
			user_punch->new_def = TRUE;
			user_punch->commands = (char *) PHRQ_malloc(sizeof(char));
			if (user_punch->commands == NULL)
				malloc_error();
			user_punch->commands[0] = '\0';
			user_punch->linebase = NULL;
			user_punch->varbase = NULL;
			user_punch->loopbase = NULL;
			user_punch->name =
				string_hsave("user defined Basic punch routine");
		case OPT_1:			/* read command */
			length = (int) strlen(user_punch->commands);
			line_length = (int) strlen(line);
			user_punch->commands =
				(char *) PHRQ_realloc(user_punch->commands,
									  (size_t) (length + line_length +
												2) * sizeof(char));
			if (user_punch->commands == NULL)
				malloc_error();
			user_punch->commands[length] = ';';
			user_punch->commands[length + 1] = '\0';
			strcat((user_punch->commands), line);
			opt_save = OPT_1;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

#if defined PHREEQ98 
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_user_graph(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate rates
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int l, length, line_length, CurveInfonr = 0;
	int return_value, opt, opt_save;
	char file_name[MAX_LENGTH];
	char token[MAX_LENGTH];
	char *next_char;
	const char *opt_list[] = {
		"start",				/* 0 */
		"end",					/* 1 */
		"heading",				/* 2 */
		"headings",				/* 3 */
		"chart_title",			/* 4 */
		"axis_titles",			/* 5 */
		"axis_scale",			/* 6 */
		"initial_solutions",	/* 7 */
		"plot_concentration_vs",	/* 8 */
		"shifts_as_points",		/* 9 */
		"grid_offset",			/* 10 */
		"connect_simulations",	/* 11 */
		"plot_csv_file"			/* 12 */
	};
	int count_opt_list = 13;
	int i;

	opt_save = OPTION_DEFAULT;
/*
 *   Read lines
 */
#ifdef PHREEQ98
	user_graph_count_headings = 0;
#endif
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT) opt = opt_save;

		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in USER_GRAPH keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* start */
			opt_save = OPTION_DEFAULT;
			break;
		case 1:				/* end */
			opt_save = OPTION_DEFAULT;
			break;
		case 2:				/* headings */
		case 3:				/* heading */
			while (copy_token(token, &next_char, &l) != EMPTY)
			{
				user_graph_headings =
					(char **) PHRQ_realloc(user_graph_headings,
								 (size_t) (user_graph_count_headings +
										   1) * sizeof(char *));
				if (user_graph_headings == NULL)
					malloc_error();
				user_graph_headings[user_graph_count_headings] =
					string_hsave(token);
				user_graph_count_headings++;
			}
			break;
/*Modifications of read_user_punch to change the chart's appearance */
		case 4:				/* chart title */
			copy_title(token, &next_char, &l);
			SetChartTitle(token);
			break;
		case 5:	/* axis titles */
			i = 0;
			while (copy_title(token, &next_char, &l) != EMPTY)
			{
				SetAxisTitles(token, i);
				i++;
			}
			break;
		case 6:	{ /* axis scales */
				char *axis = "";
				int j = 0;
				float f_min, f_max; f_min = (float) -9.999;
				prev_next_char = next_char;
				copy_token(token, &next_char, &l);
				str_tolower(token);
				if (strstr(token, "x") == token)
					axis = "x";
				else if ((strstr(token, "y") == token) && (strstr(token, "y2") != token))
					axis = "y";
				else if ((strstr(token, "s") == token) || (strstr(token, "y2") == token))
					axis = "s";
				else
				{
					input_error++;
					copy_token(token, &prev_next_char, &l);
					error_string = sformatf(
						"Found '%s', but expect axis type \'x\', \'y\', or \'sy\'.", token);
					error_msg(error_string, CONTINUE);
				}
				while ((j < 4)
				   && (i = copy_token(token, &next_char, &l)) != EMPTY)
				{
					str_tolower(token);
					if ((i == DIGIT) || (strstr(token, "a") == token))
						SetAxisScale(axis, j, token, FALSE);
					else
					{
						input_error++;
						copy_token(token, &prev_next_char, &l);
						error_string = sformatf(
						"Found '%s', but expect number or 'a(uto)'.", token);
						error_msg(error_string, CONTINUE);
					}
					if (i == DIGIT)
					{
						if (j == 0)
							f_min = (float) atof(token);
						else if (j == 1 && fabs(f_min + 9.999) > 1e-3)
						{
							f_max = (float) atof(token);
							if (f_min > f_max)
							{
								error_string = sformatf(
								"Maximum must be larger than minimum of axis_scale, interchanging both for %s axis", axis);
								warning_msg(error_string);
								SetAxisScale(axis, 0, token, FALSE);
								sprintf(token, "%e", f_min);
								SetAxisScale(axis, 1, token, FALSE);
							}
						}
					}
					j++;		/* counter for categories */
					prev_next_char = next_char;
				}
				if (j == 4)
					SetAxisScale(axis, j, 0,
							 get_true_false(next_char, FALSE)); /* anything else than false turns on log scale */
			}
			break;
		case 7:
			graph_initial_solutions = get_true_false(next_char, FALSE);
			break;
		case 8:
			prev_next_char = next_char;
			copy_token(token, &next_char, &l);
			str_tolower(token);
			if (strstr(token, "x") == token || strstr(token, "d") == token)
				chart_type = 0;
			else if (strstr(token, "t") == token)
				chart_type = 1;
			else
			{
				input_error++;
				copy_token(token, &prev_next_char, &l);
				error_string = sformatf(
						"Found '%s', but expect plot type: (\'x\' or \'dist\') for distance, (\'t\') for time.",
						token);
					error_msg(error_string, CONTINUE);
			}
			break;
		case 9:
			shifts_as_points = get_true_false(next_char, TRUE);
			if (shifts_as_points == TRUE)
				chart_type = 0;
			else
				chart_type = 1;
			break;
		case 10:
#ifdef PHREEQ98
			i = copy_token(token, &next_char, &l);
			str_tolower(token);
			if (i == DIGIT)
				sscanf(token, "%d", &RowOffset);
			i = copy_token(token, &next_char, &l);
			str_tolower(token);
			if (i == DIGIT)
				sscanf(token, "%d", &ColumnOffset);
#endif
			break;
		case 11:
			connect_simulations = get_true_false(next_char, TRUE);
			break;
		case 12:
			string_trim(next_char);
			strcpy(file_name, next_char);
			if (!OpenCSVFile(file_name))
			{
				error_string = sformatf( "Can't open file, %s. Give the full path + name, or copy the file to the working directory.", file_name);
				input_error++;
				error_msg(error_string, CONTINUE);
			}
			break;
			/* End of modifications */
		case OPTION_DEFAULT:	/* read first command */
			rate_free(user_graph);
			user_graph->new_def = TRUE;
			user_graph->commands = (char *) PHRQ_malloc(sizeof(char));
			if (user_graph->commands == NULL)
				malloc_error();
			user_graph->commands[0] = '\0';
			user_graph->linebase = NULL;
			user_graph->varbase = NULL;
			user_graph->loopbase = NULL;
			user_graph->name =
				string_hsave("user defined Basic graph routine");
		case OPT_1:			/* read command */
			length = strlen(user_graph->commands);
			line_length = strlen(line);
			user_graph->commands =
				(char *) PHRQ_realloc(user_graph->commands,
							 (size_t) (length + line_length +
									   2) * sizeof(char));
			if (user_graph->commands == NULL)
				malloc_error();
			user_graph->commands[length] = ';';
			user_graph->commands[length + 1] = '\0';
			strcat((user_graph->commands), line);
			opt_save = OPT_1;
			break;
		}

		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
#ifdef PHREEQ98
	for (i = 0; i < user_graph_count_headings; i++)
	{
		GridHeadings(user_graph_headings[i], i);
	}
#endif
	return (return_value);
}
#endif

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_solid_solutions(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads solid solution data
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int i, j, n, l;
	int count_s_s, number_s_s, count_comps;
	int n_user, n_user_end;
	char *ptr;
	char *description;
	char token[MAX_LENGTH];

	int return_value, opt;
	char *next_char;
	const char *opt_list[] = {
		"component",			/* 0 */
		"comp",					/* 1 */
		"parms",				/* 2 */
		"gugg_nondimensional",	/* 3 */
		"gugg_kj",				/* 4 */
		"activity_coefficients",	/* 5 */
		"distribution_coefficients",	/* 6 */
		"miscibility_gap",		/* 7 */
		"spinodal_gap",			/* 8 */
		"critical_point",		/* 9 */
		"alyotropic_point",		/* 10 */
		"temp",					/* 11 */
		"tempk",				/* 12 */
		"tempc",				/* 13 */
		"thompson",				/* 14 */
		"margules",				/* 15 */
		"comp1",				/* 16 */
		"comp2"					/* 17 */
	};
	int count_opt_list = 18;
/*
 *   Read ss_assemblage number
 */
	number_s_s = 0;
	ptr = line;
	read_number_description(ptr, &n_user, &n_user_end, &description);
/*
 *   Find old ss_assemblage or alloc space for new ss_assemblage
 */
	if (ss_assemblage_search(n_user, &n) != NULL)
	{
		ss_assemblage_free(&ss_assemblage[n]);
	}
	else
	{
		n = count_ss_assemblage;
		space((void **) ((void *) &ss_assemblage), count_ss_assemblage,
			  &max_ss_assemblage, sizeof(struct ss_assemblage));
		count_ss_assemblage++;
	}
/*
 *   Initialize
 */
	ss_assemblage_init(&(ss_assemblage[n]), n_user, n_user_end,
						description);
	free_check_null(description);
/*
 *   Set use data to first read
 */
	if (!use.Get_ss_assemblage_in())
	{
		use.Set_ss_assemblage_in(true);
		use.Set_n_ss_assemblage_user(n_user);
	}
/*
 *   Read solid solutions
 */
	count_s_s = 0;
	count_comps = 0;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SOLID_SOLUTIONS keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;

/*
 * New component
 */
		case 0:				/* component */
		case 1:				/* comp */
			count_comps = ss_assemblage[n].s_s[number_s_s].count_comps++;
			ss_assemblage[n].s_s[number_s_s].comps =
				(struct s_s_comp *) PHRQ_realloc(ss_assemblage[n].
												 s_s[number_s_s].comps,
												 (size_t) (count_comps +
														   1) *
												 sizeof(struct s_s_comp));
			if (ss_assemblage[n].s_s[number_s_s].comps == NULL)
				malloc_error();
			ss_assemblage[n].s_s[number_s_s].comps[count_comps].initial_moles = 0;
			ss_assemblage[n].s_s[number_s_s].comps[count_comps].delta = 0;
			/*
			 *   Read phase name of component
			 */
			ptr = next_char;
			copy_token(token, &ptr, &l);
			ss_assemblage[n].s_s[number_s_s].comps[count_comps].name =
				string_hsave(token);
			/*
			 *   Read moles of component
			 */
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
			{
				ss_assemblage[n].s_s[number_s_s].comps[count_comps].moles =
					NAN;
			}
			else
			{
				j = sscanf(token, SCANFORMAT, &dummy);
				ss_assemblage[n].s_s[number_s_s].comps[count_comps].moles =
					(LDBLE) dummy;
				if (j != 1)
				{
					error_msg("Expected moles of solid solution.", CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
				}
			}
			break;
		case 2:				/* parms */
		case 3:				/* gugg_nondimensional */
			/*
			 *   Read parameters
			 */

			ptr = next_char;
			if (copy_token(token, &ptr, &l) != EMPTY)
			{
				sscanf(token, SCANFORMAT,
					   &(ss_assemblage[n].s_s[number_s_s].p[0]));
			}
			if (copy_token(token, &ptr, &l) != EMPTY)
			{
				sscanf(token, SCANFORMAT,
					   &(ss_assemblage[n].s_s[number_s_s].p[1]));
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 0;
			break;
		case 4:				/* gugg_kj */
			ptr = next_char;
			if (copy_token(token, &ptr, &l) != EMPTY)
			{
				sscanf(token, SCANFORMAT,
					   &(ss_assemblage[n].s_s[number_s_s].p[0]));
			}
			if (copy_token(token, &ptr, &l) != EMPTY)
			{
				sscanf(token, SCANFORMAT,
					   &(ss_assemblage[n].s_s[number_s_s].p[1]));
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 7;
			break;
		case 5:				/* activity coefficients */
			ptr = next_char;
			j = 0;
			for (i = 0; i < 4; i++)
			{
				if (copy_token(token, &ptr, &l) != EMPTY)
				{
					j += sscanf(token, SCANFORMAT,
								&(ss_assemblage[n].s_s[number_s_s].p[i]));
				}
			}
			if (j != 4)
			{
				error_string = sformatf(
						"Expected 4 parameters to calculate a0 and a1 from two activity coefficients, assemblage %d, solid solution %s",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 1;
			break;
		case 6:				/* distribution coefficients */
			ptr = next_char;
			j = 0;
			for (i = 0; i < 4; i++)
			{
				if (copy_token(token, &ptr, &l) != EMPTY)
				{
					j += sscanf(token, SCANFORMAT,
								&(ss_assemblage[n].s_s[number_s_s].p[i]));
				}
			}
			if (j != 4)
			{
				error_string = sformatf(
						"Expected 4 parameters to calculate a0 and a1 from two distribution coefficients, assemblage %d, solid solution %s",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 2;
			break;
		case 7:				/* miscibility_gap */
			ptr = next_char;
			j = 0;
			for (i = 0; i < 2; i++)
			{
				if (copy_token(token, &ptr, &l) != EMPTY)
				{
					j += sscanf(token, SCANFORMAT,
								&(ss_assemblage[n].s_s[number_s_s].p[i]));
				}
			}
			if (j != 2)
			{
				error_string = sformatf(
						"Expected 2 miscibility gap fractions of component 2 to calculate a0 and a1, assemblage %d, solid solution %s",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 3;
			break;
		case 8:				/* spinodal_gap */
			ptr = next_char;
			j = 0;
			for (i = 0; i < 2; i++)
			{
				if (copy_token(token, &ptr, &l) != EMPTY)
				{
					j += sscanf(token, SCANFORMAT,
								&(ss_assemblage[n].s_s[number_s_s].p[i]));
				}
			}
			if (j != 2)
			{
				error_string = sformatf(
						"Expected 2 spinodal gap fractions of component 2 to calculate a0 and a1, assemblage %d, solid solution %s",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 4;
			break;
		case 9:				/* critical point */
			ptr = next_char;
			j = 0;
			for (i = 0; i < 2; i++)
			{
				if (copy_token(token, &ptr, &l) != EMPTY)
				{
					j += sscanf(token, SCANFORMAT,
								&(ss_assemblage[n].s_s[number_s_s].p[i]));
				}
			}
			if (j != 2)
			{
				error_string = sformatf(
						"Expected fraction of component 2 and critical temperature to calculate a0 and a1, assemblage %d, solid solution %s",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 5;
			break;
		case 10:				/* alyotropic point */
			ptr = next_char;
			j = 0;
			for (i = 0; i < 2; i++)
			{
				if (copy_token(token, &ptr, &l) != EMPTY)
				{
					j += sscanf(token, SCANFORMAT,
								&(ss_assemblage[n].s_s[number_s_s].p[i]));
				}
			}
			if (j != 2)
			{
				error_string = sformatf(
						"Expected fraction of component 2 and sigma pi at alyotropic point to calculate a0 and a1, assemblage %d, solid solution %s",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 6;
			break;
		case 12:				/* tempk */
			ptr = next_char;
			j = 0;
			if (copy_token(token, &ptr, &l) != EMPTY)
			{
				j = sscanf(token, SCANFORMAT,
						   &(ss_assemblage[n].s_s[number_s_s].tk));
			}
			if (j != 1)
			{
				error_string = sformatf(
						"Expected temperature (Kelvin) for parameters, assemblage %d, solid solution %s, using 298.15 K",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				warning_msg(error_string);
				ss_assemblage[n].s_s[number_s_s].tk = 298.15;
			}
			break;
		case 11:				/* temp */
		case 13:				/* tempc */
			ptr = next_char;
			j = 0;
			if (copy_token(token, &ptr, &l) != EMPTY)
			{
				j = sscanf(token, SCANFORMAT,
						   &(ss_assemblage[n].s_s[number_s_s].tk));
			}
			if (j != 1)
			{
				error_string = sformatf(
						"Expected temperature (Celcius) for parameters, assemblage %d, solid solution %s, using 25 C",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				warning_msg(error_string);
				ss_assemblage[n].s_s[number_s_s].tk = 25.;
			}
			ss_assemblage[n].s_s[number_s_s].tk += 273.15;
			break;
		case 14:				/* Thompson and Waldbaum */
			ptr = next_char;
			j = 0;
			for (i = 0; i < 2; i++)
			{
				if (copy_token(token, &ptr, &l) != EMPTY)
				{
					j += sscanf(token, SCANFORMAT,
								&(ss_assemblage[n].s_s[number_s_s].p[i]));
				}
			}
			if (j != 2)
			{
				error_string = sformatf(
						"Expected Wg2 and Wg1 Thompson-Waldbaum parameters to calculate a0 and a1, assemblage %d, solid solution %s",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 8;
			break;
		case 15:				/* Margules */
			ptr = next_char;
			j = 0;
			for (i = 0; i < 2; i++)
			{
				if (copy_token(token, &ptr, &l) != EMPTY)
				{
					j += sscanf(token, SCANFORMAT,
								&(ss_assemblage[n].s_s[number_s_s].p[i]));
				}
			}
			if (j != 2)
			{
				error_string = sformatf(
						"Expected alpha2 and alpha3 Margules parameters to calculate a0 and a1, assemblage %d, solid solution %s",
						ss_assemblage[n].n_user,
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
			ss_assemblage[n].s_s[number_s_s].input_case = 9;
			break;
		case 16:				/* comp1 */
			if (count_comps < 2)
			{
				ss_assemblage[n].s_s[number_s_s].count_comps = 2;
				count_comps = 2;
				ss_assemblage[n].s_s[number_s_s].comps =
					(struct s_s_comp *) PHRQ_realloc(ss_assemblage[n].
													 s_s[number_s_s].comps,
													 (size_t) (count_comps) *
													 sizeof(struct s_s_comp));
				if (ss_assemblage[n].s_s[number_s_s].comps == NULL)
					malloc_error();
			}
			/*
			 *   Read phase name of component
			 */
			ptr = next_char;
			copy_token(token, &ptr, &l);
			ss_assemblage[n].s_s[number_s_s].comps[0].name =
				string_hsave(token);
			/*
			 *   Read moles of component
			 */
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
			{
				ss_assemblage[n].s_s[number_s_s].comps[0].moles = NAN;
			}
			else
			{
				j = sscanf(token, SCANFORMAT, &dummy);
				ss_assemblage[n].s_s[number_s_s].comps[0].moles =
					(LDBLE) dummy;
				if (j != 1)
				{
					error_msg("Expected moles of solid solution.", CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
				}
			}
			break;
		case 17:				/* comp2 */
			if (count_comps < 2)
			{
				ss_assemblage[n].s_s[number_s_s].count_comps = 2;
				count_comps = 2;
				ss_assemblage[n].s_s[number_s_s].comps =
					(struct s_s_comp *) PHRQ_realloc(ss_assemblage[n].
													 s_s[number_s_s].comps,
													 (size_t) (count_comps) *
													 sizeof(struct s_s_comp));
				if (ss_assemblage[n].s_s[number_s_s].comps == NULL)
					malloc_error();
			}
			/*
			 *   Read phase name of component
			 */
			ptr = next_char;
			copy_token(token, &ptr, &l);
			ss_assemblage[n].s_s[number_s_s].comps[1].name =
				string_hsave(token);
			/*
			 *   Read moles of component
			 */
			if ((j = copy_token(token, &ptr, &l)) == EMPTY)
			{
				ss_assemblage[n].s_s[number_s_s].comps[1].moles = NAN;
			}
			else
			{
				j = sscanf(token, SCANFORMAT, &dummy);
				ss_assemblage[n].s_s[number_s_s].comps[1].moles =
					(LDBLE) dummy;
				if (j != 1)
				{
					error_msg("Expected moles of solid solution.", CONTINUE);
					error_msg(line_save, CONTINUE);
					input_error++;
				}
			}
			break;
/*
 * New solid solution
 */
		case OPTION_DEFAULT:
			number_s_s = count_s_s++;
			/*
			 *   Make space, set default
			 */

			/* realloc space for one s_s */
			ss_assemblage[n].s_s =
				(struct s_s *) PHRQ_realloc(ss_assemblage[n].s_s,
											(size_t) (count_s_s +
													  1) *
											sizeof(struct s_s));
			if (ss_assemblage[n].s_s == NULL)
				malloc_error();

			/* malloc space for one component */
			ss_assemblage[n].s_s[number_s_s].comps =
				(struct s_s_comp *) PHRQ_malloc((size_t)
												sizeof(struct s_s_comp));
			if (ss_assemblage[n].s_s[number_s_s].comps == NULL)
				malloc_error();
			count_comps = 0;
			ss_assemblage[n].s_s[number_s_s].count_comps = 0;
			ss_assemblage[n].s_s[number_s_s].total_moles = 0;
			ss_assemblage[n].s_s[number_s_s].dn = 0;
			ss_assemblage[n].s_s[number_s_s].a0 = 0.0;
			ss_assemblage[n].s_s[number_s_s].a1 = 0.0;
			ss_assemblage[n].s_s[number_s_s].ag0 = 0.0;
			ss_assemblage[n].s_s[number_s_s].ag1 = 0.0;
			ss_assemblage[n].s_s[number_s_s].s_s_in = FALSE;
			ss_assemblage[n].s_s[number_s_s].miscibility = FALSE;
			ss_assemblage[n].s_s[number_s_s].spinodal = FALSE;
			ss_assemblage[n].s_s[number_s_s].tk = 298.15;
			ss_assemblage[n].s_s[number_s_s].xb1 = 0;
			ss_assemblage[n].s_s[number_s_s].xb2 = 0;
			ss_assemblage[n].s_s[number_s_s].input_case = 0;
			ss_assemblage[n].s_s[number_s_s].p[0] = 0.0;
			ss_assemblage[n].s_s[number_s_s].p[1] = 0.0;
			ss_assemblage[n].s_s[number_s_s].p[2] = 0.0;
			ss_assemblage[n].s_s[number_s_s].p[3] = 0.0;

			ss_assemblage[n].s_s[number_s_s].comps->name = NULL;
			ss_assemblage[n].s_s[number_s_s].comps->phase = NULL;
			/*
			 *   Read solid solution name
			 */
			ptr = line;
			copy_token(token, &ptr, &l);
			ss_assemblage[n].s_s[number_s_s].name = string_hsave(token);
			ss_assemblage[n].s_s[number_s_s].total_moles = NAN;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	for (i = 0; i < ss_assemblage[n].count_s_s; i++)
	{
		if (ss_assemblage[n].s_s[i].p[0] != 0.0 ||
			ss_assemblage[n].s_s[i].p[1] != 0.0)
		{
			if (ss_assemblage[n].s_s[number_s_s].count_comps != 2)
			{
				error_string = sformatf(
						"Solid solution, %s, is nonideal. Must define exactly two components (-comp1 and -comp2).",
						ss_assemblage[n].s_s[number_s_s].name);
				error_msg(error_string, CONTINUE);
				input_error++;
			}
		}
	}
/*
 *   Sort components by name (lowercase)
 */

	ss_assemblage[n].count_s_s = count_s_s;
	qsort(ss_assemblage[n].s_s,
		  (size_t) count_s_s, (size_t) sizeof(struct s_s), s_s_compare);
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_llnl_aqueous_model_parameters(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads aqueous model parameters
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */
	int i, count_alloc;
	char token[MAX_LENGTH];

	int return_value, opt;
	char *next_char;
	const char *opt_list[] = {
		"temperatures",			/* 0 */
		"temperature",			/* 1 */
		"temp",					/* 2 */
		"adh",					/* 3 */
		"debye_huckel_a",		/* 4 */
		"dh_a",					/* 5 */
		"bdh",					/* 6 */
		"debye_huckel_b",		/* 7 */
		"dh_b",					/* 8 */
		"bdot",					/* 9 */
		"b_dot",				/* 10 */
		"c_co2",				/* 11 */
		"co2_coefs"				/* 12 */
	};
	int count_opt_list = 13;
/*
 *   Initialize
 */
/*
 *   Read aqueous model parameters
 */
	return_value = UNKNOWN;
	opt = get_option(opt_list, count_opt_list, &next_char);
	for (;;)
	{
		next_char = line;
		if (opt >= 0)
		{
			copy_token(token, &next_char, &i);
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_DEFAULT:
		case OPTION_ERROR:
			input_error++;
			error_msg
				("Unknown input in LLNL_AQUEOUS_MODEL_PARAMETERS keyword.",
				 CONTINUE);
			error_msg(line_save, CONTINUE);
			break;

/*
 * New component
 */
		case 0:				/* temperatures */
		case 1:				/* temperature */
		case 2:				/* temp */
			count_alloc = 1;
			llnl_count_temp = 0;
			i = read_lines_doubles(next_char, &(llnl_temp),
								   &(llnl_count_temp), &(count_alloc),
								   opt_list, count_opt_list, &opt);
			/*
			   ptr = next_char;
			   llnl_temp = read_list_doubles(&ptr, &count);
			   llnl_count_temp = count;
			 */
			break;
		case 3:				/* adh */
		case 4:				/* debye_huckel_a */
		case 5:				/* dh_a */
			count_alloc = 1;
			llnl_count_adh = 0;
			i = read_lines_doubles(next_char, &(llnl_adh), &(llnl_count_adh),
								   &(count_alloc), opt_list, count_opt_list,
								   &opt);
			/*
			   ptr = next_char;
			   llnl_adh = read_list_doubles(&ptr, &count);
			   llnl_count_adh = count;
			 */
			break;
		case 6:				/* bdh */
		case 7:				/* debye_huckel_b */
		case 8:				/* dh_b */
			count_alloc = 1;
			llnl_count_bdh = 0;
			i = read_lines_doubles(next_char, &(llnl_bdh), &(llnl_count_bdh),
								   &(count_alloc), opt_list, count_opt_list,
								   &opt);
			/*
			   ptr = next_char;
			   llnl_bdh = read_list_doubles(&ptr, &count);
			   llnl_count_bdh = count;
			 */
			break;
		case 9:				/* bdot */
		case 10:				/* b_dot */
			count_alloc = 1;
			llnl_count_bdot = 0;
			i = read_lines_doubles(next_char, &(llnl_bdot),
								   &(llnl_count_bdot), &(count_alloc),
								   opt_list, count_opt_list, &opt);
			/*
			   ptr = next_char;
			   llnl_bdot = read_list_doubles(&ptr, &count);
			   llnl_count_bdot = count;
			 */
			break;
		case 11:				/* c_co2 */
		case 12:				/* co2_coefs */
			count_alloc = 1;
			llnl_count_co2_coefs = 0;
			i = read_lines_doubles(next_char, &(llnl_co2_coefs),
								   &(llnl_count_co2_coefs), &(count_alloc),
								   opt_list, count_opt_list, &opt);
			/*
			   ptr = next_char;
			   llnl_co2_coefs = read_list_doubles(&ptr, &count);
			   llnl_count_co2_coefs = count;
			 */
			break;
		}
		return_value = check_line_return;
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	/* check consistency */
	if ((llnl_count_temp <= 0) ||
		(llnl_count_temp != llnl_count_adh) ||
		(llnl_count_temp != llnl_count_bdh) ||
		(llnl_count_temp != llnl_count_bdot))
	{
		error_msg
			("Must define equal number (>0) of temperatures, dh_a, dh_b, and bdot parameters\nin LLNL_AQUEOUS_MODEL",
			 CONTINUE);
		input_error++;
	}
	if (llnl_count_co2_coefs != 5)
	{
		error_msg
			("Must define 5 CO2 activity coefficient parameters in LLNL_AQUEOUS_MODEL",
			 CONTINUE);
		input_error++;
	}
	for (i = 1; i < llnl_count_temp; i++)
	{
		if (llnl_temp[i - 1] > llnl_temp[i])
		{
			error_msg
				("Temperatures must be in ascending order in LLNL_AQUEOUS_MODEL",
				 CONTINUE);
			input_error++;
		}
	}

	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_lines_doubles(char *next_char, LDBLE ** d, int *count_d,
				   int *count_alloc, const char **opt_list,
				   int count_opt_list, int *opt)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads LDBLEs on line starting at next_char
 *      and on succeeding lines. Appends to d.
 *      Stops at KEYWORD, OPTION, and EOF
 *
 *      Input Arguments:
 *	 next_char    points to line to read from
 *	 d	    points to array of LDBLEs, must be malloced
 *	 count_d      number of elements in array
 *	 count_alloc  number of elements malloced
 *
 *      Output Arguments:
 *	 d	    points to array of LDBLEs, may have been
 *			  realloced
 *	 count_d      updated number of elements in array
 *	 count_alloc  updated of elements malloced
 *
 *      Returns:
 *	 KEYWORD
 *	 OPTION
 *	 EOF
 *	 ERROR if any errors reading LDBLEs
 */

	if (read_line_doubles(next_char, d, count_d, count_alloc) == ERROR)
	{
		return (ERROR);
	}
	for (;;)
	{
		*opt = get_option(opt_list, count_opt_list, &next_char);
		if (*opt == OPTION_KEYWORD || *opt == OPTION_EOF
			|| *opt == OPTION_ERROR)
		{
			break;
		}
		else if (*opt >= 0)
		{
			break;
		}
		next_char = line;
		if (read_line_doubles(next_char, d, count_d, count_alloc) == ERROR)
		{
			return (ERROR);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_line_doubles(char *next_char, LDBLE ** d, int *count_d, int *count_alloc)
/* ---------------------------------------------------------------------- */
{
	int i, j, l, n;
	LDBLE value;
	char token[MAX_LENGTH];

	for (;;)
	{
		j = copy_token(token, &next_char, &l);
		if (j == EMPTY)
		{
			break;
		}
		if (j != DIGIT)
		{
			return (ERROR);
		}
		if (replace("*", " ", token) == TRUE)
		{
			if (sscanf(token, "%d" SCANFORMAT, &n, &value) != 2)
			{
				return (ERROR);
			}
		}
		else
		{
			sscanf(token, SCANFORMAT, &value);
			n = 1;
		}
		for (;;)
		{
			if ((*count_d) + n > (*count_alloc))
			{
				*count_alloc *= 2;
				*d = (LDBLE *) PHRQ_realloc(*d,
											(size_t) (*count_alloc) *
											sizeof(LDBLE));
				if (*d == NULL)
					malloc_error();
			}
			else
			{
				break;
			}
		}
		for (i = 0; i < n; i++)
		{
			(*d)[(*count_d) + i] = value;
		}
		*count_d += n;
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
next_keyword_or_option(const char **opt_list, int count_opt_list)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads to next keyword or option or eof
 *
 *   Returns:
 *       KEYWORD
 *       OPTION
 *       EOF
 */
	int opt;
	char *next_char;

	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_EOF)
		{						/* end of file */
			break;
		}
		else if (opt == OPTION_KEYWORD)
		{						/* keyword */
			break;
		}
		else if (opt >= 0 && opt < count_opt_list)
		{
			break;
		}
		else
		{
			error_msg("Expected a keyword or option.", CONTINUE);
			error_msg(line_save, CONTINUE);
			input_error++;
		}
	}
	return (opt);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_named_logk(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads K's that can be used to calculate K's for species
 *
 *      Arguments:
 *	 none
 *
 *      Returns:
 *	 KEYWORD if keyword encountered, input_error may be incremented if
 *		    a keyword is encountered in an unexpected position
 *	 EOF     if eof encountered while reading mass balance concentrations
 *	 ERROR   if error occurred reading data
 *
 */

	int j, l;
	int i, empty;
	struct logk *logk_ptr;
	char token[MAX_LENGTH];

	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"log_k",				/* 0 */
		"logk",					/* 1 */
		"delta_h",				/* 2 */
		"deltah",				/* 3 */
		"analytical_expression",	/* 4 */
		"a_e",					/* 5 */
		"ae",					/* 6 */
		"ln_alpha1000",			/* 7 */
		"add_logk",				/* 8 */
		"add_log_k",			/* 9 */
		"delta_v",				/* 10 */
		"deltav",				/* 11 */
		"vm"
	};
	int count_opt_list = 13;
	logk_ptr = NULL;
/*
 *   Read name followed by options
 */
	opt_save = OPTION_DEFAULT;
	return_value = UNKNOWN;
	for (;;)
	{
		opt = get_option(opt_list, count_opt_list, &next_char);
		if (opt == OPTION_DEFAULT)
		{
			opt = opt_save;
		}
		switch (opt)
		{
		case OPTION_EOF:		/* end of file */
			return_value = EOF;
			break;
		case OPTION_KEYWORD:	/* keyword */
			return_value = KEYWORD;
			break;
		case OPTION_ERROR:
			input_error++;
			error_msg("Unknown input in SPECIES keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* log_k */
		case 1:				/* logk */
			if (logk_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_log_k_only(next_char, &logk_ptr->log_k[0]);
			logk_copy2orig(logk_ptr);
			opt_save = OPTION_DEFAULT;
			break;
		case 2:				/* delta_h */
		case 3:				/* deltah */
			if (logk_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_delta_h_only(next_char, &logk_ptr->log_k[1],
							  &logk_ptr->original_units);
			logk_copy2orig(logk_ptr);
			opt_save = OPTION_DEFAULT;
			break;
		case 4:				/* analytical_expression */
		case 5:				/* a_e */
		case 6:				/* ae */
			if (logk_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_analytical_expression_only(next_char, &(logk_ptr->log_k[T_A1]));
			logk_copy2orig(logk_ptr);
			opt_save = OPTION_DEFAULT;
			break;
		case 7:				/* ln_alpha1000 */
			if (logk_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			empty = TRUE;
			for (i = T_A1; i <= T_A6; i++)
			{
				if (logk_ptr->log_k[i] != 0.0)
				{
					empty = FALSE;
					logk_ptr->log_k[i] = 0.0;
				}
			}
			if (empty == FALSE)
			{
				error_string = sformatf(
						"Analytical expression previously defined for %s in NAMED_EXPRESSIONS\nAnalytical expression will be overwritten.",
						logk_ptr->name);
				warning_msg(error_string);
			}
			read_analytical_expression_only(next_char, &(logk_ptr->log_k[T_A1]));
			for (i = T_A1; i < T_A6; i++)
			{
				logk_ptr->log_k[i] /= 1000. * LOG_10;
			}
			logk_copy2orig(logk_ptr);
			opt_save = OPTION_DEFAULT;
			break;
		case 8:				/* add_logk */
		case 9:				/* add_log_k */
			if (logk_ptr == NULL)
			{
				error_string = sformatf(
						"No reaction defined before option, %s.",
						opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			if (logk_ptr->count_add_logk == 0)
			{
				logk_ptr->add_logk =
					(struct name_coef *)
					PHRQ_malloc(sizeof(struct name_coef));
				if (logk_ptr->add_logk == NULL)
					malloc_error();
			}
			else
			{
				logk_ptr->add_logk =
					(struct name_coef *) PHRQ_realloc(logk_ptr->add_logk,
													  (size_t) ((logk_ptr->
																 count_add_logk
																 +
																 1) *
																sizeof
																(struct
																 name_coef)));
				if (logk_ptr->add_logk == NULL)
					malloc_error();
			}
			/* read name */
			if (copy_token(token, &next_char, &i) == EMPTY)
			{
				input_error++;
				error_string = sformatf(
						"Expected the name of a NAMED_EXPRESSION.");
				error_msg(error_string, CONTINUE);
				break;
			}
			logk_ptr->add_logk[logk_ptr->count_add_logk].name =
				string_hsave(token);
			/* read coef */
			i = sscanf(next_char, SCANFORMAT,
					   &logk_ptr->add_logk[logk_ptr->count_add_logk].coef);
			if (i <= 0)
			{
				logk_ptr->add_logk[logk_ptr->count_add_logk].coef = 1;
			}
			logk_ptr->count_add_logk++;
			opt_save = OPTION_DEFAULT;
			break;
		case 10:            /* delta_v */
		case 11:            /* deltav */
		case 12:            /* vm, molar volume */
			if (logk_ptr == NULL)
			{
				error_string = sformatf(
					"No reaction defined before option, %s.",
				opt_list[opt]);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			read_delta_v_only(next_char, &logk_ptr->log_k[vm0],
				&logk_ptr->original_deltav_units);
			logk_copy2orig(logk_ptr);
			opt_save = OPTION_DEFAULT;
			break;
		case OPTION_DEFAULT:
/*
 *   Get space for logk information
 */
			logk_ptr = NULL;
			j = copy_token(token, &next_char, &l);

			logk_ptr = logk_store(token, TRUE);
/*
 *   Get pointer to each species in the reaction, store new species if necessary
 */
			opt_save = OPTION_DEFAULT;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_copy(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Reads solution,
 *	 equilibrium_phases,
 *	 exchange,
 *	 surface,
 *	 solid_solution,
 *	 gas_phase,
 *	 kinetics,
 *	 mix,
 *	 reaction,
 *	 reaction_temperature
 *
 */
	int i, l, n, n_user, n_user_start, n_user_end, return_value;
	char *ptr;
	char token[MAX_LENGTH], token1[MAX_LENGTH], nonkeyword[MAX_LENGTH];
/*
 *   Read "copy"
 */
	ptr = line;
	copy_token(token, &ptr, &l);
/*
 *   Read keyword
 */
	copy_token(token, &ptr, &l);
	check_key(token);

	switch (next_keyword)
	{
	case Keywords::KEY_NONE:					/* Have not read line with keyword */
		strcpy(nonkeyword, token);
		break;
	case Keywords::KEY_SOLUTION:				/* Solution */
	case Keywords::KEY_EQUILIBRIUM_PHASES:		/* Pure phases */
	case Keywords::KEY_REACTION:				/* Reaction */
	case Keywords::KEY_MIX:						/* Mix */
	case Keywords::KEY_EXCHANGE:				/* Ex */
	case Keywords::KEY_SURFACE:					/* Surface */
	case Keywords::KEY_REACTION_TEMPERATURE:	/* Temperature */
	case Keywords::KEY_REACTION_PRESSURE:		/* Pressure */
	case Keywords::KEY_GAS_PHASE:				/* Gas */
	case Keywords::KEY_KINETICS:				/* Kinetics */
	case Keywords::KEY_SOLID_SOLUTIONS:			/* solid_solutions */
		break;
	default:
		input_error++;
		error_msg
			("Expecting keyword solution, mix, kinetics, reaction, reaction_temperature, equilibrium_phases, exchange, surface, gas_phase, or solid_solutions, or cell.",
			 CONTINUE);
		error_msg(line_save, CONTINUE);
		check_line("End of use", FALSE, TRUE, TRUE, TRUE);
		/* empty, eof, keyword, print */
		return (ERROR);
	}
/*
 *   Read source index
 */
	strcpy(token1, token);
	i = copy_token(token, &ptr, &l);
	if (i == DIGIT)
	{
		sscanf(token, "%d", &n_user);
		if (n_user < 0)
		{
			error_msg("Source index number must be a positive integer.",
					  CONTINUE);
			error_msg(line_save, CONTINUE);
			input_error++;
			return (ERROR);
		}
		if (strstr(token, "-") != NULL)
		{
			error_msg
				("COPY does not accept a range of numbers for source index",
				 CONTINUE);
			error_msg(line_save, CONTINUE);
			input_error++;
			return (ERROR);
		}
	}
	else
	{
		error_msg("Source index number must be a positive integer.",
				  CONTINUE);
		error_msg(line_save, CONTINUE);
		input_error++;
		return (ERROR);
	}
/*
 *   Read target index or range of indices
 */
	i = copy_token(token, &ptr, &l);
	if (i == DIGIT)
	{
		replace("-", " ", &token[1]);
		n = sscanf(token, "%d%d", &n_user_start, &n_user_end);
		if (n == 1)
		{
			n_user_end = n_user_start;
		}
		if (n_user_start < 0)
		{
			//error_msg("Target index number must be a positive integer.",
			//		  CONTINUE);
			//error_msg(line_save, CONTINUE);
			//input_error++;
			//return (ERROR);
		}
	}
	else
	{
		error_msg("Target index number must be a positive integer.",
				  CONTINUE);
		error_msg(line_save, CONTINUE);
		input_error++;
		return (ERROR);
	}

	switch (next_keyword)
	{
	case Keywords::KEY_NONE:
		str_tolower(nonkeyword);
		if (strstr(nonkeyword, "cell") != nonkeyword)
		{
			error_msg("Unknown input in COPY data block.", CONTINUE);
			error_msg(line_save, CONTINUE);
			input_error++;
			return (ERROR);
		}
		copier_add(&copy_solution, n_user, n_user_start, n_user_end);
		copier_add(&copy_pp_assemblage, n_user, n_user_start, n_user_end);
		copier_add(&copy_reaction, n_user, n_user_start, n_user_end);
		copier_add(&copy_mix, n_user, n_user_start, n_user_end);
		copier_add(&copy_exchange, n_user, n_user_start, n_user_end);
		copier_add(&copy_surface, n_user, n_user_start, n_user_end);
		copier_add(&copy_temperature, n_user, n_user_start, n_user_end);
		copier_add(&copy_pressure, n_user, n_user_start, n_user_end);
		copier_add(&copy_gas_phase, n_user, n_user_start, n_user_end);
		copier_add(&copy_kinetics, n_user, n_user_start, n_user_end);
		copier_add(&copy_ss_assemblage, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_SOLUTION:								/* Solution */
		copier_add(&copy_solution, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_EQUILIBRIUM_PHASES:					/* Pure phases */
		copier_add(&copy_pp_assemblage, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_REACTION:								/* Reaction */
		copier_add(&copy_reaction, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_MIX:										/* Mix */
		copier_add(&copy_mix, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_EXCHANGE:								/* Ex */
		copier_add(&copy_exchange, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_SURFACE:									/* Surface */
		copier_add(&copy_surface, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_REACTION_TEMPERATURE:					/* Temperature */
		copier_add(&copy_temperature, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_REACTION_PRESSURE:						/* Pressure */
		copier_add(&copy_pressure, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_GAS_PHASE:								/* Gas */
		copier_add(&copy_gas_phase, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_KINETICS:								/* Kinetics */
		copier_add(&copy_kinetics, n_user, n_user_start, n_user_end);
		break;
	case Keywords::KEY_SOLID_SOLUTIONS:							/* solid_solutions */
		copier_add(&copy_ss_assemblage, n_user, n_user_start, n_user_end);
		break;
	default:
		error_msg("Error in switch for READ_COPY.", STOP);
		break;
	}
	return_value = check_line("End of COPY", FALSE, TRUE, TRUE, TRUE);
	/* empty, eof, keyword, print */
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int Phreeqc::
read_reaction_pressure(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads REACTION_PRESSURE data block
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
	assert(!reading_database());

	// Make instance, set n_user, n_user_end, description
	cxxPressure atm(this->phrq_io);
	char *ptr = line;
	char *description;
	int n_user, n_user_end;
	read_number_description(ptr, &n_user, &n_user_end, &description);
	atm.Set_n_user(n_user);
	atm.Set_n_user_end(n_user);
	atm.Set_description(description);
	free_check_null(description);

	/*
	 *  Make parser
	 */
	//CParser parser(*phrq_io->get_istream(), this->phrq_io);
	CParser parser(this->phrq_io);
	if (pr.echo_input == FALSE) parser.set_echo_file(CParser::EO_NONE);
	atm.read(parser);
	if (atm.Get_base_error_count() == 0)
	{
		Rxn_pressure_map[n_user] = atm;
	}

	if (use.Get_pressure_in() == FALSE)
	{
		use.Set_pressure_in(true);
		use.Set_n_pressure_user(atm.Get_n_user());
	}

	// Make copies if necessary
	if (n_user_end > n_user)
	{
		int i;
		for (i = n_user + 1; i <= n_user_end; i++)
		{
			Utilities::Rxn_copy(Rxn_pressure_map, n_user, i);
		}
	}

	return cleanup_after_parser(parser);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_reaction_pressure_raw(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads REACTION_PRESSURE_RAW data block
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
	assert(!reading_database());

	cxxPressure atm(this->phrq_io);
	//char *ptr = line;
	//char *description;
	//int n_user, n_user_end;
	//read_number_description(ptr, &n_user, &n_user_end, &description);
	//atm.Set_n_user(n_user);
	//atm.Set_n_user_end(n_user);
	//atm.Set_description(description);
	//free_check_null(description);
	/*
	 *  Make parser
	 */
	//CParser parser(*phrq_io->get_istream(), this->phrq_io);
	CParser parser(this->phrq_io);
	if (pr.echo_input == FALSE) parser.set_echo_file(CParser::EO_NONE);
	atm.read_raw(parser);

	// Store
	if (atm.Get_base_error_count() == 0)
	{
		Rxn_pressure_map[atm.Get_n_user()] = atm;
	}

	// Make copies if necessary

	Utilities::Rxn_copies(Rxn_pressure_map, atm.Get_n_user(), atm.Get_n_user_end());


	return cleanup_after_parser(parser);
}
int Phreeqc::
cleanup_after_parser(CParser &parser)
{
	// check_key sets next_keyword
	if (parser.get_m_line_type() == PHRQ_io::LT_EOF)
	{
		strcpy(line, "");
		strcpy(line_save, "");
		next_keyword = Keywords::KEY_END;
		return(TRUE); 
	}
	int return_value = check_key(parser.line().c_str());

	// copy parser line to line and line_save
	// make sure there is enough space
	size_t l1 = strlen(parser.line().c_str()) + 1;
	size_t l2 = strlen(parser.line_save().c_str()) + 1;
	size_t l = (l1 > l2) ? l1 : l2;
	if (l >= (size_t) max_line)
	{
		max_line = l * 2;
		line_save =	(char *) PHRQ_realloc(line_save,
			(size_t) max_line * sizeof(char));
		if (line_save == NULL)
			malloc_error();
		line = (char *) PHRQ_realloc(line, (size_t) max_line * sizeof(char));
		if (line == NULL)
			malloc_error();
	}
	strcpy(line, parser.line().c_str());
	strcpy(line_save, parser.line_save().c_str());
	return return_value;
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
read_temperature(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads REACTION_TEMPERATURE data block
 *
 *      Arguments:
 *         none
 *
 *      Returns:
 *         KEYWORD if keyword encountered, input_error may be incremented if
 *                    a keyword is encountered in an unexpected position
 *         EOF     if eof encountered while reading mass balance concentrations
 *         ERROR   if error occurred reading data
 *
 */
	assert(!reading_database());

	// Make instance, set n_user, n_user_end, description
	cxxTemperature t_react(this->phrq_io);
	char *ptr = line;
	char *description;
	int n_user, n_user_end;
	read_number_description(ptr, &n_user, &n_user_end, &description);
	t_react.Set_n_user(n_user);
	t_react.Set_n_user_end(n_user);
	t_react.Set_description(description);
	free_check_null(description);

	/*
	 *  Make parser
	 */
	//CParser parser(*phrq_io->get_istream(), this->phrq_io);
	CParser parser(this->phrq_io);
	if (pr.echo_input == FALSE) parser.set_echo_file(CParser::EO_NONE);
	t_react.read(parser);
	if (t_react.Get_base_error_count() == 0)
	{
		Rxn_temperature_map[n_user] = t_react;
	}

	if (use.Get_temperature_in() == FALSE)
	{
		use.Set_temperature_in(true);
		use.Set_n_temperature_user(t_react.Get_n_user());
	}

	// Make copies if necessary
	if (n_user_end > n_user)
	{
		int i;
		for (i = n_user + 1; i <= n_user_end; i++)
		{
			Utilities::Rxn_copy(Rxn_temperature_map, n_user, i);
		}
	}

	return cleanup_after_parser(parser);
}
