#include "Phreeqc.h"
#include "phqalloc.h"
#include "phrqproto.h"

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
read_isotopes(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads master species information for isotopes
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

	int j, l;
	struct master_isotope *master_isotope_ptr, *master_isotope_ptr_major;
	char token[MAX_LENGTH];
	struct element *elt_ptr;

	int return_value, opt, opt_save;
	char *next_char;
	const char *opt_list[] = {
		"isotope",				/* 0 */
		"total_is_major"		/* 1 */
	};
	int count_opt_list = 2;
	/*
	if (svnid == NULL)
		fprintf(stderr, " ");
	*/

	master_isotope_ptr = NULL;
	elt_ptr = NULL;
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
		case 0:				/* isotope */
			if (elt_ptr == NULL)
			{
				sprintf(error_string,
						"The element of which this isotope is a minor isotope has not been defined, %s. ISOTOPES data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			/*
			 *  Save an isotope
			 */
			master_isotope_ptr = NULL;
			j = copy_token(token, &next_char, &l);
			master_isotope_ptr = master_isotope_store(token, TRUE);
			master_isotope_ptr->elt = elt_ptr;
			master_isotope_ptr->minor_isotope = TRUE;
			master_isotope_ptr->total_is_major = FALSE;
			/*
			 *  Read units
			 */
			if (copy_token(token, &next_char, &l) == EMPTY)
			{
				sprintf(error_string,
						"Expecting units for isotopic values, %s. ISOTOPES data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			master_isotope_ptr->units = string_hsave(token);
			/*
			 *  Read standard 
			 */
			if (copy_token(token, &next_char, &l) == EMPTY)
			{
				sprintf(error_string,
						"Expecting isotope ratio of standard, %s. ISOTOPES data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			sscanf(token, SCANFORMAT, &(master_isotope_ptr->standard));
			opt_save = OPTION_DEFAULT;
			break;
		case 1:				/* total_is_major_isotope */
#ifdef SKIP
			if (elt_ptr == NULL)
			{
				sprintf(error_string,
						"The element of which this isotope is a minor isotope has not been defined, %s. ISOTOPES data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			master_isotope_ptr_major->total_is_major =
				get_true_false(next_char, TRUE);
#endif
			sprintf(error_string,
					"Obsolete identifier. The total of the element must be the sum of all isotopes. ISOTOPES data block.\n%s",
					line);
			warning_msg(error_string);
			break;
		case OPTION_DEFAULT:
/*
 *   Read and element name
 */
			if (copy_token(token, &next_char, &l) == EMPTY)
			{
				sprintf(error_string,
						"Expecting an element name for isotope definition, %s. ISOTOPES data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			elt_ptr = element_store(token);
			master_isotope_ptr = master_isotope_store(token, TRUE);
			master_isotope_ptr_major = master_isotope_ptr;
			master_isotope_ptr->elt = elt_ptr;
			master_isotope_ptr->minor_isotope = FALSE;
			master_isotope_ptr->total_is_major = FALSE;
			opt_save = OPTION_DEFAULT;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
read_calculate_values(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads basic code with which to calculate calculate_value
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
	char *ptr;
	int l, length, line_length, n;
	int return_value, opt, opt_save;
	char token[MAX_LENGTH];
	struct calculate_value *calculate_value_ptr;
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
	calculate_value_ptr = NULL;
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
			error_msg("Unknown input in CALCULATE_VALUE keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case 0:				/* start */
			opt_save = OPT_1;
			break;
		case 1:				/* end */
			opt_save = OPTION_DEFAULT;
			break;
		case OPTION_DEFAULT:	/* read calculate_value name */
/*
 *   Read calculate_value name
 */
			if (copy_token(token, &next_char, &l) == EMPTY)
			{
				sprintf(error_string,
						"Expecting a name for calculate_value definition, %s. CALCULATE_VALUES data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			calculate_value_ptr = calculate_value_store(token, TRUE);
			calculate_value_ptr->new_def = TRUE;
			calculate_value_ptr->commands =
				(char *) PHRQ_malloc(sizeof(char));
			if (calculate_value_ptr->commands == NULL)
				malloc_error();
			calculate_value_ptr->commands[0] = '\0';
			calculate_value_ptr->linebase = NULL;
			calculate_value_ptr->varbase = NULL;
			calculate_value_ptr->loopbase = NULL;
			opt_save = OPT_1;
			break;

		case OPT_1:			/* read command */
			length = (int) strlen(calculate_value_ptr->commands);
			line_length = (int) strlen(line);
			calculate_value_ptr->commands =
				(char *) PHRQ_realloc(calculate_value_ptr->commands,
									  (size_t) (length + line_length +
												2) * sizeof(char));
			if (calculate_value_ptr->commands == NULL)
				malloc_error();
			calculate_value_ptr->commands[length] = ';';
			calculate_value_ptr->commands[length + 1] = '\0';
			strcat((calculate_value_ptr->commands), line);
			opt_save = OPT_1;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
/*	output_temp_msg(sformatf( "%s", calculate_value[0].commands));
 */ return (return_value);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
read_isotope_ratios(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads isotope_ratio info, ratios are calculated with
 *      Basic programs read in CALCULATE_VALUE data block
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
	char *ptr;
	int l, n;
	int return_value, opt, opt_save;
	char token[MAX_LENGTH];
	struct isotope_ratio *isotope_ratio_ptr;
	char *description;
	int n_user, n_user_end;
	char *next_char;
	const char *opt_list[] = {
		"no_options"			/* 0 */
	};
	int count_opt_list = 0;
/*
 *   Read number (not currently used)
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
	isotope_ratio_ptr = NULL;
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
			error_msg("Unknown input in ISOTOPE_RATIOS keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case OPTION_DEFAULT:	/* read isotope_ratio name */
/*
 *   Read isotope_ratio name
 */
			if (copy_token(token, &next_char, &l) == EMPTY)
			{
				sprintf(error_string,
						"Expecting a name for isotope_ratio definition, %s. ISOTOPE_RATIOS data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			isotope_ratio_ptr = isotope_ratio_store(token, TRUE);
			/*
			 *  Read isotope
			 */
			if (copy_token(token, &next_char, &l) == EMPTY)
			{
				sprintf(error_string,
						"Expecting a name of isotope for an isotope_ratio definition, %s. ISOTOPE_RATIOS data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			isotope_ratio_ptr->isotope_name = string_hsave(token);
			opt_save = OPTION_DEFAULT;
			break;
		}
		if (return_value == EOF || return_value == KEYWORD)
			break;
	}
	return (return_value);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
read_isotope_alphas(void)
/* ---------------------------------------------------------------------- */
{
/*
 *      Reads isotope_alpha info, ratios are calculated with
 *      Basic programs read in CALCULATE_VALUE data block
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
	char *ptr;
	int l, n;
	int return_value, opt, opt_save;
	char token[MAX_LENGTH];
	struct isotope_alpha *isotope_alpha_ptr;
	char *description;
	int n_user, n_user_end;
	char *next_char;
	const char *opt_list[] = {
		"no_options"			/* 0 */
	};
	int count_opt_list = 0;
/*
 *   Read number (not currently used)
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
	isotope_alpha_ptr = NULL;
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
			error_msg("Unknown input in ISOTOPE_ALPHAS keyword.", CONTINUE);
			error_msg(line_save, CONTINUE);
			break;
		case OPTION_DEFAULT:	/* read isotope_alpha name */
/*
 *   Read isotope_alpha name
 */
			if (copy_token(token, &next_char, &l) == EMPTY)
			{
				sprintf(error_string,
						"Expecting a name for isotope_alpha definition, %s. ISOTOPE_ALPHAS data block.",
						line);
				error_msg(error_string, CONTINUE);
				input_error++;
				break;
			}
			isotope_alpha_ptr = isotope_alpha_store(token, TRUE);
			isotope_alpha_ptr->name = string_hsave(token);
			if (copy_token(token, &next_char, &l) != EMPTY)
			{
				isotope_alpha_ptr->named_logk = string_hsave(token);
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
int CLASS_QUALIFIER
add_isotopes(struct solution *solution_ptr)
/* ---------------------------------------------------------------------- */
{
	int i;
	char *ptr;
	struct master_isotope *master_isotope_ptr;
	LDBLE total_moles;
	/*
	 * zero out isotopes
	 */
	for (i = 0; i < count_master_isotope; i++)
	{
		master_isotope[i]->moles = 0;
	}
	master_isotope_ptr = master_isotope_search("H");
	if (master_isotope_ptr != NULL)
	{
		total_moles = total_h_x;
		calculate_isotope_moles(master_isotope_ptr->elt, solution_ptr,
								total_moles);
	}
	master_isotope_ptr = master_isotope_search("O");
	if (master_isotope_ptr != NULL)
	{
		total_moles = total_o_x;
		calculate_isotope_moles(master_isotope_ptr->elt, solution_ptr,
								total_moles);
	}
	for (i = 0; solution_ptr->totals[i].description != NULL; i++)
	{
		ptr = solution_ptr->totals[i].description;
		master_isotope_ptr = master_isotope_search(ptr);
		if (master_isotope_ptr == NULL)
			continue;
		if (master_isotope_ptr->minor_isotope == FALSE)
		{
			total_moles = total(master_isotope_ptr->name);
			calculate_isotope_moles(master_isotope_ptr->elt, solution_ptr,
									total_moles);
		}
	}
	/*
	 * Set isotopes flag
	 */
	initial_solution_isotopes = FALSE;
	for (i = 0; i < count_master_isotope; i++)
	{
		if (master_isotope[i]->minor_isotope == TRUE
			&& master_isotope[i]->moles > 0)
		{
			initial_solution_isotopes = TRUE;
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
calculate_isotope_moles(struct element *elt_ptr,
						struct solution *solution_ptr, LDBLE total_moles)
/* ---------------------------------------------------------------------- */
{
	int i, j, l_iter;
	int count_isotopes, total_is_major;
	char *ptr;
	struct master_isotope *master_isotope_ptr, *master_isotope_ptr1;
	struct master_isotope list[MAX_ELTS];
	LDBLE m_major, tot;
	/*
	 *  Get total concentration of elt_ptr
	 */
	if (total_moles <= 0)
	{
		sprintf(error_string,
				"Can not calculate molality of isotopes, molality of element is zero, %s",
				elt_ptr->name);
		warning_msg(error_string);
		return (ERROR);
	}
	m_major = total_moles;
	/*
	 *  Make a list of isotopes
	 */
	count_isotopes = 0;
	total_is_major = FALSE;
	master_isotope_ptr = master_isotope_search("H");
	if ((master_isotope_ptr != NULL) && (master_isotope_ptr->elt == elt_ptr))
	{
		memcpy(&(list[count_isotopes]), master_isotope_ptr,
			   sizeof(struct master_isotope));
		list[count_isotopes].ratio = 1.0;
		if (list[count_isotopes].minor_isotope == FALSE)
		{
			total_is_major = list[count_isotopes].total_is_major;
		}
		count_isotopes++;
	}
	master_isotope_ptr = master_isotope_search("O");
	if ((master_isotope_ptr != NULL) && (master_isotope_ptr->elt == elt_ptr))
	{
		memcpy(&(list[count_isotopes]), master_isotope_ptr,
			   sizeof(struct master_isotope));
		list[count_isotopes].ratio = 1.0;
		if (list[count_isotopes].minor_isotope == FALSE)
		{
			total_is_major = list[count_isotopes].total_is_major;
		}
		count_isotopes++;
	}
	for (i = 0; solution_ptr->totals[i].description != NULL; i++)
	{
		ptr = solution_ptr->totals[i].description;
		master_isotope_ptr = master_isotope_search(ptr);
		if (master_isotope_ptr == NULL)
			continue;
		if (master_isotope_ptr->elt != elt_ptr)
			continue;
		memcpy(&(list[count_isotopes]), master_isotope_ptr,
			   sizeof(struct master_isotope));
		if (list[count_isotopes].minor_isotope == FALSE)
		{
			total_is_major = list[count_isotopes].total_is_major;
		}
		count_isotopes++;
	}
	/*
	 *   Loop to calculate isotope molalities
	 */
	for (l_iter = 0; l_iter < itmax; l_iter++)
	{
		tot = 0;
		for (i = 0; i < count_isotopes; i++)
		{
			if (list[i].minor_isotope == FALSE)
			{
				list[i].moles = m_major;
				tot += m_major;
				continue;
			}
			if (strcmp_nocase(list[i].units, "permil") == 0)
			{
				from_permil(&(list[i]), m_major);
				tot += list[i].moles;
				continue;
			}
			if (strcmp_nocase(list[i].units, "pct") == 0)
			{
				from_pct(&(list[i]), total_moles);
				tot += list[i].moles;
				continue;
			}
			if (strcmp_nocase(list[i].units, "pmc") == 0)
			{
				from_pct(&(list[i]), total_moles);
				tot += list[i].moles;
				continue;
			}
			if (strcmp_nocase(list[i].units, "tu") == 0)
			{
				from_tu(&(list[i]));
				tot += list[i].moles;
				continue;
			}
			if (strcmp_nocase(list[i].units, "pci/l") == 0)
			{
				from_pcil(&(list[i]));
				tot += list[i].moles;
				continue;
			}
			sprintf(error_string, "Isotope units not recognized, %s",
					list[i].units);
			input_error++;
			error_msg(error_string, CONTINUE);
		}
		if (total_is_major == TRUE)
			break;
		if (fabs(total_moles - tot) < convergence_tolerance * total_moles)
		{
			break;
		}
		else
		{
			m_major = m_major * total_moles / tot;
		}
	}
	if (l_iter >= itmax)
	{
		error_msg("Failed to converge in CALCULATE_ISOTOPE_MOLES.", STOP);
	}
	/*
	 *  Update master_isotope
	 */
	for (j = 0; j < count_master_isotope; j++)
	{
		for (i = 0; i < count_isotopes; i++)
		{
			if (list[i].name == master_isotope[j]->name)
			{
				memcpy(master_isotope[j], &(list[i]),
					   sizeof(struct master_isotope));
			}
		}
	}
	/*
	 * Update solution
	 */
	master_isotope_ptr1 = master_isotope_search("H");
	if (master_isotope_ptr1 != NULL && master_isotope_ptr1->elt == elt_ptr)
	{
		total_h_x = m_major;
	}
	master_isotope_ptr1 = master_isotope_search("O");
	if (master_isotope_ptr1 != NULL && master_isotope_ptr1->elt == elt_ptr)
	{
		total_o_x = m_major;
	}
	for (i = 0; solution_ptr->totals[i].description != NULL; i++)
	{
		ptr = solution_ptr->totals[i].description;
		master_isotope_ptr = master_isotope_search(ptr);
		if (master_isotope_ptr == NULL)
			continue;
		if (master_isotope_ptr->elt != elt_ptr)
			continue;
		solution_ptr->totals[i].moles = master_isotope_ptr->moles;
		solution_ptr->totals[i].input_conc = master_isotope_ptr->moles;
	}

	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
from_permil(struct master_isotope *master_isotope_ptr, LDBLE major_total)
/* ---------------------------------------------------------------------- */
{
	LDBLE r;

	r = (master_isotope_ptr->ratio / 1000. +
		 1.0) * master_isotope_ptr->standard;
	master_isotope_ptr->moles = major_total * r;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
from_pct(struct master_isotope *master_isotope_ptr, LDBLE total_moles)
/* ---------------------------------------------------------------------- */
{
	master_isotope_ptr->moles =
		master_isotope_ptr->ratio / 100 * master_isotope_ptr->standard *
		total_moles;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
from_tu(struct master_isotope *master_isotope_ptr)
/* ---------------------------------------------------------------------- */
{
	master_isotope_ptr->moles =
		master_isotope_ptr->ratio * master_isotope_ptr->standard *
		mass_water_aq_x / gfw_water;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
from_pcil(struct master_isotope *master_isotope_ptr)
/* ---------------------------------------------------------------------- */
{
	master_isotope_ptr->moles =
		master_isotope_ptr->ratio * master_isotope_ptr->standard *
		mass_water_aq_x / gfw_water;
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
print_initial_solution_isotopes(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print isotopes for initial solution
 */
	int i, j;
	int print_isotope;

	if (pr.initial_isotopes == FALSE || pr.all == FALSE)
		return (OK);
	if (state != INITIAL_SOLUTION)
		return (OK);
	if (initial_solution_isotopes == FALSE)
		return (OK);
/*
 *   Print heading
 */
	print_centered("Isotopes");
	output_temp_msg(sformatf( "%10s\t%12s\t%12s\t%12s\t%12s\n\n", "Isotope",
			   "Molality", "Moles", "Ratio", "Units"));
	for (i = 0; i < count_master_isotope; i++)
	{
		if (master_isotope[i]->minor_isotope == FALSE)
		{
			print_isotope = FALSE;
			for (j = 0; j < count_master_isotope; j++)
			{
				if ((master_isotope[j]->elt == master_isotope[i]->elt) &&
					(master_isotope[j]->minor_isotope == TRUE) &&
					(master_isotope[j]->moles > 0))
				{
					print_isotope = TRUE;
					break;
				}
			}
			if (print_isotope == FALSE)
				continue;
			/*
			 *  Print isotope values
			 */
			output_temp_msg(sformatf( "%10s\t%12.5e\t%12.5e\n",
					   master_isotope[i]->name,
					   (double) (master_isotope[i]->moles / mass_water_aq_x),
					   (double) master_isotope[i]->moles));
			for (j = 0; j < count_master_isotope; j++)
			{
				if (i == j)
					continue;
				if ((master_isotope[j]->elt == master_isotope[i]->elt) &&
					(master_isotope[j]->minor_isotope == TRUE))
				{
					output_temp_msg(sformatf(
							   "%10s\t%12.5e\t%12.5e\t%12.5e\t%12s\n",
							   master_isotope[j]->name,
							   (double) (master_isotope[j]->moles /
										 mass_water_aq_x),
							   (double) master_isotope[j]->moles,
							   (double) master_isotope[j]->ratio,
							   master_isotope[j]->units));
				}
			}
			output_temp_msg(sformatf( "\n"));
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
punch_isotopes(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Punch isotope ratios relative to standards
 */
	int i;
	LDBLE iso;
	struct isotope_ratio *isotope_ratio_ptr;
	struct master_isotope *master_isotope_ptr;

	if (punch.in == FALSE || punch.isotopes == FALSE)
		return (OK);
	if (punch.count_isotopes == 0)
		return (OK);
	for (i = 0; i < punch.count_isotopes; i++)
	{
		iso = MISSING;
		if (state == INITIAL_SOLUTION)
		{
			isotope_ratio_ptr = isotope_ratio_search(punch.isotopes[i].name);
			if (isotope_ratio_ptr != NULL)
			{
				master_isotope_ptr =
					master_isotope_search(isotope_ratio_ptr->isotope_name);
				if (master_isotope_ptr != NULL
					&& master_isotope_ptr->minor_isotope == TRUE)
				{
					iso = master_isotope_ptr->ratio;
				}
			}
		}
		else
		{
			isotope_ratio_ptr = isotope_ratio_search(punch.isotopes[i].name);
			if (isotope_ratio_ptr != NULL)
			{
				iso = isotope_ratio_ptr->converted_ratio;
			}
		}
		if (punch.high_precision == FALSE)
		{
			fpunchf(sformatf("I_%s", punch.isotopes[i].name), "%12.4e\t",
				(double) iso);
		}
		else
		{
			fpunchf(sformatf("I_%s", punch.isotopes[i].name), "%20.12e\t",
				(double) iso);
		}

	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
punch_calculate_values(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Punch calculate values
 */
	int i;
	LDBLE result;
	struct calculate_value *calculate_value_ptr;

	if (punch.in == FALSE || punch.calculate_values == FALSE)
		return (OK);
	if (punch.count_calculate_values == 0)
		return (OK);
	for (i = 0; i < punch.count_calculate_values; i++)
	{
		result = MISSING;
		calculate_value_ptr =
			calculate_value_search(punch.calculate_values[i].name);
		if (calculate_value_ptr != NULL)
		{
			result = calculate_value_ptr->value;
		}
		if (punch.high_precision == FALSE)
		{
			fpunchf(sformatf("V_%s", punch.calculate_values[i].name),
				"%12.4e\t", (double) result);
		}
		else
		{
			fpunchf(sformatf("V_%s", punch.calculate_values[i].name),
				"%20.12e\t", (double) result);
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
print_isotope_ratios(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print isotopes for initial solution
 */
	int i, j;
	int print_isotope;
	struct master *master_ptr;
	struct master_isotope *master_isotope_ptr;
	struct isotope_ratio *isotope_ratio_ptr;
	char token[MAX_LENGTH];


	if (pr.isotope_ratios == FALSE || pr.all == FALSE)
		return (OK);
	if (state == INITIAL_SOLUTION)
		return (OK);
/*
 *   Print heading
 */
	print_isotope = FALSE;
	for (i = 0; i < count_master_isotope; i++)
	{
		if (master_isotope[i]->minor_isotope == FALSE)
			continue;
		master_ptr = master_bsearch(master_isotope[i]->name);
		if (master_ptr == NULL)
			continue;
		if (master_ptr->total > 0 || master_ptr->s->moles > 0)
		{
			print_isotope = TRUE;
			break;
		}
	}
	if (print_isotope == FALSE)
		return (OK);

	print_centered("Isotope Ratios");
	output_temp_msg(sformatf( "%25s\t%12s\t%15s\n\n", "Isotope Ratio",
			   "Ratio", "Input Units"));

	for (j = 0; j < count_isotope_ratio; j++)
	{
		if (isotope_ratio[j]->ratio == MISSING)
			continue;
		isotope_ratio_ptr = isotope_ratio[j];
		master_isotope_ptr =
			master_isotope_search(isotope_ratio[j]->isotope_name);
		/*
		 *  Print isotope ratio
		 */
		strcpy(token, isotope_ratio[j]->name);
		while (replace("_", " ", token) == TRUE);
		output_temp_msg(sformatf( "     %-20s\t%12.5e\t%15.5g  %-10s\n",
				   token, (double) isotope_ratio[j]->ratio,
				   (double) isotope_ratio[j]->converted_ratio,
				   master_isotope_ptr->units));
	}
	output_temp_msg(sformatf( "\n"));
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
print_isotope_alphas(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Print isotopes for initial solution
 */
	int i, j;
	int print_isotope;
	struct master *master_ptr;
	struct isotope_alpha *isotope_alpha_ptr;
	char token[MAX_LENGTH];
	LDBLE log_alpha;

	if (pr.isotope_alphas == FALSE || pr.all == FALSE)
		return (OK);
	if (state == INITIAL_SOLUTION)
		return (OK);
/*
 *   Print heading
 */
	print_isotope = FALSE;
	for (i = 0; i < count_master_isotope; i++)
	{
		if (master_isotope[i]->minor_isotope == FALSE)
			continue;
		master_ptr = master_bsearch(master_isotope[i]->name);
		if (master_ptr == NULL)
			continue;
		if (master_ptr->total > 0 || master_ptr->s->moles > 0)
		{
			print_isotope = TRUE;
			break;
		}
	}
	if (print_isotope == FALSE)
		return (OK);

	print_centered("Isotope Alphas");
	output_temp_msg(sformatf( "%75s\n", "1000ln(Alpha)"));
	output_temp_msg(sformatf( "%79s\n", "----------------------"));
	output_temp_msg(sformatf( "%-37s%14s%14s%12.1f C\n\n",
			   "     Isotope Ratio", "Solution alpha", "Solution",
			   (double) tc_x));

	for (j = 0; j < count_isotope_alpha; j++)
	{
		if (isotope_alpha[j]->value == MISSING)
			continue;
		isotope_alpha_ptr = isotope_alpha[j];
		/*
		 *  Print isotope ratio
		 */
		strcpy(token, isotope_alpha[j]->name);
		while (replace("_", " ", token) == TRUE);
		if (isotope_alpha[j]->named_logk != NULL)
		{
			if (isotope_alpha[j]->value <= 0)
			{
				log_alpha = -999.999;
			}
			else
			{
				log_alpha = 1000 * log(isotope_alpha[j]->value);
			}
			output_temp_msg(sformatf( "%-37s%14.5g%14.5g%14.5g\n", token,
					   (double) isotope_alpha[j]->value, (double) log_alpha,
					   (double) (1000 *
								 calc_logk_n(isotope_alpha[j]->named_logk) *
								 LOG_10)));
		}
		else
		{
			output_temp_msg(sformatf( "%-37s%14.5g%14.5g\n", token,
					   (double) isotope_alpha[j]->value,
					   (double) (1000 * log(isotope_alpha[j]->value))));
		}

	}
	output_temp_msg(sformatf( "\n"));
	return (OK);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
calculate_values(void)
/* ---------------------------------------------------------------------- */
{
	int j;
	struct calculate_value *calculate_value_ptr;
	struct isotope_ratio *isotope_ratio_ptr;
	struct isotope_alpha *isotope_alpha_ptr;
	struct master_isotope *master_isotope_ptr;
	char l_command[] = "run";


	/*
	 * initialize ratios as missing
	 */
	for (j = 0; j < count_calculate_value; j++)
	{
		calculate_value[j]->calculated = FALSE;
		calculate_value[j]->value = MISSING;
	}

	for (j = 0; j < count_calculate_value; j++)
	{
		calculate_value_ptr = calculate_value[j];
		rate_moles = NAN;
		if (calculate_value_ptr->new_def == TRUE)
		{
			if (basic_compile
				(calculate_value[j]->commands, &calculate_value[j]->linebase,
				 &calculate_value[j]->varbase,
				 &calculate_value[j]->loopbase) != 0)
			{
				sprintf(error_string,
						"Fatal Basic error in CALCULATE_VALUES %s.",
						calculate_value[j]->name);
				error_msg(error_string, STOP);
			}
			calculate_value_ptr->new_def = FALSE;
		}
		if (basic_run
			(l_command, calculate_value[j]->linebase,
			 calculate_value[j]->varbase, calculate_value[j]->loopbase) != 0)
		{
			sprintf(error_string, "Fatal Basic error in calculate_value %s.",
					calculate_value[j]->name);
			error_msg(error_string, STOP);
		}
		if (rate_moles == NAN)
		{
			sprintf(error_string, "Calculated value not SAVE'd for %s.",
					calculate_value[j]->name);
			error_msg(error_string, STOP);
		}
		else
		{
			calculate_value[j]->calculated = TRUE;
			calculate_value[j]->value = rate_moles;
		}
	}
	for (j = 0; j < count_isotope_ratio; j++)
	{
		isotope_ratio_ptr = isotope_ratio[j];
		master_isotope_ptr =
			master_isotope_search(isotope_ratio_ptr->isotope_name);
		calculate_value_ptr = calculate_value_search(isotope_ratio_ptr->name);
		/*
		 *  Calculate converted isotope ratio
		 */
		if (calculate_value_ptr->value == MISSING)
		{
			isotope_ratio_ptr->ratio = MISSING;
			isotope_ratio_ptr->converted_ratio = MISSING;
		}
		else
		{
			isotope_ratio_ptr->ratio = calculate_value_ptr->value;
			isotope_ratio_ptr->converted_ratio =
				convert_isotope(master_isotope_ptr,
								calculate_value_ptr->value);
		}
	}
	for (j = 0; j < count_isotope_alpha; j++)
	{
		isotope_alpha_ptr = isotope_alpha[j];
		calculate_value_ptr = calculate_value_search(isotope_alpha_ptr->name);
		/*
		 *  Calculate converted isotope ratio
		 */
		if (calculate_value_ptr->value == MISSING)
		{
			isotope_alpha_ptr->value = MISSING;
		}
		else
		{
			isotope_alpha_ptr->value = calculate_value_ptr->value;
		}
	}
	return (OK);
}

/* ---------------------------------------------------------------------- */
LDBLE CLASS_QUALIFIER
convert_isotope(struct master_isotope * master_isotope_ptr, LDBLE ratio)
/* ---------------------------------------------------------------------- */
{
	char *units;
	units = master_isotope_ptr->units;

	if (strcmp_nocase(units, "permil") == 0)
	{
		return ((ratio / master_isotope_ptr->standard - 1) * 1000);
	}
	if (strcmp_nocase(units, "pct") == 0)
	{
		return (ratio / master_isotope_ptr->standard * 100.);
	}
	if (strcmp_nocase(units, "pmc") == 0)
	{
		return (ratio / master_isotope_ptr->standard * 100.);
	}
	if (strcmp_nocase(units, "tu") == 0)
	{
		return (ratio / master_isotope_ptr->standard);
	}
	if (strcmp_nocase(units, "pci/l") == 0)
	{
		return (ratio / master_isotope_ptr->standard);
	}
	sprintf(error_string,
			"Did not recognize isotope units in convert_isotope, %s", units);
	error_msg(error_string, STOP);
	return (-99.0);
}

/*
 *  Utility routines for master_isotope
 */

/* ---------------------------------------------------------------------- */
struct master_isotope * CLASS_QUALIFIER
master_isotope_store(const char *name, int replace_if_found)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for master_isotope.
 *
 *   Pointer to a master_isotope structure is always returned.
 *
 *   If the string is not found, a new entry is made in the hash table. Pointer to 
 *      the new structure is returned.
 *   If "name" is found and replace is true, pointers in old master_isotope structure
 *      are freed and replaced with additional input.
 *   If "name" is found and replace is false, the old master_isotope structure is not
 *      modified and a pointer to it is returned.
 *
 *   Arguments:
 *      name    input, character string to be found in "master_isotope".
 *      replace_if_found input, TRUE means reinitialize master_isotope structure if found
 *                     FALSE means just return pointer if found.
 *
 *   Returns:
 *      pointer to master_isotope structure "master_isotope" where "name" can be found.
 */
	int n;
	struct master_isotope *master_isotope_ptr;
	ENTRY item, *found_item;
	char token[MAX_LENGTH];
/*
 *   Search list
 */
	strcpy(token, name);

	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(master_isotope_hash_table, item, FIND);

	if (found_item != NULL && replace_if_found == FALSE)
	{
		master_isotope_ptr = (struct master_isotope *) (found_item->data);
		return (master_isotope_ptr);
	}
	else if (found_item != NULL && replace_if_found == TRUE)
	{
		master_isotope_ptr = (struct master_isotope *) (found_item->data);
		master_isotope_init(master_isotope_ptr);
	}
	else
	{
		n = count_master_isotope++;
		/* make sure there is space in s */
		if (count_master_isotope >= max_master_isotope)
		{
			space((void **) ((void *) &master_isotope), count_master_isotope,
				  &max_master_isotope, sizeof(struct master_isotope *));
		}
		/* Make new master_isotope structure */
		master_isotope[n] = master_isotope_alloc();
		master_isotope_ptr = master_isotope[n];
	}
	/* set name and z in pointer in master_isotope structure */
	master_isotope_ptr->name = string_hsave(token);
/*
 *   Update hash table
 */
	item.key = master_isotope_ptr->name;
	item.data = (void *) master_isotope_ptr;
	found_item = hsearch_multi(master_isotope_hash_table, item, ENTER);
	if (found_item == NULL)
	{
		sprintf(error_string, "Hash table error in master_isotope_store.");
		error_msg(error_string, CONTINUE);
	}

	return (master_isotope_ptr);
}

/* ---------------------------------------------------------------------- */
struct master_isotope * CLASS_QUALIFIER
master_isotope_alloc(void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a master_isotope structure, initializes
 *      arguments: void
 *      return: pointer to a master_isotope structure
 */
{
	struct master_isotope *master_isotope_ptr;
	master_isotope_ptr =
		(struct master_isotope *) PHRQ_malloc(sizeof(struct master_isotope));
	if (master_isotope_ptr == NULL)
		malloc_error();
/*
 *   set pointers in structure to NULL, variables to zero
 */
	master_isotope_init(master_isotope_ptr);

	return (master_isotope_ptr);
}

/* ---------------------------------------------------------------------- */
STATIC int CLASS_QUALIFIER
master_isotope_init(struct master_isotope *master_isotope_ptr)
/* ---------------------------------------------------------------------- */
/*
 *      return: pointer to a master_isotope structure
 */
{
/*
 *   set pointers in structure to NULL
 */
	master_isotope_ptr->name = NULL;
	master_isotope_ptr->units = NULL;
	master_isotope_ptr->standard = 0;
	return (OK);
}

/* ---------------------------------------------------------------------- */
struct master_isotope * CLASS_QUALIFIER
master_isotope_search(const char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for master_isotope.
 *
 *   Arguments:
 *      name    input, character string to be found in "master_isotope".
 *
 *   Returns:
 *      pointer to master_isotope structure "master_isotope" where "name" can be found.
 *      or NULL if not found.
 */
	struct master_isotope *master_isotope_ptr;
	ENTRY item, *found_item;
	char token[MAX_LENGTH];
/*
 *   Search list
 */
	strcpy(token, name);

	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(master_isotope_hash_table, item, FIND);

	if (found_item != NULL)
	{
		master_isotope_ptr = (struct master_isotope *) (found_item->data);
		return (master_isotope_ptr);
	}
	return (NULL);
}

/*
 *  Utility routines for calculate_value
 */

/* ---------------------------------------------------------------------- */
struct calculate_value * CLASS_QUALIFIER
calculate_value_store(const char *name, int replace_if_found)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for calculate_value.
 *
 *   Pointer to a calculate_value structure is always returned.
 *
 *   If the string is not found, a new entry is made in the hash table. Pointer to 
 *      the new structure is returned.
 *   If "name" is found and replace is true, pointers in old calculate_value structure
 *      are freed and replaced with additional input.
 *   If "name" is found and replace is false, the old calculate_value structure is not
 *      modified and a pointer to it is returned.
 *
 *   Arguments:
 *      name    input, character string to be found in "calculate_value".
 *      replace_if_found input, TRUE means reinitialize calculate_value structure if found
 *                     FALSE means just return pointer if found.
 *
 *   Returns:
 *      pointer to calculate_value structure "calculate_value" where "name" can be found.
 */
	int n;
	struct calculate_value *calculate_value_ptr;
	char token[MAX_LENGTH];
	ENTRY item, *found_item;
/*
 *   Search list
 */
	strcpy(token, name);
	str_tolower(token);
	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(calculate_value_hash_table, item, FIND);

	if (found_item != NULL && replace_if_found == FALSE)
	{
		calculate_value_ptr = (struct calculate_value *) (found_item->data);
		return (calculate_value_ptr);
	}
	else if (found_item != NULL && replace_if_found == TRUE)
	{
		calculate_value_ptr = (struct calculate_value *) (found_item->data);
		calculate_value_free(calculate_value_ptr);
		calculate_value_init(calculate_value_ptr);
	}
	else
	{
		n = count_calculate_value++;
		/* make sure there is space in s */
		if (count_calculate_value >= max_calculate_value)
		{
			space((void **) ((void *) &calculate_value),
				  count_calculate_value, &max_calculate_value,
				  sizeof(struct calculate_value *));
		}
		/* Make new calculate_value structure */
		calculate_value[n] = calculate_value_alloc();
		calculate_value_ptr = calculate_value[n];
	}
	/* set name and z in pointer in calculate_value structure */
	calculate_value_ptr->name = string_hsave(name);
/*
 *   Update hash table
 */
	item.key = string_hsave(token);
	item.data = (void *) calculate_value_ptr;
	found_item = hsearch_multi(calculate_value_hash_table, item, ENTER);
	if (found_item == NULL)
	{
		sprintf(error_string, "Hash table error in calculate_value_store.");
		error_msg(error_string, CONTINUE);
	}

	return (calculate_value_ptr);
}

/* ---------------------------------------------------------------------- */
struct calculate_value * CLASS_QUALIFIER
calculate_value_alloc(void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a calculate_value structure, initializes
 *      arguments: void
 *      return: pointer to a calculate_value structure
 */
{
	struct calculate_value *calculate_value_ptr;
	calculate_value_ptr =
		(struct calculate_value *)
		PHRQ_malloc(sizeof(struct calculate_value));
	if (calculate_value_ptr == NULL)
		malloc_error();
/*
 *   set pointers in structure to NULL, variables to zero
 */
	calculate_value_init(calculate_value_ptr);

	return (calculate_value_ptr);
}

/* ---------------------------------------------------------------------- */
STATIC int CLASS_QUALIFIER
calculate_value_init(struct calculate_value *calculate_value_ptr)
/* ---------------------------------------------------------------------- */
/*
 *      return: pointer to a calculate_value structure
 */
{
/*
 *   set pointers in structure to NULL
 */
	calculate_value_ptr->name = NULL;
	calculate_value_ptr->commands = NULL;
	calculate_value_ptr->linebase = NULL;
	calculate_value_ptr->varbase = NULL;
	calculate_value_ptr->loopbase = NULL;

	return (OK);
}

/* ---------------------------------------------------------------------- */
struct calculate_value * CLASS_QUALIFIER
calculate_value_search(const char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for calculate_value.
 *
 *   Arguments:
 *      name    input, character string to be found in "calculate_value".
 *
 *   Returns:
 *      pointer to calculate_value structure "calculate_value" where "name" can be found.
 *      or NULL if not found.
 */
	struct calculate_value *calculate_value_ptr;
	char token[MAX_LENGTH];
	ENTRY item, *found_item;
/*
 *   Search list
 */
	strcpy(token, name);
	str_tolower(token);
	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(calculate_value_hash_table, item, FIND);

	if (found_item != NULL)
	{
		calculate_value_ptr = (struct calculate_value *) (found_item->data);
		return (calculate_value_ptr);
	}
	return (NULL);
}

/* ---------------------------------------------------------------------- */
int CLASS_QUALIFIER
calculate_value_free(struct calculate_value *calculate_value_ptr)
/* ---------------------------------------------------------------------- */
{
/*
 *   Frees memory allocated within calculate_value[i], does not free calculate_value structure
 *   Input: i, number of calculate_value
 *   Return: OK
 */
	char cmd[] = "new; quit";

	if (calculate_value_ptr == NULL)
		return (ERROR);
	calculate_value_ptr->commands =
		(char *) free_check_null(calculate_value_ptr->commands);
	basic_run(cmd, calculate_value_ptr->linebase,
			  calculate_value_ptr->varbase, calculate_value_ptr->loopbase);
	calculate_value_ptr->linebase = NULL;
	calculate_value_ptr->varbase = NULL;
	calculate_value_ptr->loopbase = NULL;
	return (OK);
}

/*
 *  Utility routines for isotope_ratio
 */

/* ---------------------------------------------------------------------- */
struct isotope_ratio * CLASS_QUALIFIER
isotope_ratio_store(const char *name, int replace_if_found)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for isotope_ratio.
 *
 *   Pointer to a isotope_ratio structure is always returned.
 *
 *   If the string is not found, a new entry is made in the hash table. Pointer to 
 *      the new structure is returned.
 *   If "name" is found and replace is true, pointers in old isotope_ratio structure
 *      are freed and replaced with additional input.
 *   If "name" is found and replace is false, the old isotope_ratio structure is not
 *      modified and a pointer to it is returned.
 *
 *   Arguments:
 *      name    input, character string to be found in "isotope_ratio".
 *      replace_if_found input, TRUE means reinitialize isotope_ratio structure if found
 *                     FALSE means just return pointer if found.
 *
 *   Returns:
 *      pointer to isotope_ratio structure "isotope_ratio" where "name" can be found.
 */
	int n;
	struct isotope_ratio *isotope_ratio_ptr;
	char token[MAX_LENGTH];
	ENTRY item, *found_item;
/*
 *   Search list
 */
	strcpy(token, name);
	str_tolower(token);
	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(isotope_ratio_hash_table, item, FIND);

	if (found_item != NULL && replace_if_found == FALSE)
	{
		isotope_ratio_ptr = (struct isotope_ratio *) (found_item->data);
		return (isotope_ratio_ptr);
	}
	else if (found_item != NULL && replace_if_found == TRUE)
	{
		isotope_ratio_ptr = (struct isotope_ratio *) (found_item->data);
		isotope_ratio_init(isotope_ratio_ptr);
	}
	else
	{
		n = count_isotope_ratio++;
		/* make sure there is space in s */
		if (count_isotope_ratio >= max_isotope_ratio)
		{
			space((void **) ((void *) &isotope_ratio), count_isotope_ratio,
				  &max_isotope_ratio, sizeof(struct isotope_ratio *));
		}
		/* Make new isotope_ratio structure */
		isotope_ratio[n] = isotope_ratio_alloc();
		isotope_ratio_ptr = isotope_ratio[n];
	}
	/* set name and z in pointer in isotope_ratio structure */
	isotope_ratio_ptr->name = string_hsave(name);
/*
 *   Update hash table
 */
	item.key = string_hsave(token);
	item.data = (void *) isotope_ratio_ptr;
	found_item = hsearch_multi(isotope_ratio_hash_table, item, ENTER);
	if (found_item == NULL)
	{
		sprintf(error_string, "Hash table error in isotope_ratio_store.");
		error_msg(error_string, CONTINUE);
	}

	return (isotope_ratio_ptr);
}

/* ---------------------------------------------------------------------- */
struct isotope_ratio * CLASS_QUALIFIER
isotope_ratio_alloc(void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a isotope_ratio structure, initializes
 *      arguments: void
 *      return: pointer to a isotope_ratio structure
 */
{
	struct isotope_ratio *isotope_ratio_ptr;
	isotope_ratio_ptr =
		(struct isotope_ratio *) PHRQ_malloc(sizeof(struct isotope_ratio));
	if (isotope_ratio_ptr == NULL)
		malloc_error();
/*
 *   set pointers in structure to NULL, variables to zero
 */
	isotope_ratio_init(isotope_ratio_ptr);

	return (isotope_ratio_ptr);
}

/* ---------------------------------------------------------------------- */
STATIC int CLASS_QUALIFIER
isotope_ratio_init(struct isotope_ratio *isotope_ratio_ptr)
/* ---------------------------------------------------------------------- */
/*
 *      return: pointer to a isotope_ratio structure
 */
{
/*
 *   set pointers in structure to NULL
 */
	isotope_ratio_ptr->name = NULL;
	isotope_ratio_ptr->isotope_name = NULL;
	isotope_ratio_ptr->ratio = MISSING;
	isotope_ratio_ptr->converted_ratio = MISSING;

	return (OK);
}

/* ---------------------------------------------------------------------- */
struct isotope_ratio * CLASS_QUALIFIER
isotope_ratio_search(const char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for isotope_ratio.
 *
 *   Arguments:
 *      name    input, character string to be found in "isotope_ratio".
 *
 *   Returns:
 *      pointer to isotope_ratio structure "isotope_ratio" where "name" can be found.
 *      or NULL if not found.
 */
	struct isotope_ratio *isotope_ratio_ptr;
	char token[MAX_LENGTH];
	ENTRY item, *found_item;
/*
 *   Search list
 */
	strcpy(token, name);
	str_tolower(token);
	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(isotope_ratio_hash_table, item, FIND);

	if (found_item != NULL)
	{
		isotope_ratio_ptr = (struct isotope_ratio *) (found_item->data);
		return (isotope_ratio_ptr);
	}
	return (NULL);
}

/*
 *  Utility routines for isotope_alpha
 */

/* ---------------------------------------------------------------------- */
struct isotope_alpha * CLASS_QUALIFIER
isotope_alpha_store(const char *name, int replace_if_found)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for isotope_alpha.
 *
 *   Pointer to a isotope_alpha structure is always returned.
 *
 *   If the string is not found, a new entry is made in the hash table. Pointer to 
 *      the new structure is returned.
 *   If "name" is found and replace is true, pointers in old isotope_alpha structure
 *      are freed and replaced with additional input.
 *   If "name" is found and replace is false, the old isotope_alpha structure is not
 *      modified and a pointer to it is returned.
 *
 *   Arguments:
 *      name    input, character string to be found in "isotope_alpha".
 *      replace_if_found input, TRUE means reinitialize isotope_alpha structure if found
 *                     FALSE means just return pointer if found.
 *
 *   Returns:
 *      pointer to isotope_alpha structure "isotope_alpha" where "name" can be found.
 */
	int n;
	struct isotope_alpha *isotope_alpha_ptr;
	char token[MAX_LENGTH];
	ENTRY item, *found_item;
/*
 *   Search list
 */
	strcpy(token, name);
	str_tolower(token);
	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(isotope_alpha_hash_table, item, FIND);

	if (found_item != NULL && replace_if_found == FALSE)
	{
		isotope_alpha_ptr = (struct isotope_alpha *) (found_item->data);
		return (isotope_alpha_ptr);
	}
	else if (found_item != NULL && replace_if_found == TRUE)
	{
		isotope_alpha_ptr = (struct isotope_alpha *) (found_item->data);
		isotope_alpha_init(isotope_alpha_ptr);
	}
	else
	{
		n = count_isotope_alpha++;
		/* make sure there is space in s */
		if (count_isotope_alpha >= max_isotope_alpha)
		{
			space((void **) ((void *) &isotope_alpha), count_isotope_alpha,
				  &max_isotope_alpha, sizeof(struct isotope_alpha *));
		}
		/* Make new isotope_alpha structure */
		isotope_alpha[n] = isotope_alpha_alloc();
		isotope_alpha_ptr = isotope_alpha[n];
	}
	/* set name and z in pointer in isotope_alpha structure */
	isotope_alpha_ptr->name = string_hsave(name);
/*
 *   Update hash table
 */
	item.key = string_hsave(token);
	item.data = (void *) isotope_alpha_ptr;
	found_item = hsearch_multi(isotope_alpha_hash_table, item, ENTER);
	if (found_item == NULL)
	{
		sprintf(error_string, "Hash table error in isotope_alpha_store.");
		error_msg(error_string, CONTINUE);
	}

	return (isotope_alpha_ptr);
}

/* ---------------------------------------------------------------------- */
struct isotope_alpha * CLASS_QUALIFIER
isotope_alpha_alloc(void)
/* ---------------------------------------------------------------------- */
/*
 *   Allocates space to a isotope_alpha structure, initializes
 *      arguments: void
 *      return: pointer to a isotope_alpha structure
 */
{
	struct isotope_alpha *isotope_alpha_ptr;
	isotope_alpha_ptr =
		(struct isotope_alpha *) PHRQ_malloc(sizeof(struct isotope_alpha));
	if (isotope_alpha_ptr == NULL)
		malloc_error();
/*
 *   set pointers in structure to NULL, variables to zero
 */
	isotope_alpha_init(isotope_alpha_ptr);

	return (isotope_alpha_ptr);
}

/* ---------------------------------------------------------------------- */
STATIC int CLASS_QUALIFIER
isotope_alpha_init(struct isotope_alpha *isotope_alpha_ptr)
/* ---------------------------------------------------------------------- */
/*
 *      return: pointer to a isotope_alpha structure
 */
{
/*
 *   set pointers in structure to NULL
 */
	isotope_alpha_ptr->name = NULL;
	isotope_alpha_ptr->named_logk = NULL;
	isotope_alpha_ptr->value = MISSING;

	return (OK);
}

/* ---------------------------------------------------------------------- */
struct isotope_alpha * CLASS_QUALIFIER
isotope_alpha_search(const char *name)
/* ---------------------------------------------------------------------- */
{
/*
 *   Function locates the string "name" in the hash table for isotope_alpha.
 *
 *   Arguments:
 *      name    input, character string to be found in "isotope_alpha".
 *
 *   Returns:
 *      pointer to isotope_alpha structure "isotope_alpha" where "name" can be found.
 *      or NULL if not found.
 */
	struct isotope_alpha *isotope_alpha_ptr;
	char token[MAX_LENGTH];
	ENTRY item, *found_item;
/*
 *   Search list
 */
	strcpy(token, name);
	str_tolower(token);
	item.key = token;
	item.data = NULL;
	found_item = hsearch_multi(isotope_alpha_hash_table, item, FIND);

	if (found_item != NULL)
	{
		isotope_alpha_ptr = (struct isotope_alpha *) (found_item->data);
		return (isotope_alpha_ptr);
	}
	return (NULL);
}
