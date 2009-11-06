#include "Phreeqc.h"
//class Phreeqc
//{
//	Phreeqc(void);
//	~Phreeqc(void);
//};
#include <algorithm>			// std::replace

Phreeqc::Phreeqc(void)
{
	struct const_key keyword_temp[] = {
	{"eof", 0},
	{"end", 0},
	{"solution_species", 0},
	{"solution_master_species", 0},
	{"solution", 0},
	{"phases", 0},
	{"pure_phases", 0},
	{"reaction", 0},
	{"mix", 0},
	{"use", 0},
	{"save", 0},
	{"exchange_species", 0},
	{"exchange_master_species", 0},
	{"exchange", 0},
	{"surface_species", 0},
	{"surface_master_species", 0},
	{"surface", 0},
	{"reaction_temperature", 0},
	{"inverse_modeling", 0},
	{"gas_phase", 0},
	{"transport", 0},
	{"debug", 0},
	{"selected_output", 0},
	{"select_output", 0},
	{"knobs", 0},
	{"print", 0},
	{"equilibrium_phases", 0},
	{"equilibria", 0},
	{"equilibrium", 0},
	{"pure", 0},
	{"title", 0},
	{"comment", 0},
	{"advection", 0},
	{"kinetics", 0},
	{"incremental_reactions", 0},
	{"incremental", 0},
	{"rates", 0},
	{"solution_s", 0},
	{"user_print", 0},
	{"user_punch", 0},
	{"solid_solutions", 0},
	{"solid_solution", 0},
	{"solution_spread", 0},
	{"spread_solution", 0},
	{"selected_out", 0},
	{"select_out", 0},
	{"user_graph", 0},
	{"llnl_aqueous_model_parameters", 0},
	{"llnl_aqueous_model", 0},
	{"database", 0},
	{"named_analytical_expression", 0},
	{"named_analytical_expressions", 0},
	{"named_expressions", 0},
	{"named_log_k", 0},
	{"isotopes", 0},
	{"calculate_values", 0},
	{"isotope_ratios", 0},
	{"isotope_alphas", 0},
	{"copy", 0},
	{"pitzer", 0},
	{"sit", 0}
#ifdef PHREEQC_CPP
	,
	{"solution_raw", 0},
	{"exchange_raw", 0},
	{"surface_raw", 0},
	{"equilibrium_phases_raw", 0},
	{"kinetics_raw", 0},
	{"solid_solutions_raw", 0},
	{"gas_phase_raw", 0},
	{"reaction_raw", 0},
	{"mix_raw", 0},
	{"reaction_temperature_raw", 0},
	{"dump", 0},
	{"solution_modify", 0},
	{"equilibrium_phases_modify", 0},
	{"exchange_modify", 0},
	{"surface_modify", 0},
	{"solid_solutions_modify", 0},
	{"gas_phase_modify", 0},
	{"kinetics_modify", 0},
	{"delete", 0},
	{"run_cells", 0}
#endif
	};
	NKEYS = (sizeof(keyword_temp) / sizeof(struct const_key));	/* Number of valid keywords */

	int i;

	//keyword = (struct const_key *) PHRQ_malloc((size_t) (NKEYS * sizeof(const_key)));
	keyword = new const_key[NKEYS];
	for (i = 0; i < NKEYS; i++)
	{
		keyword[i].name = string_duplicate(keyword_temp[i].name);
		keyword[i].keycount = 0;
	}
	//cl1.c
	x_arg = NULL, res_arg = NULL, scratch = NULL;
	x_arg_max = 0, res_arg_max = 0, scratch_max = 0;

	// dw.c
	GASCON = 0.461522e0;
	TZ = 647.073e0;
	AA = 1.e0;
	G1 = 11.e0; 
	G2 = 44.333333333333e0;
	GF = 3.5e0;

	// model.c
	min_value = 1e-10;
	normal = NULL;
	ineq_array = NULL;
	res = NULL;
	cu = NULL;
	zero =	NULL;
	delta1 = NULL;
	iu = NULL;
	*is = NULL;
	*back_eq = NULL;
	normal_max = 0;
	ineq_array_max = 0;
	res_max = 0;
	cu_max = 0;
	zero_max = 0;
	delta1_max = 0;
	iu_max = 0;
	is_max = 0;
	back_eq_max = 0;

	// output.c
	output_callbacks = new Phreeqc::output_callback[MAX_CALLBACKS];
	count_output_callback = 0;
	forward_output_to_log = 0;

	// phqalloc.c
	s_pTail = NULL;

	// phreeqc_files.c

	default_data_base = string_duplicate("phreeqc.dat");

//FILE *input_file = NULL;
//FILE *database_file = NULL;
//FILE *output = NULL;		/* OUTPUT_MESSAGE */
//FILE *log_file = NULL;	/* OUTPUT_LOG */
//FILE *punch_file = NULL;	/* OUTPUT_PUNCH */
//FILE *error_file = NULL;	/* OUTPUT_ERROR */
//FILE *dump_file = NULL;	/* OUTPUT_DUMP */
}

Phreeqc::~Phreeqc(void)
{

	int i;
	for (i = 0; i < NKEYS; i++)
	{
		keyword[i].name = (char *) free_check_null((void *) keyword[i].name);
	}
	delete[] keyword;

	free_check_null(default_data_base);

	//this->clean_up();
}
