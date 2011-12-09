#if defined(REVISED_GASES)
#include "Phreeqc.h"
/* ---------------------------------------------------------------------- */
int Phreeqc::
setup_gas_phase(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Fill in data for gas phase unknown (sum of partial pressures)
 *   in unknown structure
 */
	int i;
	if (use.gas_phase_ptr == NULL)
		return (OK);
/*
 *   One for each gas component
 */
	gas_unknowns.clear();
	gas_unknown = NULL;
	use.gas_phase_ptr->total_moles = 0;
	for (i = 0; i < use.gas_phase_ptr->count_comps; i++)
	{
		x[count_unknowns]->type = GAS_MOLES;
		x[count_unknowns]->description = string_hsave("gas moles");
		x[count_unknowns]->phase = use.gas_phase_ptr->comps[i].phase;
		x[count_unknowns]->moles = use.gas_phase_ptr->comps[i].moles;
		if (x[count_unknowns]->moles <= 0)
		{
			x[count_unknowns]->moles = MIN_TOTAL;
		}
		x[count_unknowns]->ln_moles = log(x[count_unknowns]->moles);
		x[count_unknowns]->gas_phase = use.gas_phase_ptr;
		gas_unknowns.push_back(x[count_unknowns]);
		use.gas_phase_ptr->total_moles += x[count_unknowns]->moles;
		x[count_unknowns]->phase->moles_x = x[count_unknowns]->moles;
		count_unknowns++;
	}
	if (gas_unknowns.size() > 0)
	{
		gas_unknown = gas_unknowns[0];
	}

#ifdef SKIP
	if (use.gas_phase_ptr->type == VOLUME)
	{
		x[count_unknowns]->type = GAS_PRESSURE;
		x[count_unknowns]->description = string_hsave("gas pressure");
		x[count_unknowns]->moles = 0.0;
		x[count_unknowns]->gas_phase = use.gas_phase_ptr;
		gas_pressure_unknown = x[count_unknowns];
		x[count_unknowns]->pressure = patm_x;
		count_unknowns++;
	}
#endif
	return (OK);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
build_gas_phase(void)
/* ---------------------------------------------------------------------- */
{
/*
 *   Put coefficients into lists to sum iaps to test for equilibrium
 *   Put coefficients into lists to build jacobian for 
 *      sum of partial pressures equation and
 *      mass balance equations for elements contained in gases
 */
	int i, j;
	int row, col;
	struct master *master_ptr;
	struct rxn_token *rxn_ptr;
	struct gas_comp *gas_comp_ptr;
	struct phase *phase_ptr;
	struct unknown *unknown_ptr;
	LDBLE coef, coef_elt;

	if (gas_unknown == NULL)
		return (OK);
	for (i = 0; i < use.gas_phase_ptr->count_comps; i++)
	{
/*
 *   Determine elements in gas component
 */
		count_elts = 0;
		paren_count = 0;
		gas_comp_ptr = &(use.gas_phase_ptr->comps[i]);
		phase_ptr = gas_comp_ptr->phase;
		if (phase_ptr->rxn_x == NULL)
			continue;
		add_elt_list(phase_ptr->next_elt, 1.0);
#ifdef COMBINE
		change_hydrogen_in_elt_list(0);
#endif
/*
 *   Build mass balance sums for each element in gas
 */
		if (debug_prep == TRUE)
		{
			output_msg(sformatf( "\n\tMass balance summations %s.\n\n",
					   gas_comp_ptr->phase->name));
		}

		/* All elements in gas */
		for (j = 0; j < count_elts; j++)
		{
			unknown_ptr = NULL;
			if (strcmp(elt_list[j].elt->name, "H") == 0)
			{
				unknown_ptr = mass_hydrogen_unknown;
			}
			else if (strcmp(elt_list[j].elt->name, "O") == 0)
			{
				unknown_ptr = mass_oxygen_unknown;
			}
			else
			{
				if (elt_list[j].elt->primary->in == TRUE)
				{
					unknown_ptr = elt_list[j].elt->primary->unknown;
				}
				else if (elt_list[j].elt->primary->s->secondary != NULL)
				{
					unknown_ptr =
						elt_list[j].elt->primary->s->secondary->unknown;
				}
			}
			if (unknown_ptr != NULL)
			{
				coef = elt_list[j].coef;
				//store_mb(&(gas_comp_ptr->phase->moles_x), &(unknown_ptr->f),
				//		 coef);
				store_mb(&(gas_unknowns[i]->moles), &(unknown_ptr->f), coef);
				if (debug_prep == TRUE)
				{
					output_msg(sformatf( "\t\t%-24s%10.3f\n",
							   unknown_ptr->description, (double) coef));
				}
			}
		}
		if (use.gas_phase_ptr->type == PRESSURE)
		{
			/* Total pressure of gases */
			store_mb(&(gas_comp_ptr->phase->p_soln_x), &(gas_unknown->f),
					 1.0);
		}
/*
 *   Build jacobian sums for mass balance equations
 */
		if (debug_prep == TRUE)
		{
			output_msg(sformatf( "\n\tJacobian summations %s.\n\n",
					   phase_ptr->name));
		}
		for (j = 0; j < count_elts; j++)
		{
			unknown_ptr = NULL;
			if (strcmp(elt_list[j].elt->name, "H") == 0)
			{
				unknown_ptr = mass_hydrogen_unknown;
			}
			else if (strcmp(elt_list[j].elt->name, "O") == 0)
			{
				unknown_ptr = mass_oxygen_unknown;
			}
			else
			{
				if (elt_list[j].elt->primary->in == TRUE)
				{
					unknown_ptr = elt_list[j].elt->primary->unknown;
				}
				else if (elt_list[j].elt->primary->s->secondary != NULL)
				{
					unknown_ptr =
						elt_list[j].elt->primary->s->secondary->unknown;
				}
			}
			if (unknown_ptr == NULL)
			{
				continue;
			}
			if (debug_prep == TRUE)
			{
				output_msg(sformatf( "\n\t%s.\n",
						   unknown_ptr->description));
			}
			row = unknown_ptr->number * (count_unknowns + 1);
			coef_elt = elt_list[j].coef;
			for (rxn_ptr = phase_ptr->rxn_x->token + 1;
				 rxn_ptr->s != NULL; rxn_ptr++)
			{

				if (rxn_ptr->s->secondary != NULL
					&& rxn_ptr->s->secondary->in == TRUE)
				{
					master_ptr = rxn_ptr->s->secondary;
				}
				else
				{
					master_ptr = rxn_ptr->s->primary;
				}
				if (master_ptr == NULL)
				{
					sprintf(error_string,
							"Element needed for gas component, %s, is not in model.",
							phase_ptr->name);
					warning_msg(error_string);
					//error_msg(error_string, CONTINUE);
					//input_error++;
					continue;
				}
				if (debug_prep == TRUE)
				{
					output_msg(sformatf( "\t\t%s\n",
							   master_ptr->s->name));
				}
				if (master_ptr->unknown == NULL)
				{
					continue;
				}
				if (master_ptr->in == FALSE)
				{
					sprintf(error_string,
							"Element, %s, in phase, %s, is not in model.",
							master_ptr->elt->name, phase_ptr->name);
					error_msg(error_string, CONTINUE);
					input_error++;
				}
				col = master_ptr->unknown->number;
				coef = coef_elt * rxn_ptr->coef;
				store_jacob(&(gas_comp_ptr->phase->moles_x),
							&(array[row + col]), coef);
				if (debug_prep == TRUE)
				{
					output_msg(sformatf( "\t\t%-24s%10.3f\t%d\t%d\n",
							   master_ptr->s->name, (double) coef,
							   row / (count_unknowns + 1), col));
				}
			}
			if (use.gas_phase_ptr->type == PRESSURE)
			{
				/* derivative wrt total moles of gas */
				store_jacob(&(gas_comp_ptr->phase->fraction_x),
							&(array[row + gas_unknown->number]), coef_elt);
				if (debug_prep == TRUE)
				{
					output_msg(sformatf( "\t\t%-24s%10.3f\t%d\t%d\n",
							   "gas moles", (double) elt_list[j].coef,
							   row / (count_unknowns + 1),
							   gas_unknown->number));
				}
			}
		}
/*
 *   Build jacobian sums for sum of partial pressures equation
 */
		if (use.gas_phase_ptr->type != PRESSURE)
			continue;
		if (debug_prep == TRUE)
		{
			output_msg(sformatf( "\n\tPartial pressure eqn %s.\n\n",
					   phase_ptr->name));
		}
		unknown_ptr = gas_unknown;
		row = unknown_ptr->number * (count_unknowns + 1);
		for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL; rxn_ptr++)
		{
			if (rxn_ptr->s != s_eminus && rxn_ptr->s->in == FALSE)
			{
				sprintf(error_string,
					"Element in species, %s, in phase, %s, is not in model.",
					rxn_ptr->s->name, phase_ptr->name);
				warning_msg(error_string);
				//error_msg(error_string, CONTINUE);
				//input_error++;
			}
			else
			{
				if (rxn_ptr->s->secondary != NULL
					&& rxn_ptr->s->secondary->in == TRUE)
				{
					master_ptr = rxn_ptr->s->secondary;
				}
				else 
				{
					master_ptr = rxn_ptr->s->primary;
				}

				if (master_ptr == NULL)
				{
					sprintf(error_string,
						"Master species for %s, in phase, %s, is not in model.",
						rxn_ptr->s->name, phase_ptr->name);
					error_msg(error_string, CONTINUE);
					input_error++;
				}
				else
				{
					if (debug_prep == TRUE)
					{
						output_msg(sformatf( "\t\t%s\n", master_ptr->s->name));
					}
					if (master_ptr->unknown == NULL)
					{
						assert(false);
						continue;
					}
					if (master_ptr->in == FALSE)
					{
						sprintf(error_string,
							"Element, %s, in phase, %s, is not in model.",
							master_ptr->elt->name, phase_ptr->name);
						warning_msg(error_string);
						//error_msg(error_string, CONTINUE);
						//input_error++;
					}
					col = master_ptr->unknown->number;
					coef = rxn_ptr->coef;
					store_jacob(&(gas_comp_ptr->phase->p_soln_x), &(array[row + col]), coef);
					if (debug_prep == TRUE)
					{
						output_msg(sformatf( "\t\t%-24s%10.3f\t%d\t%d\n",
							master_ptr->s->name, (double) coef,
							row / (count_unknowns + 1), col));
					}
				}
			}
		}
	}
	return (OK);
}
/* ---------------------------------------------------------------------- */
LDBLE Phreeqc::
calc_PR(void)
/* ---------------------------------------------------------------------- */
/*  Calculate fugacity and fugacity coefficient for gas pressures if critical T and P
    are defined.
  1) Solve molar volume V_m or total pressure P from Peng-Robinson's EOS:
  P = R * T / (V_m - b) - a * aa / (V_m^2 + 2 * b * V_m - b^2)
     a = 0.457235 * (R * T_c)^2 / P_c
     b = 0.077796 * R * T_c / P_c
     aa = (1 + kk * (1 - T_r^0.5))^2
     kk = 0.37464 + 1.54226 * omega - 0.26992 * omega^2
     T_r = T / T_c
  multicomponent gas phase:
     use: b_sum = Sum(x_i * b), x_i is mole-fraction
          a_aa_sum = Sum_i( Sum_j(x_i * x_j * (a_i * aa_i * a_j * aa_j)^0.5) )
  2) Find the fugacity coefficient phi for gas i:
  log(phi_i) = B_ratio * (z - 1) - log(z - B) + A / (2.8284 * B) * (B_ratio - 2 / a_aa_sum * a_aa_sum2) *\
           log((z + 2.4142 * B) / (z - 0.4142 * B))
     B_ratio = b_i / b_sum
     A = a_aa_sum * P / R_TK^2
     B = b_sum * P / R_TK
     a_aa_sum2 = Sum_j(x_j * (a_aa_i * a_aa_j)^0.5
  3) correct the solubility of gas i with:
  pr_si_f = log10(phi_i) -  Delta_V_i * (P - 1) / (2.303 * R * TK);
*/
{
	//int i, i1;
	LDBLE T_c, P_c;
	LDBLE A, B, B_r, b2, kk, oo, a_aa, T_r;
	LDBLE m_sum, b_sum, a_aa_sum, a_aa_sum2;
	LDBLE phi;
	LDBLE R_TK, R = R_LITER_ATM; /* L atm / (K mol) */
	LDBLE r3[4], r3_12, rp, rp3, rq, rz, ri, ri1, one_3 = 0.33333333333333333;
	LDBLE disct, vinit, v1, ddp, dp_dv, dp_dv2;
	int it;
	struct phase *phase_ptr;
	LDBLE V_m = 0, P = 0;

	LDBLE TK = tk_x;
	R_TK = R * TK;
	m_sum = b_sum = a_aa_sum = 0.0;
	size_t i;

	if (gas_unknowns.size() == 0)
	{
		error_msg("No gas unknowns.", STOP);
	}
	struct gas_phase * gas_phase_ptr = use.gas_phase_ptr;
	
	for (i = 0; i < gas_unknowns.size(); i++)
	{
		m_sum += gas_unknowns[i]->moles;
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->t_c == 0.0 || phase_ptr->p_c == 0.0)
			continue;
		if (!phase_ptr->pr_a)
		{
			T_c = phase_ptr->t_c;
			P_c = phase_ptr->p_c;
			phase_ptr->pr_a = 0.457235 * R * R * T_c * T_c / P_c;
			phase_ptr->pr_b = 0.077796 * R * T_c / P_c;
			T_r = TK / T_c;
			oo = phase_ptr->omega;
			kk = 0.37464 + oo * (1.54226 - 0.26992 * oo);
			phase_ptr->pr_alpha = pow(1 + kk * (1 - sqrt(T_r)), 2);
			phase_ptr->pr_tk = TK;
			phase_ptr->pr_in = true;
		}
		if (phase_ptr->pr_tk != TK)
		{
			T_r = TK / phase_ptr->t_c;
			oo = phase_ptr->omega;
			kk = 0.37464 + oo * (1.54226 - 0.26992 * oo);
			phase_ptr->pr_alpha = pow(1 + kk * (1 - sqrt(T_r)), 2);
			phase_ptr->pr_tk = TK;
			phase_ptr->pr_in = true;
		}
	}
	if (m_sum == 0)
			return (OK);
	gas_phase_ptr->v_m = use.gas_phase_ptr->volume / m_sum;
	for (i = 0; i < gas_unknowns.size(); i++)
	{
		phase_ptr = gas_unknowns[i]->phase;
		phase_ptr->fraction_x = gas_unknowns[i]->moles / m_sum;							// phase_ptr->fraction_x updated
	}

	for (i = 0; i < gas_unknowns.size(); i++)
	{
		a_aa_sum2 = 0.0;
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->t_c == 0.0 || phase_ptr->p_c == 0.0)
			continue;
		b_sum += phase_ptr->fraction_x * phase_ptr->pr_b;
		size_t i1;
		struct phase *phase_ptr1;
		for (i1 = 0; i1 <  gas_unknowns.size(); i1++)
		{
			phase_ptr1 = gas_unknowns[i1]->phase;
			if (phase_ptr1->t_c == 0.0 || phase_ptr1->p_c == 0.0)
				continue;
			if (phase_ptr1->fraction_x == 0)
				continue;
			a_aa = sqrt(phase_ptr->pr_a * phase_ptr->pr_alpha *
				        phase_ptr1->pr_a * phase_ptr1->pr_alpha);
			a_aa_sum += phase_ptr->fraction_x * phase_ptr1->fraction_x * a_aa;
			a_aa_sum2 += phase_ptr1->fraction_x * a_aa;
		}
		phase_ptr->pr_aa_sum2 = a_aa_sum2;
	}
	b2 = b_sum * b_sum;

	if (gas_phase_ptr->type == VOLUME)
	{
		V_m = gas_phase_ptr->volume / m_sum;
		P = R_TK / (V_m - b_sum) - a_aa_sum / (V_m * (V_m + 2 * b_sum) - b2);
		if (P < 150)
		{
			// check for 3-roots...
			r3[1] = b_sum - R_TK / P;
			r3[2] = -3.0 * b2 + (a_aa_sum - R_TK * 2.0 * b_sum) / P;
			r3[3] = b2 * b_sum + (R_TK * b2 - b_sum * a_aa_sum) / P;
			// the discriminant of the cubic eqn...
			disct = 18. * r3[1] * r3[2] * r3[3] -
				4. * pow(r3[1], 3) * r3[3] + 
				r3[1] * r3[1] * r3[2] * r3[2] -
				4. * pow(r3[2], 3) - 
				27. * r3[3] * r3[3];
			if (disct > 0)
			{
				// 3-roots, find the largest P...
				it = 0;
				ddp = 1e-9;
				v1 = vinit = 0.4;
				dp_dv = -R_TK / ((v1 - b_sum) * (v1 - b_sum)) +
					a_aa_sum * (2 * v1 + 2 * b_sum) /
					pow((v1 * v1 + 2. * b_sum * v1 - b2), 2);
				while (fabs(dp_dv) > 1e-11 && it < 40)
				{
					it +=1;
					v1 -= ddp;
					dp_dv2 = -R_TK / ((v1 - b_sum) * (v1 - b_sum)) +
						a_aa_sum * (2 * v1 + 2 * b_sum) /
						pow((v1 * v1 + 2. * b_sum * v1 - b2), 2);
					v1 -= (dp_dv * ddp / (dp_dv - dp_dv2) - ddp);
					if (v1 > vinit || v1 < 0.03)
					{
						//if (vinit < 0.1)
						//	vinit -= 0.01;
						//else 
							vinit -= 0.05;
						if (vinit < 0.06) // 0.01
							it = 40;
						v1 = vinit;
					}
					dp_dv = -R_TK / ((v1 - b_sum) * (v1 - b_sum)) +
						a_aa_sum * (2 * v1 + 2 * b_sum) /
						pow((v1 * v1 + 2. * b_sum * v1 - b2), 2);
				}
				if (it == 40)
				{
// accept a (possible) whobble in the curve...
//					error_msg("No convergence when calculating P in Peng-Robinson.", STOP);
				}
				else if (V_m < v1)
					P = R_TK / (v1 - b_sum) - a_aa_sum / (v1 * (v1 + 2 * b_sum) - b2);
			}
		}
		gas_phase_ptr->total_p = P;												// phase_ptr->total_p updated
	}
	else
	{
		P = gas_phase_ptr->total_p;
		r3[1] = b_sum - R_TK / P;
		r3_12 = r3[1] * r3[1];
		r3[2] = -3.0 * b2 + (a_aa_sum - R_TK * 2.0 * b_sum) / P;
		r3[3] = b2 * b_sum + (R_TK * b2 - b_sum * a_aa_sum) / P;
		// solve t^3 + rp*t + rq = 0.
		// molar volume V_m = t - r3[1] / 3... 
		rp = r3[2] - r3_12 / 3;
		rp3 = rp * rp * rp;
		rq = (2.0 * r3_12 * r3[1] - 9.0 * r3[1] * r3[2]) / 27 + r3[3];
		rz = rq * rq / 4 + rp3 / 27;
		if (rz >= 0) // Cardono's method...
		{
			ri = sqrt(rz);
			if (ri + rq / 2 <= 0)
			{
				V_m = pow(ri - rq / 2, one_3) + pow(- ri - rq / 2, one_3) - r3[1] / 3;
			} else
			{
				ri = - pow(ri + rq / 2, one_3);
				V_m = ri - rp / (3.0 * ri) - r3[1] / 3;
			}
		} else // use complex plane...
		{
			ri = sqrt(- rp3 / 27); // rp < 0
			ri1 = acos(- rq / 2 / ri);
			V_m = 2.0 * pow(ri, one_3) * cos(ri1 / 3) - r3[1] / 3;
		}
		gas_phase_ptr->v_m = V_m;												 // phase_ptr->fraction_x updated
	}
 // calculate the fugacity coefficients...
	for (i = 0; i < gas_unknowns.size(); i++)
	{
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->fraction_x == 0.0)
		{
			phase_ptr->pr_p = 0;
			phase_ptr->pr_phi = 1;
			phase_ptr->pr_si_f = 0.0;
			continue;
		}	
		phase_ptr->pr_p = phase_ptr->fraction_x * P;

		if (phase_ptr->t_c == 0.0 || phase_ptr->p_c == 0.0)
		{
			phase_ptr->pr_phi = 1;
			continue;
		}
		rz = P * V_m / R_TK;
		A = a_aa_sum * P / (R_TK * R_TK);
		B = b_sum * P / R_TK;
		B_r = phase_ptr->pr_b / b_sum;
		if (rz > B)
			phi = B_r * (rz - 1) - log(rz - B) + A / (2.828427 * B) * (B_r - 2.0 * phase_ptr->pr_aa_sum2 / a_aa_sum) *
				  log((rz + 2.41421356 * B) / (rz - 0.41421356 * B));
		else
			phi = -3.0; // fugacity coefficient > 0.05
		phase_ptr->pr_phi = exp(phi);
		phase_ptr->pr_si_f = phi / LOG_10;										// pr_si_f updated
		phase_ptr->pr_in = true;
	}
	return (V_m);
}
/* ---------------------------------------------------------------------- */
int Phreeqc::
calc_gas_pressures(void)
/* ---------------------------------------------------------------------- */
{
	int n_g = 0;
	LDBLE lp, V_m = 0;
	struct rxn_token *rxn_ptr;
	struct phase *phase_ptr;
	bool PR = false, pr_done = false;
	size_t i;
/*
 *   moles and partial pressures for gases
 */
	if (use.gas_phase_ptr == NULL)
		return (OK);

	use.gas_phase_ptr->total_moles = 0;

	for (i = 0; i < gas_unknowns.size(); i++)
	{
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->in == TRUE)
		{
			if (!PR && phase_ptr->t_c > 0 && phase_ptr->p_c > 0)
				PR = true;
			n_g++;
		}
		use.gas_phase_ptr->total_moles += gas_unknowns[i]->moles;
	}
	if (use.gas_phase_ptr->type == PRESSURE)
	{
			//calc_PR(phase_ptrs, n_g, gas_unknown->gas_phase->total_p, tk_x, 0);
		calc_PR();
	} else
	{
		if (PR && use.gas_phase_ptr->total_moles > 0)
		{
			//V_m = use.gas_phase_ptr->volume / use.gas_phase_ptr->total_moles;
			//V_m = use.gas_phase_ptr->volume / gas_unknown->moles;
			//if (V_m < 0.035)
			//{
			//	V_m = 0.035;
			//} else if (V_m > 1e4)
			//{
			//	V_m = 1e4;
			//}
			/* need to warn about minimal V_m and maximal P... */
			//if (use.gas_phase_ptr->v_m > 0.035)
			//if (V_m > 0.035)
			//{
			//	if (V_m < 0.07)
			//		V_m = (3. * use.gas_phase_ptr->v_m + V_m) / 4;
			//	else
			//		V_m = (1. * use.gas_phase_ptr->v_m + V_m) / 2;
			//}
			//if (iterations > 100)
			//{
			//	V_m *= 1.0; /* debug */
			//}
			//calc_PR(phase_ptrs, n_g, 0, tk_x, V_m);
			calc_PR();
			pr_done = true;
			use.gas_phase_ptr->total_moles = 0;
			//if (fabs(use.gas_phase_ptr->total_p - patm_x) > 0.01)
			//{
			//	same_pressure = FALSE;
			//	patm_x = use.gas_phase_ptr->total_p;
			//	k_temp(tc_x, patm_x);
			//}
		} else
		{
			use.gas_phase_ptr->total_p = 0;
			use.gas_phase_ptr->total_moles = 0;
		}
	}

	for (i = 0; i < gas_unknowns.size(); i++)
	{
		//gas_comp_ptr = &(use.gas_phase_ptr->comps[i]);
		phase_ptr = gas_unknowns[i]->phase;
		if (phase_ptr->in == TRUE)
		{
			lp = -phase_ptr->lk;
			for (rxn_ptr = phase_ptr->rxn_x->token + 1; rxn_ptr->s != NULL;
				 rxn_ptr++)
			{
				lp += rxn_ptr->s->la * rxn_ptr->coef;
			}
			phase_ptr->p_soln_x = exp(LOG_10 * (lp - phase_ptr->pr_si_f));

			if (use.gas_phase_ptr->type == PRESSURE)
			{
				phase_ptr->moles_x = phase_ptr->p_soln_x *
					use.gas_phase_ptr->total_moles / use.gas_phase_ptr->total_p;
				phase_ptr->fraction_x =	phase_ptr->moles_x / use.gas_phase_ptr->total_moles;
			}
			//else
			{
				if (pr_done)
				{
					//lp = phase_ptr->p_soln_x / use.gas_phase_ptr->total_p *
					//	use.gas_phase_ptr->volume / V_m;
					lp = phase_ptr->p_soln_x / use.gas_phase_ptr->total_p *
						use.gas_phase_ptr->volume / use.gas_phase_ptr->v_m;
					phase_ptr->moles_x = lp;
				}
				else
				{
					phase_ptr->moles_x = phase_ptr->p_soln_x *
						use.gas_phase_ptr->volume / (R_LITER_ATM * tk_x);
					use.gas_phase_ptr->total_p += phase_ptr->p_soln_x;
				}
				use.gas_phase_ptr->total_moles += phase_ptr->moles_x;
			}
		}
		else
		{
			phase_ptr->moles_x = 0;
			phase_ptr->fraction_x = 0;
		}
	}

	if (use.gas_phase_ptr->type == VOLUME && !PR)
	{
		/*
		 * Fixed-volume gas phase reacting with a solution
		 * Change pressure used in logK to pressure of gas phase
		 */
		//if (fabs(use.gas_phase_ptr->total_p - patm_x) > 0.01)
		//{
		//	same_pressure = FALSE;
		//	patm_x = use.gas_phase_ptr->total_p;
		//	k_temp(tc_x, patm_x);
		//}
	}


	//delete phase_ptrs;
	return (OK);
}
#endif