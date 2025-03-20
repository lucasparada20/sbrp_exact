
#ifndef EXACT_SBRP_H
#define EXACT_SBRP_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <vector>
#include "NodeSBRP.h"
#include "DriverSBRP.h"
#include "ProblemDefinition.h"
#include "ExactSbrpGraph.h"
#include "ExactSbrpSep.h"
#include "ExactSbrpCallBacks.h"
#include "RecourseLowerBound.h"


class ExactSbrpO
{
	public:
		ExactSbrpO()
		{
			max_time = 300;
			lazy_call = NULL;
			usercut_call = NULL;
		}
		
		void Solve(Prob * prob);
		void Init(IloEnv env);
		void SolveProblem(IloEnv env);
		void Clear();
		void AddInfSetInq(IloEnv env);

		int max_time;
		double time_taken;
		double start_time;
		int status;
		double lb;
		double ub;
		double ub_recourse;
		double ub_distance;
		int nb_inf_sets;
		int nb_inf_paths;
		int nb_sub_tours;
		int nb_sub_tour_frac;
		int nb_l_cuts;
		int nb_p_cuts;
		int nb_frac_l_cuts;
		int nb_sorted_l_cuts;
		int nb_benders_cuts;

		double cplex_distances;
		double cplex_recourse;
		double cplex_relative_gap;
		bool solvedAtRoot;

	private:

		Prob * prob;
		ExactSbrpGraphO * graph;

		IloModel model;
		IloObjective obj_func;
		IloCplex cplex;
		IloNumVarArray x;

		IloArray<IloNumVarArray> wp; //w_plus
		IloArray<IloNumVarArray> wm; //w_minus

		IloNumVar theta; // Second Stage Recourse Cost
		IloNumVarArray thetas; // Variables for the dissagregation of recourse
		IloNumVar z; //Nb of vehicles
		

		ExactSbrpSepO * _sep;
		ExactSbrpLazyCallBackO * lazy_call;
		ExactSbrpUserCutCallBackO * usercut_call;
		
		IloRangeArray constraints; //Container of constraints for Post-Mortem 
};

#endif
