#ifndef EXACT_SBRP_ORIENTED_H_MULTICUT
#define EXACT_SBRP_ORIENTED_H_MULTICUT

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
#include "ExactSbrpSepMulticut.h"
#include "ExactSbrpCallBackMulticut.h"

#include "RecourseLowerBound.h"
#include "FirstStageDuals.h"
#include <ctime>

class ExactSbrpMulticut
{
	public:
		ExactSbrpMulticut()
		{
			max_time = 300;
			lazy_call = NULL;
			usercut_call = NULL;
		}
		~ExactSbrpMulticut(){}


		void Solve(Prob * prob);

		int max_time;
		double time_taken;
		clock_t start_time;
		int status;
		double lb;
		double ub;
		double ub_recourse;
		double ub_distance;
		int nb_sub_tours;
		int nb_sub_tours_from_frac;
		int nb_benders_cuts;
		int nb_benders_feasibility_cuts;


    	//To print to re_.txt file
		double cplex_distances;
		double cplex_recourse;
		double cplex_relative_gap;

	private:
		void Init(IloEnv env);
		void SolveProblem(IloEnv env);
		void Clear();

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
		
		IloNumVarArray thetas_multicut;

		ExactSbrpSepMulticut * _sep;
		ExactSbrpMulticutLazyCallBack * lazy_call;
		ExactSbrpMulticutUserCutCallBack * usercut_call;
		
		IloRangeArray constraints; //Container of constraints for Post-Mortem analysis ..

};




#endif