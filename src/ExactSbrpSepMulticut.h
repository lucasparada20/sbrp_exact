#ifndef EXACT_SBRP_SEP_H_MULTICUT
#define EXACT_SBRP_SEP_H_MULTICUT

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
#include "RouteFeasibility.h"
#include "FirstStageDuals.h"

class ExactSbrpSepMulticut
{
	public:
		ExactSbrpSepMulticut(IloEnv env, ExactSbrpGraphO * graph, IloNumVarArray x, IloNumVar theta, IloNumVarArray thetas, FirstStageDuals *Duals);
		ExactSbrpSepMulticut(IloEnv env, ExactSbrpGraphO * graph, IloNumVarArray x, IloNumVar theta, IloNumVarArray thetas);

		void SeperateFrac(IloRangeArray array);
		void SeperateInt(IloRangeArray array);

		bool SeperateUserCapRecursive(IloRangeArray array);    
		void SeperateBendersCut(IloRangeArray array);

		double best_sol;
		double best_sol_recourse;
		double best_sol_distance;
		std::vector< std::vector<Node*> > best_solution;
		int nb_sub_tours;
		int nb_sub_tours_from_frac;
		int nb_benders_cuts;
		int nb_benders_feasibility_cuts;
    
		std::vector<int> _component;

	private:
 
		void SeperateSubTourInequalities(IloRangeArray array);

		void ResearchConnectedComponent(Node* n, int comp);
		void ResearchDepotComponent(Node* n, int comp);

		//return the number of added constraints
		int TestAndAddFeasBendersCut(IloRangeArray array);

		FirstStageDuals * _Dual; //For Benders
	
		IloEnv _env;
		ExactSbrpGraphO * _graph;
		IloNumVarArray _x;
		IloNumVar _theta;
		IloNumVarArray _thetas;

		Prob* _prob;

		int nbseps;

		std::vector<Node*> _path;
		
     
};

#endif