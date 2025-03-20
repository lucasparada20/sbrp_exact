
#ifndef EXACT_SBRP_SEP_H
#define EXACT_SBRP_SEP_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <vector>
#include <map>

#include "NodeSBRP.h"
#include "DriverSBRP.h"
#include "ProblemDefinition.h"
#include "ExactSbrpGraph.h"
#include "RouteFeasibility.h"
#include "Parameters.h"
#include "RecourseLowerBound.h"
#include "FirstStageDuals.h"


class SubPathSbrp
{	public:
	SubPathSbrp(){}
	double cost, sort, thetas;
	int from, to, nb;
	bool operator < (const SubPathSbrp & str) const {return (sort > str.sort);}
};

class ExactSbrpSepO
{
	public:
		
		ExactSbrpSepO(IloEnv env, ExactSbrpGraphO * graph, IloNumVarArray x, IloNumVar theta, IloNumVarArray thetas, FirstStageDuals *Duals);		
		ExactSbrpSepO(IloEnv env, ExactSbrpGraphO * graph, IloNumVarArray x, IloNumVar theta, IloNumVarArray thetas);

		void SeparateFrac(IloRangeArray array);
		void SeparateInt(IloRangeArray array);
		bool SeparateUserCapRecursive(IloRangeArray array);
		void SeparateDissagregatedRecourse(IloRangeArray array);		
		void SeparateBendersCut(IloRangeArray array);
		void SeparateSubTourInequalities(IloRangeArray array);
		void ResearchConnectedComponent(Node* n, int comp);
		void ResearchDepotComponent(Node* n, int comp);

		//These return the number of added constraints
		int TestAndAddInfeasiblePath(std::vector<Node*> & path, IloRangeArray array);
		int TestAndAddAllInfeasiblePath(std::vector<Node*> & path, IloRangeArray array);
		int TestAndAddSortedLCuts(std::vector<Node*> & path,IloRangeArray array);
		int TestAndAddFeasBendersCut(IloRangeArray array);		

		double best_sol;
		double best_sol_recourse;
		double best_sol_distance;
		std::vector< std::vector<Node*> > best_solution;
		
		//Counters for cuts
		int nb_opt_cuts;
		int nb_inf_sets;
		int nb_inf_paths;
		int nb_sub_tours;
		int nb_sub_tour_frac;
		int nb_l_cuts;
		int nb_p_cuts;
		int nb_frac_l_cuts;
		int nb_sorted_l_cuts;
		int nb_benders_cuts;
    
		std::vector<int> _component;

	private:
	 
		FirstStageDuals * _Dual; //For Benders
		
		IloEnv _env;
		ExactSbrpGraphO * _graph;
		IloNumVarArray _x;
		IloNumVar _theta;
		IloNumVarArray _thetas;

		Prob* _prob;

};

#endif