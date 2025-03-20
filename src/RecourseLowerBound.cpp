#include "RecourseLowerBound.h"
#include "RouteFeasibility.h"
#include "Parameters.h"
#include <algorithm>    // std::reverse
#include <numeric> // std::accumulate

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <cmath>

double RecourseLowerBound::Calculate(Prob* prob)
{	
	std::vector<Node*> stations;
	stations.push_back(prob->GetNode(prob->GetNodeCount() - 2));
	for (int i = 0; i < prob->GetCustomerCount(); i++)
		stations.push_back(prob->GetCustomer(i));
	stations.push_back(prob->GetNode(prob->GetNodeCount() - 1));

	double lb = Calculate(prob, stations);

	//printf("L1:%.3lf Cost:%.2lf time:%.3lf\n",l1,lb, time_taken_l1);

	return lb;
}

double RecourseLowerBound::Calculate(Prob* prob, std::vector<Node*>& stations)
{
	feasible = true;
	time_taken_l1 = time_taken_l2 = 0;
	l1 = l2 = 0;
	inf_scenario = -1;

	l1 = 0;

	for (int e = 0;e < prob->GetScenarioCount(); e++)
	{
		double l = Calculate(prob, stations, e);
    
		//printf("scenario:%d cost:%.1lf\n", e,l);
		if(l >= 9999999999)
		{
			l1 = 9999999999;
			feasible = false;
			inf_scenario = e;
			break;
		}
		l1 += l;
   
	}
	l1 /= prob->GetScenarioCount();
	double cost = l1 * Parameters::GetCminEpsilon();

	return cost;

}

double RecourseLowerBound::Calculate(Prob* prob, std::vector<Node*>& stations, int scenario)
{
	int nb_iter = 1;
	int Q = prob->GetDriver(0)->capacity;
	int dmds = 0;
	int lb_left = 0, ub_left = 0, lb_right = 0, ub_right = 0;
	for (size_t l = 0; l < stations.size(); l++)
	{
		if(stations[l]->type != NODE_TYPE_CUSTOMER) continue;

		Node * n = stations[l];
		int dmd = n->demands[scenario];
		dmds += dmd;

		if(dmd > 0) { lb_left += dmd; ub_left += dmd; }
		ub_left += n->w_minus;

		if(dmd < 0) { lb_right -= dmd; ub_right -= dmd; }
		ub_right += n->w_plus;
	}

	while(ub_left+Q < lb_right || lb_left > ub_right+Q)
	{
		Q += prob->GetDriver(0)->capacity;
		nb_iter++;
	}

	//printf("ub_left:%d < lb_right:%d || lb_left:%d > ub_right:%d\n", ub_left, lb_right,lb_left , ub_right);
	if(ub_left+Q  < lb_right || lb_left > ub_right+Q)
		return 9999999999;
	else
		return std::max(0, std::abs(dmds) - Q);
}

double RecourseLowerBound::CalculateWithMinDriverCount(Prob* prob)
{
	std::vector<Node*> stations;
	for (int i = 0; i < prob->GetCustomerCount(); i++)
		stations.push_back(prob->GetCustomer(i));
	return CalculateWithMinDriverCount(prob, stations);
}
double RecourseLowerBound::CalculateWithMinDriverCount(Prob* prob, std::vector<Node*>& stations)
{
	//first sum the demands
	std::vector<int> dmds(prob->GetScenarioCount(), 0);
	std::vector<int> lb_left(prob->GetScenarioCount(), 0);
	std::vector<int> ub_left(prob->GetScenarioCount(), 0);
	std::vector<int> lb_right(prob->GetScenarioCount(), 0);
	std::vector<int> ub_right(prob->GetScenarioCount(), 0);
	for (int e = 0;e < prob->GetScenarioCount(); e++)
		for (size_t l = 0; l < stations.size(); l++)
		{
			if(stations[l]->type != NODE_TYPE_CUSTOMER) continue;

			Node * n = stations[l];
			int dmd = n->demands[e];
			dmds[e] += dmd;

			if(dmd > 0) { lb_left[e] += dmd; ub_left[e] += dmd; }
			ub_left[e] += n->w_minus;

			if(dmd < 0) { lb_right[e] -= dmd; ub_right[e] -= dmd; }
			ub_right[e] += n->w_plus;
		}

	//get the minimum number of vehicles
	int Q = prob->GetDriver(0)->capacity;
	int nb_drivers = 1;
	for (int e = 0;e < prob->GetScenarioCount(); e++)
		while(ub_left[e]+Q*nb_drivers < lb_right[e] || lb_left[e] > ub_right[e]+Q*nb_drivers)
			nb_drivers++;

	double lb = 0;
	for (int e = 0;e < prob->GetScenarioCount(); e++)
		lb += std::max(0, std::abs(dmds[e]) - Q*nb_drivers);
	lb /= prob->GetScenarioCount();
	lb *= Parameters::GetCminEpsilon();

	return lb;
}


int RecourseLowerBound::GetDriverCount(Prob* prob)
{	
	std::vector<Node*> stations;
	for (int i = 0; i < prob->GetCustomerCount(); i++)
		stations.push_back(prob->GetCustomer(i));
	return GetDriverCount(prob, stations);
}

int RecourseLowerBound::GetDriverCount(Prob* prob, std::vector<Node*>& stations)
{
	int nb_drivers=1;
	for (int e=0;e<prob->GetScenarioCount();e++)
	{
		int Q = prob->GetDriver(0)->capacity * nb_drivers;
		int dmds = 0;
		int lb_left = 0, ub_left = 0, lb_right = 0, ub_right = 0;
		for (size_t l = 0; l < stations.size(); l++)
		{
			if(stations[l]->type != NODE_TYPE_CUSTOMER) continue;
			Node * n = stations[l];
			
			//n->Show();
			
			int dmd = n->demands[e];
			dmds += dmd;

			if(dmd > 0) { lb_left += dmd; ub_left += dmd; }
			ub_left += n->w_minus;

			if(dmd < 0) { lb_right -= dmd; ub_right -= dmd; }
			ub_right += n->w_plus;
		}

		while(ub_left+Q < lb_right || lb_left > ub_right+Q)
		{
			Q += prob->GetDriver(0)->capacity;
			nb_drivers++;
		}
	}
	return nb_drivers;
}

void RecourseLowerBound::SetWorstScenario(Prob* prob){
	std::vector<Node*> stations;
	stations.push_back(prob->GetNode(prob->GetNodeCount() - 2));
	for (int i = 0; i < prob->GetCustomerCount(); i++)
		stations.push_back(prob->GetCustomer(i));
	stations.push_back(prob->GetNode(prob->GetNodeCount() - 1));

	double worst_cost = -1;
	int worst_scenario = -1;
	for (int e=0;e<prob->GetScenarioCount();e++)
	{
		double cost = Calculate(prob, stations, e);
		if(cost > worst_cost)
		{
			worst_cost = cost;
			worst_scenario = e;
		}
	}
	Parameters::SetWorstScenario(worst_scenario);
	printf("worst_scenario:%d cost:%.1lf\n",worst_scenario,worst_cost);

}

void RecourseLowerBound::SortFromWorstScenarios(Prob* prob)
{
	std::vector<Node*> stations;
	stations.push_back(prob->GetNode(prob->GetNodeCount() - 2));
	for (int i = 0; i < prob->GetCustomerCount(); i++)
		stations.push_back(prob->GetCustomer(i));
	stations.push_back(prob->GetNode(prob->GetNodeCount() - 1));

	RecourseLowerBound rec;

	std::vector<LBandIndex> Container;
	std::vector<double> ScenarioL(prob->GetScenarioCount(), -1);
	for (size_t e = 0; e < ScenarioL.size(); e++)
	{
		LBandIndex element;
		element.Index = e;
		element.LB = ScenarioL[e] = rec.Calculate(prob, stations, e);
		Container.push_back(element);
	}
	
	std::vector<double> aux(prob->GetScenarioCount(), -1);
	for (size_t e = 0; e < Container.size(); e++)
		aux[e] = Container[e].LB;

	std::sort(aux.begin(), aux.end());
	std::reverse(aux.begin(), aux.end());
	
	for (size_t e = 0; e < aux.size(); e++)
	{
		for (size_t i = 0; i < Container.size(); i++)
		{
			if (aux[e] == Container[i].LB)
				if (Container[i].AlreadyAdded == false)
				{
					prob->AddSortedScenarios(Container[i].Index);
					Container[i].AlreadyAdded = true;
					break;
				}
		}

	}
	
}

double RecourseLowerBound::CalculateLP(Prob* prob, std::vector<Node*> & stations) //LP for a route
{
	IloEnv env;
	IloModel model(env);

	IloArray<IloNumVarArray> flow(env, prob->GetScenarioCount());
	IloArray<IloNumVarArray> zeta(env, stations.size()); //2D Variable for z, for fixed scenario.
	IloArray<IloNumVarArray> w_minus(env, prob->GetScenarioCount());
	IloArray<IloNumVarArray> w_plus(env, prob->GetScenarioCount());

	//V A R I A B L E S

	// Definition of z_{ik} variables
	for (size_t i = 0; i < stations.size(); i++)
	{
		zeta[i] = IloNumVarArray(env, stations.size(), 0, 1, ILOFLOAT); //z_{eik}=1 if station i, visited in position k and scenario e.

		for (size_t k = 0; k < stations.size(); k++)
		{
			char name[40];
			sprintf(name, "z%d_%d", (int)i, (int)k);
			zeta[i][k].setName(name);
		}
	}

	for (int e = 0; e < prob->GetScenarioCount(); e++)
	{
		w_minus[e] = IloNumVarArray(env, stations.size(), 0, IloInfinity, ILOFLOAT);
		w_plus[e] = IloNumVarArray(env, stations.size(), 0, IloInfinity, ILOFLOAT);

		for (size_t k = 0; k < stations.size(); k++)
		{
			char name[40];

			sprintf(name, "wm%d_%d", e, (int)k);
			w_minus[e][k].setName(name);

			sprintf(name, "wp%d_%d", e, (int)k);
			w_plus[e][k].setName(name);
		}
	}


	//Definition of Flow Variables
	for (int e = 0; e < prob->GetScenarioCount(); e++)
	{
		flow[e] = IloNumVarArray(env, stations.size(), 0, prob->GetDriver(0)->capacity, ILOFLOAT);
		for (size_t k = 0; k < stations.size(); k++)
		{
			char name[40];
			sprintf(name, "f%d_%d", e, (int)k);
			flow[e][k].setName(name);
		}
	}


	//C O N S T R A I N T S

	//Upper Bounds on W variables
	// w^{+}_{k} \leq \sum^{n}_{k=1} \bar w^{+}_{i}*z_{ik}, for each scenario.
	// w^{-}_{k} \leq \sum^{n}_{k=1} \bar w^{-}_{i}*z_{ik}, for each scenario.

	for (int e = 0; e < prob->GetScenarioCount(); e++)
	{
		for (size_t k = 1; k + 1 < stations.size(); k++)
		{
			IloExpr expr1(env);
			expr1 += w_plus[e][k];
			for (size_t i = 1; i + 1 < stations.size(); i++)
				expr1 -= stations[i]->w_plus * zeta[i][k];
			model.add(expr1 <= 0);
			expr1.end();

			IloExpr expr(env);
			expr += w_minus[e][k];
			for (size_t i = 1; i + 1 < stations.size(); i++)
				expr -= stations[i]->w_minus * zeta[i][k];
			model.add(expr <= 0);
			expr.end();
		}
	}



	//\sum^{n}_{k=1} z_{ik} = 1 , \forall i \in \{1,...,n\}
	//One z for each station.
	for (size_t i = 1; i + 1 < stations.size(); i++)
	{
		IloExpr expr(env);
		for (size_t k = 1; k + 1 < stations.size(); k++)
			expr += zeta[i][k];
		model.add(expr == 1);
		expr.end();
	}


	//\sum^{n}_{k=1} z_{ik} = 1 , \forall k \in \{1,...,n\}
	//One z for each position of the sequence.
	for (size_t k = 1; k + 1 < stations.size(); k++)
	{
		IloExpr expr(env);
		for (size_t i = 1; i + 1 < stations.size(); i++)
			expr += zeta[i][k];
		model.add(expr == 1);
		expr.end();
	}


	//Flow Conservation on the sequence
	for (int e = 0; e < prob->GetScenarioCount(); e++)
	{
		for (size_t k = 1; k + 1 < stations.size(); k++)
		{
			IloExpr expr(env);
			expr += flow[e][k];
			expr -= flow[e][k - 1];
			for (size_t i = 1; i + 1 < stations.size(); i++)
				expr -= stations[i]->demands[e] * zeta[i][k];
			expr += w_plus[e][k];
			expr -= w_minus[e][k];
			model.add(expr == 0);
			expr.end();
		}
	}



	//std::cout << "C_{min} is: " << Parameters::GetCmin() << std::endl;

	IloExpr recourse(env);
	for (int e = 0; e < prob->GetScenarioCount(); e++)
		for (size_t i = 1; i + 1 < stations.size(); i++)
			recourse += (w_minus[e][i] + w_plus[e][i]);
	model.add(IloMinimize(env, recourse));
	recourse.end();


	//std::cout << model << std::endl;

	IloCplex cplex(model);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	cplex.setParam(IloCplex::Param::Threads, 1);
	//cplex.exportModel("recourse.lp");
	//cplex.setParam(cplex.TiLim, 300);

	bool re = cplex.solve();

	double sol = re ? cplex.getObjValue() : 9999999999;
	printf("Recourse Cost:%.lf\n", sol);

	//Print the Sequence
	/*for (int i = 1; i + 1 < stations.size(); i++)
	{
		for (int k = 1; k + 1 < stations.size(); k++)
			printf("%.2lf ", (double)cplex.getValue(zeta[i][k]) + 0.000001);
		printf("\n");
	}*/


	env.end();
	return sol;
}

