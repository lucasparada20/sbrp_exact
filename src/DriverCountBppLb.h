



#ifndef DRIVERCOUNTBPPLB
#define DRIVERCOUNTBPPLB

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ilcplex/ilocplex.h>
#include "Solution.h"
#include "ProblemDefinition.h"
#include "NodeSBRP.h"
#include "DriverSBRP.h"
#include "time.h"
#include <vector>

class DriverCountBppColumn
{
public:
	int id;
	std::vector<int> nodes;
	double value;

	void Add(int n)
	{
		nodes.push_back(n);
	}

	bool Contains(int n)
	{
		for(size_t i=0;i<nodes.size();i++)
			if(nodes[i] == n)
				return true;
		return false;
	}
	void Show()
	{
		printf("Column:%d stations:", id);
		for(size_t i=0;i<nodes.size();i++)
			printf("%d-", nodes[i]);
		printf("\n");
	}
};

class DriverCountBppLbLazyCallBack : public IloCplex::LazyConstraintCallbackI
{
	public:
		DriverCountBppLbLazyCallBack(IloEnv env, Prob * prob, IloNumVarArray z);
		~DriverCountBppLbLazyCallBack(){}

		void main();
		IloCplex::CallbackI *duplicateCallback() const{return new (getEnv()) DriverCountBppLbLazyCallBack(*this);}
    private:
    	Prob * _prob;
    	IloNumVarArray _z;
};

class DriverCountBppLb
{
	public:
		DriverCountBppLb(){}
		~DriverCountBppLb(){}

		int GetLB(Sol * sol);
		void InitMasterProblem(IloEnv env);
		void InitPricingProblem(IloEnv env);
		bool ExactPricing();

		double taken_time;
		int lb;
	private:
		Sol * _sol;
		Prob * _prob;
		std::vector<Node*> _nodes;
		int n;

		std::vector<DriverCountBppColumn> columns;

		double pricing_cost;

		//variables for the master problem
		IloModel model;
		IloObjective obj;
		IloCplex cplex;
		IloNumVarArray x;
		IloRangeArray consts_vars;		//constraints for the items
		IloNumArray duals_vars;

		//variables for the pricing problem
		IloModel model2;
		IloObjective obj2;
		IloCplex cplex2;
		IloNumVarArray z;
};



#endif
