#ifndef EXACT_SBRP_MULTICUT_CALL_BACK_H
#define EXACT_SBRP_MUTICUT_CALL_BACK_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <vector>

#include "ExactSbrpGraph.h"
#include "ExactSbrpSepMulticut.h"


class ExactSbrpMulticutLazyCallBack : public IloCplex::LazyConstraintCallbackI
{
	public:
		ExactSbrpMulticutLazyCallBack(IloEnv env, ExactSbrpGraphO * graph, IloNumVarArray x, IloNumVarArray thetas, ExactSbrpSepMulticut * sep);
		~ExactSbrpMulticutLazyCallBack();

		IloCplex::CallbackI *duplicateCallback() const
		{
        	return new (getEnv()) ExactSbrpMulticutLazyCallBack(*this);
        }

        void main();

    private:
    	ExactSbrpGraphO * _graph;
    	IloNumVarArray _x;
		IloNumVarArray _thetas;
    	ExactSbrpSepMulticut * _sep;
    public:
    	IloRangeArray added_constraints;
		bool add_constraints;
};

class ExactSbrpMulticutUserCutCallBack : public IloCplex::UserCutCallbackI
{
	public:
		ExactSbrpMulticutUserCutCallBack(IloEnv env, ExactSbrpGraphO * graph, IloNumVarArray x, IloNumVarArray thetas, ExactSbrpSepMulticut * sep);
		~ExactSbrpMulticutUserCutCallBack();

		IloCplex::CallbackI *duplicateCallback() const
		{
        	return new (getEnv()) ExactSbrpMulticutUserCutCallBack(*this);
        }

        void main();

	private:
    	ExactSbrpGraphO * _graph;
    	IloNumVarArray _x;
		IloNumVarArray _thetas;
    	ExactSbrpSepMulticut * _sep;
    public:
    	IloRangeArray added_constraints;
		bool add_constraints;
};

#endif
