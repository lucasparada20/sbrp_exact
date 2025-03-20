
#ifndef EXACT_SBRP_CALL_BACK_H
#define EXACT_SBRP_CALL_BACK_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <cstring>
#include <ilcplex/ilocplex.h>
#include <algorithm>
#include <vector>

#include "ExactSbrpGraph.h"
#include "ExactSbrpSep.h"


class ExactSbrpLazyCallBackO : public IloCplex::LazyConstraintCallbackI
{
	public:
		ExactSbrpLazyCallBackO(IloEnv env, ExactSbrpGraphO * graph, IloNumVarArray x, IloNumVarArray thetas, ExactSbrpSepO * sep);
		~ExactSbrpLazyCallBackO();

		IloCplex::CallbackI *duplicateCallback() const
		{
        	return new (getEnv()) ExactSbrpLazyCallBackO(*this);
        }

        void main();

    private:
    	ExactSbrpGraphO * _graph;
    	IloNumVarArray _x;
		IloNumVarArray _thetas;
    	ExactSbrpSepO * _sep;
    public:
    	IloRangeArray added_constraints;
		bool add_constraints;
};

class ExactSbrpUserCutCallBackO : public IloCplex::UserCutCallbackI
{
	public:
		ExactSbrpUserCutCallBackO(IloEnv env, ExactSbrpGraphO * graph, IloNumVarArray x, IloNumVarArray thetas, ExactSbrpSepO * sep);
		~ExactSbrpUserCutCallBackO();

		IloCplex::CallbackI *duplicateCallback() const
		{
        	return new (getEnv()) ExactSbrpUserCutCallBackO(*this);
        }

        void main();

	private:
    	ExactSbrpGraphO * _graph;
    	IloNumVarArray _x;
		IloNumVarArray _thetas;
    	ExactSbrpSepO * _sep;
    public:
    	IloRangeArray added_constraints;
		bool add_constraints;
};

#endif
