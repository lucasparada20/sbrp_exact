
#include "DriverCountBppLb.h"
#include <ilcplex/ilocplex.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>

#include "Parameters.h"
#include "RecourseLowerBound.h"

ILOSTLBEGIN

typedef IloArray<IloNumVarArray> NumVarMatrix;

int DriverCountBppLb::GetLB(Sol * sol)
{
	_sol = sol;
	_prob = sol->GetProb();
	n = _prob->GetCustomerCount();
	if(n == 0) return 0;
	//#ifdef WITH_OUTPUT
	//printf("DriverCountBppLb::GetLB() nb:%d\n", n);
	//#endif

	_nodes.clear();
	for(int i=0;i<n;i++)
		_nodes.push_back( _prob->GetCustomer(i) );

	clock_t begin = clock();
	IloEnv env;
	IloEnv env2;
  
	InitMasterProblem(env);
	InitPricingProblem(env2);
  
	//cplex.exportModel("DriverCountBppLb.lp");


	pricing_cost = 9999999999;
	bool re = false;
	int cont = 0;
	int cplex_status = 0;
	int iter = 0;
	int iter_max = 1000;
	bool has_timed_out = false;
	int best_lb2 = 0;
	do
	{
		cont = 0;
		re = cplex.solve();
		cplex_status = (int)cplex.getCplexStatus();
		double sol = cplex.getObjValue();
		double lb_cost= sol/pricing_cost;
		int lb1 = (int)ceil(sol-0.0001);
		int lb2 = (int)ceil(lb_cost-0.0001);
		best_lb2 = std::max(best_lb2, lb2);
		printf("Solved:%d status:%d iter:%d value:%.2lf ",re,cplex_status,iter,sol);
		printf("price:%lf lb_cost:%lf lb1:%d lb2:%d\n", pricing_cost, lb_cost,lb1,lb2);

		if(lb1 <= _prob->GetDriverCountLB()) break;
		if(lb2 >= lb1) break;
		
		clock_t now_time = clock();
		
		
		if(iter >= iter_max || (double)(now_time - begin)/CLOCKS_PER_SEC >= 1200)
		{
			has_timed_out = true;
			break;
		}

		cplex.getDuals(duals_vars, consts_vars);
		if(ExactPricing()) cont = 1;
		iter++;
	}
	while(cont == 1);


	lb = (int)ceil(cplex.getObjValue()-0.0001);
	//for(int i=0;i<x.getSize();i++)
	//	if(cplex.getValue(x[i]) > EPSILON)
		{
			//printf("Pattern:%d nb:%d v:%.3lf\t",i,nb++,(double)cplex.getValue(x[i]));
			//for(int j=0;j<patterns[i].size();j++)
			//	printf("%d ", patterns[i][j]);
			//printf("\n");
		}
	if(has_timed_out)
		lb = best_lb2;
	
	clock_t end_time = clock();
	
	taken_time = (double)(end_time - begin) / CLOCKS_PER_SEC;


	env.end();
	env2.end();

	//#ifdef WITH_OUTPUT
	printf("final bound:%d iter:%d time:%.4lf\n", lb,iter, taken_time);
	//#endif
	return lb;
}

void DriverCountBppLb::InitMasterProblem(IloEnv env)
{
	int i;
	model = IloModel(env);
	consts_vars = IloRangeArray(env);
	duals_vars = IloNumArray(env);

	int nb_to_add = 0;
	for(int i=0;i<_sol->GetDriverCount();i++)
		if(_sol->GetRouteLength(i) >= 2)
			nb_to_add++;
	x = IloNumVarArray(env,n+nb_to_add,0,1,ILOFLOAT);

	columns.resize(n);
	for(i=0;i<n;i++)
	{
		columns[i].id = i;
		columns[i].Add(i);
	}
	for(int i=0;i<_sol->GetDriverCount();i++)
		if(_sol->GetRouteLength(i) >= 2)
		{
			DriverCountBppColumn col;
			col.id = columns.size();

			Node * cur = _sol->GetNode( _sol->GetDriver(i)->StartNodeID );
			while(cur->type != NODE_TYPE_END_DEPOT)
			{
				if(cur->type == NODE_TYPE_CUSTOMER)
					col.Add(cur->no - 1);
				cur = _sol->Next[cur->id];
			}
			columns.push_back(col);
		}

	IloExpr obj1(env);
	for(i=0;i<x.getSize();i++)
		obj1 += x[i];

	obj = IloMinimize(env, obj1);
	model.add(obj);
	obj1.end();

	for(i=0;i<n;i++)
	{
		IloExpr expr(env);
		//expr += x[i];
		for(size_t j=0;j<columns.size();j++)
			if(columns[j].Contains(i))
				expr += x[j];
		consts_vars.add(expr >= 1);
	}

	model.add(consts_vars);

	cplex = IloCplex(model);
	cplex.setParam(IloCplex::Param::Threads,1);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
}
void DriverCountBppLb::InitPricingProblem(IloEnv env)
{
	model2 = IloModel(env);
	z = IloNumVarArray(env,n,0,1,ILOINT);

	IloExpr obj1(env);
	obj2 = IloMaximize(env, obj1);
	model2.add(obj2);
	obj1.end();

	for(int i=0;i<n;i++)
	{
		char name[20];
		sprintf(name,"z%d",i+1);
		z[i].setName(name);
	}

	for(int e = 0; e < _prob->GetScenarioCount(); e++)
	{
		IloExpr expr1(env);
		IloExpr expr2(env);
    
		for(int i = 0; i < n; i++)
		{
			//Delivery Knapsack
			if (_nodes[i]->demands[e] < 0)
				expr1 += std::abs(_nodes[i]->demands[e])* z[i];
			if (_nodes[i]->demands[e] > 0)
				expr1 -= _nodes[i]->demands[e] * z[i];
			expr1 -= _nodes[i]->w_minus * z[i];

			//Pickup Knapsack
			if (_nodes[i]->demands[e] > 0)
				expr2 += _nodes[i]->demands[e] * z[i];
			if (_nodes[i]->demands[e] < 0)
				expr2 -= std::abs(_nodes[i]->demands[e]) * z[i];
			expr2 -= _nodes[i]->w_plus * z[i];
		}

		model2.add(expr1 <= _prob->GetDriver(0)->capacity);
		model2.add(expr2 <= _prob->GetDriver(0)->capacity);

		expr1.end();
		expr2.end();
	}

	cplex2 = IloCplex(model2);
	//cplex2.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
	cplex2.setOut(env.getNullStream());
	cplex2.setWarning(env.getNullStream());
	cplex2.setParam(IloCplex::Param::Threads,1);

	//DriverCountBppLbLazyCallBack * lazycallback = new (env) DriverCountBppLbLazyCallBack(env, _prob, z);
	//cplex2.use(lazycallback);
	//cplex2.exportModel("knapsack_subproblem.lp");
}
bool DriverCountBppLb::ExactPricing()
{
	bool found_new = false;
	//for(int i=0;i<n;i++)
		//if(duals_vars[i] > 0)
			//printf("item:%d dual:%.3lf\n", i, duals_vars[i]);

	int mul_const = 100000;
	//for(int i=0;i<count;i++)
	//	printf("item:%d x:%d w:%d\n", i, items[i]->x, items[i]->w);

	for(int i=0;i<n;i++)
	{
		char name[20];
		sprintf(name,"z%d",i);z[i].setName(name);
		obj2.setLinearCoef(z[i], ((int)mul_const * duals_vars[i]));
	}

	//cplex2.setParam(IloCplex::ObjLLim, mul_const);
	cplex2.solve();

	double dsol = (int)cplex2.getObjValue();
	pricing_cost = dsol / mul_const;
	int sol = (int)dsol;
	printf("sol:%d\n", sol);
	if(sol <= mul_const) return false;

	/*std::vector<Node*> nodes;
	for(int i=0;i<n;i++)
		if(cplex2.getValue(z[i]) > .9)
			nodes.push_back( _prob->GetCustomer( i ) );
	if(RecourseLowerBound::CalculateLP( _prob, nodes) >= 9999999999)
	{
		IloExpr expr(z.getEnv());
		for(size_t i=0;i<nodes.size();i++)
			expr += z[ nodes[i]->no - 1];
		model2.add(expr <= (int)(nodes.size()-1));
		expr.end();
		continue;
	}*/
	{
		found_new = true;
		IloNumVar newx = IloNumVar(x.getEnv(),0,IloInfinity, ILOFLOAT);
		x.add(newx);
		obj.setLinearCoef(newx, 1);
		//for(int i=0;i<n;i++)
			//printf("z[%d] = %lf\n",i,cplex2.getValue(z[i]));

		DriverCountBppColumn col;
		col.id = columns.size();

		for(int i=0;i<n;i++)
			if(cplex2.getValue(z[i]) > .9)
				col.Add(i);
		col.Show();
		columns.push_back(col);

		for(size_t i=0;i<col.nodes.size();i++)
		{
			consts_vars[ col.nodes[i] ].setLinearCoef(newx, 1);
		}
			

	}

	return found_new;
}

DriverCountBppLbLazyCallBack::DriverCountBppLbLazyCallBack(IloEnv env, Prob * prob, IloNumVarArray z) :
											 IloCplex::LazyConstraintCallbackI(env), _prob(prob), _z(z)
{

}

void DriverCountBppLbLazyCallBack::main()
{

	IloNumArray values(getEnv());
	getValues(values,_z);
	std::vector<Node*> nodes;
	for(int i = 0 ; i < _prob->GetCustomerCount();i++)
		if(values[i] >= 0.9)
			nodes.push_back( _prob->GetCustomer(i) );

	//printf("DriverCountBppLbLazyCallBack obj:%.2lf nodes:", (double)getObjValue());
	//for(size_t i=0;i<nodes.size();i++)
	//	printf("%d-", nodes[i]->no-1);
	//printf("\n");

	if(RecourseLowerBound::CalculateLP( _prob, nodes) >= 9999999999)
	{
		IloExpr expr(_z.getEnv());
		for(size_t i=0;i<nodes.size();i++)
			expr += _z[ nodes[i]->no - 1];
		add(expr <= (int)(nodes.size()-1));
		expr.end();
	}
}
