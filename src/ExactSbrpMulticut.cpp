#include "ExactSbrpMulticut.h"
#include "Constants.h"
#include "Parameters.h"

void ExactSbrpMulticut::Solve(Prob* pprob)
{
	start_time = clock();
	prob = pprob;
	graph = new ExactSbrpGraphO(prob);

	IloEnv env;
	Init(env);
	SolveProblem(env);
	Clear();
	env.end();
}

void ExactSbrpMulticut::SolveProblem(IloEnv env)
{
	try {

		bool re = cplex.solve();

		if ( ( cplex.getStatus() == IloAlgorithm::Infeasible ) ||
		( cplex.getStatus() == IloAlgorithm::InfeasibleOrUnbounded  ) ) {
			std::cout <<  std::endl << "No solution - starting Conflict refinement" <<  std::endl;

			IloConstraintArray infeas(env);
			IloNumArray preferences(env);

			infeas.add(lazy_call->added_constraints);
			infeas.add(usercut_call->added_constraints);
			if ( lazy_call->added_constraints.getSize() || usercut_call->added_constraints.getSize() ) {
				std::cout << "Lazy Constraints and User Cuts ignored" << std::endl;
			}
			for(IloInt i = 0; i<x.getSize(); i++)
			{
				infeas.add(IloBound(x[i], IloBound::Lower));
				infeas.add(IloBound(x[i], IloBound::Upper));			 
			}
			infeas.add(IloBound(z, IloBound::Lower));
			infeas.add(IloBound(z, IloBound::Upper));			 
			for(IloInt i = 0; i<thetas_multicut.getSize(); i++)
			{
				infeas.add(IloBound(thetas_multicut[i], IloBound::Lower));
				infeas.add(IloBound(thetas_multicut[i], IloBound::Upper));			 
			}
			infeas.add(IloBound(theta, IloBound::Lower));
			infeas.add(IloBound(theta, IloBound::Upper));	 

			for (IloInt i = 0; i<infeas.getSize(); i++) {
				preferences.add(1.0);  // user may wish to assign unique preferences
			}

			if ( cplex.refineConflict(infeas, preferences) ) {
				IloCplex::ConflictStatusArray conflict = cplex.getConflict(infeas);
				env.getImpl()->useDetailedDisplay(IloTrue);
				std::cout << "Conflict :" <<  std::endl;
				for (IloInt i = 0; i<infeas.getSize(); i++) {
				if ( conflict[i] == IloCplex::ConflictMember)
					std::cout << "Proved  : " << infeas[i] <<  std::endl;
				if ( conflict[i] == IloCplex::ConflictPossibleMember)
					std::cout << "Possible: " << infeas[i] <<  std::endl;
				}
			}
		else
		std::cout << "Conflict could not be refined" <<  std::endl;
		std::cout <<  std::endl;
    }
   
	int cplex_status = (int)cplex.getCplexStatus();
	double sol = re?cplex.getObjValue():0;
	
	clock_t end_time = clock();	
	time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;
	printf("re:%d sol:%.2lf status:%d nbnodes:%d time:%.3lf\n", (int)re, sol, cplex_status, (int)cplex.getNnodes(), time_taken);

	if (re && cplex_status != 11)
	{
		status = STATUS_SOLVED;
		lb = ub = _sep->best_sol;
		ub_distance = _sep->best_sol_distance;
		ub_recourse = _sep->best_sol_recourse;
		//build the solution
		for (int i = 0; i < graph->GetArcCount(); i++)
			graph->GetArc(i)->value = cplex.getValue(x[i]);
		graph->AssignPositiveValues();
		graph->MakePaths();
		//graph->ShowPosValueArcs();
		graph->ShowPaths();
		std::cout << "Theta: " << cplex.getValue(theta) << std::endl;

		double sum_d = 0;
		for (int i = 0; i < graph->GetPathCount(); i++)
		{
			std::vector<Node*>& path = graph->GetPath(i);
			for (size_t j = 1; j < path.size(); j++)
			{
				ExSbrpArcO* arc = graph->GetArc(path[j - 1]->no, path[j]->no);
				sum_d += arc->cost;
			}
		}
		std::cout << "Distances: " << sum_d << std::endl;

	}
	else //timed out
	{
		ub = _sep->best_sol;
		ub_distance = _sep->best_sol_distance;
		ub_recourse = _sep->best_sol_recourse;
		lb = cplex.getBestObjValue();
		status = STATUS_UNSOLVED;

	}

	nb_sub_tours = _sep->nb_sub_tours;
	nb_sub_tours_from_frac = _sep->nb_sub_tours_from_frac;
	nb_benders_cuts = _sep->nb_benders_cuts;
	nb_benders_feasibility_cuts = _sep->nb_benders_feasibility_cuts;

	printf("nb_sub_tours:%d\n", _sep->nb_sub_tours);
	printf("nb_sub_tours_from_frac:%d\n", _sep->nb_sub_tours_from_frac);
	printf("nb_benders_cuts:%d\n",_sep->nb_benders_cuts);
	printf("nb_benders_feasibility_cuts:%d\n",_sep->nb_benders_feasibility_cuts);


	} catch (IloException &ex) {
	   std::cerr << ex << std::endl;
	}

}

void ExactSbrpMulticut::Init(IloEnv env)
{
	model = IloModel(env);
	
	constraints = IloRangeArray(env);

	std::vector<Node*> stations;
	stations.resize(graph->GetNodeCount(), NULL);
	for (int i = 0; i < graph->GetNodeCount(); i++)
		stations[i] = graph->GetNode(i);

	x = IloNumVarArray(env, graph->GetArcCount(), 0, 1, ILOINT);
	z = IloNumVar(env, 1, prob->GetDriverCount(), ILOINT);
	theta = IloNumVar(env, 0, IloInfinity, ILOFLOAT);

	z.setName("z");
	theta.setName("t");

	//thetas_i variables for dissagregated recourse
	thetas = IloNumVarArray(env, prob->GetCustomerCount()+1, 0, IloInfinity, ILOFLOAT);

	//IF MULTICUT
	thetas_multicut = IloNumVarArray(env,prob->GetScenarioCount(),0,IloInfinity,ILOFLOAT);
	for (int i = 0; i < thetas_multicut.getSize(); i++)
	{
		char name[40];
		sprintf(name, "t_multi%d", i);
		thetas_multicut[i].setName(name);
	}	
	////////////////////////////////

	for (int i = 0; i < thetas.getSize(); i++)
	{
		char name[40];
		sprintf(name, "t%d", i);
		thetas[i].setName(name);
	}

	for (int i = 0; i < graph->GetArcCount(); i++)
	{
		ExSbrpArcO* ar = graph->GetArc(i);
		//printf("arc:%d index:%d cost:%lf from:%d to:%d\n",i,ar->index,ar->cost,ar->from->no,ar->to->no);
		char name[40];
		sprintf(name, "x%d_%d", ar->from->no, ar->to->no);
		x[i].setName(name);
		if (ar->from->no != 0)
			x[i].setBounds(0, 1);
	}

	IloExpr obj1(env);
	for (int i = 0; i < graph->GetArcCount(); i++)
		obj1 += graph->GetArc(i)->cost * x[i];
	obj1 += theta;
	obj_func = IloMinimize(env, obj1);
	model.add(obj_func);
	obj1.end();

	/* FIRST STAGE CONSTRAINTS*/

	//\sum_{ i \in C} x_{0i} = 2K
	{
		IloExpr expr(env);
		for (int i = 0; i < graph->GetArcsInOfCount(0); i++)
			expr += x[graph->GetArcInOf(0, i)->index];
		model.add(expr == z);
		constraints.add(expr - z == 0);
		expr.end();
	}

	{
		IloExpr expr(env);
		for (int i = 0; i < graph->GetArcsOutOfCount(0); i++)
			expr += x[graph->GetArcOutOf(0, i)->index];
		model.add(expr == z);
		constraints.add(expr - z == 0);
		expr.end();
	}

	//station's in and out degree constrains
	{
		for (int i = 1; i < graph->GetNodeCount(); i++)
		{
			IloExpr expr(env);
			for (int j = 0; j < graph->GetArcsInOfCount(i); j++)
				expr += x[graph->GetArcInOf(i, j)->index];
			model.add(expr == 1);
			constraints.add(expr == 1);
			expr.end();
		}

		for (int i = 1; i < graph->GetNodeCount(); i++)
		{
			IloExpr expr(env);
			for (int j = 0; j < graph->GetArcsOutOfCount(i); j++)
				expr += x[graph->GetArcOutOf(i, j)->index];
			model.add(expr == 1);
			constraints.add(expr == 1);
			expr.end();
		}
	}

	//x_{ij} + x_{ji} \leq 1
	//Necessary for preventing 2-cycles
	{
		for (int i = 1; i < graph->GetNodeCount(); i++)
			for (int j = i + 1; j < graph->GetNodeCount(); j++)
				if (i != j)
				{
					ExSbrpArcO* ar = graph->GetArc(i, j);
					ExSbrpArcO* arr = graph->GetArc(j, i);
					if(ar == NULL || arr == NULL) continue;

					IloExpr expr(env);
					expr += x[ar->index];
					expr += x[arr->index];
					model.add(expr <= 1);
					constraints.add(expr <= 1);
					expr.end();
				}
	}

	//Infeasible Arc Inequalities (to outgoing station 'h')
	std::vector<Node*> path(3, NULL);
	for (int i = 0; i < graph->GetArcCount(); i++)
	{
		ExSbrpArcO* ar = graph->GetArc(i);
		path[0] = ar->from;
		path[1] = ar->to;
		if(ar->to->no == 0) continue;

		IloExpr expr(env);
		expr += x[ar->index];

		int nb=0;
		for (int j = 0; j < graph->GetArcsOutOfCount(ar->to->no); j++)
		{
			ExSbrpArcO* arr = graph->GetArcOutOf(ar->to->no, j);
			path[2] = arr->to;

			if(!RouteFeasibility::IsFeasible(prob, path))
			{
				expr += x[arr->index];
				nb++;
			}
		}

		if(nb >= 1)
		{
			//std::cout << expr << "<= 1" << std::endl;
			model.add(expr <= 1);
			constraints.add(expr <= 1);
		}
		expr.end();
	}

	//Infeasible Arc Inequalities (from incoming station 'h')
	for (int i = 0; i < graph->GetArcCount(); i++)
	{
		ExSbrpArcO* ar = graph->GetArc(i);
		path[1] = ar->from;
		path[2] = ar->to;
		if(ar->from->no == 0) continue;

		IloExpr expr(env);
		expr += x[ar->index];

		int nb=0;
		for (int j = 0; j < graph->GetArcsInOfCount(ar->from->no); j++)
		{
			ExSbrpArcO* arr = graph->GetArcInOf(ar->from->no, j);
			path[0] = arr->from;

			if(!RouteFeasibility::IsFeasible(prob, path))
			{
				expr += x[arr->index];
				nb++;
			}
		}

		if(nb >= 1)
		{
			//std::cout << expr << "<= 1" << std::endl;
			model.add(expr <= 1);
			constraints.add(expr <= 1);
		}
		expr.end();
	}

	// theta >= sum of all thetas_i
	{
		IloExpr expr(env);
		expr += theta;

		for (int i = 0; i < thetas_multicut.getSize(); i++)
			expr -= thetas_multicut[i];

		//std::cout << expr << " >= 0" << std::endl;
		model.add(expr >= 0);
		constraints.add(expr >= 0);
		expr.end();
	}
	
	cplex = IloCplex(model);

	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::TimeLimit,max_time);

	_sep = new ExactSbrpSepMulticut(env, graph, x, theta, thetas_multicut);
	lazy_call = new (env) ExactSbrpMulticutLazyCallBack(env, graph, x, thetas_multicut, _sep);
	usercut_call = new (env) ExactSbrpMulticutUserCutCallBack(env, graph,x, thetas_multicut,_sep);

	cplex.setWarning(env.getNullStream());
	cplex.use(lazy_call);
	cplex.use(usercut_call);
	
	cplex.setParam(IloCplex::Param::MIP::Limits::Nodes,9999999999);	
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
}

void ExactSbrpMulticut::Clear()
{
	delete graph;
	delete _sep;
	if (lazy_call != NULL)
	{
		cplex.remove(lazy_call);
		delete lazy_call;

	}
	if (usercut_call != NULL)
	{
		delete usercut_call;
		cplex.remove(usercut_call);

	}

}

