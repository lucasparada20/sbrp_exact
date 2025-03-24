#include "ExactSbrp.h"
#include "Constants.h"
#include "Parameters.h"

void ExactSbrpO::Solve(Prob* pprob)
{
	start_time = clock();
	prob = pprob;
	graph = new ExactSbrpGraphO(prob);

	IloEnv env;
	Init(env);
	SolveProblem(env);
	
	Clear(); //Releases resource from cplex callback objects and user defined objects graph and sep
	
	env.end();
}

void ExactSbrpO::SolveProblem(IloEnv env)
{
	try
	{
		bool re = cplex.solve();
		clock_t end_time = clock();
		time_taken = (double)(end_time - start_time) / CLOCKS_PER_SEC;
		
		int cplex_status = (int)cplex.getCplexStatus();
		double sol = re?cplex.getObjValue():0;
		
		printf("re:%d sol:%.3lf status:%d nbnodes:%d time:%.3lf\n", (int)re, sol, cplex_status, (int)cplex.getNnodes(), time_taken);

		if (re && cplex_status != 11)
		{
			status = STATUS_SOLVED;
			if(!solvedAtRoot)
			{
				lb = ub = _sep->best_sol;
				ub_distance = _sep->best_sol_distance;
				ub_recourse = _sep->best_sol_recourse;
			} else {
				lb = ub = cplex.getObjValue();
				ub_recourse = cplex.getValue(theta);
				ub_distance = cplex.getObjValue() - cplex.getValue(theta);
			}

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
				for (int j = 1; j < path.size(); j++)
				{
					ExSbrpArcO* arc = graph->GetArc(path[j - 1]->no, path[j]->no);
					sum_d += arc->cost;
				}
			}
			std::cout << "Distances: " << sum_d << std::endl;
			//graph->PrintGraph((char*)"sol.dot");
		}
		else //timed out
		{
			
			ub = _sep->best_sol;
			ub_distance = _sep->best_sol_distance;
			ub_recourse = _sep->best_sol_recourse;
			lb = cplex.getBestObjValue();
			status = STATUS_UNSOLVED;

		}

		nb_inf_sets = _sep->nb_inf_sets;
		nb_inf_paths = _sep->nb_inf_paths;
		nb_sub_tours = _sep->nb_sub_tours;
		nb_sub_tour_frac = _sep->nb_sub_tour_frac;
		nb_l_cuts = _sep->nb_l_cuts;
		nb_p_cuts = _sep->nb_p_cuts;
		nb_frac_l_cuts = _sep->nb_frac_l_cuts;
		nb_sorted_l_cuts = _sep->nb_sorted_l_cuts;
		nb_benders_cuts = _sep->nb_benders_cuts;
		
		int cuts_type = Parameters::GetTypeOfOptimalityCuts();
		printf("Type of cuts: %s\n",cuts_type == 1 ? "P&L" : cuts_type == 2 ? "Benders" : "Hybrid");
		printf("nb_inf_sets:%d\n", _sep->nb_inf_sets);
		printf("nb_inf_paths:%d\n",_sep->nb_inf_paths);
		printf("nb_sub_tours:%d\n", _sep->nb_sub_tours);
		printf("nb_sub_tours_from_frac:%d\n", _sep->nb_sub_tour_frac);
		printf("nb_p_cuts:%d\n", _sep->nb_p_cuts);
		printf("nb_l_cuts:%d\n", _sep->nb_l_cuts);
		printf("nb_frac_l_cuts:%d\n", _sep->nb_frac_l_cuts);
		printf("nb_sorted_l_cuts:%d\n",_sep->nb_sorted_l_cuts); 
		printf("nb_benders_cuts:%d\n",_sep->nb_benders_cuts);

	} catch (IloException &ex) {
	   std::cerr << ex << std::endl;
	}

}

void ExactSbrpO::Init(IloEnv env)
{
	int Q = prob->GetDriver(0)->capacity;
	model = IloModel(env);

	std::vector<Node*> stations;
	stations.resize(graph->GetNodeCount(), NULL);
	for (int i = 0; i < graph->GetNodeCount(); i++)
		stations[i] = graph->GetNode(i);

	x = IloNumVarArray(env, graph->GetArcCount(), 0, 1, ILOINT);
	z = IloNumVar(env, prob->GetDriverCountLB(), prob->GetDriverCount(), ILOINT);
	theta = IloNumVar(env, 0, IloInfinity, ILOFLOAT);

	z.setName("z");
	theta.setName("t");

	//thetas_i variables for dissagregated recourse
	thetas = IloNumVarArray(env, prob->GetCustomerCount()+1, 0, IloInfinity, ILOFLOAT);

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


	//\sum_{ i \in C} x_{0i} = 2K
	{
		IloExpr expr(env);
		for (int i = 0; i < graph->GetArcsInOfCount(0); i++)
			expr += x[graph->GetArcInOf(0, i)->index];
		model.add(expr == z);
		expr.end();
	}

	{
		IloExpr expr(env);
		for (int i = 0; i < graph->GetArcsOutOfCount(0); i++)
			expr += x[graph->GetArcOutOf(0, i)->index];
		model.add(expr == z);
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
			expr.end();
		}

		for (int i = 1; i < graph->GetNodeCount(); i++)
		{
			IloExpr expr(env);
			for (int j = 0; j < graph->GetArcsOutOfCount(i); j++)
				expr += x[graph->GetArcOutOf(i, j)->index];
			model.add(expr == 1);
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
		}
		expr.end();
	}

	if(Parameters::GetTypeOfOptimalityCuts() == OPT_CUT_TYPE_PL || Parameters::GetTypeOfOptimalityCuts() == OPT_CUT_TYPE_HYBRID)
	{
		std::map<int, bool> are_in;
		for (int i = 0; i < graph->GetArcCount(); i++)
		{
			ExSbrpArcO* ar = graph->GetArc(i);
			if(ar->from->no == 0 || ar->to->no == 0) continue;

			int key = std::max(ar->from->no, ar->to->no) * 1000 + std::min(ar->from->no, ar->to->no);
			if(are_in[key]) continue;
			are_in[key] = true;

			double cost = 0;
			int nb = 0;
			for(int e = 0;e<prob->GetScenarioCount();e++)
				if(std::abs(ar->from->demands[e] + ar->to->demands[e]) > Q)
				{
					cost += std::abs(ar->from->demands[e] + ar->to->demands[e]) - Q;
					nb++;
				}
			cost = (cost / prob->GetScenarioCount()) * Parameters::GetCminEpsilon();
			//"Only add this cut if the total violation across all scenarios is at least 2 units per scenario on average."
			if(nb >= 1 && cost >= 2*Parameters::GetCminEpsilon())
			{
				IloExpr expr(env);
				expr -= cost * x[ar->index];
				expr += thetas[ ar->from->no ];
				expr += thetas[ ar->to->no ];
				ExSbrpArcO* arc = graph->GetArc(ar->to->no , ar->from->no);
				if (arc != NULL)
					expr -= cost * x[arc->index];
				//std::cout << expr << " <= 0" << std::endl;
				model.add(expr >= 0);
				
			}
		}
	}

	//theta >= sum over all arcs (i,j) of sigma_ij x_ij
	//where sigma_ij is the recourse cost of traversing arc (i,j)
	{
		IloExpr expr(env);
		expr += theta;
		int nb1=0;
		for (int i = 0; i < graph->GetArcCount(); i++)
		{
			ExSbrpArcO* ar = graph->GetArc(i);
			if(ar->from->no == 0 || ar->to->no == 0) continue;

			double cost = 0;
			int nb = 0;
			for(int e = 0;e<prob->GetScenarioCount();e++)
				if(std::abs(ar->from->demands[e] + ar->to->demands[e]) > Q)
				{
					cost += std::abs(ar->from->demands[e] + ar->to->demands[e]) - Q;
					nb++;
				}
			cost = (cost / prob->GetScenarioCount()) * Parameters::GetCminEpsilon();
			if(nb >= 1)
			{
				expr -= cost * x[ar->index];
				nb1++;
			}
		}

		//std::cout << expr << " >= 0" << std::endl;
		if(nb1 >= 1)
			model.add(expr >= 0);
		expr.end();
	}

	//Bounding theta >= sum_{i \in arcs} L * x_i 
	{
		IloExpr expr(env);
		expr += theta;
	
		RecourseLowerBound lb;
    	double l = lb.Calculate(prob);

		for (int i = 0; i < graph->GetArcCount(); i++)
		{
			ExSbrpArcO* ar = graph->GetArc(i);
			if(ar->from->no == 0 || ar->to->no == 0) continue;
			expr -= l * x[ar->index];
		}
		model.add(expr >= (l* (3 - graph->GetNodeCount())) );
		//std::cout << expr << " >= " << (l* (3 - graph->GetNodeCount())) << std::endl;
		expr.end();
	}

	//L_1 (\underscore(K)_1) inequality
	{
		IloExpr expr(env);

		for(int i=0;i<prob->GetCustomerCount();i++)
			expr += thetas[ prob->GetCustomer(i)->no ];

		RecourseLowerBound rec;
		double lb_recourse = rec.CalculateWithMinDriverCount(prob);

		double rhs = lb_recourse + Q * Parameters::GetCminEpsilon() * prob->GetDriverCountLB();
		expr += z * Q * Parameters::GetCminEpsilon();

		//std::cout << expr << " >= " << rhs << std::endl;
		model.add(expr >= rhs);
		expr.end();
	}
  
	// theta >= sum of all thetas_i
	{
		IloExpr expr(env);
		expr += theta;

		for (size_t i = 0; i < thetas.getSize(); i++)
			expr -= thetas[i];

		//std::cout << expr << " >= 0" << std::endl;
		model.add(expr >= 0);
		expr.end();
	}

	AddInfSetInq(env);
	
	//First, the root node
	cplex = IloCplex(model);

	cplex.setParam(IloCplex::Param::Threads, 1);
	cplex.setParam(IloCplex::Param::TimeLimit,max_time);

	_sep = new ExactSbrpSepO(env, graph, x, theta, thetas);
	lazy_call = new (env) ExactSbrpLazyCallBackO(env, graph, x, thetas, _sep);
	usercut_call = new (env) ExactSbrpUserCutCallBackO(env, graph,x, thetas,_sep);

	cplex.setWarning(env.getNullStream());
	cplex.use(lazy_call);
	cplex.use(usercut_call);
	
	lazy_call->add_constraints = usercut_call->add_constraints = true;
	cplex.setParam(IloCplex::Param::MIP::Limits::Nodes,0);

	if(cplex.solve())
	{
		printf("Solved the problem in the Root Node!\n"); solvedAtRoot=true;
	}
		
	double obj = cplex.getBestObjValue();
	printf("Best Obj:%.3lf\n", obj);
	model.add(lazy_call->added_constraints);
	model.add(usercut_call->added_constraints);

	//Second, add root cuts and solve the model
	lazy_call->add_constraints = usercut_call->add_constraints = false;
	_sep->best_sol = 9999999999; _sep->best_sol_recourse = 9999999999; _sep->best_sol_distance = 9999999999;
	
	cplex.setParam(IloCplex::Param::MIP::Limits::Nodes,9999999999);	
	cplex.setParam(IloCplex::Param::MIP::Tolerances::MIPGap, 0);
	cplex.setParam(IloCplex::Param::MIP::Tolerances::UpperCutoff, prob->GetUpperBound() + EPSILON + 0.01);

}

void ExactSbrpO::Clear()
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
		cplex.remove(usercut_call);
		delete usercut_call;
	}
}

void ExactSbrpO::AddInfSetInq(IloEnv env)
{
	std::vector<Node*> path(3, NULL);
	for(int i=1;i<graph->GetNodeCount();i++)
		for(int j=i+1;j<graph->GetNodeCount();j++)
			for(int k=j+1;k<graph->GetNodeCount();k++)
			{
				path[0]=graph->GetNode(i);
				path[1]=graph->GetNode(j);
				path[2]=graph->GetNode(k);
				
				int nb_veh = RecourseLowerBound::GetDriverCount(prob,path);
				if(nb_veh==1) 
				{
					continue;
				}
				
				int cntr=0;
				IloExpr expr(env);
					for(size_t u = 0; u < path.size(); u++)
						for(size_t v = 0; v < path.size(); v++)
							if (u != v)
							{
								ExSbrpArcO* arc = graph->GetArc(path[u]->no, path[v]->no);
								if (arc == NULL) continue;
								expr += x[arc->index];
								cntr++;
							}
				if(cntr > 3)
				{
					model.add(expr <= 3 - nb_veh);
					//std::cout << "root inf set:" << expr << " <= " << 3-nb_veh << std::endl;				
				}
				expr.end();						
			}
}
