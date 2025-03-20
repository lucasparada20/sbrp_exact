#include "ExactSbrpSep.h"

ExactSbrpSepO::ExactSbrpSepO(IloEnv env, ExactSbrpGraphO* graph, IloNumVarArray x, IloNumVar theta, IloNumVarArray thetas) : _env(env), _graph(graph), _x(x), _theta(theta), _prob(graph->GetProblem()),_thetas(thetas)
{
	nb_inf_sets = 0;
	nb_inf_paths = 0;
	nb_sub_tours = 0;
	nb_sub_tour_frac = 0;  
	nb_l_cuts = 0;
	nb_p_cuts = 0;
	nb_frac_l_cuts = 0;
	nb_sorted_l_cuts = 0;
	nb_benders_cuts = 0;
	best_sol = best_sol_recourse = best_sol_distance = 9999999999;

	_component.resize(_graph->GetNodeCount(), -1);
}

//Separates inequalities from fractional solutions
void ExactSbrpSepO::SeparateFrac(IloRangeArray array)
{	
	
	for(int i = 0; i < _graph->GetNodeCount(); i++)//stations + depot
		_component[i] = -1;

	int comp = -1;
	for(int i=1;i<_graph->GetNodeCount();i++)
		if(_component[_graph->GetNode(i)->no] == -1)
			ResearchConnectedComponent(_graph->GetNode(i), ++comp);
	comp++;

	std::vector< std::vector<Node*> > components(comp);
	for(int i=1;i<_graph->GetNodeCount();i++)
	{
		Node * n = _graph->GetNode(i);
		components[ _component[n->no] ].push_back(n);
	}
	
	for(int i=0;i<comp;i++)
	{
		if(components[i].size()<=2) continue;//Added to the root
		
		double sumv = 0;
		int nb=components[i].size();
		RecourseLowerBound rec;
		int nb_veh = rec.GetDriverCount(_prob, components[i]);
		IloExpr expr(array.getEnv());
		for(int m = 0; m < nb; m++)
			for(int n = 0; n < nb; n++)
				if (m != n)
				{
					ExSbrpArcO* arc = _graph->GetArc(components[i][m]->no, components[i][n]->no);
					if (arc == NULL) continue;
					sumv += arc->value;
					expr += _x[arc->index];
				}
			
		if(sumv > nb - nb_veh + 0.2)
		{
			array.add(expr <= nb - nb_veh);
			nb_sub_tour_frac++;

		}
		expr.end();
		
		double L = rec.Calculate(_prob, components[i]);
		double sumt=0.0;
		for(int m = 0; m < components[i].size(); m++) 
				sumt += _graph->GetTheta(components[i][m]->no);
		double lhs = sumt - L*sumv;
		double rhs = L * (nb_veh + 1 - nb);

		if(L-0.6 > sumt && lhs+0.6<rhs && L > 0.2)
		{			
			IloExpr expr(array.getEnv());
			for(int m = 0; m < components[i].size(); m++) 
				expr += _thetas[ components[i][m]->no ];
			
			for(int m = 0; m < components[i].size(); m++)
				for(int n = 0; n < components[i].size(); n++)
					if (m != n)
					{
						ExSbrpArcO* arc = _graph->GetArc(components[i][m]->no, components[i][n]->no);
						if (arc != NULL) expr -= L * _x[arc->index];
					}
			//std::cout << " frac L-cut:" << expr << " >= " << rhs << std::endl;
			array.add(expr >= rhs);
			nb_frac_l_cuts++;
			expr.end();
		}	

		
	}		
	
	//Look for tours and see if they are infeasible
	for(int i = 0; i < _graph->GetNodeCount(); i++)//stations + depot
		_component[i] = -1;
	ResearchDepotComponent(_graph->GetNode(0), 0);

	std::vector<Node*> tour;
	for(int i=1;i<_component.size();i++)
		if(_component[i]==0)
			tour.push_back(_graph->GetNode(i));
	double sumv=0.0;
	for(int m = 1; m < tour.size(); m++)
	{
		ExSbrpArcO* arc = _graph->GetArc(tour[m-1]->no, tour[m]->no);
		if (arc != NULL) sumv += arc->value;
	}
	int nb = tour.size();
	tour.insert(tour.begin(),_graph->GetNode(0)); //Inf set method needs depots...
	tour.push_back(_graph->GetNode(0));
	
	if(tour.size()<=5) return; // All 2,3-size inf sets are added at the root as inf arcs.
	//Additionally, all single routes are feasible.
	//So an inf set needs to be at least: 0-i-j-k-l-0, size>=6

	RecourseLowerBound rec;
	int nb_veh = rec.GetDriverCount(_prob, tour);
	if( sumv > nb - 2 - nb_veh + 0.2 )
		TestAndAddInfeasiblePath(tour,array);
}

//Separates inequalities from integer solutions
void ExactSbrpSepO::SeparateInt(IloRangeArray array)
{
	if(!SeparateUserCapRecursive(array))
	{
		_graph->MakePaths();
		for(int i=0;i<_graph->GetPathCount();i++)
		{
			std::vector<Node*> path = _graph->GetPath(i);
			if(path.size()<=5) continue; // All 2,3-size inf sets are added at the root as inf arcs.
			//Additionally, all single routes are feasible.
			//So an inf set needs to be at least: 0-i-j-k-l-0, size>=6
			TestAndAddInfeasiblePath(path,array);
		}
	}	

}

bool ExactSbrpSepO::SeparateUserCapRecursive(IloRangeArray array)
{
	bool added_inq=false;
	
	for(int i = 0; i < _graph->GetNodeCount(); i++)
		_component[i] = -1;

	int comp = -1;
	for(int i = 0; i < _graph->GetNodeCount(); i++)
		if (_component[_graph->GetNode(i)->no] == -1)
			ResearchConnectedComponent(_graph->GetNode(i), ++comp);
	
	//Could be done in the same fashion as the fractional components
	for(int co = 1; co <= comp; co++)
	{
		//get the ordered cycle
		std::vector<Node*> tour;
		for(int i=1;i<_graph->GetNodeCount();i++)
			if(_component[_graph->GetNode(i)->no] == co)
			{
				ExSbrpArcO* ar = _graph->GetArcsOutPos(_graph->GetNode(i)->no, 0);
				tour.push_back(ar->from);

				Node* from = ar->from;
				Node* to = ar->to;
				while (1)
				{
					tour.push_back(to);
					from = to;
					ar = _graph->GetArcsOutPos(to->no, 0);
					to = ar->to;
					if(to->no == _graph->GetNode(i)->no) break;
				}
				break;
			}
		if (tour.size() == 0) continue;

		int r_S = RecourseLowerBound::GetDriverCount(_prob, tour);
		
		double sumv=0.0;
		IloExpr expr(array.getEnv());
		for(int m = 0; m < tour.size(); m++)
			for(int n = 0; n < tour.size(); n++)
				if (m != n)
				{
					ExSbrpArcO* arc = _graph->GetArc(tour[m]->no, tour[n]->no);
					if (arc != NULL) expr += _x[arc->index];
					if (arc != NULL) sumv += arc->value;
				}
		double rhs = (tour.size()) - r_S;

		array.add(expr <= rhs);
		nb_sub_tours++;
		expr.clear();
		added_inq=true;
		
		nb_l_cuts =+ TestAndAddSortedLCuts(tour,array);

		expr.end();

	}// End connected component
	
	if(!added_inq) return added_inq; // Easier to find inf sets later
	
	int cur_path = 0;
	for(int i = 0; i < _graph->GetArcsOutPosCount(0); i++)
	{
		std::vector<Node*> tour;
		int nb_arcs=0;
		ExSbrpArcO* ar = _graph->GetArcsOutPos(0, i);
		tour.push_back(ar->from);
		nb_arcs++;

		Node* from = ar->from;
		Node* to = ar->to;
		while (1)
		{
			tour.push_back(to);
			from = to;
			if (from->no == 0) break;
			ar = _graph->GetArcsOutPos(to->no, 0);
			to = ar->to;
			nb_arcs++;
		}

		if(tour.size()<=5) return added_inq; // All 2,3-size inf sets are added at the root as inf arcs.
		//Additionally, all single routes are feasible.
		//So an inf set needs to be at least: 0-i-j-k-l-0, size>=6

		int inf_set_inq_found = TestAndAddInfeasiblePath(tour, array);
	}
	return added_inq;

}

void ExactSbrpSepO::ResearchConnectedComponent(Node* n, int comp)
{
	_component[n->no] = comp;
	for(int i = 0; i < _graph->GetArcsOutPosCount(n->no); i++)
	{
		ExSbrpArcO* a = _graph->GetArcsOutPos(n->no, i);
		if (a->to->no!=0 && _component[a->to->no] == -1 && a->value >= 0.45)//more subtours than inf set
		{
			_component[a->to->no] = comp;
			ResearchConnectedComponent(_graph->GetNode(a->to->no),comp);
		}
	}
}

void ExactSbrpSepO::ResearchDepotComponent(Node* n, int comp)
{
	_component[n->no] = comp;
	for(int i = 0; i < _graph->GetArcsOutPosCount(n->no); i++)
	{
		ExSbrpArcO* a = _graph->GetArcsOutPos(n->no, i);
		if (_component[a->to->no] == -1 && a->value >= 0.60)
		{
			_component[a->to->no] = comp;
			ResearchDepotComponent(_graph->GetNode(a->to->no),comp);
		}
	}
}

// Test if a the path is infeasible
// In such case, it searches for the smallest infeasible set
// The smallest infeasible path
// Note that the first and last nodes are the depot
int ExactSbrpSepO::TestAndAddInfeasiblePath(std::vector<Node*> & path, IloRangeArray array)
{
	int nb_inq = 0;
	if(RouteFeasibility::IsFeasible(_prob, path)) return 0;
	if(path.size() <= 3) return 0;

	std::vector<Node*> subpath;

	//First, check for infeasible sets
	bool found_inf = false;
	for(int k=1;k+2<path.size();k++)
	{
		for(int j=1;j+2<path.size() && j+k+1<path.size();j++)
		{
			subpath.clear();
			subpath.push_back(path[0]);

			for(int l=j;l<=j+k && l+1 < path.size();l++)
				subpath.push_back(path[l]);
			subpath.push_back(path.back());

			int nb_drivers = RecourseLowerBound::GetDriverCount(_prob, subpath);
			if(nb_drivers <= 1) continue;

			IloExpr expr(array.getEnv());
			for(int l=1;l+1<subpath.size();l++)
				for(int k=1;k+1<subpath.size();k++)
					if(l != k)
					{
						ExSbrpArcO* arc = _graph->GetArc( subpath[l]->no, subpath[k]->no);
						if(arc != NULL)
							expr += _x[arc->index];
					}

			array.add(expr <= ((int)subpath.size()) - 2 - nb_drivers);
			expr.end();
			nb_inf_sets++;
			nb_inq++;
			found_inf = true;
			break;
		}//end for j
		if(found_inf) break;

	}//end for k
	if(found_inf) return nb_inq;

	//Second, check for infeasible path (since a set bounds more solutions, from here onwards the code is not executed)
	{
		IloExpr expr(array.getEnv());
		int nb = 0;
		for(int j = 1; j < path.size(); j++)
		{
			if (path[j - 1]->type != NODE_TYPE_CUSTOMER || path[j]->type != NODE_TYPE_CUSTOMER) continue;

			ExSbrpArcO* arc = _graph->GetArc(path[j - 1]->no, path[j]->no);
			expr += _x[arc->index];
			nb++;
		}

		printf("infeasible path\n");
		std::cout << expr << " <= " << nb-1 << std::endl;
		array.add(expr <= nb - 1);
		nb_inf_paths++;
		nb_inq++;
		expr.end();
	}

	if(path.size() <= 4) return nb_inq;


	for(int k=1;k+2<path.size();k++)
	{
		bool found_inf = false;
		for(int j=1;j+2<path.size() && j+k+1<path.size();j++)
		{
			subpath.clear();
			subpath.push_back(path[0]);

			for(int l=j;l<=j+k && l+1 < path.size();l++)
				subpath.push_back(path[l]);
			subpath.push_back(path.back());

			if(RouteFeasibility::IsFeasible(_prob, subpath)) continue;

			IloExpr expr(array.getEnv());
			int nb = 0; double sum_v=0.0;
			for(int l = 1; l < subpath.size(); l++)
			{
				if(subpath[l - 1]->type != NODE_TYPE_CUSTOMER || subpath[l]->type != NODE_TYPE_CUSTOMER) continue;

				ExSbrpArcO* arc = _graph->GetArc(subpath[l - 1]->no, subpath[l]->no);
				if(arc != NULL)
				{
					expr += _x[arc->index];
					sum_v += arc->value;
					nb++;
				}
			}
			if(sum_v > nb - 1 + 0.2 )
			{
				array.add(expr <= nb-1);
				nb_inf_paths++;
				nb_inq++;
				expr.end();
				found_inf = true;				
			}

		}//end for j

		if(found_inf) break;
	}//end for k

	return nb_inq;
}

//Path is without depots
int ExactSbrpSepO::TestAndAddSortedLCuts(std::vector<Node*>& path,IloRangeArray array)
{
	int added_sorted_l_cuts=0;
	if (path.size()<=1) return 0;
	
	path.insert(path.begin(), _prob->GetNode(_prob->GetNodeCount() - 1));
	path.push_back(_prob->GetNode(_prob->GetNodeCount() - 1));
	
	std::vector< SubPathSbrp > l_paths;
	int path_size = path.size();
	
	std::vector<Node*> subpath;

	for(int k = path_size-1; k >= std::max(1,path_size-5); k--)
	{
		for(int j = 1; j + 2 < path.size() && j + k + 1 < path.size(); j++)
		{
			subpath.clear();
			subpath.push_back(path[0]);
			for(int l = j; l <= j + k && l + 1 < path.size(); l++)
				subpath.push_back(path[l]);
			subpath.push_back(path.back());

			RecourseLowerBound rec;
			double L = rec.Calculate(_prob, subpath);

			SubPathSbrp pp;
			pp.nb = subpath.size()-2;
			pp.cost = L;
			pp.from = j;
			pp.to = j + pp.nb;
			pp.thetas = 0;
			for(size_t l=1;l+1<subpath.size();l++)
				pp.thetas += _graph->GetTheta( subpath[l]->no );
			pp.sort = (L-pp.thetas)/pp.nb;
			
			double sumv = 0;
			for(int m = 0; m < subpath.size(); m++)
				for(int n = 0; n < subpath.size(); n++)
					if (m != n)
					{
						ExSbrpArcO* arc = _graph->GetArc(path[m]->no, path[n]->no);
						if (arc == NULL) continue;
						sumv += arc->value;
					}
			
			double lhs = pp.thetas - L*sumv;
			int nb_veh = rec.GetDriverCount(_prob,subpath);
			double rhs = L * (nb_veh + 1 - pp.nb);
			
			if(L-pp.thetas > 0.6 && lhs+0.6<rhs && L>0.2)
				l_paths.push_back(pp); //Pool of L-Cuts in this function
				
		}
		//exit(1);
		std::sort(l_paths.begin(), l_paths.end());
		
		for(size_t k=0;k<std::min(5, (int)l_paths.size());k++)
		{
			SubPathSbrp & pp = l_paths[k];
			double L = pp.cost;

			IloExpr expr(array.getEnv());
			for(int l = pp.from; l <pp.to; l++)
				expr += _thetas[ path[l]->no ];

			for(int m=pp.from;m<pp.to;m++)
				for(int n=m+1;n<pp.to;n++)
					if( path[m]->no < path[n]->no )
					{
						ExSbrpArcO * arc = _graph->GetArc(path[m]->no,path[n]->no);
						if (arc != NULL) expr -= L*_x[arc->index];
					}		

			subpath.clear();
			for(int l = pp.from; l <pp.to; l++)
				subpath.push_back(path[l]);
			
			RecourseLowerBound rec;
			int nb_veh = rec.GetDriverCount(_prob,subpath);
			double rhs = L * (nb_veh + 1 - pp.nb);

			array.add(expr >= rhs);
			expr.end();
			added_sorted_l_cuts++;
			nb_sorted_l_cuts++;
		}
	}
	return added_sorted_l_cuts;
}

//Separates P & L cuts 
void ExactSbrpSepO::SeparateDissagregatedRecourse(IloRangeArray array) 
{
	std::vector<Node*> subpath;
	double sum_path_recourse = 0;

	for(int i = 0; i < _graph->GetPathCount(); i++)
	{
		std::vector<Node*>& path = _graph->GetPath(i);

		RecourseLowerBound rec;
		double path_recourse = RouteFeasibility::RecourseCost(_prob,path);

		sum_path_recourse += path_recourse;
		double sum_thetas = 0;
		for(int k = 1; k + 1 < path.size(); k++)
			sum_thetas += _graph->GetTheta(path[k]->no);
		if(path_recourse - sum_thetas <= 0.0001) continue;

		std::vector< SubPathSbrp > l_paths;
		std::vector< SubPathSbrp > p_paths;
		int path_size = path.size();
		
		for(int k = path_size-1; k >= 1; k--)
		{
			double max_rec_sub_path = 0;
			for(int j = 1; j + 2 < path.size() && j + k + 1 < path.size(); j++)
			{
				subpath.clear();
				subpath.push_back(path[0]);
				for(int l = j; l <= j + k && l + 1 < path.size(); l++)
					subpath.push_back(path[l]);
				subpath.push_back(path.back());

				double L = rec.Calculate(_prob, subpath);

				SubPathSbrp pp;
				pp.nb = subpath.size()-2;
				pp.cost = L;
				pp.from = j;
				pp.to = j + pp.nb;
				pp.thetas = 0;
				for(size_t l=1;l+1<subpath.size();l++)
					pp.thetas += _graph->GetTheta( subpath[l]->no );
				pp.sort = (L-pp.thetas)/pp.nb;
				
				if(L-pp.thetas >= 0.2)
					l_paths.push_back(pp);

				L = (subpath.size() < path.size())?RouteFeasibility::RecourseCost(_prob, subpath):path_recourse;
				max_rec_sub_path = std::max(max_rec_sub_path, L);
				if(L > path_recourse+0.1)
				{
					printf("Recourse function is not monotonic\n");
					printf("Recourse:%.3lf Path ",path_recourse);
					for(int k = 0; k < path.size(); k++)
						printf("%d-", path[k]->no);
					printf("\n");
					printf("Recourse:%.3lf SubPath ",L);
					for(int k = 0; k < subpath.size(); k++)
						printf("%d-", subpath[k]->no);
					printf("\n");
					exit(1);
				}

				pp.cost = L;
				pp.sort = (L-pp.thetas)/pp.nb;

				if(pp.nb == path.size()-2)
					pp.sort += 10000;

				if(L-pp.thetas >= 0.2 || (subpath.size() == path.size() && L-pp.thetas >= 0.0001))
					p_paths.push_back(pp);
				
			}
			if(path_recourse < 0.01 * max_rec_sub_path) break;
		}

		std::sort(l_paths.begin(), l_paths.end());
		std::sort(p_paths.begin(), p_paths.end());
		
		for(int k=0;k<std::min(3, (int)l_paths.size());k++)	
		{
			SubPathSbrp & pp = l_paths[k];
			double L = pp.cost;

			IloExpr expr(array.getEnv());
			for(int l = pp.from; l <pp.to; l++)
				expr += _thetas[ path[l]->no ];

			for(int m=pp.from;m<pp.to;m++)
				for(int n=m+1;n<pp.to;n++)
					if( path[m]->no < path[n]->no )
					{
						ExSbrpArcO * arc = _graph->GetArc(path[m]->no,path[n]->no);
						if (arc != NULL) expr -= L*_x[arc->index];
					}
			double rhs = L * (2 - pp.nb);

			array.add(expr >= rhs);
			nb_l_cuts++;
			expr.end();
		}
		
		for(int k=0;k<std::min(3, (int)p_paths.size());k++)
		{
			SubPathSbrp & pp = p_paths[k];
			double L = pp.cost;

			IloExpr expr(array.getEnv());
			for(int l = pp.from; l < pp.to; l++)
				expr += _thetas[ path[l]->no ];

			for(int m = pp.from+1; m < pp.to; m++)
			{
				ExSbrpArcO* arc = _graph->GetArc(path[m - 1]->no, path[m]->no);
				expr -= L *_x[arc->index];
			}
			double rhs = L * (2-pp.nb);

			array.add(expr >= rhs);
			nb_p_cuts++;
			expr.end();
		}

	}//end for each path

	double sol_distance = _graph->GetCost();
	if(sol_distance + sum_path_recourse < best_sol)
	{
		best_sol_distance = sol_distance;
		best_sol_recourse = sum_path_recourse;
		best_sol = sol_distance + sum_path_recourse;
		printf("New solution:%.1lf dist:%.1lf rec:%.1lf paths:%d\n", best_sol, best_sol_distance, best_sol_recourse,_graph->GetPathCount());
		best_solution.clear();
		for(int i = 0; i < _graph->GetPathCount(); i++)
			best_solution.push_back( _graph->GetPath(i) );
	}
}

//Separates Benders optimality cuts
void ExactSbrpSepO::SeparateBendersCut(IloRangeArray array)
{
	double Q = _prob->GetDriver(0)->capacity;
	_Dual = new FirstStageDuals(_graph);
	
	for(int i = 0;i<_graph->GetArcCount();i++)
		_graph->UnWalkOnArc(i);
	for(int i=0;i<_graph->GetPosArcCount();i++)
	{
		ExSbrpArcO* arc = _graph->GetPosArc(i);
		_graph->WalkOnArc(arc->index);
	}
	
	IloExpr expr_dual(array.getEnv());	
	for(int e=0;e<_prob->GetScenarioCount();e++)
	{
		if(e==0)
		{ 
		  _Dual->InitModelDual(0);
		  _Dual->SolveDual();		  
		}
		else
		{
		  _Dual->UpdateObjectiveDual(e);
		  _Dual->SolveDual();
		} 
		
		for (int j = 1; j < _graph->GetNodeCount(); j++)
		{
			Node* n = _graph->GetNode(j);	
			expr_dual -= (_Dual->v[j-1]*n->w_plus + _Dual->pi[j-1]*n->w_minus + _Dual->phi[j-1]*n->demands[e]);
		}
		for (int i = 0; i < _graph->GetArcCount(); i++)
		{
			ExSbrpArcO* arc = _graph->GetArc(i);
			expr_dual -= _Dual->ksi[i]*Q*_x[arc->index];
		}
	}//End scenario

	expr_dual /= _prob->GetScenarioCount();
	expr_dual *= Parameters::GetCminEpsilon(); 
	expr_dual += _theta;	
	array.add(expr_dual>=0);
	nb_benders_cuts++;
	expr_dual.end();

	double sum_path_recourse=0.0;
	for(int i=0;i<_graph->GetPathCount();i++)
		sum_path_recourse += RouteFeasibility::RecourseCost(_prob,_graph->GetPath(i));
	
	double sol_distance = _graph->GetCost();
	if(sol_distance + sum_path_recourse < best_sol)
	{
		best_sol_distance = sol_distance;
		best_sol_recourse = sum_path_recourse;
		best_sol = sol_distance + sum_path_recourse;
		best_solution.clear();
		printf("new_sol dist:%.1lf rec:%.1lf\n",sol_distance,sum_path_recourse);
		for(int i = 0; i < _graph->GetPathCount(); i++)
			best_solution.push_back( _graph->GetPath(i) );
	}
	
	delete _Dual;
}