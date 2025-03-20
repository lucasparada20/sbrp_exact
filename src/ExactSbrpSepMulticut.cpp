#include "ExactSbrpSepMulticut.h"
#include "Parameters.h"
#include "RecourseLowerBound.h"
#include <vector>
#include <algorithm>    // std::reverse
#include <map>

ExactSbrpSepMulticut::ExactSbrpSepMulticut(IloEnv env, ExactSbrpGraphO* graph, IloNumVarArray x, IloNumVar theta, IloNumVarArray thetas) : _env(env), _graph(graph), _x(x), _theta(theta), _prob(graph->GetProblem()),_thetas(thetas)
{
	nb_sub_tours = 0;
	nb_sub_tours_from_frac = 0;  
	nb_benders_cuts = 0;
	nb_benders_feasibility_cuts=0;
	best_sol = best_sol_recourse = best_sol_distance = 9999999999;

  _component.resize(_graph->GetNodeCount(), -1);
}

void ExactSbrpSepMulticut::SeperateFrac(IloRangeArray array) //With components
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
		
		if (sumv > nb - 1 + 0.2)
		{
			array.add(expr <= nb - 1);
			nb_sub_tours_from_frac++;		
		}
	}		
}

void ExactSbrpSepMulticut::SeperateInt(IloRangeArray array)
{
	if(!SeperateUserCapRecursive(array))
	{
		_graph->MakePaths();
		TestAndAddFeasBendersCut(array);
	}	

}

bool ExactSbrpSepMulticut::SeperateUserCapRecursive(IloRangeArray array)
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

		{
			
			
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
			double rhs = (tour.size()) - 1;
			//std::cout << "Subtour: " << expr << " <= " << rhs << std::endl;
			array.add(expr <= rhs);
			nb_sub_tours++;
			expr.clear();
			added_inq=true;
			

			expr.end();
                
		}
	}// End connected co
	return added_inq;
}

void ExactSbrpSepMulticut::ResearchConnectedComponent(Node* n, int comp)
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

void ExactSbrpSepMulticut::ResearchDepotComponent(Node* n, int comp)
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

int ExactSbrpSepMulticut::TestAndAddFeasBendersCut(IloRangeArray array)
{
	int nb_inq=0;
	double Q = _prob->GetDriver(0)->capacity;

	for(int i = 0;i<_graph->GetArcCount();i++)
		_graph->UnWalkOnArc(i);
	for(int i=0;i<_graph->GetPosArcCount();i++)
	{
		ExSbrpArcO* arc = _graph->GetPosArc(i);
		_graph->WalkOnArc(arc->index);
	}

	_Dual = new FirstStageDuals(_graph);
	_Dual->Benders_feas=true;//Changes the optimizer to primal simplex
		
	for(int e=0;e<_prob->GetScenarioCount();e++)
	{
		if(e==0)
		{
			_Dual->InitModelDual(0);
			_Dual->SolveDual();	
		}else
		{
			_Dual->UpdateObjectiveDual(e);
			_Dual->SolveDual();	
		}
		if(_Dual->_status == 2)//Unbounded dual -> infeasible primal
		{
			IloExpr expr_dual(array.getEnv());
			double sum=0.0;
			for (int j = 1; j < _graph->GetNodeCount(); j++)
			{
				Node* n = _graph->GetNode(j);

				expr_dual -= (_Dual->v[j-1]*n->w_plus + _Dual->pi[j-1]*n->w_minus + _Dual->phi[j-1]*n->demands[e]);

				double prod1 = _Dual->v[j-1]*n->w_plus + _Dual->pi[j-1]*n->w_minus + _Dual->phi[j-1]*n->demands[e];		
				
				sum+=prod1;
			}
			for (int i = 0; i < _graph->GetArcCount(); i++)
			{
				ExSbrpArcO* arc = _graph->GetArc(i);
				expr_dual -= _Dual->ksi[i]*Q*_x[arc->index];
			}		

			array.add(expr_dual>=0);
			nb_inq++;
			nb_benders_feasibility_cuts++;			
		}

	}

	delete _Dual;
	
	return nb_inq;
}

void ExactSbrpSepMulticut::SeperateBendersCut(IloRangeArray array)
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
		
	for(int e=0;e<_prob->GetScenarioCount();e++)
	{
		IloExpr expr_dual(array.getEnv());
		IloExpr expr_primal(array.getEnv());
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
			
		double sum=0.0;
		
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

		expr_dual /= _prob->GetScenarioCount();
		expr_dual *= Parameters::GetCminEpsilon(); 
		expr_dual += _thetas[e];	
		
		array.add(expr_dual>=0);

		nb_benders_cuts++;
		
		expr_dual.end();
		expr_primal.end();
		
	}
	
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

