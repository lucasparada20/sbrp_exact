#ifndef _SOLUTION_H
#define _SOLUTION_H

#include <stdio.h>
#include <stdlib.h>
#include "Constants.h"
#include "ProblemDefinition.h"
#include "CostFunctionSBRP.h"
#include <vector>
#include <sstream>
#include <iostream>
#include "NodeSBRP.h"
#include "DriverSBRP.h"

class Sol
{
	public:
		Sol():show_output(true){}
		Sol(Prob * prob, CostFunctionSBRP * cost_func):
										  Next(0),Prev(0),AssignTo(0),
										  RoutesLength(prob->GetDriverCount()),
										  _prob(prob),UnassignedCount(0),
										  Unassigneds(0),UnassignedIndex(0),
										  _cost_func(cost_func),show_output(true),_last_cost(0),_is_feasible(true), _total_distances(0)

		{
			for(int i = 0 ; i < _prob->GetNodeCount() ; i++)
			{
				Next.push_back(NULL);
				Prev.push_back(NULL);
				AssignTo.push_back(NULL);
				Unassigneds.push_back(NULL);
				UnassignedIndex.push_back(-1);
			}

			for(int i = 0 ; i < _prob->GetDriverCount() ; i++)
			{
				RoutesLength[i] = 0;
				Driver * d = GetDriver(i);
				Node * n1 = GetNode(d->StartNodeID);
				Node * n2 = GetNode(d->EndNodeID);

				Next[n1->id] = n2;
				Prev[n2->id] = n1;
				AssignTo[n1->id] = d;
				AssignTo[n2->id] = d;
			}
		}

		void AddNode()
		{
			Next.push_back(NULL);
			Prev.push_back(NULL);
			AssignTo.push_back(NULL);
			Unassigneds.push_back(NULL);
			UnassignedIndex.push_back(-1);
		}

		void AddDriver()
		{
			RoutesLength.push_back(0);
		}

		void InsertAfter(Node * n, Node * prev)
		{
			RoutesLength[ AssignTo[prev->id]->id ]++;
			AssignTo[n->id] = AssignTo[prev->id];
			Next[n->id] = Next[	prev->id ];
			Prev[n->id] = prev;
			if(Next[prev->id] != NULL)
				Prev[ Next[prev->id]->id ] = n;

			Next[prev->id] = n;
		}

		void Remove(Node * n)
		{
			RoutesLength[ AssignTo[n->id]->id ]--;
			if(Next[n->id] != NULL) Prev[ Next[n->id]->id ] = Prev[n->id];
			if(Prev[n->id] != NULL) Next[ Prev[n->id]->id ] = Next[n->id];
			AssignTo[n->id] = NULL;
		}

		void AddToUnassigneds(Node * n)
		{
			UnassignedIndex[n->id] = UnassignedCount;
			Unassigneds[UnassignedCount] = n;
			UnassignedCount++;
		}

		void RemoveFromUnassigneds(Node * n)
		{
			int ind = UnassignedIndex[n->id];
			//printf("RemoveFromUnassigneds:%d index:%d count:%d\n", n->id, ind,UnassignedCount);
			if(UnassignedCount > 1)
			{
				Node * rep = Unassigneds[UnassignedCount - 1];
				Unassigneds[ind] = rep;
				UnassignedIndex[ rep->id ] = ind;
			}

			Unassigneds[UnassignedCount - 1] = NULL;
			UnassignedCount--;
			UnassignedIndex[n->id] = -1;
		}
		void RemoveAndUnassign(Node * n)
		{
			Remove(n);
			AddToUnassigneds(n);
		}

		void PutAllNodesToUnassigned()
		{
			for(int i=0;i<GetCustomerCount();i++)
				AddToUnassigneds( GetCustomer(i));
		}

		void MakePath(int driver, std::vector<Node*> & path)
		{
			Driver * d = GetDriver(driver);
			Node * prev = GetNode(d->StartNodeID);
			for(size_t j=0;j<path.size();j++)
				//if((path[j]->type & NODE_TYPE_CUSTOMER) == NODE_TYPE_CUSTOMER)
				if(path[j]->type == NODE_TYPE_CUSTOMER)
				{
					if(AssignTo[ path[j]->id ] != NULL)
						Remove( path[j] );
					if( UnassignedIndex[ path[j]->id ] != -1)
						RemoveFromUnassigneds( path[j] );
					InsertAfter(path[j], prev);
					prev = path[j];
				}
		}

		void GetPath(Driver * d, std::vector<Node*> & path)
		{
			path.clear();
			Node * n = GetNode(d->StartNodeID);
			while(n != NULL)
			{
				path.push_back(n);
				n = Next[n->id];
			}
		}

		void Show()
		{
			int nbnonempty = 0;
			for(int i=0;i<GetDriverCount();i++)
				if(RoutesLength[i] >= 1)
					nbnonempty++;
			printf("Solution non-empty routes:%d routes:%d cost:%.4lf\n", nbnonempty,GetDriverCount(), _cost_func->GetCost(*this) );
			for(int i=0;i<GetDriverCount();i++)
				if(show_output && RoutesLength[i] >= 1)
					Show(GetDriver(i));

			if(show_output && GetUnassignedCount() >= 1)
			{
				printf("Unassigneds:");
				for(int i=0;i<GetUnassignedCount();i++)
					printf("%d ", GetUnassigned(i)->no);
				printf("\n");
			}
		}

		void Show(Driver * d)
		{
			if(_cost_func != NULL)
				_cost_func->Show(this, d);
		}

		int GetCustomerCount(){ return _prob->GetCustomerCount();}
		Node * GetCustomer(int i){ return _prob->GetCustomer(i);}

		int GetNodeCount(){ return _prob->GetNodeCount();}
		Node * GetNode(int i){ return _prob->GetNode(i);}

		int GetDriverCount(){ return _prob->GetDriverCount();}
		Driver * GetDriver(int i){ return _prob->GetDriver(i);}

		int GetUsedDriverCount()
		{
			int nb = 0;
			for(int i=0;i<GetDriverCount();i++)
				if( RoutesLength[ GetDriver(i)->id ] >= 1)
					nb++;
			return nb;
		}
		int GetUnassignedCount(){ return UnassignedCount;}
		Node * GetUnassigned(int i){ return Unassigneds[i];}

		double GetCost(){ _last_cost = _cost_func->GetCost(*this); return _last_cost;}
		double GetCost(Driver * d){ return _cost_func->GetCost(*this,d);}
		double GetCost(int i){ return _cost_func->GetCost(*this,GetDriver(i));}
		double GetLastCalculatedCost(){ return _last_cost; }
		CostFunctionSBRP * GetCostFunction(){return _cost_func;}

		bool IsFeasible(){return _is_feasible;}
		void SetIsFeasible(bool f){_is_feasible = f;}

		void Update(){ _cost_func->Update(*this);}
		void Update(Driver * d){ _cost_func->Update(*this,d);}

		double GetTotalDistances(){return _total_distances;}
		void SetTotalDistances(double d){_total_distances = d;}

		double GetTotalRecourse(){return _total_recourse;}
		void SetTotalRecourse(double d){_total_recourse = d;}

		bool IsUnassigned(Node * n){return UnassignedIndex[n->id] != -1;}

		Prob * GetProblemDefinition(){return _prob;}
		Prob * GetProb(){return _prob;}

		double ** GetDistances(){ return _prob->GetDistances();}
		double GetDist(Node * n1, Node * n2){return _prob->GetDistances()[n1->distID][n2->distID];}
		int GetRouteLength(int i){return RoutesLength[i];}
		int GetRouteLength(Driver * d){return RoutesLength[d->id];}

		std::vector<Node*> Next;
		std::vector<Node*> Prev;
		std::vector<Driver*> AssignTo;
		std::vector<int> RoutesLength;//# of customers in each route

		Node * GetNext(Node * n){return Next[n->id];}
		Node * GetPrev(Node * n){return Prev[n->id];}
    
	private:
		Prob * _prob;
		int UnassignedCount;
		std::vector<Node*> Unassigneds;
		std::vector<int> UnassignedIndex;
		CostFunctionSBRP * _cost_func;
		
		bool show_output;
		double _last_cost;
		bool _is_feasible;
		double _total_distances;
		double _total_recourse;
};

#endif
