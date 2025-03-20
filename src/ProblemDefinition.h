#ifndef PROBLEM_DEF_H
#define PROBLEM_DEF_H

#include <stddef.h>
#include <vector>
#include <cstring>
#include "Constants.h"
#include "NodeSBRP.h"
#include "DriverSBRP.h"


class Prob
{
	public:
		Prob(): _nodes(0),_customers(0),_drivers(0),_distances(NULL),_dimension(0),_driver_count_lb(1), _upper_bound(9999999999.9),_delete_matrices(true),_scenarios(0)
		{
			_nodes.reserve(4000);		
			_customers.reserve(4000);
			_drivers.reserve(300);
			_scenarios.reserve(4000);
		}

		~Prob()
		{
			if(_delete_matrices)
			{
				if(	_distances != NULL)
				{
					for(int i = 0 ; i < _dimension ; i++)
						delete [] _distances[i];
					delete [] _distances;
				}
				_distances = NULL;
			}
		}

		void AddCustomer(Node * n){ _customers.push_back(n->id); }
		void AddNode(Node & n){ _nodes.push_back(n); }

		int GetCustomerCount(){ return (int)_customers.size();}
		int GetNodeCount(){ return (int)_nodes.size();}

		Node* GetCustomer(int i){ return &_nodes[ _customers[i] ]; }
		Node* GetNode(int i){ return &_nodes[i]; }

		void AddDriver(Driver & d){ _drivers.push_back(d);}
		int GetDriverCount(){ return (int)_drivers.size();}
		Driver* GetDriver(int i){ return &_drivers[i];}

		/* SECOND STAGE */
		void AddScenario(int i) { _scenarios.push_back(i);}
		int GetScenario(int i) { return (int) _scenarios[i]; }
		int GetScenarioCount() { return (int)_scenarios.size(); }

		void AddSortedScenarios(int i) { _sortedscenarios.push_back(i); }
		int GetSortedScenario(int i) { return (int)_sortedscenarios[i]; }
		std::vector<int> & GetVectorOfSortedScenarios() { return _sortedscenarios; }

		void SetMatrices(double ** d, int dim){ _distances = d; _dimension = dim;}
		double ** GetDistances(){ return _distances;}
		double GetDistance(Node * i, Node * j){ return _distances[i->distID][j->distID];}
		double GetDist(Node * i, Node * j){ return _distances[i->distID][j->distID];}

		void ShowNodes()
		{
			for(size_t i=0;i<_nodes.size();i++)
				_nodes[i].Show();
		}
		void ShowNodes(int scenario)
		{
			for(size_t i=0;i<_nodes.size();i++)
				_nodes[i].Show(scenario);
		}

		double GetUpperBound(){return _upper_bound;}
		void SetUpperBound(double ub){_upper_bound = ub;}

		int GetDriverCountLB(){return _driver_count_lb;}
		void SetDriverCountLB(int d){_driver_count_lb = d;}
	
	private:
		std::vector<Node> _nodes;			//list of nodes
		std::vector<int> _customers;		//list of nodes that are customers
		std::vector<Driver> _drivers;		//list of drivers

		double ** _distances;
		int _dimension;

		int _driver_count_lb;
		double _upper_bound;
		bool _delete_matrices;		//if the problem definition comes from a copy it is false,
											//if it is original it is true

		std::vector<int> _scenarios;
		std::vector<int> _sortedscenarios;
  
};

#endif
