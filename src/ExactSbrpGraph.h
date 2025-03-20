

#ifndef EXACT_SBRP_GRAPH_H
#define EXACT_SBRP_GRAPH_H

#include <algorithm>
#include <vector>
#include <map>
#include "NodeSBRP.h"
#include "DriverSBRP.h"
#include "ProblemDefinition.h"

#define EXACT_SBRP_MIN_VIOLATION 0.2

class ExSbrpArcO
{
public:
	Node* from;
	Node* to;
	int index;
	double cost;
	double time;
	double value;
	char walked_on;
};

class ExSbrpSortPath
{
public:
	int path_no;
	int nb_customers;
};

class ExSbrpSorterPath
{
public:
	bool operator()(const ExSbrpSortPath& a, const ExSbrpSortPath& b)
	{
		if (a.nb_customers - b.nb_customers != 0)
			return a.nb_customers < b.nb_customers;
		else
			return a.path_no < b.path_no;
	}
};

class ExactSbrpGraphO
{
public:
	ExactSbrpGraphO(Prob* prob);

	int GetNodeCount() { return (int)_nodes.size(); }
	Node* GetNode(int index) { return _nodes[index]; }
	Node* GetDepot() { return _nodes[0]; }

	int GetArcCount() { return (int)_arcs.size(); }
	ExSbrpArcO* GetArc(int index) { return &_arcs[index]; }
	ExSbrpArcO* GetArc(int from, int to);

	//Originally set in the constructor
	int GetDemandMax() { return demandMax; }

	int GetPosArcCount() { return (int)_pos_arcs.size(); }
	int GetArcsOfCount(int index) { return (int)_arcs_of[index].size(); }
	int GetArcsInOfCount(int index) { return (int)_arcs_in_of[index].size(); }
	int GetArcsOutOfCount(int index) { return (int)_arcs_out_of[index].size(); }
	int GetArcsOfPosCount(int i) { return (int)_pos_arcs_of[i].size(); }
	int GetArcsInPosCount(int i) { return (int)_pos_arcs_in[i].size(); }
	int GetArcsOutPosCount(int i) { return (int)_pos_arcs_out[i].size(); }

	ExSbrpArcO* GetArcsOfPos(int i, int j) { return _pos_arcs_of[i][j]; }
	ExSbrpArcO* GetArcOf(int index, int j) { return _arcs_of[index][j]; }
	ExSbrpArcO* GetArcInOf(int index, int j) { return _arcs_in_of[index][j]; }
	ExSbrpArcO* GetArcOutOf(int index, int j) { return _arcs_out_of[index][j]; }
	ExSbrpArcO* GetPosArc(int i) { return _pos_arcs[i]; }
	ExSbrpArcO* GetArcsInPos(int j, int i) { return _pos_arcs_in[j][i]; }
	ExSbrpArcO* GetArcsOutPos(int i, int j) { return _pos_arcs_out[i][j]; }
	ExSbrpArcO* GetArcsPos(int from, int to);
	void UnWalkOnPosArc(int i){_pos_arcs[i]->walked_on = 0;}
	void WalkOnPosArc(std::vector<ExSbrpArcO*> CycleArcs)
	{
		for(size_t i = 0; i < CycleArcs.size(); i++)
			for (size_t j = 0; j < _pos_arcs.size(); j++)
			{
			  if(CycleArcs[i]->index == _pos_arcs[j]->index)
			  {
				_pos_arcs[j]->walked_on = 1;
				break;
			  }
			}
	}
  
	void UnWalkOnArc(int i){_arcs[i].walked_on = 0;}
	void WalkOnArc(int i){_arcs[i].walked_on = 1;}
	void WalkOnArc(std::vector<ExSbrpArcO*> CycleArcs){
		for(size_t i = 0; i < CycleArcs.size(); i++)
			for (size_t j = 0; j < _arcs.size(); j++)
			{
			  if(CycleArcs[i]->index == _arcs[j].index)
			  {
				_arcs[j].walked_on = 1;
				break;
			  }
			}
	}
  
	void WalkOnPosArc(int i){_pos_arcs[i]->walked_on = 1;}

	int GetPathCount() { return (int)_paths.size(); }
	std::vector<Node*>& GetPath(int i) { return _paths[i]; }
	std::vector<Node*>& GetSortedPath(int i) { return _paths[_sorted_paths[i].path_no]; }

	double GetSumArcValue(int node)
	{
		double v = 0;
		for (int i = 0; i < GetArcsOfPosCount(node); i++)
			v += GetArcsOfPos(node, i)->value;
		return v;
	}

	double GetTheta(int i){return thetas[i];}
	void SetTheta(int i, double v){thetas[i] = v;}
	void SetThetasSize(int size){thetas.resize(size);}

	double GetCost() { return _cost; }
	bool IsInteger() { return _is_integer; }
	void AssignPositiveValues();
	void ShowPosValueArcs();
	void ShowValueArcs();
	void MakePaths();
	void ShowPaths();
	void PrintGraph(char* filename);

	void SetNodeID(int nodeid) { _nodeid = nodeid; }
	int GetNodeID() { return _nodeid; }

	Prob* GetProblem() { return _prob; }
	int removed_arcs;
	
private:
	
	int _nodeid;
	Prob* _prob;
	std::vector<Node*> _nodes;	//nodes [0, ..., n] 0 is the depot, others are customers
	std::vector< ExSbrpArcO > _arcs;
	std::vector< std::vector<ExSbrpArcO*> > _arcs_of;
	std::vector< std::vector<ExSbrpArcO*> > _arcs_in_of;
	std::vector< std::vector<ExSbrpArcO*> > _arcs_out_of;

	std::vector< ExSbrpArcO* > _pos_arcs; //arcs with a positive value
	std::vector< std::vector<ExSbrpArcO*> > _pos_arcs_of; //arcs with a positive value
	std::vector< std::vector<ExSbrpArcO*> > _pos_arcs_out; //arcs with a positive value
	std::vector< std::vector<ExSbrpArcO*> > _pos_arcs_in; //arcs with a positive value

	std::vector< std::vector<Node*> > _paths;
	std::vector< ExSbrpSortPath > _sorted_paths;
	std::map<int, ExSbrpArcO*> _map_arcs;
	std::vector<std::vector<ExSbrpArcO*>> _matrix_of_arcs;
	std::map<int, ExSbrpArcO*> _map_pos_arcs;
	int demandMax;
	bool _is_integer;
	double _cost;

	std::vector<double> thetas;
	double recourse_theta;
};


#endif
