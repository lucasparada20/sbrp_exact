#include "ExactSbrpGraph.h"
#include "Network.h"
#include "RouteFeasibility.h"

ExactSbrpGraphO::ExactSbrpGraphO(Prob* prob)
{
	_prob = prob;
	demandMax = prob->GetDriver(0)->capacity;

	_nodes.clear(); _nodes.reserve( prob->GetCustomerCount()+1 );
	_nodes.push_back(prob->GetNode(prob->GetDriver(0)->StartNodeID));
	for (int i = 0; i < prob->GetCustomerCount(); i++)
		_nodes.push_back(prob->GetCustomer(i));

	std::vector<Node*> path;
	removed_arcs = 0;
	int k = 0;
	for (size_t i = 0; i < _nodes.size(); i++)
		for (size_t j = 0; j < _nodes.size(); j++)
		{
			if (i == j) continue;

			if(_nodes[i]->type == NODE_TYPE_CUSTOMER && _nodes[j]->type == NODE_TYPE_CUSTOMER)
			{
				path.clear(); path.reserve(4);
				path.push_back(_nodes[0]);
				path.push_back(_nodes[i]);
				path.push_back(_nodes[j]);
				path.push_back(_nodes[0]);

				if(!RouteFeasibility::IsFeasible(prob, path))
				{
					//printf("From i:%d to j:%d infeasible\n", _nodes[i]->no, _nodes[j]->no);
					removed_arcs++;
					continue;
				}
			}

			ExSbrpArcO a;
			a.from = _nodes[i];
			a.to = _nodes[j];
			a.index = k;
			a.cost = prob->GetDist(_nodes[i], _nodes[j]);
			_arcs.push_back(a);

			k++;
		}
	printf("arcs:%d removed:%d\n", (int)_arcs.size(), removed_arcs);

	_arcs_of.resize(_nodes.size());
	_arcs_in_of.resize(_nodes.size());
	_arcs_out_of.resize(_nodes.size());
	
	for (size_t i = 0; i < _arcs.size(); i++)
	{
		ExSbrpArcO* ar = &_arcs[i];
		_arcs_of[ar->from->no].push_back(ar);
		_arcs_of[ar->to->no].push_back(ar);

		_arcs_out_of[ar->from->no].push_back(ar);
		_arcs_in_of[ar->to->no].push_back(ar);
	}
	
	_paths.resize(prob->GetDriverCount());
	_pos_arcs_of.resize(_nodes.size());
	_pos_arcs_in.resize(_nodes.size());
	_pos_arcs_out.resize(_nodes.size());
 
	_matrix_of_arcs.resize(GetNodeCount());
	for (int i = 0; i < GetNodeCount(); i++)
		_matrix_of_arcs[i].resize(GetNodeCount());

	for (size_t i = 0; i < _arcs.size(); i++)
	{
		ExSbrpArcO* ar = &_arcs[i];
		_matrix_of_arcs[ar->from->no][ar->to->no] = ar;
	}
}

ExSbrpArcO* ExactSbrpGraphO::GetArc(int from, int to)
{
	return _matrix_of_arcs[from][to];
}

ExSbrpArcO* ExactSbrpGraphO::GetArcsPos(int from, int to)
{
	if (from == to) return NULL;
	return _map_pos_arcs[from * 10000 + to];
}

void ExactSbrpGraphO::AssignPositiveValues()
{
	_pos_arcs.clear();
	_map_pos_arcs.clear();
	for (int i = 0; i < GetNodeCount(); i++)
	{
		_pos_arcs_of[i].clear();
		_pos_arcs_out[i].clear();
		_pos_arcs_in[i].clear();
	}

	_is_integer = true;
	_cost = 0;
	for (int i = 0; i < GetArcCount(); i++)
	{
		ExSbrpArcO* ar = GetArc(i);
		_cost += ar->value * ar->cost;
		if (ar->value > EPSILON)
		{
			_pos_arcs.push_back(ar);
			_pos_arcs_of[ar->from->no].push_back(ar);
			_pos_arcs_of[ar->to->no].push_back(ar);

			_pos_arcs_in[ar->to->no].push_back(ar);
			_pos_arcs_out[ar->from->no].push_back(ar);

			int no1 = ar->from->no * 10000 + ar->to->no;
			_map_pos_arcs[no1] = ar;

			if (ar->value <= 1 - EPSILON) _is_integer = false;
		}
	}
}

void ExactSbrpGraphO::ShowValueArcs()
{
	for (int i = 0; i < GetArcCount(); i++)
	{
		ExSbrpArcO* ar = GetArc(i);
		printf("i:%d from:%d to:%d va:%.3lf\n", (int)i, ar->from->no, ar->to->no, ar->value);
	}
}

void ExactSbrpGraphO::ShowPosValueArcs()
{
	for (size_t i = 0; i < _pos_arcs.size(); i++)
	{
		ExSbrpArcO* ar = _pos_arcs[i];
		printf("i:%d from:%d to:%d cost:%lf va:%.3lf\n", (int)i, ar->from->no, ar->to->no, ar->cost, ar->value);
	}

}

void ExactSbrpGraphO::PrintGraph(char* filename)
{
	Network* n = new Network(GetNodeCount());
	for (size_t i = 0; i < _pos_arcs.size(); i++)
		n->AddArc(_pos_arcs[i]->from->no, _pos_arcs[i]->to->no, _pos_arcs[i]->index);
	for (size_t i = 0; i < _pos_arcs.size(); i++)
		n->SetArc(_pos_arcs[i]->index, 1, _pos_arcs[i]->value);
	n->PrintGraphViz(filename);
	delete n;
}

void ExactSbrpGraphO::ShowPaths()
{
	for (size_t i = 0; i < _paths.size(); i++)
	{
		printf("path:%d nb:%d nodes:", (int)i, (int)_paths[i].size());
		for (size_t j = 0; j < _paths[i].size(); j++)
			printf("%d-", _paths[i][j]->no);
		printf("\n");
	}
}

void ExactSbrpGraphO::MakePaths()
{
	//printf("MakePaths()\n");
	for (size_t i = 0; i < _pos_arcs.size(); i++)
		_pos_arcs[i]->walked_on = 0;

	int nb_paths = GetArcsOutPosCount(0);

	_paths.resize(nb_paths);
	_paths.clear();
	//printf("nboutarcs:%d size:%d\n", GetArcsOutOfCount(0),nb_paths);

	int cur_path = 0;
	for (int i = 0; i < GetArcsOutPosCount(0); i++)
	{
		ExSbrpArcO* ar = GetArcsOutPos(0, i);
		if (ar->walked_on == 1) continue;

		std::vector<Node*> vv;
		_paths.push_back(vv);

		_paths[cur_path].push_back(ar->from);

		ar->walked_on = 1;

		Node* from = ar->from;
		Node* to = ar->to;
		while (1)
		{
			_paths[cur_path].push_back(to);
			from = to;
			if (from->no == 0) break;
			ar = GetArcsOutPos(to->no, 0);
			to = ar->to;
			ar->walked_on = 1;
		}

		//printf("%d \n", to->no);
		cur_path++;
	}

}
