
#include "InsRmvMethodSBRP.h"
#include "RouteFeasibility.h"
#include "RecourseLowerBound.h"
#include "Parameters.h"
#include <algorithm>

void InsRmvMethodSBRP::FillInsertionList(Sol & s, std::vector<Node*> & list)
{
	list.clear(); list.reserve( _prob->GetCustomerCount() );
	for(int i=0;i<s.GetUnassignedCount();i++)
		list.push_back(s.GetUnassigned(i));
}

void InsRmvMethodSBRP::InsertCost(Sol & s, Node * n, Driver * d, Move & mo)
{
	//printf("InsertCost n:%d d:%d\n", n->id, d->id);
	mo.IsFeasible = false;
	mo.DeltaCost = INFINITE;
	mo.n = n;
	mo.to = d;

	path.clear();
	Node * prev = s.GetNode( d->StartNodeID );
	path.push_back(prev);
	path.push_back(n);
	while(prev->type != NODE_TYPE_END_DEPOT)
	{
		Node * next = s.Next[ prev->id ];
		path.push_back(next);
		prev = next;
	}
	if(RecourseLowerBound::GetDriverCount(s.GetProb(), path) >= 2) return;

	//s.Show(d);

	int pos = 0;
	prev = s.GetNode( d->StartNodeID );
	while(prev->type != NODE_TYPE_END_DEPOT)
	{
		Node * next = s.Next[ prev->id ];

		double deltaDist = s.GetDist(prev,n) + s.GetDist(n,next) - s.GetDist(prev,next);
		double newcost = deltaDist;

		//printf("Prev:%d next:%d\n", prev->id, next->id);
		if(prev->type != NODE_TYPE_START_DEPOT)
			iter_swap(path.begin() + pos, path.begin() + pos + 1);

		//for(size_t i=0;i<path.size();i++)
		//	printf("%d-", path[i]->id);
		//printf("\n");


		if(newcost < mo.DeltaCost && RouteFeasibility::IsFeasible(s.GetProb(), path))
		{
			double rec = RouteFeasibility::RecourseCost(s.GetProb(), path);
			if(newcost + rec < mo.DeltaCost)
			{
				mo.DeltaDistance = deltaDist;
				mo.DeltaCost = newcost + rec;
				mo.IsFeasible = true;
				mo.move.prev = prev;
			}
		}

		prev = next;
		pos++;
	}
}

void InsRmvMethodSBRP::ApplyInsertMove(Sol & s, Move & m)
{
	if(m.from != NULL)
	{
		s.Remove(m.n);
	}
	else if(s.IsUnassigned(m.n))
		s.RemoveFromUnassigneds(m.n);

	s.InsertAfter(m.n, m.move.prev);
}
