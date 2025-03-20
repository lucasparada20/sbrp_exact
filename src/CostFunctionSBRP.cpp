#include "CostFunctionSBRP.h"
#include "Constants.h"
#include "RouteFeasibility.h"
#include "Solution.h"

double CostFunctionSBRP::GetCost(Sol & s)
{
	double cost = 0;
	for(int i = 0 ; i < s.GetDriverCount() ; i++)
		cost += GetCost(s, s.GetDriver(i));
	return cost + s.GetUnassignedCount()*UNASSIGNED_COST;
}

double CostFunctionSBRP::GetCost(Sol & s,Driver * d)
{
	double dist = 0;
	Node * prev = s.GetNode( d->StartNodeID );
	path.clear();
	while(prev->type != NODE_TYPE_END_DEPOT)
	{
		Node * next = s.Next[ prev->id ];
		dist += s.GetDist(prev,next);
		path.push_back(prev);
		prev = next;
	}
	path.push_back(prev);

	if(!RouteFeasibility::IsFeasible(s.GetProb(), path))
		return INF_ROUTE_COST;

	return dist + RouteFeasibility::RecourseCost(s.GetProb(), path);
}


void CostFunctionSBRP::Update(Sol & s)
{
	double sumd = 0, sumr = 0;
	for(int i = 0 ; i < s.GetDriverCount() ; i++)
	{
		Driver * d = s.GetDriver(i);
		Update(s, d);
		sumd += d->curDistance;
		sumr += d->curRecourse;
	}
	s.SetTotalDistances(sumd);
	s.SetTotalRecourse(sumr);
}

void CostFunctionSBRP::Update(Sol & s, Driver * d)
{
	d->curDistance = 0;

	Node * prev = s.GetNode( d->StartNodeID );
	path.clear();
	while(prev->type != NODE_TYPE_END_DEPOT)
	{
		Node * next = s.Next[ prev->id ];
		path.push_back(prev);
		d->curDistance += s.GetDist(prev,next);
		prev = next;
	}
	path.push_back(prev);
	d->curRecourse = RouteFeasibility::RecourseCost(s.GetProb(), path);
}

void CostFunctionSBRP::Show(Sol * s, Driver * d)
{
	GetCost(*s,d);
	printf("Route:%d cost:%.2lf rec:%.2lf len:%d:", d->id, d->curDistance, d->curRecourse, s->RoutesLength[d->id]);
	Node * cur = s->GetNode( d->StartNodeID );
	while(cur != NULL)
	{
		printf("%d-", cur->id);
		cur = s->Next[ cur->id ];
	}
	printf("\n");
}
