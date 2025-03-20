

#ifndef MOVE_H_
#define MOVE_H_
#include <stdlib.h>
#include "NodeSBRP.h"
#include "DriverSBRP.h"

class MoveSBRP
{
	public:
	Node * prev;
};

class Move
{
	public:
		Node * n;
		Driver * to;
		Driver * from;
		double DeltaCost;
		double DeltaDistance;
		bool IsFeasible;
		MoveSBRP move;

		Move(): n(NULL), to(NULL),from(NULL), DeltaCost(INFINITE), IsFeasible(false){}


		bool operator < (const Move & m) const
		{
			return (DeltaCost < m.DeltaCost);
		}

		bool operator < (const Move * m) const
		{
			return (DeltaCost < m->DeltaCost);
		}
};



#endif
