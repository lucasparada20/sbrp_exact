#ifndef INSRMV_SBRP
#define INSRMV_SBRP

#include "NodeSBRP.h"
#include "DriverSBRP.h"
#include "Move.h"
#include "ProblemDefinition.h"
#include "Solution.h"
#include "Constants.h"

class InsRmvMethodSBRP
{
	public:
		InsRmvMethodSBRP(Prob & prob) { _prob = &prob;}  ;

		void InsertCost(Sol & s, Node * n, Driver * d, Move & m);
		void ApplyInsertMove(Sol & s, Move & m);
		void RemoveCost(Sol & s, Node * n, Move & m);
		void CheckMove(Sol & s, Move & m){};
		void FillInsertionList(Sol & s, std::vector<Node*> & list);

	private:
		std::vector<Node*> path;
		Prob * _prob;
};


#endif
