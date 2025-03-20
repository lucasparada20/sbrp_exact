#ifndef SEQUENTIAL_INSERTION
#define SEQUENTIAL_INSERTION

#include "Solution.h"
#include "Move.h"
#include "InsRmvMethodSBRP.h"

class SequentialInsertionSBRP 
{
	public:
		SequentialInsertionSBRP(InsRmvMethodSBRP & insrmv): _insrmv(insrmv){}

	void Insert(Sol & s);

	private:
		InsRmvMethodSBRP & _insrmv;

};


#endif