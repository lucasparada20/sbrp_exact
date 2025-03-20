#ifndef COST_FUNCTION_SBRP
#define COST_FUNCTION_SBRP

class Sol;  // Forward declaration of Sol

#include "NodeSBRP.h"
#include "DriverSBRP.h"

class CostFunctionSBRP
{
	public:
		CostFunctionSBRP(){}
		~CostFunctionSBRP(){}

		double GetCost(Sol & s);
		double GetCost(Sol & s, Driver * d);
		void Update(Sol & s);
		void Update(Sol & s, Driver * d);
		void Show(Sol * s, Driver * d);
	private:
		std::vector<Node*> path;
};


#endif
