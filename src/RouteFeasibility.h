
#ifndef ROUTE_FEASIBILITY
#define ROUTE_FEASIBILITY

#include <vector>
#include "NodeSBRP.h"
#include "DriverSBRP.h"
#include "ProblemDefinition.h"


class RouteFeasibility
{
public:

	static bool IsFeasible(Prob* prob, std::vector<Node*>& path);
  static bool IsFeasible(Prob* prob, std::vector<Node*>& path, int e); //e is the scenario

	static double RecourseCost(Prob* prob, std::vector<Node*>& path);
	static double RecourseCost(Prob* prob, std::vector<Node*>& path, int scenario);
	static std::vector<double> MissingBicycles(Prob* prob, std::vector<Node*>& path);

private:
	
};

#endif
