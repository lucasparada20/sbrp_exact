#include "RouteFeasibility.h"
#include "Parameters.h"
#include <algorithm>    // std::reverse
#include <algorithm>

bool RouteFeasibility::IsFeasible(Prob* prob, std::vector<Node*>& path)
{
 
 	std::vector<int> SortedScenarios = prob->GetVectorOfSortedScenarios();

	for (int e = 0; e < prob->GetScenarioCount(); e++)
		if (!IsFeasible(prob, path, SortedScenarios[e]))
			return false;
	return true;
}

bool RouteFeasibility::IsFeasible(Prob* prob, std::vector<Node*> & path, int e) //e is the scenario
{
	int Q = prob->GetDriver(0)->capacity;
	int min_lambda = 0;
	int sum_lambda = 0;
	int max_mu = 0;
	int sum_mu = 0;
	for (size_t i=0;i<path.size();i++)
	{
		Node* n = path[i];
		if(n->type != NODE_TYPE_CUSTOMER) continue;//ignore depots

		int lambda_i = std::max(-Q, n->demands[e] - n->w_plus);
		int mu_i = std::min(Q, n->demands[e] + n->w_minus);
		sum_lambda += lambda_i;
		min_lambda = std::min(sum_lambda, min_lambda);
		sum_mu += mu_i;
		max_mu = std::max(sum_mu, max_mu);

		int lb = sum_lambda - min_lambda;
		int ub = sum_mu + Q - max_mu;
		if(lb > ub)
			return false;
	}
	return true;
}

double RouteFeasibility::RecourseCost(Prob* prob, std::vector<Node*>& path)
{
	double cost = 0;
	int Q = prob->GetDriver(0)->capacity;
	int BigM = 9999;
	int * prev = new int[Q+1];
	int * next = new int[Q+1];
	for(int e = 0; e < prob->GetScenarioCount(); e++)
	{
		for(int q=0;q<=Q;q++)
			next[q] = prev[q] = 0;

		for(int i=path.size()-2;i>=1;i--)
		{
			Node * n = path[i];
			for(int q=0;q<=Q;q++)
			{
				int best = BigM;
				for(int u=0;u<=n->w_plus;u++)
					if(q+n->demands[e]-u >= 0 && q+n->demands[e]-u<=Q)
						best = std::min(best, u+prev[q+n->demands[e]-u]);

				for(int u=0;u<=n->w_minus;u++)
					if(q+n->demands[e]+u >= 0 && q+n->demands[e]+u<=Q)
						best = std::min(best, u+prev[q+n->demands[e]+u]);

				next[q] = best;
			}
			int * temp = next;
			next = prev;
			prev = temp;
		}

		int cost2 = BigM;
		for(int q=0;q<=Q;q++)
			cost2 = std::min(cost2 , prev[q]);
		cost += cost2;

		//printf("scenario:%d cost:%d\n", e, cost2);
	}
 
  cost = (cost/prob->GetScenarioCount()) * Parameters::GetCminEpsilon();
	delete [] prev;
	delete [] next;
	return cost;
}

double RouteFeasibility::RecourseCost(Prob* prob, std::vector<Node*>& path, int sce)
{
  int Q = prob->GetDriver(0)->capacity;
	int BigM = 9999;
	std::vector< std::vector<int> > memo(path.size());
	for(size_t i=0;i<path.size();i++)
		memo[i].resize(Q+1, 0);

	for(int i=path.size()-2;i>=1;i--)
	{
		Node * n = path[i];
		for(int q=0;q<=Q;q++)
		{
			int best = BigM;
			for(int u=0;u<=n->w_plus;u++)
				if(q+n->demands[sce]-u >= 0 && q+n->demands[sce]-u<=Q)
					best = std::min(best, u+memo[i+1][q+n->demands[sce]-u]);

			for(int u=0;u<=n->w_minus;u++)
				if(q+n->demands[sce]+u >= 0 && q+n->demands[sce]+u<=Q)
					best = std::min(best, u+memo[i+1][q+n->demands[sce]+u]);

			memo[i][q] = best;      
		}
	}

	int cost = BigM;
	for(int q=0;q<=Q;q++)
		cost = std::min(cost , memo[1][q]);
	return cost;
}
