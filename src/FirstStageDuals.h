#ifndef FIRST_STAGE_DUALS
#define FIRST_STAGE_DUALS

#include "ExactSbrpGraph.h"
#include "NodeSBRP.h"
#include "DriverSBRP.h"
#include "ProblemDefinition.h"
#include <ilcplex/ilocplex.h>
#include <map>

class FirstStageDuals
{
  public:
	FirstStageDuals(ExactSbrpGraphO * graph);
	~FirstStageDuals();

	void InitModelPrimal(int scenario);
	void InitModelDual(int scenario);

	void UpdateBalanceRHSPrimal(int scenario);
	void UpdateObjectiveDual(int scenario);

	void UpdateFlowBounds(int scenario);

	void SolvePrimal();
	void SolveDual();

	//Containers for the dual variables
	std::vector<double> ksi;
	std::vector<double> phi;
	std::vector<double> pi;
	std::vector<double> v;

	int _scenario=-1;
	int _status;
	double _rec;
	bool Benders_feas=false;
	bool found_inf=false;

	std::vector<double> ray;

	struct dual_class_node
	{
		int no; //no from the graph
		double phi_val=0.0;
		double pi_val=0.0;
		double v_val=0.0;
		int phi_id=-1;
		int pi_id=-1;
		int v_id=-1;
	};
	struct dual_class_arc
	{
		int from=-1;
		int to=-1;
		double ksi_val=0.0;
		int ksi_id=-1;
		bool isPos = false;
		int index = -1;
	};

	std::vector<dual_class_node> nodes;
	std::vector<dual_class_arc> arcs;

	double _sol=0.0;

private:
	
	Prob* _prob;
	ExactSbrpGraphO * _graph;

	//----------------------------------------------------------------
	//PRIMAL MODEL
	//----------------------------------------------------------------
	IloEnv env;
	IloModel model;
	IloCplex cplex;
	IloNumVarArray flow;
	IloNumVarArray w_minus;
	IloNumVarArray w_plus;
	IloRangeArray flow_consts_vars;		//constraints for the flows
	IloRangeArray balance_consts_vars;		//constraints for the flows 
	IloRangeArray wp_consts_vars;		//constraints for the wp
	IloRangeArray wm_consts_vars;		//constraints for the wp
	//----------------------------------------------------------------
	//DUAL MODEL
	//----------------------------------------------------------------
	IloEnv env_dual;
	IloModel model_dual;
	IloObjective obj;
	IloCplex cplex_dual;
	IloNumVarArray phi_dual;
	IloNumVarArray pi_dual;
	IloNumVarArray v_dual;
	IloNumVarArray ksi_dual;

	IloRangeArray dual_flow_consts;		//constraints for the flows
	IloRangeArray dual_balance_consts;		//constraints for the flows 
	IloRangeArray dual_wp_consts;		//constraints for the wp
	IloRangeArray dual_wm_consts;		//constraints for the wp  
	//----------------------------------------------------------------

	//To store the values if needed
	IloNumArray ksi_values;
	IloNumArray phi_values;
	IloNumArray pi_values;
	IloNumArray v_values;

	std::vector<Node*> _path;

};

#endif