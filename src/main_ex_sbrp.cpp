#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <csignal>
#include <string.h> //strcm

#include "NodeSBRP.h"
#include "DriverSBRP.h"
#include "LoadSBRP.h"

#include "Parameters.h"
#include "ProblemDefinition.h"
#include "ExactSbrp.h"
#include "RecourseLowerBound.h"
#include "RouteFeasibility.h"
#include "DriverCountBppLb.h"
#include "Solution.h"
#include "CostFunctionSBRP.h"
#include "SequentialInsertionSBRP.h"
#include "ExactSbrpMulticut.h"
#include "UbManager.h"

#include <ctime>
#include <fstream>
#include <iostream>

int main(int arg, char ** argv)
{
	
	if(arg == 1)
	{
		printf("usage: executable, instance_file, epsilon, delta, cuts_type, instance_type, algorithm optional:initial_solution_file \n");
		printf("Instance_type: dins or pcg\n");
		printf("Cuts: P&L=3 Benders=6 Hybrid=7\n");
		printf("Algorithm: DL-shaped=dl Multicut=m(Only accepts Benders Opt Cuts=6)\n");
		printf("exiting.\n");
		exit(1);
	}

	Parameters param;
	param.Read(arg,argv);
	
	if(std::strcmp(Parameters::GetAlgorithm(),"dl")!=0 && std::strcmp(Parameters::GetAlgorithm(),"m")!=0)
	{
		printf("No algorithm specificied. Exiting ...\n");
		printf("Algorithm: DL-shaped=dl Multicut=m\n");
		exit(1);			
	}		
		
	if(std::strcmp(Parameters::GetAlgorithm(),"m")==0 && Parameters::GetTypeOfOptimalityCuts() != 6)
	{
		printf("Multicut requires Benders optimality cuts. Exiting ...\n");
		printf("Reminder for optimality cuts: P&L=3 Benders=6 Hybrid=7\n");
		exit(1);
	}

	Prob pr;
	LoadSBRP Load;

	if(std::strcmp(Parameters::GetInstanceType(),"dins")==0) 
		Load.Load_dins(pr, Parameters::GetInstanceFileName());
	else if(std::strcmp(Parameters::GetInstanceType(),"pcg")==0) 
		Load.Load_pcg(pr, Parameters::GetInstanceFileName());
	else {
		printf("Wrong file type. Exiting ... \n"); exit(1);
	}
 
	if(std::strcmp(Parameters::GetInstanceType(),"dins")==0 && 
		std::strcmp(Parameters::GetAlgorithm(),"dl")==0)
		{
			UbManager ub_manager;
			ub_manager.Load((char*)"instances_dins/all_upper_bounds.txt");
			double ub = ub_manager.GetUpperBound(Parameters::GetInstanceFileName(), Parameters::GetEpsilon(), Parameters::GetDelta());
			printf("UbFile:%.1lf\n", ub);
			pr.SetUpperBound( ub );
		}
 
	CostFunctionSBRP cost_func;
	
	Sol sol(&pr,&cost_func);
	sol.PutAllNodesToUnassigned();
	
	InsRmvMethodSBRP method(pr);
	SequentialInsertionSBRP seq(method);
	
	RecourseLowerBound::SetWorstScenario(&pr);
	RecourseLowerBound::SortFromWorstScenarios(&pr);
	
	seq.Insert(sol); //Sequential insertion of nodes to build an initial solution and store in sol
	
	int min_drv_cnt1 = RecourseLowerBound::GetDriverCount(&pr);
	DriverCountBppLb driver_lb;
	int min_drv_cnt2 = driver_lb.GetLB(&sol); //The set covering model uses the initial solution
	pr.SetDriverCountLB( std::max( min_drv_cnt1, min_drv_cnt2) );
	
	RecourseLowerBound rec;	
	double lb_recourse = rec.CalculateWithMinDriverCount(&pr);
	printf("Recourse LB with min drivers:%.2lf MinDriverCount1:%d MinDriverCount2:%d\n", lb_recourse, min_drv_cnt1,min_drv_cnt2);
  
	//The Dissagregated Integer L-Shaped Method
	if(std::strcmp(Parameters::GetAlgorithm(),"dl")==0) 
	{
		ExactSbrpO ex;
		ex.max_time = 300;
		ex.Solve(&pr);
	}
	//The L-Shaped Multicut Method
	if(std::strcmp(Parameters::GetAlgorithm(),"m")==0)
	{
		ExactSbrpMulticut exm;
		exm.max_time = 300;
		exm.Solve(&pr);		
	}

	return 0;
}


