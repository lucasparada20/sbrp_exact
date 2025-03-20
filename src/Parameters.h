#ifndef _PARAMETERS_H
#define _PARAMETERS_H

#include <time.h>
#include <math.h>

#define OPT_CUT_TYPE_PL 3
#define OPT_CUT_TYPE_BENDERS 6
#define OPT_CUT_TYPE_HYBRID 7

class Parameters
{

public:
	void Read(int arg, char ** argv);

	static double GetEpsilon() { return Epsilon; }
	static double GetDelta() { return delta; }
	static double GetCmin() { return Cmin; }
	static double GetCminEpsilon(){return ceil(Cmin*Epsilon);}
	static double GetL() { return L; }
	static int GetWorstScenario() { return WorstScenario; }
	static int GetOppositeScenario() { return OppositeScenario; }
	static int GetTypeOfOptimalityCuts(){return opt_cut_type;}
	static char* GetAlgorithm(){return algorithm;}
	static char* GetInstanceFileName(){return instance_file;}
	static char* GetInstanceType(){return instance_type;}
	static char* GetInitialSolFileName(){return initial_solution_file;}
	static char* GetReFileName(){return re_file;}
	static char* GetOutputFileName(){return output_file;}


	static void SetEpsilon(double e) { Epsilon = e; }
	static void SetDelta(double g) { delta = g; }
	static void SetCmin(double c) { Cmin = c; }
	static void SetL(double ele) { L = ele; }
	static void SetWorstScenario(int e) { WorstScenario = e; }
	static void SetOppositeScenario(int e) { OppositeScenario = e; }

	static void SetTypeOfOptimalityCuts(int t){opt_cut_type = t;}
private:

	static double delta;
	static double Epsilon;
	static double Cmin;
	static double L;
	static int WorstScenario;
	static int OppositeScenario;
	static int opt_cut_type;

	static char * output_file;
	static char * re_file;
	static char * solution_file;
	static char * instance_file;
	static char * algorithm;
	static char * instance_type;
	static char * initial_solution_file;	
};


#endif
