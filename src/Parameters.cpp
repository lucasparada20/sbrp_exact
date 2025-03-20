#include "Parameters.h"
#include "Constants.h"
#include <stdio.h>
#include <stdlib.h>
#include <ctime>
#include <string>
#include <string.h>
#include <vector>
#include <limits>
#include<algorithm>

double Parameters::delta = 0;
double Parameters::Epsilon = 0;
double Parameters::Cmin = 0;
double Parameters::L = 0;
int Parameters::WorstScenario = -1;
int Parameters::OppositeScenario = -1;
int Parameters::opt_cut_type = OPT_CUT_TYPE_PL;

char output_file_temp[100];
char re_file_temp[100];
char solution_file_temp[100];
char instance_file_temp[100];
char algorithm_temp[100];
char instance_type_temp[100];
char initial_solution_file_temp[100];

char * Parameters::output_file = NULL;
char * Parameters::re_file = NULL;
char * Parameters::solution_file = NULL;
char * Parameters::instance_file = NULL;
char * Parameters::algorithm = NULL;
char * Parameters::instance_type = NULL;
char * Parameters::initial_solution_file = NULL;


void Parameters::Read(int arg, char ** argv)
{
	printf("Reading parameters\n");
	for(int i=0;i<arg;i++)
	{
		char * first = strtok (argv[i]," ;=");
		char * second = strtok (NULL, " ;=");
		printf ("Parameter:%s value:%s\n",first,second);
		if(second == NULL) continue;

		if(strcmp(first, "instance_file") == 0)
		{
			strcpy(instance_file_temp, second);
			instance_file = instance_file_temp;
		}
		else if(strcmp(first, "re_file") == 0)
		{
			strcpy(re_file_temp, second);
			re_file = re_file_temp;
		}
		else if(strcmp(first, "epsilon") == 0)
		{
			sscanf(second, "%lf", &Epsilon);
		}
		else if(strcmp(first, "delta") == 0)
		{
			sscanf(second, "%lf", &delta);
		}
		else if(strcmp(first, "opt_cuts") == 0)
		{
			sscanf(second, "%d", &opt_cut_type);
		}
		else if(strcmp(first, "output_file") == 0)
		{
			strcpy(output_file_temp, second);
			output_file = output_file_temp;
		}
		else if(strcmp(first, "instance_type") == 0)
		{
			strcpy(instance_type_temp, second);
			instance_type = instance_type_temp;
		}				
		else if(strcmp(first, "algorithm") == 0)
		{
			strcpy(algorithm_temp, second);
			algorithm = algorithm_temp;
		}
		else if(strcmp(first, "initial_solution_file") == 0)
		{
			strcpy(initial_solution_file_temp, second);
			initial_solution_file = initial_solution_file_temp;
		}	
		
	}
}
