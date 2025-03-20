
#ifndef UB_MANAGER_H
#define UB_MANAGER_H

#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <string>
#include <string.h>

/*
	Class that load a list of upper bounds for different problems

	Just pass the problem name (or path with name)
	to get the upper bound of the problem
*/
class UbManager_Item
{
public:
	std::string name;
	double delta;
	double epsilon;
	double ub;
};

class UbManager
{
	public:
		UbManager(){}
		~UbManager(){}


		void Load(char * filename);
		double GetUpperBound(char * instance, double epsilon, double delta);

	private:
		std::vector<UbManager_Item> _bounds;
};




#endif
