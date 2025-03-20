#ifndef NODE_CVRP
#define NODE_CVRP

#include "Constants.h"
#include "Parameters.h"
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

class Node
{
	public:
		Node() : id(-1), origin_id(-1), stationcapacity(0), occupancy(0), distID(-1), no(0), type(0), demands(0), w_plus(0), w_minus(0), IsNotInCycle(true), arc_index(-1) {}

		int id;			//from 0 to n-1
		int origin_id;  //ID at the begining, it is never modified, only customers have a valid unique ID
		int stationcapacity;   // Station capacity
		int occupancy;   //Initial occupational level of the station
		int distID;		//indicate which line and column that contains the distance/time in the matrices
		int no;			//personnal identifier, put what you want
		char type;		//type of the node
		std::vector<int> demands; //scenario demands

		int w_plus;
		int w_minus;
    
		bool IsNotInCycle;

		int arc_index; //used for the induced graph
    
		void Show()
		{
			if(type == NODE_TYPE_CUSTOMER)
			{
				printf("Node:%d type:Cust cap:%d occ:%d wp:%d wm:%d demands:%d dmd:", id, stationcapacity, occupancy, w_plus, w_minus, (int)demands.size());
				for(size_t e=0;e<demands.size();e++)
					printf("%d ", demands[e]);
				printf("\n");
			}
			else if(type == NODE_TYPE_START_DEPOT)
			{
				printf("Node:%d type:stadepot\n", id);
			}
			else if(type == NODE_TYPE_END_DEPOT)
			{
				printf("Node:%d type:enddepot\n", id);
			}
		}
		void Show(const int scenario)
		{
			if(type == NODE_TYPE_CUSTOMER)
			{
				printf("Node:%d type:Cust cap:%d occ:%d wp:%d wm:%d e:%d dmd:", id, stationcapacity, occupancy, w_plus, w_minus, scenario);
				//for(size_t e=0;e<demands.size();e++)
					printf("%d ", demands[scenario]);
				printf("\n");
			}
			else if(type == NODE_TYPE_START_DEPOT)
			{
				printf("Node:%d type:stadepot\n", id);
			}
			else if(type == NODE_TYPE_END_DEPOT)
			{
				printf("Node:%d type:enddepot\n", id);
			}			
		}

		void UpdateW()
		{
			int cap = (int)ceil(Parameters::GetDelta() * stationcapacity);
			w_minus = std::min(cap, occupancy);
			w_plus = std::min(cap, stationcapacity - occupancy);
		}

		bool IsCustomer()
		{
			return (type & NODE_TYPE_CUSTOMER) == NODE_TYPE_CUSTOMER;
		}
};

#endif
