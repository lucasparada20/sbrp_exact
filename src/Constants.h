#ifndef CONSTANTS_H
#define CONSTANTS_H

#define EPSILON 10E-5
#define INF   99999999
#define INFINITE 99999999.9

#define NODE_TYPE_CUSTOMER 0x01
#define NODE_TYPE_START_DEPOT 0x02
#define NODE_TYPE_END_DEPOT 0x04
#define NODE_TYPE_PICKUP 0x03
#define NODE_TYPE_DROP 0x05

#define STATUS_UNDETERMINED 2		//we tried to solve but we couldn't
#define STATUS_SOLVED 1				//it has been solved!
#define STATUS_UNSOLVED 0			//not yet solved
#define STATUS_PACK_TIMEOUT 3		//got a packing timeout
#define STATUS_ROUTING_TIMEOUT 4	//got a timeout from a routing routine
#define STATUS_TIMEOUT 7			//got a timeout
#define STATUS_FEASIBLE 5			//the problem is feasible
#define STATUS_INFEASIBLE 6			//the problem is infeasible

#define INF_ROUTE_COST 1000000
#define UNASSIGNED_COST 1000000


#endif
