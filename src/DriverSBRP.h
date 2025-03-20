#ifndef DRIVER_SBRP
#define DRIVER_SBRP

class Driver
{
	public:
		int id;				//from 0 to nb drivers - 1
		int StartNodeID;	//
		int EndNodeID;		//
		int capacity;
		int sum_demand;
		double curDistance;
		double curRecourse;
		bool is_feasible;

		int GetStartNodeId(){return StartNodeID;}
		int GetEndNodeId(){return EndNodeID;}
		
		void Show(){ printf("Driver:%d cap:%d\n",id,capacity); }
};

#endif
