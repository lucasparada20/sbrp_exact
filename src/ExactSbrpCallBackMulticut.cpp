#include "ExactSbrpCallBackMulticut.h"

ExactSbrpMulticutUserCutCallBack::ExactSbrpMulticutUserCutCallBack(IloEnv env, ExactSbrpGraphO * graph,
													IloNumVarArray x, IloNumVarArray thetas, ExactSbrpSepMulticut * sep) :
													IloCplex::UserCutCallbackI(env), _graph(graph), _x(x),_thetas(thetas),_sep(sep),
													added_constraints(env),add_constraints(false)
{

}
ExactSbrpMulticutUserCutCallBack::~ExactSbrpMulticutUserCutCallBack()
{
	added_constraints.end();
}
void ExactSbrpMulticutUserCutCallBack::main()
{
	//printf("ExactSbrpCallBackMulticut obj:%.2lf\n", (double)getObjValue());
	if(!isAfterCutLoop()) return;

	IloNumArray values(getEnv());
	getValues(values,_x);
	for(int i = 0 ; i < _graph->GetArcCount();i++)
	{
		ExSbrpArcO * arc = _graph->GetArc(i);
		arc->value = (double)values[arc->index];
	}
	values.clear();
	getValues(values,_thetas);

	_graph->SetThetasSize( _thetas.getSize() );
	for(int i=0;i<_thetas.getSize();i++)
		_graph->SetTheta(i, (double)values[i]);

	values.end();
	_graph->AssignPositiveValues();

	IloRangeArray inq(getEnv());
	_sep->SeperateFrac(inq);
	
	if(inq.getSize() >= 1 && add_constraints)
		added_constraints.add(inq);
		
	for(int i=0;i<inq.getSize();i++)
	{
		add(inq[i]);
	}
	
	inq.end();
}



ExactSbrpMulticutLazyCallBack::ExactSbrpMulticutLazyCallBack(IloEnv env, ExactSbrpGraphO * graph,
											 IloNumVarArray x, IloNumVarArray thetas, ExactSbrpSepMulticut * sep) :
											 IloCplex::LazyConstraintCallbackI(env), _graph(graph), _x(x),_thetas(thetas),_sep(sep),
											 added_constraints(env),add_constraints(false)
{

}
ExactSbrpMulticutLazyCallBack::~ExactSbrpMulticutLazyCallBack()
{
	added_constraints.end();
}
void ExactSbrpMulticutLazyCallBack::main()
{
	//printf("ExactSbrpLazyCallBackO obj:%.2lf\n", (double)getObjValue());
	IloNumArray values(getEnv());
	getValues(values,_x);
	for(int i = 0 ; i < _graph->GetArcCount();i++)
	{
		ExSbrpArcO * arc = _graph->GetArc(i);
		arc->value = (double)values[arc->index];
	}
	values.clear();
 	getValues(values,_thetas);

	_graph->SetThetasSize( _thetas.getSize() );
	for(int i=0;i<_thetas.getSize();i++)
		_graph->SetTheta(i, (double)values[i]);

	values.end();
	_graph->AssignPositiveValues();
	//_graph->ShowPosValueArcs();

	IloRangeArray inq(getEnv());
	_sep->SeperateInt(inq);
	if (inq.getSize() == 0)
	{
		_graph->MakePaths();

		if (inq.getSize() == 0)
			switch(Parameters::GetTypeOfOptimalityCuts())
			{

				case OPT_CUT_TYPE_BENDERS: // Number 6
					_sep->SeperateBendersCut(inq);
					break;
			}
	}
	
	if(inq.getSize() >= 1 && add_constraints)
		added_constraints.add(inq);
	for(int i=0;i<inq.getSize();i++)
		add(inq[i]);	
	inq.end();
}
