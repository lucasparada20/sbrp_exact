#include "ExactSbrpCallBacks.h"


ExactSbrpUserCutCallBackO::ExactSbrpUserCutCallBackO(IloEnv env, ExactSbrpGraphO * graph,
													IloNumVarArray x, IloNumVarArray thetas, ExactSbrpSepO * sep) :
													IloCplex::UserCutCallbackI(env), _graph(graph), _x(x),_thetas(thetas),_sep(sep),
													added_constraints(env),add_constraints(false)
{

}
ExactSbrpUserCutCallBackO::~ExactSbrpUserCutCallBackO()
{
	added_constraints.end();
}
void ExactSbrpUserCutCallBackO::main()
{
	//printf("ExactSbrpUserCutCallBackO obj:%.2lf\n", (double)getObjValue());
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
	_sep->SeparateFrac(inq);
	
	if(inq.getSize() >= 1 && add_constraints)
		added_constraints.add(inq);
		
	for(int i=0;i<inq.getSize();i++)
	{
		add(inq[i]);
	}
	inq.end();
}



ExactSbrpLazyCallBackO::ExactSbrpLazyCallBackO(IloEnv env, ExactSbrpGraphO * graph,
											 IloNumVarArray x, IloNumVarArray thetas, ExactSbrpSepO * sep) :
											 IloCplex::LazyConstraintCallbackI(env), _graph(graph), _x(x),_thetas(thetas),_sep(sep),
											 added_constraints(env),add_constraints(false)
{

}
ExactSbrpLazyCallBackO::~ExactSbrpLazyCallBackO()
{
	added_constraints.end();
}
void ExactSbrpLazyCallBackO::main()
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
	_sep->SeparateInt(inq);
	if (inq.getSize() == 0)
	{
		_graph->MakePaths();

		if (inq.getSize() == 0)
			switch(Parameters::GetTypeOfOptimalityCuts())
			{
				case OPT_CUT_TYPE_PL: // Number 3
					_sep->SeparateDissagregatedRecourse(inq);
					break;
				case OPT_CUT_TYPE_BENDERS: // Number 6
					_sep->SeparateBendersCut(inq);
					break;
				case OPT_CUT_TYPE_HYBRID: // Number 7
					if(_graph->GetPathCount()==1) _sep->SeparateBendersCut(inq);
					else _sep->SeparateDissagregatedRecourse(inq);
					break;
			}
	}
	
	if(inq.getSize() >= 1 && add_constraints)
		added_constraints.add(inq);
	for(int i=0;i<inq.getSize();i++)
		add(inq[i]);	
	inq.end();
}

