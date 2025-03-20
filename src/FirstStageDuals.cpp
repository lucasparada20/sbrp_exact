#include "FirstStageDuals.h"
#include "RouteFeasibility.h"

FirstStageDuals::FirstStageDuals(ExactSbrpGraphO* graph):_graph(graph)
{
  _prob = graph->GetProblem();
}

FirstStageDuals::~FirstStageDuals()
{
	env.end();
	env_dual.end();
}

void FirstStageDuals::InitModelDual(int scenario)
{
	nodes.clear();
	arcs.clear();
	_scenario=scenario;
	double Q = _prob->GetDriver(0)->capacity;
	
	model_dual = IloModel(env_dual);
	phi_dual = IloNumVarArray(env_dual,_graph->GetNodeCount()-1,-IloInfinity, IloInfinity, ILOFLOAT);
	pi_dual = IloNumVarArray(env_dual,_graph->GetNodeCount()-1,-IloInfinity, 0, ILOFLOAT);
	v_dual =IloNumVarArray(env_dual,_graph->GetNodeCount()-1,-IloInfinity, 0, ILOFLOAT);
	//ONLY pos_acrs arc <=0. The rest are free variables!
	ksi_dual =IloNumVarArray(env_dual,_graph->GetArcCount(),-IloInfinity, 0, ILOFLOAT);
		
	dual_flow_consts = IloRangeArray(env_dual);
	dual_balance_consts = IloRangeArray(env_dual);
	dual_wp_consts = IloRangeArray(env_dual);
	dual_wm_consts = IloRangeArray(env_dual);
	
	//Dual for the nodes: Each non-depot node has 3 dual variables.
	std::map<int,int> map_nodes;
	for (int i = 1; i < _graph->GetNodeCount(); i++)
	{
		Node* n = _graph->GetNode(i);
		if(n->type != NODE_TYPE_CUSTOMER) continue;	
		char name[40];
		sprintf(name, "phi%d", n->no);
		phi_dual[i-1].setName(name);
		sprintf(name, "pi%d", n->no);
		pi_dual[i-1].setName(name);
		sprintf(name, "v%d", n->no);
		v_dual[i-1].setName(name);
		
		map_nodes.insert(std::pair<int,int>(n->no,i-1));
		
		dual_class_node node;
		node.no=n->no;
		node.phi_id=phi_dual[i-1].getId();
		node.pi_id=pi_dual[i-1].getId();
		node.v_id=v_dual[i-1].getId();
		nodes.push_back(node);
		
		//Wm and Wp constraints
		IloExpr expr(env_dual);
		expr += (-phi_dual[i-1]+pi_dual[i-1]);
		dual_wm_consts.add(expr<=1);
		model_dual.add(dual_wm_consts);
		expr.clear();
		
		expr += (phi_dual[i-1]+v_dual[i-1]);
		dual_wp_consts.add(expr<=1);
		model_dual.add(dual_wp_consts);
		expr.end();		
    }
	
	//Dual for the arcs: Each arc has 1 dual variable.
	std::map<int,int> map_arcs;
	for(int i=0;i<_graph->GetArcCount();i++)
	{
		ExSbrpArcO* arc = _graph->GetArc(i);
		char name[40];
		sprintf(name, "ksi%d_%d",arc->from->no,arc->to->no);
		ksi_dual[i].setName(name);
		dual_class_arc a;
		a.from=arc->from->no;
		a.to=arc->to->no;
		a.ksi_id=(int)ksi_dual[i].getId();
		if(arc->value > 0.0) 
		{
			a.isPos = true;
		}
		a.index = arc->index;
		arcs.push_back(a);
		
		map_arcs.insert(std::pair<int,int>(arc->index,i));
	}

	//Flow on non-depot node constraints
	for(int i=0;i<_graph->GetArcCount();i++)
	{
		ExSbrpArcO* arc = _graph->GetArc(i);
		if(arc->from->no ==0 || arc->to->no == 0) continue;
		auto it_node_from = map_nodes.find(arc->from->no);
		auto it_node_to = map_nodes.find(arc->to->no);
		auto it_arc = map_arcs.find(arc->index);
		
		IloExpr expr(env_dual);
		expr += ( -phi_dual[it_node_from->second] + phi_dual[it_node_to->second] + ksi_dual[it_arc->second] );

		dual_flow_consts.add(expr<=0);
		model_dual.add(dual_flow_consts);
		expr.end();
	}
	
	//depot outgoing, incoming arc constraints
	for(int i=0;i<_graph->GetArcCount();i++)
	{
		ExSbrpArcO* arc = _graph->GetArc(i);
		if(arc->from->no ==0)
		{
			auto it_arc = map_arcs.find(arc->index);
			auto it_node = map_nodes.find(arc->to->no);
			IloExpr expr(env_dual);
			expr += ( phi_dual[it_node->second] + ksi_dual[it_arc->second]);
			dual_flow_consts.add(expr<=0);
			model_dual.add(dual_flow_consts);
			//model_dual.add(expr<=0);
			expr.end();
		}else if(arc->to->no==0)
		{
			auto it_arc = map_arcs.find(arc->index);
			auto it_node = map_nodes.find(arc->from->no);
			IloExpr expr(env_dual);
			expr += ( -phi_dual[it_node->second] + ksi_dual[it_arc->second]);
			dual_flow_consts.add(expr<=0);
			model_dual.add(dual_flow_consts);
			//model_dual.add(expr<=0);
			expr.end();
		}
	}
	
	IloExpr expr(env_dual);
	for(int i = 1; i < _graph->GetNodeCount(); i++)
	{
		Node* n = _graph->GetNode(i);
		expr += (phi_dual[i-1]*n->demands[_scenario] + pi_dual[i-1]*n->w_minus + v_dual[i-1]*n->w_plus);
	}
	for(int i=0;i<_graph->GetArcCount();i++)
	{
		ExSbrpArcO* arc = _graph->GetArc(i);

		if(arc->value > 0.0)
		{
			auto it = map_arcs.find(arc->index);
			expr+=Q*ksi_dual[ it->second ];
		}
	}
	obj = IloObjective(env_dual,expr,IloObjective::Maximize,"Dual_obj");
	model_dual.add(obj);
	expr.end();
}

void FirstStageDuals::UpdateObjectiveDual(int scenario)
{
	_scenario=scenario;
	IloNumArray demands(env_dual,_graph->GetNodeCount()-1);
	for(int i=1;i<_graph->GetNodeCount();i++)
	{
		Node * n = _graph->GetNode(i);
		demands[i-1]=(double)n->demands[scenario];
	}
	obj.setLinearCoefs(phi_dual,demands);
	
	demands.end();
}

void FirstStageDuals::SolveDual() 
{
	IloCplex cplex_dual(model_dual);
	if(Benders_feas)
	{	
		cplex_dual.setParam(IloCplex::Param::Preprocessing::Presolve, 0);
		cplex_dual.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Primal);
		cplex_dual.setParam(IloCplex::Param::NodeAlgorithm, IloCplex::Primal);
		cplex_dual.setParam(IloCplex::Param::Threads,1);
	}
	
	cplex_dual.setParam(IloCplex::Param::Threads,1);
	cplex_dual.setOut(env_dual.getNullStream());
	cplex_dual.setWarning(env_dual.getNullStream());
	
	bool re = cplex_dual.solve();
	
	if(cplex_dual.getStatus() == IloAlgorithm::Unbounded)
	{

		int cplex_status = (int)cplex_dual.getCplexStatus();
		_status=cplex_status;
		_sol=9999999999;
        
		ksi.clear();
		phi.clear();
		pi.clear();
		v.clear();
		
		//printf("dual sce:%d status:%d sol:%.0lf\n",_scenario,cplex_status,_sol);
		
		IloNumVarArray var(env_dual);
        IloNumArray val(env_dual);
        cplex_dual.getRay(val, var);
		
		int cnt=0;
		for(size_t i=0;i<nodes.size();i++)
		{
			int _id=nodes[i].phi_id;
			for(int j=0;j<var.getSize();j++)
				if(_id==(int)var[j].getId())
				{
					nodes[i].phi_val = (double)val[j];
					break;
				}
					
			_id=nodes[i].pi_id;
			for(int j=0;j<var.getSize();j++)
				if(_id==(int)var[j].getId())
				{
					nodes[i].pi_val = (double)val[j];
					break;
				}
					
			_id=nodes[i].v_id;
			for(int j=0;j<var.getSize();j++)
				if(_id==(int)var[j].getId()) 
				{
					nodes[i].v_val = (double)val[j];
					break;
				}
			
			if(nodes[i].phi_val<0.0)cnt++;
			if(nodes[i].pi_val<0.0)cnt++;
			if(nodes[i].v_val<0.0)cnt++;
			phi.push_back(nodes[i].phi_val);
			pi.push_back(nodes[i].pi_val);
			v.push_back(nodes[i].v_val);
		}
		
		cnt=0;
		for(size_t i=0;i<arcs.size();i++)
		{
			int _id=arcs[i].ksi_id;
			for(int j=0;j<var.getSize();j++)
			{
				if(_id==(int)var[j].getId())
				{
					arcs[i].ksi_val = (double)val[j];
					cnt++;
					break;
				}					
										
			}
		}
		
		//printf("Added value to %d arcs from ray...\n",cnt);
		
		for(int i=0;i<_graph->GetArcCount();i++)
		{
			ExSbrpArcO* arc = _graph->GetArc(i);
			for(size_t j=0;j<arcs.size();j++)
				if(arc->index == arcs[j].index)
				{
					ksi.push_back(arcs[i].ksi_val);
					break;
				}
					
		}

		var.end();
		val.end();
		
	}else if( cplex_dual.getStatus() == IloAlgorithm::Optimal )
	{

		int cplex_status = (int)cplex_dual.getCplexStatus();
		_status=cplex_status;
		double sol = re ? cplex_dual.getObjValue() : 9999999999;
		_sol=sol;


		ksi_values = IloNumArray(env_dual,(int)ksi_dual.getSize());
		phi_values = IloNumArray(env_dual,(int)phi_dual.getSize());
		pi_values = IloNumArray(env_dual,(int)pi_dual.getSize());
		v_values = IloNumArray(env_dual,(int)v_dual.getSize());
	
		ksi.clear(); ksi.reserve( (int)ksi_dual.getSize() );
		phi.clear(); phi.reserve( (int)phi_dual.getSize() );
		pi.clear(); pi.reserve( (int)pi_dual.getSize() );
		v.clear(); v.reserve( (int)v_dual.getSize() );

		for(int j = 0; j < ksi_dual.getSize(); j++)
			ksi.push_back((double)cplex_dual.getValue(ksi_dual[j]));
		  
		for(int j = 0; j < phi_dual.getSize(); j++)
			phi.push_back((double)cplex_dual.getValue(phi_dual[j]));
		  
		for(int j = 0; j < v_dual.getSize(); j++)
			v.push_back((double)cplex_dual.getValue(v_dual[j]));
		  
		for(int j = 0; j < pi_dual.getSize(); j++)	
			pi.push_back((double)cplex_dual.getValue(pi_dual[j]));
	}
}

void FirstStageDuals::InitModelPrimal(int scenario)
{
    _scenario=scenario;
	double Q = _prob->GetDriver(0)->capacity;
  
    model = IloModel(env);
    
	flow = IloNumVarArray(env,_graph->GetArcCount(),0, IloInfinity, ILOFLOAT);
	w_minus = IloNumVarArray(env, _graph->GetNodeCount(), 0, IloInfinity, ILOFLOAT);
	w_plus = IloNumVarArray(env, _graph->GetNodeCount(), 0, IloInfinity, ILOFLOAT);
    
    flow_consts_vars = IloRangeArray(env);
    balance_consts_vars = IloRangeArray(env);
    wp_consts_vars = IloRangeArray(env);
    wm_consts_vars = IloRangeArray(env);
    
    for (int i = 1; i < _graph->GetNodeCount(); i++)
	{
		Node* n = _graph->GetNode(i);

		char name[40];
		sprintf(name, "wm%d", n->no);
		w_minus[n->no].setName(name);

		sprintf(name, "wp%d", n->no);
		w_plus[n->no].setName(name);
   
 		IloExpr expr(env);
		expr += w_minus[n->no];
		expr -= n->w_minus; 
		wm_consts_vars.add(expr <= 0);
    
		IloExpr expr1(env);
		expr1 += w_plus[n->no];
		expr1 -= n->w_plus; 
		wp_consts_vars.add(expr1 <= 0);
				
		model.add(wp_consts_vars);
		model.add(wm_consts_vars);
		expr1.end();
		expr.end();  
    }
 
	for(int i = 0; i < _graph->GetArcCount(); i++)
	{
		ExSbrpArcO* arc = _graph->GetArc(i);
		if(arc == NULL) continue;
		char name[40];
		sprintf(name, "f%d_%d", arc->from->no,arc->to->no);
		flow[i].setName(name);

		IloExpr expr(env);
		expr += flow[i];
		if(arc->walked_on == 1)
			flow_consts_vars.add(expr <= Q);
		else
			flow_consts_vars.add(expr == 0);
		model.add(flow_consts_vars);
		expr.end();
	}

	//Flow Conservation on AllStations
	for (int i = 1; i < _graph->GetNodeCount(); i++)
	{
		IloExpr expr(env);

		for (int j = 0; j < _graph->GetArcsInOfCount(i); j++)
			expr -= flow[_graph->GetArcInOf(i, j)->index];

		for (int j = 0; j < _graph->GetArcsOutOfCount(i); j++)
			expr += flow[_graph->GetArcOutOf(i, j)->index];

		Node* n = _graph->GetNode(i);  

		expr += w_plus[n->no];
		expr -= w_minus[n->no];
		balance_consts_vars.add(expr == n->demands[scenario]);
		model.add(balance_consts_vars);
		expr.end();
	}

	IloExpr recourse(env);
	for (int i = 1; i < _graph->GetNodeCount(); i++)
	{
		Node* n = _graph->GetNode(i);
		recourse += (w_minus[n->no] + w_plus[n->no]); 
	}	
	model.add(IloMinimize(env, recourse));
	recourse.end();
    
}

void FirstStageDuals::UpdateBalanceRHSPrimal(int scenario)
{
    
    _scenario=scenario;
	for (int i = 1; i < _graph->GetNodeCount(); i++)
    {
	  Node* n = _graph->GetNode(i);  
      double q = n->demands[scenario];
      balance_consts_vars[i-1].setBounds(q,q);     
    }
      
}

void FirstStageDuals::SolvePrimal() 
{

  cplex = IloCplex(model);
	
	if(Benders_feas)
	{
		cplex.setParam(IloCplex::Param::Preprocessing::Presolve, IloFalse);
		cplex.setParam(IloCplex::Param::Preprocessing::Reduce, 0);
		cplex.setParam(IloCplex::Param::RootAlgorithm, IloCplex::Dual);
		cplex.setParam(IloCplex::Param::NodeAlgorithm, IloCplex::Dual); 		

	}
    
	cplex_dual.setParam(IloCplex::Param::Threads,1);
	cplex.setOut(env.getNullStream());
	cplex.setWarning(env.getNullStream());
	  
	IloNumArray flow_duals_vars(env);
	IloNumArray balance_duals_vars(env);
	IloNumArray wp_duals_vars(env);
	IloNumArray wm_duals_vars(env); 

	bool re = cplex.solve();

	int cplex_status = (int)cplex.getCplexStatus();

	_status=cplex_status;

	double sol = re ? cplex.getObjValue() : 9999999999;
	//printf("primal sce:%d status:%d sol:%.1lf\n",_scenario,cplex_status,sol);
	_sol=sol;

   
	if(cplex_status == 1)
	{
		cplex.getDuals(flow_duals_vars, flow_consts_vars);
		cplex.getDuals(balance_duals_vars, balance_consts_vars);
		cplex.getDuals(wp_duals_vars, wp_consts_vars);
		cplex.getDuals(wm_duals_vars, wm_consts_vars);

		//printf("Primal sol ksi_size:%d phi_size:%d pi_size:%d v_size:%d\n",flow_duals_vars.getSize(),balance_duals_vars.getSize(),wm_duals_vars.getSize(),wp_duals_vars.getSize());
		  
		ksi.clear(); ksi.reserve( (int)ksi_dual.getSize() );
		phi.clear(); phi.reserve( (int)phi_dual.getSize() );
		pi.clear(); pi.reserve( (int)pi_dual.getSize() );
		v.clear(); v.reserve( (int)v_dual.getSize() );
		
		for(int j = 0; j < flow_duals_vars.getSize(); j++)
		ksi.push_back((double)flow_duals_vars[j]); 
		  
		for(int j = 0; j < balance_duals_vars.getSize(); j++)
		phi.push_back((double)balance_duals_vars[j]);
		  
		for(int j = 0; j < wp_duals_vars.getSize(); j++)
		v.push_back((double)wp_duals_vars[j]);
		  
		for(int j = 0; j < wm_duals_vars.getSize(); j++)
		pi.push_back((double)wm_duals_vars[j]);

		flow_duals_vars.end();		
		balance_duals_vars.end();		
		wp_duals_vars.end();		
		wm_duals_vars.end();

		cplex.end();	  
	}

}
