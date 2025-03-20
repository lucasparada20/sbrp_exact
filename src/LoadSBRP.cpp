
#include "LoadSBRP.h"
#include "Parameters.h"
#include <iostream>     
#include <fstream>      
#include <cstdio>       
#include <cstdlib>      
#include <vector>       
#include <string>       
#include <sstream>      
#include <iterator>     
#include <limits>       
#include <algorithm>
#include <cfloat>    

void LoadSBRP::Load_pcg(Prob& pr, const char* filename) {
    FILE* ff = fopen(filename, "r");
    if (!ff) {
        fprintf(stderr, "Error in the input filename: %s\n", filename);
        exit(1);
    }

    int Q, Stations, Scenarios = 0;
    fscanf(ff, "%d", &Q);
    fscanf(ff, "%d", &Stations);
    printf("Q: %d N: %d\n", Q, Stations);

    std::vector<int> Station_cap(Stations), Station_occ(Stations);
    for (int& cap : Station_cap) fscanf(ff, "%d", &cap);
    for (int& occ : Station_occ) fscanf(ff, "%d", &occ);

    fscanf(ff, "%d", &Scenarios);
    printf("Scenarios: %d\n", Scenarios);

    std::vector<std::vector<int>> Demands(Stations);
    for (auto& demand : Demands) {
        demand.resize(Scenarios);
        for (int& d : demand) fscanf(ff, "%d", &d);
    }

    std::vector<std::vector<double>> Distances(Stations, std::vector<double>(Stations));
    for (auto& row : Distances)
        for (double& d : row) fscanf(ff, "%lf", &d);

    double c_min = std::numeric_limits<double>::max();
    int from = -1, to = -1;
    for (int i = 0; i < Stations; i++) {
        for (int j = 0; j < Stations; j++) {
            if (i != j && Distances[i][j] < c_min) {
                c_min = Distances[i][j];
                from = i;
                to = j;
            }
        }
    }

    if (c_min < 1.0) {
        fprintf(stderr, "c_min: %.1lf At least two stations are stacked on top. From:%d to:%d Stopping...\n", c_min,from,to);
        exit(1);
    }

    Parameters::SetCmin(c_min);

    for (int i = 0; i < Stations - 1; i++) {
		Node n;
		n.id = i;
		n.no = i+1;
		n.distID = i+1;
		n.type = NODE_TYPE_CUSTOMER;
		n.occupancy = Station_occ[i+1];
		n.stationcapacity = Station_cap[i+1];

        n.demands.reserve( Scenarios );
		for(int e=0; e<Scenarios;e++) 
			n.demands.push_back(Demands[i+1][e]);
        
		n.UpdateW();
        pr.AddNode(n);
    }

    for (int i = 0; i < Scenarios; i++) {
        pr.AddScenario(i);
    }

    for (int i = 0; i < Stations; i++) {
		Node dep1;
		dep1.id = Stations - 1 + i*2;
		dep1.no = 0;
		dep1.distID = 0;
		dep1.type = NODE_TYPE_START_DEPOT;
		dep1.occupancy = Station_occ[0];
		dep1.stationcapacity = Station_cap[0];

		Node dep2(dep1);
		dep2.id = Stations + i*2;
		dep2.type = NODE_TYPE_END_DEPOT;
		dep2.occupancy = Station_occ[0];
		dep2.stationcapacity = Station_cap[0];

		Driver d;
		d.capacity = Q;
		d.StartNodeID = dep1.id;
		d.EndNodeID = dep2.id;
		d.id = i;

		pr.AddNode(dep1);
		pr.AddNode(dep2);
		pr.AddDriver(d);
    }

    for (int i = 0; i < Stations - 1; i++) {
        pr.AddCustomer(pr.GetNode(i));
    }

    int dim = Stations;
    auto d = new double*[dim];
    for (int i = 0; i < dim; i++) {
        Node* n1 = pr.GetNode(i);
        d[n1->distID] = new double[dim];
        for (int j = 0; j < dim; j++) {
            Node* n2 = pr.GetNode(j);
            d[n1->distID][n2->distID] = (i == j) ? 0.0 : Distances[n1->distID][n2->distID];
        }
    }

    pr.SetMatrices(d, dim);
    fclose(ff);
}

void LoadSBRP::Load_dins(Prob& pr, const char* filename) {
    pr = Prob();

    std::ifstream ff(filename);
    if (!ff) {
        fprintf(stderr, "Error in the input filename: %s\n", filename);
        exit(1);
    }

    int Stations;
    ff >> Stations;

    std::vector<int> StationCapacity(Stations), Occupancy(Stations);
    for (int& cap : StationCapacity) ff >> cap;
    for (int& occ : Occupancy) ff >> occ;

    int Scenarios;
    ff >> Scenarios;

    std::vector<double> Probability(Scenarios);
    for (double& prob : Probability) ff >> prob;

	std::vector< std::vector<int> > Demands(Scenarios, std::vector<int>(Stations, 0));
	for (int i = 0; i < Scenarios; i++)
		for (int j = 0; j < Stations; j++)
			ff >> Demands[i][j];


    int FleetSize, VehicleCap;
    ff >> FleetSize >> VehicleCap;

    std::vector<std::vector<double>> Distance(Stations, std::vector<double>(Stations));
    for (auto& row : Distance)
        for (double& d : row) ff >> d;

    double c_min = DBL_MAX;
    for (const auto& row : Distance)
        for (double value : row)
            if (value > 0) c_min = std::min(c_min, value);

    Parameters::SetCmin(c_min);

    for (int i = 0; i < Stations - 1; i++) {
		Node n;
		n.id = i;
		n.no = i+1;
		n.distID = i+1;
		n.type = NODE_TYPE_CUSTOMER;
		n.occupancy = Occupancy[i+1];
		n.stationcapacity = StationCapacity[i+1];

		n.demands.reserve( Scenarios );
        for(int e=0; e<Scenarios;e++)
			n.demands.push_back(Demands[e][i+1]);
		
        n.UpdateW();
        pr.AddNode(n);
		
    }

    for (int i = 0; i < Scenarios; i++) {
        pr.AddScenario(i);
    }

    for (int i = 0; i < FleetSize; i++) {
		Node dep1;
		dep1.id = Stations - 1 + i*2;
		dep1.no = 0;
		dep1.distID = 0;
		dep1.type = NODE_TYPE_START_DEPOT;
		dep1.occupancy = Occupancy[0];
		dep1.stationcapacity = StationCapacity[0];

		Node dep2(dep1);
		dep2.id = Stations + i*2;
		dep2.type = NODE_TYPE_END_DEPOT;
		dep2.occupancy = Occupancy[0];
		dep2.stationcapacity = StationCapacity[0];

		Driver d;
		d.capacity = VehicleCap;
		d.StartNodeID = dep1.id;
		d.EndNodeID = dep2.id;
		d.id = i;

		pr.AddNode(dep1);
		pr.AddNode(dep2);
		pr.AddDriver(d);
    }
		
    for (int i = 0; i < Stations - 1; i++) {
        pr.AddCustomer(pr.GetNode(i));
    }

    int dim = Stations;
    auto d = new double*[dim];
    for (int i = 0; i < dim; i++) {
        Node* n1 = pr.GetNode(i);
        d[n1->distID] = new double[dim];
        for (int j = 0; j < dim; j++) {
            Node* n2 = pr.GetNode(j);
            d[n1->distID][n2->distID] = (i == j) ? 0.0 : Distance[n1->distID][n2->distID];
        }
    }

    pr.SetMatrices(d, dim);
}

