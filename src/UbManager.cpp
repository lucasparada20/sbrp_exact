
#include "UbManager.h"
#include <map>
#include <string>

void UbManager::Load(char * filename)
{
	FILE * f = fopen(filename, "r");

	if(f == NULL)
	{
		printf("UbManager::Load\nError at the opening of the file:%s\n", filename);
		return;
	}

	int nb;
	fscanf(f, "%d\n",&nb);


	for(int i=0;i<nb;i++)
	{
		char name[40];
		double delta, epsilon, bound;

		fscanf(f,"%s %lf %lf %lf\n", name, &epsilon, &delta, &bound);
		printf("%s %lf %lf %lf\n", name, epsilon, delta, bound);
		UbManager_Item item;
		item.name = std::string(name);
		item.delta = delta;
		item.epsilon = epsilon;
		item.ub = bound;

		_bounds.push_back(item);
	}
   //std::cout << "nb bounds loaded: " << _bounds.size() << std::endl; 
}

double UbManager::GetUpperBound(char * ins, double epsilon, double delta)
{
    int len = strlen(ins);
    char name[40]; 
    memset(name, 0, 40 * sizeof(char));
    
    int i, p1 = -1, p2 = -1;

    // Find last '/' and '.'
    for (i = 0; i < len; i++) {
        if (ins[i] == '/') p1 = i;
        if (ins[i] == '.') p2 = i;
    }

    // Ensure valid positions
    if (p1 != -1 && p2 != -1 && p1 < p2) {
        for (i = p1 + 1; i < p2; i++)
            name[i - p1 - 1] = ins[i];
        
        name[i - p1 - 1] = '\0'; // Proper null termination
    }

    std::string st(name);
    //std::cout << "Extracted filename: " << st << std::endl; 
	//std::cout << "ins: " << ins << std::endl;

	for(size_t i=0;i<_bounds.size();i++)
	{
		//std::cout << "St: " << st << std::endl;
		//std::cout << "Name from '_bounds[i]': " << _bounds[i].name << std::endl;
		if(_bounds[i].name == st && _bounds[i].delta == delta && _bounds[i].epsilon == epsilon)
				return _bounds[i].ub;
	}
		
	return 999999999.9;
}
