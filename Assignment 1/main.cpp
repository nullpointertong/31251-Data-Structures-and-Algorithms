#include <iostream>
#include <vector> 
#include "weighted_graph.hpp" 
#include "test.h"

int main()
{
	weighted_graph<int> g; 
	
	g.add_vertex(0);
	g.add_vertex(1);
	g.add_vertex(2);
	g.add_vertex(3);
	
	g.add_edge(0,1,1);
	g.add_edge(1,2,2);
	g.add_edge(2,3,3);
	g.add_edge(3,0,4);
	
	if(g.are_adjacent(3,0))
	{
		std::cout<<"Test"<<std::endl;
	}
	
	
	
	g.mst(); 
}