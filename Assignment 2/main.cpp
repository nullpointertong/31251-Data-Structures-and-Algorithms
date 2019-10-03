#include <iostream>
#include "weighted_graph.hpp"
#include "graph_algorithms.cpp"

int main(){
 weighted_graph<int> g; 
	g.add_vertex(0);
	g.add_vertex(1);
	g.add_vertex(2);
	g.add_vertex(3);
	g.add_vertex(4);
	
	g.add_edge(0,1,1);
	g.add_edge(1,2,2);
	g.add_edge(2,3,3);
	g.add_edge(3,4,4);
	g.add_edge(4,0,4);
	
	
	for(auto l : dijkstras(g, 0))
	{
	  std::cout<< l.first << ", " << l.second <<std::endl;	
	}
	
	
	auto result = articulation_points(g);
	
	for (auto i : result) std::cout << i << " ";
	std::cout << std::endl;
			
}
	
	