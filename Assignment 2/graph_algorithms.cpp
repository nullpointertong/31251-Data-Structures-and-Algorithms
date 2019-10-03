#ifndef GRAPH_ALGS
#define GRAPH_ALGS

#include <map>
#include <vector>
#include <queue>
#include <stack>
#include <list>
#include <deque>
#include <set>
#include <map>
#include <unordered_set>
#include <unordered_map>
#include <iostream>
#include <utility>
#include <algorithm>
#include <limits>
#include "weighted_graph.hpp"
#include "easy_weighted_graph_algorithms.cpp"

//Returns true if the graph is connected, false otherwise.

template <typename vertex>
bool is_connected(const weighted_graph<vertex>& g) {
    if (g.num_vertices() == 0 || g.num_vertices() == depth_first(g, g.get_vertices()[0]).size()) {
        return true; //Launches a Depth First Traversal and compares the values of the dft 
    } //versus a vector containing all the vertices
    return false; //Returns true if all values of vertices are found in the dft otherwise false     
}

//Returns a vector of weighted graphs, where each weighted graph is a connected
//component of the input graph.

template <typename vertex>
std::vector<weighted_graph<vertex>> connected_components(const weighted_graph<vertex>& g) {
    std::vector<weighted_graph < vertex>> connected_comp; //Sets up a vector containing all vertices within g and a graph vector to

    if (g.num_vertices() == 0) { //Checks for empty graph
        return connected_comp;
    } else {
        bool visited[g.num_vertices()]; //NOTE: Arrays are used over Maps here intentionally for speed reasons                      
        for (auto i = g.begin(); i != g.end(); i++) { //also note the indexing function which converts Vertices into position int values
            visited[g.index(*i)] = false;
        }

        for (auto v : g.get_vertices()) {
            while (!visited[g.index(v)]) { //Finds unvisted Vertex
                weighted_graph<vertex> sub_graph_g; //Creates a new graph to store the sub graph
                std::vector<vertex> dftDis = depth_first(g, v); //Stores the vertices of the DFT
                for (auto u : dftDis) {
                    sub_graph_g.add_vertex(u); //Adds found vertices into sub graph and marks added vertices as visited
                    visited[g.index(u)] = true;
                }
                for (auto u : dftDis) {
                    for (auto w : dftDis) {
                        if (g.are_adjacent(u, w)) { //Adds any previous edges from g into the sub graph
                            sub_graph_g.add_edge(u, w, g.get_edge_weight(u, w));
                        }
                    }
                }
                connected_comp.push_back(sub_graph_g); //pushes sub graph back into the graph vector
            }
        }
    }
    return connected_comp; //returns vector of sub graphs 
}

//Returns a map of the vertices of the weighted graph g and their distances from
//the given starting vertex v.

template <typename vertex>
std::map<vertex, int> dijkstras(const weighted_graph<vertex>& g, const vertex& v) {

    std::map<vertex, int> vertex_distance_map; //Intialisation of a vertex/distance map, visted array, current vertex and a vector of all vertices in g

    if (g.num_vertices() == 0) { //Checks for empty graph
        return vertex_distance_map;
    } else {
        bool visited[g.num_vertices()]; //NOTE: Arrays are used over Maps here intentionally for speed reasons
        vertex current; //also note the indexing function which converts vertices into position int values

        for (auto i = g.begin(); i != g.end(); i++) {
            vertex_distance_map[*i] = __INT_MAX__; //Intialisation of all non-v values to Int Max and bool array
            visited[g.index(*i)] = false;
        }
        vertex_distance_map[v] = 0; //Sets tentative distance of all vertices to "infinity"(or int_max here) besides the current
        current = v; //which is set to 0, v is set as current 
        while (!visited[g.index(current)]) {
            for (auto w : g.get_vertices()) { //Checks if current is visited and searchs for edges incident to current
                if (g.are_adjacent(current, w)) {
                    if (vertex_distance_map[w] > vertex_distance_map[current] + g.get_edge_weight(current, w)) {
                        vertex_distance_map[w] = vertex_distance_map[current] + g.get_edge_weight(current, w);
                    }//Compares tentative distance of incident vertex to the tentative distance of current + edge distance
                }
            }
            visited[g.index(current)] = true; //Sets current to visited
            int min = __INT_MAX__;
            for (auto l : g.get_vertices()) {
                if (vertex_distance_map[l] < min && !visited[g.index(l)]) {
                    min = vertex_distance_map[l];
                    current = l; //Choses next current vertex based on the lowest tentative distance
                }
            }
        }
    }
    return vertex_distance_map; //Returns map of vertices and shortest distance
}

//Recursive helper method used to update variables(via Parameters) needed in order 
//to find Articulation points 
template <typename vertex>
void articulation(vertex v, int d, const weighted_graph<vertex>& g, std::vector<vertex>& articu_points,
        std::map<vertex, bool>& visited, std::map<vertex, int>& depth, std::map<vertex, int>& low,
        std::map<vertex, vertex>& parent) { 

    int child_count = 0; //Sets up number of children for 1 parent vertex             
    bool is_artic = false; //Returns true if vertex is articulation point
    visited[v] = true; //Sets current vertex as visited
    depth[v] = d; //Records how deep a vertex is within a graph                             
    low[v] = d; //Records the lowest depth value within graph

    for (auto u : g.get_vertices()) {
        if (g.are_adjacent(u, v)) {
            if (!visited[u]) //Iterates over every neighbour of v where not visited
            {
                parent[u] = v; //Sets u as child to v 
                articulation(u, d + 1, g, articu_points, visited, depth, low, parent); //Recursively increases depth and moves to neighbour of v
                child_count++; //Increase the child count 
                low[v] = std::min(low[v], low[u]); //Sets the value of low to the lowest depth between v and neighbour of v(u)
                if (low[u] >= depth[v]) //if the depth of selected vertex v is the same or smaller than the lowest depth of neighbour u 
                { //v is an aritculate point 
                    is_artic = true;
                }
            } else if (u != parent[v])//if the neighbour of v(u) is not a parent of v depth of u is compared to the lowest depth v                                   
            {
                low[v] = std::min(low[v], depth[u]);
            }
        }
    }
    if (parent.count(v) > 0 && is_artic) // if v has parent vertices and is marked as an articulation point the vertex will be pushed back
    {
        articu_points.push_back(v);
    } else if (parent.count(v) == 0 && child_count > 1) // if v has no parent vertices and it has more than 1 child vertex it is an articulation
    { //Point is pushed back
        articu_points.push_back(v);
    }
}

//Returns a vector containing all the articulation points of the
//input weighted graph g.

template <typename vertex>
std::vector<vertex> articulation_points(const weighted_graph<vertex>& g) {
    std::vector<vertex> articu_points; //Stores data gathered from recursive helper method 
    std::map<vertex, bool> visited;
    std::map<vertex, int> depth;
    std::map<vertex, int> low;
    std::map<vertex, vertex> parent;

    if (g.num_vertices() > 0) { //Inputs variables into helper function to retrieve Articulation points
        articulation(*g.begin(), 0, g, articu_points, visited, depth, low, parent); //Calls helper method
    }
    return articu_points; //Returns Articulation points 
}

#endif
