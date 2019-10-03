#ifndef WEIGHTED_GRAPH_H
#define WEIGHTED_GRAPH_H 
//A large selection of data structures from the standard
//library. You need not feel compelled to use them all,
//but as you can't add any, they're all here just in case.
#include <iostream>
#include <vector>
#include <queue>
#include <stack>
#include <iostream>
#include <unordered_set>
#include <unordered_map>
#include <set>
#include <array>
#include <list>
#include <forward_list>
#include <deque>
#include <map>

template <typename vertex>
class weighted_graph {
private:
    std::vector<std::vector<int> > adj_matrix; //Matrix that stores all edges found
    std::vector<vertex> vertices; // Vector containing all vertices


    int num_of_Edges = 0; //Counts all edges, vertices and keeps track of total weight of all edges
    int num_of_Vertices = 0;
    int totalWeight = 0;


    //You will need to add some data members here
    //to actually represent the graph internally,
    //and keep track of whatever you need to.

    //The graph_iterator class provides an iterator
    //over the vertices of the graph.
    //This is one of the harder parts, so if you're
    //not too comfortable with C++ leave this for last.
    //If you are, there are many ways of doing this,
    //as long as it passes the tests, it's okay.

    class graph_iterator {
    private:
        //You may need data members here. 
        weighted_graph<vertex> owner; // Class construct for Weighted graph
        int pos = 0; //Position of Iterator


    public: //Iterates over the vertices of a graph
        graph_iterator(const weighted_graph &); //All comments in function
        graph_iterator(const weighted_graph &, size_t);
        ~graph_iterator();
        graph_iterator operator=(const graph_iterator&);
        bool operator==(const graph_iterator&) const;
        bool operator!=(const graph_iterator&) const;
        graph_iterator operator++();
        graph_iterator operator++(int);
        const vertex operator*(); //Built look like a pointer and returns what it is pointing at in this case a vertex
        const vertex* operator->();
    };

    //The neighbour_iterator class provides an iterator
    //over the neighbours of a given vertex. This is
    //probably harder (conceptually) than the graph_iterator.
    //Unless you know how iterators work.

    class neighbour_iterator {
    private: //Iterators over the neighbours of a chosen vertice //All comments in function

        //You may need data members here.
        weighted_graph<vertex> owner;
        int pos = 0;
        vertex target_vertex;

    public:
        neighbour_iterator(const neighbour_iterator&);
        neighbour_iterator(const weighted_graph &, const vertex&);
        neighbour_iterator(const weighted_graph &, const vertex&, size_t);
        ~neighbour_iterator();
        neighbour_iterator operator=(const neighbour_iterator& it);
        bool operator==(const neighbour_iterator&) const;
        bool operator!=(const neighbour_iterator&) const;
        neighbour_iterator operator++();
        neighbour_iterator operator++(int);
        const std::pair<vertex, int> operator*(); //Returns pair as edge 
        const std::pair<vertex, int>* operator->();
    };

public:

    weighted_graph(); //A constructor for weighted_graph. It should start empty.
    ~weighted_graph(); //A destructor. Depending on how you do things, this may
    //not be necessary.

    bool are_adjacent(const vertex&, const vertex&) const; //Returns true if the two vertices are
    //adjacent, false otherwise.
    bool has_vertex(const vertex&) const; //Returns true if the passed in vertex is 
    //a vertex of the graph, false otherwise.
    void add_vertex(const vertex&); //Adds the passed in vertex to the graph (with no edges).
    void add_edge(const vertex&, const vertex&, const int&); //Adds an edge between the two vertices
    //with the given weight (as an int).
    void remove_vertex(const vertex&); //Removes the given vertex. Should also clear any incident edges. 
    int index(const vertex&) const; // Gets the position of the vectex within the Vectices vector
    void remove_edge(const vertex&, const vertex&); //Removes the edge between the two vertices, if it exists.
    void set_edge_weight(const vertex&, const vertex&, const int&); //Changes the edge weight between the two
    //vertices to the new weight (the int).
    int get_edge_weight(const vertex&, const vertex&) const; //Returns the weight on the edge between the two vertices.
    int degree(const vertex&) const; //Returns the degree of the vertex.
    int weighted_degree(const vertex&); //Returns the sum of the weights on all the edges incident to the vertex.
    int num_vertices() const; //Returns the total number of vertices in the graph.
    int num_edges() const; //Returns the total number of edges in the graph (just the count, not the weight).
    int total_weight(); //Returns the sum of all the edge weights in the graph.

    std::vector<vertex> get_vertices(); //Returns a vector containing all the vertices.
    std::vector<vertex> get_neighbours(const vertex&); //Returns a vector containing the neighbours of the given vertex.

    graph_iterator begin(); //Returns a graph_iterator pointing to the start of the vertex set.
    graph_iterator end(); //Returns a graph_iterator pointing to one-past-the-end of the vertex set.

    neighbour_iterator neighbours_begin(const vertex&); //Returns a neighbour_iterator pointing to the start
    //of the neighbour set for the given vertex.
    neighbour_iterator neighbours_end(const vertex&); //Returns a neighbour_iterator pointing to one-past-the-end
    //of the neighbour set for the given vertex.

    std::vector<vertex> depth_first(const vertex&); //Returns the vertices of the graph in the order they
    //are visited in by a depth-first traversal starting at
    //the given vertex.
    std::vector<vertex> breadth_first(const vertex&); //Returns the vertices of the graph in the order they
    //are visisted in by a breadth-first traversal starting
    //at the given vertex.

    weighted_graph<vertex> mst(); //Returns a minimum spanning tree of the graph.

};

//Define all your methods down here (or move them up into the header, but be careful you don't double up).
//Although these are just the same names copied from above, you may find a few more clues in the full
//method headers. Note also that C++ is sensitive to the order you declare and define things in - you
//have to have it available before you use it.

template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g) {
    owner = g; //Intialises the Class Object owner as well as set's it to position 0 	
    pos = 0;
}

template <typename vertex> weighted_graph<vertex>::graph_iterator::graph_iterator(const weighted_graph & g, size_t start_pos) {
    owner = g; //Sets position to start 
    pos = start_pos;
}

template <typename vertex> weighted_graph<vertex>::graph_iterator::~graph_iterator() {
}

template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator=(const graph_iterator& it) { //Attaches iterators to the class object and the position
    owner = it.owner; //Attaches iterators to variables 
    pos = it.pos;
    return *this;
}

template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator==(const graph_iterator& it) const {
    return pos == it.pos; //Checks if the iterator is the same as the position
}

template <typename vertex> bool weighted_graph<vertex>::graph_iterator::operator!=(const graph_iterator& it) const {
    return pos != it.pos; //Checks if iterator is different from the position 
}

template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++() {
    ++pos;
    return *this; //Moves the position
}

template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::graph_iterator::operator++(int) {
    pos++;
    return *this; // Moves the position
}

template <typename vertex> const vertex weighted_graph<vertex>::graph_iterator::operator*() {
    return owner.get_vertices()[pos]; //Returns iterated vertices
}

template <typename vertex> const vertex* weighted_graph<vertex>::graph_iterator::operator->() {
    return (owner.get_vertices()[pos]); //Returns iterated vertices 
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u) {
    owner = g;
    pos = 0; //Sets parameters

    if (g.has_vertex(u)) {
        target_vertex = u;
    }
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const weighted_graph & g, const vertex& u, size_t start_pos) {
    owner = g;
    pos = start_pos; //Sets starting position 
    if (g.has_vertex(u)) {
        target_vertex = u;
    }
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::neighbour_iterator(const neighbour_iterator& it) {
    pos = it.pos;
    owner = it.owner; //Ataches iterators to all parameters
    target_vertex = it.target_vertex;
}

template <typename vertex> weighted_graph<vertex>::neighbour_iterator::~neighbour_iterator() {
}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator=(const neighbour_iterator& it) {
    pos = it.pos;
    owner = it.owner; //Copy constructor 
    target_vertex = it.target_vertex;
    return *this;
}

template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator==(const neighbour_iterator& it) const {
    return pos == it.pos; //Checks position of iterator
}

template <typename vertex> bool weighted_graph<vertex>::neighbour_iterator::operator!=(const neighbour_iterator& it) const {
    return pos != it.pos; //Checks position of iterator
}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++() {
    ++pos; //Increments iterator along the list
    return *this;
}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbour_iterator::operator++(int) {
    ++pos; //Increments iterator along the list
    return *this;
}

template <typename vertex> const std::pair<vertex, int> weighted_graph<vertex>::neighbour_iterator::operator*() {
    std::vector<vertex> neighboursIt = owner.get_neighbours(target_vertex); //Returns iteraterated neighbour 
    neighboursIt[pos];
    return neighboursIt;
}

template <typename vertex> const std::pair<vertex, int>* weighted_graph<vertex>::neighbour_iterator::operator->() {
    //Returns iteraterated neighbour 
    auto l = std::pair<vertex, int>(owner.get_neighbours(target_vertex)[pos], owner.get_edge_weight(target_vertex, owner.get_neighbours(target_vertex)[pos]));
    return &l;
}

template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::begin() {
    return graph_iterator(*this); //Defines starting point for the iterator
}

template <typename vertex> typename weighted_graph<vertex>::graph_iterator weighted_graph<vertex>::end() {
    return graph_iterator(*this, this->num_of_Vertices); //Defines the ending point for the iterator
}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_begin(const vertex& u) {
    return neighbour_iterator(*this, u); //Defines starting pointer for the iterator
}

template <typename vertex> typename weighted_graph<vertex>::neighbour_iterator weighted_graph<vertex>::neighbours_end(const vertex& u) {
    return neighbour_iterator(*this, u, this->degree(u)); //Defines the ending pointer for the iterator
}

template <typename vertex> weighted_graph<vertex>::weighted_graph() {
} //Constructor 

template <typename vertex> weighted_graph<vertex>::~weighted_graph() {
} //Destructor 

template <typename vertex> bool weighted_graph<vertex>::has_vertex(const vertex& u) const {
    if (index(u) == -1) //Checks if Vertex is in the graph or not
    { //Note: If Vertex is not found in the index function the index function will equal =-1
        return false; //Note2: This function is used almost throughout the graph to validate Vertices to prevent crashes
    }//from invalid vertices 
    else
        return true;
}

template <typename vertex> bool weighted_graph<vertex>::are_adjacent(const vertex& u, const vertex& v) const {
    if (has_vertex(u) && has_vertex(v)) { //Checks if vertices are valid
        if (adj_matrix[index(u)][index(v)] != 0) {
            return true;
        } else //Checks if an edge is between a pair of vertices by seeing if edgeweight is a non zero value
            return false;
    } else
        return false;
}

template <typename vertex> void weighted_graph<vertex>::add_vertex(const vertex& v) {
    //Adds Vertex to graph 
    ++num_of_Vertices;
    //Increments Vertice Counter
    vertices.push_back(v); //Pushs the vertex into the vector

    adj_matrix.resize(num_of_Vertices); //Resizes the adj matrix to make space for vertex
    for (unsigned i = 0; i < num_of_Vertices; i++) {
        adj_matrix[i].resize(num_of_Vertices, 0);
    }
}

template <typename vertex> void weighted_graph<vertex>::add_edge(const vertex& u, const vertex& v, const int& weight) {
    //Add edge between vertices
    if (u != v && has_vertex(u) && has_vertex(v)) { //Checks for self loops and for valid vertexs 
        adj_matrix[index(u)][index(v)] = weight; //Sets an edge between 2 vertices a weight.
        adj_matrix[index(v)][index(u)] = weight; //As graph is undirected weighted is added for all combinations of vertices

        ++num_of_Edges; //Increments the edge counter; 
        totalWeight = totalWeight + weight; // Adds weight to total weight
    }
}

template <typename vertex> void weighted_graph<vertex>::remove_vertex(const vertex& u) {
    if (has_vertex(u)) //Checks for vertex 
    {
        int i = index(u); //Finds position of vertex
        vertices.erase(vertices.begin() + i); //Removes vertex from vector/ reduces the num_of_Vertices
        --num_of_Vertices;

        adj_matrix.resize(num_of_Vertices); //Resizes the adj matrix 
        for (unsigned i = 0; i < num_of_Vertices; i++) {
            adj_matrix[i].resize(num_of_Vertices, 0);
        }
    }
}

template <typename vertex> void weighted_graph<vertex>::remove_edge(const vertex& u, const vertex& v) {
    //Removes edge
    if (has_vertex(u) && has_vertex(v) && adj_matrix[u][v] != 0) { //Checks if vertex and edge exists 
        if ((u >= 0) && (v >= 0) && (u != v)) {
            adj_matrix[index(u)][index(v)] = 0;
            adj_matrix[index(v)][index(u)] = 0; //Sets edge weight to 0 and reduces the edge counter by 1
            --num_of_Edges;
        }
    }
}

template <typename vertex> void weighted_graph<vertex>::set_edge_weight(const vertex& u, const vertex& v, const int& weight) {
    if (has_vertex(u) && has_vertex(v) && adj_matrix[u][v] != 0) { //Checks if vertex and edge exists 
        adj_matrix[index(u)][index(v)] = weight; //Modifies the weight of an existing edge
        adj_matrix[index(v)][index(u)] = weight;
    }
}

template <typename vertex> int weighted_graph<vertex>::get_edge_weight(const vertex& u, const vertex& v) const {
    if (has_vertex(u) && has_vertex(v)) {
        if (adj_matrix[index(u)][index(v)] != 0)
            return adj_matrix[index(u)][index(v)]; //Returns the weight of an existing edge(if any) between 2 vertices
    }
}

template <typename vertex> int weighted_graph<vertex>::degree(const vertex& u) const {
    int degree = 0;
    for (int i = 0; i < vertices.size(); i++) //Finds the number of neighbours a vertex has and returns it
    {
        if (has_vertex(u)) {
            if (adj_matrix[index(u)][i] != 0) { //Finds edges and increments the degree variable and returns it
                ++degree;
            }
        }
    }
    return degree;
}

template <typename vertex> int weighted_graph<vertex>::weighted_degree(const vertex& u) {
    int weightedDegree = 0;
    for (int i = 0; i < vertices.size(); i++) //Adds up the weight of the neighbour edges and returns it
    {
        if (has_vertex(u)) {
            if (adj_matrix[index(u)][i] != 0) { //Finds edges and adds edgesweight to the weightedDegree variable and returns it
                weightedDegree = weightedDegree + adj_matrix[index(u)][i];
            }
        }
    }
    return weightedDegree;
}

template <typename vertex> int weighted_graph<vertex>::num_vertices() const {
    return num_of_Vertices; //Returns number of vertices
}

template <typename vertex> int weighted_graph<vertex>::num_edges() const {
    return num_of_Edges; //Returns number of edges 
}

template <typename vertex> int weighted_graph<vertex>::total_weight() {
    return totalWeight; //Returns weight of all edges added together
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::get_vertices() {
    return vertices; //Returns the vector containing the vertices
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::get_neighbours(const vertex& u) {
    std::vector<vertex> neighbours;
    for (int i = 0; i < vertices.size(); i++) //Finds all the neighbour vertices by finding the position which creates a non-zero edge-weight
    {
        if (has_vertex(u)) {
            if (adj_matrix[index(u)][i] != 0) {
                neighbours.push_back(vertices[i]); //Pushes vertex into neighbour when an edge is found
            }
        }
    }
    return neighbours;
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::depth_first(const vertex& start_vertex) {
    bool visited[num_of_Vertices]; //Traveres the graph depth wise using a stack verses a queue
    for (unsigned i = 0; i < num_of_Vertices; i++) {
        visited[i] = false; //intailizes a boolean array allowing for vertices to be marked visted or not visted
    }

    std::stack<vertex> unprocessed; //Intialises a stack for unprocessed vertices, a vector for the ordered vertices as well as pushing a starting vertex
    unprocessed.push(start_vertex);
    std::vector<vertex> ordered;
   
    while (!unprocessed.empty()) { //When starting vertex is pushed into the vector the while statement is triggered
        vertex n = unprocessed.top(); //n takes top value of the unprocessed stack and afterwards the top value is removed
        unprocessed.pop();
        if (!visited[index(n)] && has_vertex(n)) { //Checks if vertice has been visted of not and marks it is visted if it hasn't been visted yet
            visited[index(n)] = true; //Marks the vertex position as visted
            ordered.push_back(index(n)); //Marks the position of vertex as visted, pushes vertex into the ordered vector
            for (unsigned i = num_of_Vertices; i != 0; i--) { //Finds the next neighbour and pushes it into unprocessed.
                if (adj_matrix[index(n)][i - 1]) {
                    unprocessed.push(vertices[i - 1]); //Note that due to the use of a stack it is traversed depth wise 
                }
            }
        }
    }
    return ordered;
}

template <typename vertex> std::vector<vertex> weighted_graph<vertex>::breadth_first(const vertex& start_vertex) {
    bool visited[num_of_Vertices]; //Traveres the graph breadth wise
    for (unsigned i = 0; i < num_of_Vertices; i++) {
        visited[i] = false; //intailizes a boolean array allowing for vertices to be marked visted or not visted
    }

    std::queue<vertex> unprocessed; //Intialises a queue for unprocessed vertices, a vector for the ordered vertices as well as pushing a starting vertex
    unprocessed.push(start_vertex);
    std::vector<vertex> ordered;

    while (!unprocessed.empty()) { //When starting vertex is pushed into the vector the while statement is triggered
        vertex n = unprocessed.front(); //n takes the front value and afterwards the front value is removed
        unprocessed.pop();
        if (!visited[index(n)] && has_vertex(n)) { //Checks if vertices has been visted or not and marks as visted if it hasn't been visted yet
            visited[index(n)] = true; //Marks the vertex position as visted
            ordered.push_back(n); //Pushes vertex into ordered vector  
            for (unsigned i = 0; i < num_of_Vertices; i++) { //Finds the next neighbour and pushes it into unprocessed.
                if (adj_matrix[index(n)][i] != 0) {
                    unprocessed.push(vertices[i]); //Checks for neighbour to push into unprocessed
                }
            }
        }
    }
    return ordered; //Note due to the use of a queue it is traversed breadth wise 
}

template <typename vertex> weighted_graph<vertex> weighted_graph<vertex>::mst() {
    weighted_graph<vertex> mst; //Creates a spanning tree with the minimal edge weight required to reach all vertices
    mst.add_vertex(vertices[0]); //Pushes vertex in first position
    while (mst.num_vertices() < num_vertices()) {
        int weight_min = 9999; //Weight Min set to large value so that it can be compared to the edge weights
        vertex min_v; //Intialises variables min_v and min_u to add to mst 
        vertex min_u;

        for (vertex v : mst.get_vertices()) { //Using a for each loop I can loop through every vertice within the mst and it's weighted graph neighour
            for (vertex u : get_neighbours(v)) {
                if (!mst.has_vertex(u) && adj_matrix[index(u)][index(v)] < weight_min) { //Stops cycles from happening by checking if u is in weighted graph
                    min_v = v; //Checks for the smallest weight in weighted graph with u and v
                    min_u = u; //Using these values I take v and u as min_v  and min_u and weight as weight_min
                    weight_min = adj_matrix[index(u)][index(v)];
                }
            }
        }
        mst.add_vertex(min_u);
        mst.add_edge(min_u, min_v, weight_min); //Adds a vertex (min_u) and an edge to the mst and then returns to looping due to the for each loop 
    }
    return mst;
}

template <typename vertex> int weighted_graph<vertex>::index(const vertex& v) const {
    for (int i = 0; i < vertices.size(); i++) //Find Position of Vertex within the vertices vector 
    {
        if (vertices[i] == v) { //Returns iterative element when vertice is found
            return i;
        }
    }
    return -1; //Returns -1 if not found allowing for this function to be reused in has_vertex
}

#endif