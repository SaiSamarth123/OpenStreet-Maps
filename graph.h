// graph.h
// Completed by: Sai Samarth
//
// Basic graph class using adjacency list representation.
// Currently there is no limit to the size of adjacency list.
//
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
//

#pragma once

#include <iostream>
#include <stdexcept>
#include <vector>
#include <set>
#include <unordered_map>
#include <unordered_set>
using namespace std;

template<typename VertexT, typename WeightT>
class graph {
    private:
    unordered_map <VertexT, unordered_map<VertexT, WeightT>> AdjList;
    vector<VertexT>  Vertices;
    
    public:
    // default constructor:
    // Constructs an empty graph object.
    // There is no size limit as we're using an adjacency list.
    graph() {}
    //
    // NumVertices
    //
    // Returns the # of vertices currently in the graph.
    //
    int NumVertices() const {
        // return the size of AdjList which is the number of vertices.
        return static_cast<int>(this->AdjList.size());
    }
    //
    // NumEdges
    //
    // Returns the # of edges currently in the graph.
    //
    int NumEdges() const {
      // Return the sum of sizes of the inner maps of every vertex.
      int count = 0;
      // Loop through the AdjList and get size of nested maps
      for (auto edge : AdjList) {
          count += edge.second.size();
      }
      return count;
    }
    //
    // addVertex
    //
    // Adds the vertex v to the graph if not present, and if so returns true.  If the vertex already
    // exists in the graph, then false is returned.
    //
    bool addVertex(VertexT v) {
        // is the vertex already in the graph?
        // If so, we do not insert again
        // We check if atleast one occurrence of the vertex is in AdjList.
        if (AdjList.count(v) > 0) {
            return false;
        }
        // if we get here, vertex does not exist so insert.
        // We create an empty map for the edges
        // and add the vertex to the AdjList.
        unordered_map <VertexT, WeightT> edges;
        AdjList[v] = edges;
        return true;
    }
    //
    // addEdge
    //
    // Adds the edge (from, to, weight) to the graph, and returns
    // true.  If the vertices do not exist, false is returned.
    //
    // NOTE: if the edge already exists, the existing edge weight
    // is overwritten with the new edge weight.
    //
    bool addEdge(VertexT from, VertexT to, WeightT weight) {
        // We check if the starting vertex is in the AdjList.
        int fromCount = AdjList.count(from);
        if (fromCount == 0) {  // return false if not found:
            return false;
        }
        // We check if the destination vertex is in the AdjList
        int toCount = AdjList.count(to);
        if (toCount == 0) {  // not found:
            return false;
        }
        // Get the nested map of from vertex.
        auto edges = AdjList.at(from);
        // Add the edge to the inner map.
        edges[to] = weight;
        // Add the inner map back to AdjList.
        AdjList[from] = edges;
        return true;
    }
    //
    // getWeight
    //
    // Returns the weight associated with a given edge.  If
    // the edge exists, the weight is returned via the reference
    // parameter and true is returned.  If the edge does not
    // exist, the weight parameter is unchanged and false is
    // returned.
    //
    bool getWeight(VertexT from, VertexT to, WeightT& weight) const {
        int fromCount = AdjList.count(from);
        int toCount = AdjList.count(to);
        if (fromCount == 0 || toCount == 0) {  // return false if not found
            return false;
        }
        // Get the nested edges if present
        auto edges = AdjList.at(from);
        // Check if from vertex has any edges and the given edge
        if (edges.empty() || edges.count(to) == 0) {
            return false;
        }
        // Get the weight of edge
        weight = edges.at(to);
        return true;
    }
    //
    // neighbors
    //
    // Returns a set containing the neighbors of v, i.e. all
    // vertices that can be reached from v along one edge.
    // Since a set is returned, the neighbors are returned in
    // sorted order; use foreach to iterate through the set.
    //
    set<VertexT> neighbors(VertexT v) const {
        set<VertexT>  S;
        // we need to check if v is in the AdjList
        int vertexCount = AdjList.count(v);
        if (vertexCount == 0) {  // not found:
            return S;
        }
        // Get the nested edges of vertex from AdjList
        auto edges = AdjList.at(v);
        // Return empty set if vertex has no edges.
        if (edges.empty()) {  // not found:
            return S;
        }
        // Loop through each edge and add the vertex of the edge to set
        for (auto neighbor : edges) {
            S.insert(neighbor.first);
        }
        return S;
    }
    //
    // getVertices
    //
    // Returns a vector containing all the vertices currently in
    // the graph.
    //
    vector<VertexT> getVertices() const {
        vector <VertexT> vertices;
        // Loop through the AdjList and add all the vertices to the vector
        for (auto &vertex : AdjList) {
            vertices.push_back(vertex.first);
        }
        return vertices;
    }
    //
    // dump
    //
    // Dumps the internal state of the graph for debugging purposes.
    //
    // Example:
    //    graph<string,int>  G(26);
    //    ...
    //    G.dump(cout);  // dump to console
    //
    void dump(ostream& output) const {
        output << "***************************************************" << endl;
        output << "********************* GRAPH ***********************" << endl;
        // Number of vertices and edges
        output << "**Num vertices: " << this->NumVertices() << endl;
        output << "**Num edges: " << this->NumEdges() << endl;
        output << endl;
        // Output all the vertices of the graph.
        output << "**Vertices:" << endl;
        vector<VertexT> vertices = this->getVertices();
        for (int i = 0; i < this->NumVertices(); ++i) {
            output << " " << i << ". " << vertices[i] << endl;
        }
        output << endl;
        // Output the edges of the graph.
        output << "**Edges:" << endl;
        for (auto v : vertices) {
            set<VertexT> neighbors = this->neighbors(v);
            for (auto n : neighbors) {
                WeightT weight;
                if (this->getWeight(v, n, weight)) {
                    output << "(" << v << "," << n << "," << weight << ") ";
                } else {
                    output << "(" << v << "," << n << "," << "???" << ") ";
                }
            }
        }
        output << "**************************************************" << endl;
    }
};
