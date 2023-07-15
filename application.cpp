// application.cpp
// Completed by: Sai Samarth
// This is the application part that writes a console-based C++ program to input a campus map (e.g. UIC’s East campus) and navigate between buildings via footways
// Gives the shortest path between two buildings if it exists.
// University of Illinois at Chicago
// CS 251: Spring 2021
// Project #7 - Openstreet Maps
//
// References:
// TinyXML: https://github.com/leethomason/tinyxml2
// OpenStreetMap: https://www.openstreetmap.org
// OpenStreetMap docs:
//   https://wiki.openstreetmap.org/wiki/Main_Page
//   https://wiki.openstreetmap.org/wiki/Map_Features
//   https://wiki.openstreetmap.org/wiki/Node
//   https://wiki.openstreetmap.org/wiki/Way
//   https://wiki.openstreetmap.org/wiki/Relation
//

#include <iostream>
#include <iomanip> /*setprecision*/
#include <string>
#include <queue>
#include <map>
#include <stack>
#include <limits>
#include <cstdlib>
#include <cstring>
#include <cassert>
#include "graph.h"
#include "tinyxml2.h"
#include "dist.h"
#include "osm.h"
using namespace std;
using namespace tinyxml2;
// Priotize class used for the priority queue
class prioritize {
public:
    // Function that is responsible for sorting the priority queue
    bool operator()(pair<long long, double> const &p1,
    pair<long long, double> const &p2) const {
        if (p1.second == p2.second) {
            return p1.first > p2.first;
        }
        return p1.second > p2.second;
    }
};
// findShorterDist
// Dijkstra HELPER FUNCTION
// Helper function finds the shortest distance between the neighbors of a point
// Function is called from the dijkstra's function.
void findShorterDist(long long curV, double weight, map<long long,
double> &distances, priority_queue<pair<long long, double>,
vector<pair<long long, double> >, prioritize> &pq,
map<long long, long long> &pred, graph<long long, double> &G) {
    // Getting the neighbors of curV
    set<long long> neighbors = G.neighbors(curV);
    // Loop through the neighbors set
    for (auto adjV : neighbors) {
        double edgeWeight = 0;
        G.getWeight(curV, adjV, edgeWeight);
        // calculate the alternate path distance
        double altPathDist = weight + edgeWeight;
        // Update the distance if alternate distance less than current distance
        if (altPathDist < distances.at(adjV)) {
            distances[adjV] = altPathDist;
            pq.emplace(adjV, altPathDist);
            pred[adjV] = curV;
        }
    }
}
//
// Dijkstra:
//
// Performs Dijkstra's shortest weighted path algorithm from
// the given start vertex.  Returns a map of (string,int)
// pairs where the string is a vertex V and the int is the
// distance from the start vertex to V; if no such path exists,
// the distance is INF.
// Also returns the predecessor array through the parameter.
void Dijkstra(graph<long long, double> &G, long long startV,
map<long long, double> &distances, map<long long, long long> &pred) {
    double INF = numeric_limits<double>::max();  // Declaring INF
    set<long long> visitedSet;  // Set to check if visited
    vector<long long> vertices = G.getVertices();  // Vector of all vertices
    // Priority queue to go through the algorithm
    priority_queue<pair<long long, double>, vector<pair<long long, double> >, prioritize> pq;
    // Set distances to INF
    for (auto v : vertices) {
        distances[v] = INF;
        pq.push(make_pair(v, INF));
        pred[v] = 0;
    }
    // Set start to 0
    distances.at(startV) = 0;
    pq.emplace(startV, 0);
    // Loop until the queue is empty
    while (!pq.empty()) {
        auto curV = pq.top();
        pq.pop();
        // Check if the curV is visited
        if (curV.second == INF) {
            break;
        } else if (visitedSet.count(curV.first) == 1) {
            continue;
        } else {
            visitedSet.insert(curV.first);
        }
        // Call the helper function that finds the shorter distances
        findShorterDist(curV.first, curV.second, distances, pq, pred, G);
    }
}
//
// BuildGraph
// runApplication HELPER FUNCTION
// Function that builds the graph using the Footways and Nodes
void BuildGraph(graph<long long, double> &G, vector<FootwayInfo> &Footways,
map<long long, Coordinates> &Nodes) {
    // Loop through the Nodes and add vertices
    for (auto node : Nodes) {
        G.addVertex(node.first);
    }
    // Loop through the footways and add the edges
    for (auto Footway : Footways) {
        vector<long long> nodes = Footway.Nodes;
        for (unsigned int i = 0; i < nodes.size(); i++) {
            if (i + 1 < nodes.size()) {
                long long v1 = nodes.at(i);
                long long v2 = nodes.at(i + 1);
                Coordinates point1 = Nodes.at(v1);
                Coordinates point2 = Nodes.at(v2);
                // Get the distance between the 2 vertices
                double weight = distBetween2Points(point1.Lat, point1.Lon, point2.Lat, point2.Lon);
                // Add an edge in both directions
                G.addEdge(v1, v2, weight);
                G.addEdge(v2, v1, weight);
            }
        }
    }
}
//
// findBuildings
// runApplication HELPER FUNCTION
// Function that finds the buildings that the
// user enters from the Buildings vector
// Finds both the starting and destination
// buildings and returns them through the parameters.
void findBuildings(vector<BuildingInfo> &Buildings, bool &start,
bool &dest, Coordinates &startB, Coordinates &destB, BuildingInfo &b,
BuildingInfo &b1, string &destBuilding, string &startBuilding) {
    // Search for the starting building and update the parameters
    for (auto building : Buildings) {  // Loop through the Buildings
        size_t foundStart = building.Fullname.find(startBuilding);
        // Look for the abbreviation and fullname
        if (startBuilding == building.Abbrev || (foundStart < startBuilding.size() && foundStart >= 0)) {
            start = true;
            startB = building.Coords;
            b = building;
            break;
        }
    }
    // Search for the destination building and update the parameters
    for (auto building : Buildings) {  // Loop through the buildings
        size_t found = building.Fullname.find(destBuilding);
        // Look for the abbreviation and fullname
        if (destBuilding == building.Abbrev || (found < destBuilding.size() && found >= 0)) {
            dest = true;
            destB = building.Coords;
            b1 = building;
            break;
        }
    }
}
//
// findStartNode
// findPath HELPER FUNCTION
// Function that finds only the startNode
// Function is called only when both the destination
// and starting buildings are the same.
// Finds the node that is closest to the buildings by finding the min distance
void findStartNode(vector<FootwayInfo> &Footways, Coordinates &startNode,
map<long long, Coordinates> &Nodes, Coordinates &startB) {
    double closestSNode = 0;
    for (auto footway : Footways) {  // Loop through the footways
        vector<long long> nodes = footway.Nodes;
        for (auto point : nodes) {  // Loop through the nodes of a footway
            Coordinates p = Nodes.at(point);
            double tempStart = distBetween2Points(startB.Lat, startB.Lon, p.Lat, p.Lon);
            // Check if this is the closestNode
            if (tempStart < closestSNode || closestSNode == 0) {
                if (closestSNode != tempStart) {
                    // Only update if the distance is not same
                    closestSNode = tempStart;
                    startNode = p;
                }
            }
        }
    }
}
//
// findBothNodes
// findPath HELPER FUNCTION
// Function that finds both the starting and
// destination nodes for the given buildings
// Function is called when the start and destination are different
void findBothNodes(vector<FootwayInfo> &Footways, Coordinates &startNode,
Coordinates &destNode, map<long long, Coordinates> &Nodes,
Coordinates &startB, Coordinates &destB) {
    double closestSNode = 0;
    double closestDNode = 0;
    // Loop through the footways
    for (auto footway : Footways) {
        // Loop through every node in a footway
        vector<long long> nodes = footway.Nodes;
        for (auto point : nodes) {
            Coordinates p = Nodes.at(point);
            double tempStart = distBetween2Points(startB.Lat, startB.Lon, p.Lat, p.Lon);
            double tempDest = distBetween2Points(destB.Lat, destB.Lon, p.Lat, p.Lon);
            // Check the distance between start and node
            if (tempStart < closestSNode || closestSNode == 0) {
                if (closestSNode != tempStart) {
                    // Don't change the node if the distance is same
                    closestSNode = tempStart;
                    startNode = p;
                }
            }
            // Check the distance between destination and node
            if (tempDest < closestDNode || closestDNode == 0) {
                if (closestDNode != tempDest) {
                    // Don't change the node if the distance is same
                    closestDNode = tempDest;
                    destNode = p;
                }
            }
        }
    }
}
//
// showPath
// findPath HELPER FUNCTION
// Function that prints the path from starting building to destination
void showPath(Coordinates &destNode, map<long long, long long> &pred,
map<long long, double> &distances, Coordinates &startNode) {
    // Creating variables to find the path
    long long path = destNode.ID;
    unsigned int length = 0;
    stack<long long> printPath;
    cout << endl;
    cout << "Navigating with Dijkstra..." << endl;
    if (pred.at(destNode.ID) != 0) {  // Checks if there is a path
        cout << "Distance to dest: " << distances.at(destNode.ID);
        cout << " miles" << endl;
        cout << "Path: ";
        while (pred.at(path) != 0) {
            // Add the predecessor from destination till start to the stack
            printPath.push(path);
            path = pred.at(path);
        }
        printPath.push(startNode.ID);  // Add the start to the stack
        while (!printPath.empty()) {
            // Loop till stack is empty and print the path
            cout << printPath.top();
            printPath.pop();
            if (printPath.size() != 0) {
                cout << "->";
            }
        }
    } else {  // Prints this if there is no path
        cout << "Sorry, destination unreachable" << endl;
    }
}
//
// printNodes
// findPath HELPER FUNCTION
// Function that prints the starting and destination nodes
void printNodes(Coordinates startNode, Coordinates destNode) {
    // Prints the information of startNode and destNode
    cout << endl;
    cout << "Nearest start node:" << endl;
    cout << " " << startNode.ID << endl;
    cout << " (" << startNode.Lat << ", " << startNode.Lon << ")" << endl;
    cout << "Nearest destination node:" << endl;
    cout << " " << destNode.ID << endl;
    cout << " (" << destNode.Lat << ", " << destNode.Lon << ")" << endl;
}
//
// printBuildings
// findPath HELPER FUNCTION
// Function that prints the starting and destination buildings
void printBuildings(Coordinates startB, Coordinates destB,
BuildingInfo b, BuildingInfo b1) {
    // Prints the information of the buildings including coordinates
    cout << "Starting point:\n" << " " << b.Fullname << endl;
    cout << " (" << startB.Lat << ", " << startB.Lon << ")" << endl;
    cout << "Destination point:\n" << " " << b1.Fullname << endl;
    cout << " (" << destB.Lat << ", " << destB.Lon << ")" << endl;
}
//
// printDetails
// runApplication HELPER FUNCTION
// Function that prints the details of the map
void printDetails(map<long long, Coordinates> Nodes,
vector<FootwayInfo> Footways, vector<BuildingInfo> Buildings) {
    // Prints the sizes of the Nodes, Buildings and Footways data structures
    cout << endl;
    cout << "# of nodes: " << Nodes.size() << endl;
    cout << "# of footways: " << Footways.size() << endl;
    cout << "# of buildings: " << Buildings.size() << endl;
}
//
// printGraphDetails
// runApplication HELPER FUNCTION
// Function prints the details of the graph
void printGraphDetails(graph<long long, double> G) {
    // Prints the size of the graph and it's edges
    cout << "# of vertices: " << G.NumVertices() << endl;
    cout << "# of edges: " << G.NumEdges() << endl;
    cout << endl;
}
//
// findPath
// runApplication HELPER FUNCTION
// Function finds the path between starting and destination Nodes
// Calls the dijkstra's function which does most of the work
// Function just sets up to call other helper functions and prints outputs
void findPath(bool &start, bool &dest, Coordinates &startB,
Coordinates &destB, BuildingInfo b, BuildingInfo b1,
vector<FootwayInfo> Footways, map<long long, Coordinates> Nodes,
graph<long long, double> &G) {
    if (start && dest) {  // Checks if both buildings exist
        // Calls function to print the buildings
        printBuildings(startB, destB, b, b1);
        Coordinates startNode; Coordinates destNode;
        // Checks if both the start and destination are the same building
        if (startB.ID == destB.ID) {
            // Call the function to find starting node
            findStartNode(Footways, startNode, Nodes, startB);
            destNode = startNode;
        } else {
            // Call the function to find both nodes if they're different
            findBothNodes(Footways, startNode, destNode, Nodes, startB, destB);
        }
        printNodes(startNode, destNode);  // Calls function to print the nodes
        if (destNode.ID == startNode.ID) {
            // If the destination and start are the same
            // we just output the distance and path
            cout << endl;
            cout << "Navigating with Dijkstra..." << endl;
            // Distance is 0 since it is same
            cout << "Distance to dest: 0 miles" << endl;
            // Path is same as the destination node
            cout << "Path: " << destNode.ID;
        } else {  // Else we call the Dijkstra's function and show the path
            map<long long, double> distances;
            long long startV = startNode.ID;
            map<long long, long long> pred;
            Dijkstra(G, startV, distances, pred);
            showPath(destNode, pred, distances, startNode);
        }
        cout << endl;
    } else {
        // Print error statements if either one of the buildings are not found
        if (!start) {
            cout << "Start building not found" << endl;
        } else if (!dest) {
            cout << "Destination building not found" << endl;
        }
    }
}
//
// runApplication
//
// Function is responsible for running the application
// and asking input from user
// Contains the loop that controls the application
// Calls all the other helper functions to find paths
// between two points on the map.
void runApplication(map<long long, Coordinates> Nodes,
vector<FootwayInfo> Footways, vector<BuildingInfo> Buildings,
graph<long long, double> &G) {
    printDetails(Nodes, Footways, Buildings);
    BuildGraph(G, Footways, Nodes);
    printGraphDetails(G);
    // Navigation from building to building
    string startBuilding, destBuilding;
    cout << "Enter start (partial name or abbreviation), or #> ";
    getline(cin, startBuilding);
    // Loop runs until user exits
    while (startBuilding != "#") {
        cout << "Enter destination (partial name or abbreviation)> ";
        getline(cin, destBuilding);
        bool start = false; bool dest = false;
        Coordinates startB; Coordinates destB;
        BuildingInfo b; BuildingInfo b1;
        // Calls function to find buildings.
        findBuildings(Buildings, start, dest, startB, destB, b, b1, destBuilding, startBuilding);
        // Calls function to find the path between the buildings.
        findPath(start, dest, startB, destB, b, b1, Footways, Nodes, G);
        // another navigation?
        cout << endl;
        cout << "Enter start (partial name or abbreviation), or #> ";
        getline(cin, startBuilding);
    }
}
//
// main
//
// Contains information to create the application
int main() {
    // maps a Node ID to it's coordinates (lat, lon)
    map<long long, Coordinates> Nodes;
    // info about each footway, in no particular order
    vector<FootwayInfo> Footways;
    // info about each building, in no particular order
    vector<BuildingInfo> Buildings;
    XMLDocument xmldoc;
    graph<long long, double> G;
    cout << "** Navigating UIC open street map **" << endl; cout << endl;
    cout << std::setprecision(8);
    string def_filename = "map.osm"; string filename;
    cout << "Enter map filename> ";
    getline(cin, filename);
    if (filename == "") {
        filename = def_filename;
    }
    // Load XML-based map file
    if (!LoadOpenStreetMap(filename, xmldoc)) {
        cout << "**Error: unable to load open street map." << endl;
        cout << endl;
        return 0;
    }
    // Read the nodes, which are the various known positions on the map:
    int nodeCount = ReadMapNodes(xmldoc, Nodes);
    // Read the footways, which are the walking paths:
    int footwayCount = ReadFootways(xmldoc, Footways);
    // Read the university buildings:
    int buildingCount = ReadUniversityBuildings(xmldoc, Nodes, Buildings);
    // Stats
    assert(nodeCount == (int)Nodes.size());
    assert(footwayCount == (int)Footways.size());
    assert(buildingCount == (int)Buildings.size());
    // Calls the function that runs the application
    // This function calls several other functions which
    // calls other functions that help in running the application.
    runApplication(Nodes, Footways, Buildings, G);
    // done:
    cout << "** Done **" << endl;
    return 0;
}
