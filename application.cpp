
#include "application.h"

#include <cassert>
#include <cstdlib>
#include <cstring>
#include <iomanip> 
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <numeric>


#include "dist.h"
#include "graph.h"
#include "osm.h"
#include "tinyxml2.h"

using namespace std;
using namespace tinyxml2;

double INF = numeric_limits<double>::max();

struct prioritize {
    bool operator()(const pair<long long, double>& p1, const pair<long long, double>& p2) const {
        return p1.second > p2.second;
    }
};

graph<long long, double> buildGraph(
    const map<long long, Coordinates>& Nodes,
    const vector<FootwayInfo>& Footways,
    const vector<BuildingInfo>& Buildings) {
    graph<long long, double> G;

for (auto it = Nodes.begin(); it != Nodes.end(); ++it) {
    const auto& node = *it;
    G.addVertex(node.first);
}

for (const auto& footway : Footways) {
    for (auto it = footway.Nodes.begin(); it != footway.Nodes.end() - 1; ++it) {
        long long fromID = *it;
        long long toID = *(it + 1);
        if (Nodes.count(fromID) && Nodes.count(toID)) {
            double lat1 = Nodes.at(fromID).Lat;
            double lon1 = Nodes.at(fromID).Lon;
            double lat2 = Nodes.at(toID).Lat;
            double lon2 = Nodes.at(toID).Lon;
            double dist = distBetween2Points(lat1, lon1, lat2, lon2);
            G.addEdge(fromID, toID, dist);
            G.addEdge(toID, fromID, dist);  
            }
        }
    }

for (auto it = Buildings.begin(); it != Buildings.end(); ++it) {
    const auto& building = *it;
    auto [buildingNUM, buildingSide, buildingVer] = std::make_tuple(building.Coords.ID, building.Coords.Lat, building.Coords.Lon);
    G.addVertex(buildingNUM);

    for (auto node_it = Nodes.begin(); node_it != Nodes.end(); ++node_it) {
        const auto& node = *node_it;
        if (node.second.OnFootway) { 
        // Gets latitude and longitude of the current node
        double latitude = node.second.Lat;
        double longitude = node.second.Lon;
        // Calculates the distance between the building and the current node
        double dist = distBetween2Points(buildingSide, buildingVer, latitude, longitude);
            double maxDistSquared = 0.041 * 0.041; 
            if (dist * dist <= maxDistSquared) { 
                G.addEdge(node.first, buildingNUM, dist); 
                G.addEdge(buildingNUM, node.first, dist); 
            }
        }
    }
}
    return G;
}


vector<long long> dijkstra(
    const graph<long long, double>& G,
    long long start,
    long long target,
    const set<long long>& ignoreNodes) {
    
    priority_queue<pair<long long, double>, vector<pair<long long, double>>, prioritize> worklist;
    map<long long, double> shortestDistances;
    map<long long, long long> previousNode;
    double INF = numeric_limits<double>::max();

    // Initialize distances and priority queue
    for (const auto& vertex : G.getVertices()) {
        shortestDistances[vertex] = (vertex == start ? 0 : INF);
        worklist.emplace(vertex, shortestDistances[vertex]);
    }

while (!worklist.empty()) {
    auto [currentNode, currentDist] = worklist.top();
    worklist.pop();

    if (currentNode == target) break; // If target is reached, exit early

    if (currentNode != start && currentNode != target && ignoreNodes.count(currentNode)) {
        continue; // Skip processing for ignored nodes, except for start or target
    }

    // Process each adjacent vertex
    for (const auto& [adjNode, weight] : G.getAdjacent(currentNode)) {
        if (ignoreNodes.count(adjNode) && adjNode != target) continue; // Skip ignored nodes unless it's the target

        double newDist = currentDist + weight;
        if (newDist < shortestDistances[adjNode]) {
            shortestDistances[adjNode] = newDist;
            previousNode[adjNode] = currentNode;
            worklist.push(make_pair(adjNode, newDist)); // Update priority queue
        }
    }
}


    // Reconstruct the path from end to start using the previous node map
    vector<long long> path;
    if (shortestDistances[target] != INF) {
        for (long long at = target; at != start; at = previousNode[at]) {
            path.push_back(at);
        }
        path.push_back(start);
        reverse(path.begin(), path.end());
    }

    return path;
}




double pathLength(const graph<long long, double>& G, const vector<long long>& path) {
    double length = 0.0;
    double weight;
    for (size_t i = 0; i + 1 < path.size(); i++) {
        bool res = G.getWeight(path.at(i), path.at(i + 1), weight);
        assert(res);
        length += weight;
    }
    return length;
}

void outputPath(const vector<long long>& path) {
    for (size_t i = 0; i < path.size(); i++) {
        cout << path.at(i);
        if (i != path.size() - 1) {
            cout << "->";
        }
    }
    cout << endl;
}

void application(
    const vector<BuildingInfo>& Buildings,
    const graph<long long, double>& G) {
    string person1Building, person2Building;

    set<long long> buildingNodes;
    for (const auto& building : Buildings) {
        buildingNodes.insert(building.Coords.ID);
    }

    cout << endl;
    cout << "Enter person 1's building (partial name or abbreviation), or #> ";
    getline(cin, person1Building);

    while (person1Building != "#") {
        cout << "Enter person 2's building (partial name or abbreviation)> ";
        getline(cin, person2Building);

        //
        // find the building coordinates
        //
        bool foundP1 = false;
        bool foundP2 = false;
        Coordinates P1Coords, P2Coords;
        string P1Name, P2Name;

        for (const BuildingInfo& building : Buildings) {
            if (building.Abbrev == person1Building) {
                foundP1 = true;
                P1Name = building.Fullname;
                P1Coords = building.Coords;
            }
            if (building.Abbrev == person2Building) {
                foundP2 = true;
                P2Name = building.Fullname;
                P2Coords = building.Coords;
            }
        }

        for (const BuildingInfo& building : Buildings) {
            if (!foundP1 &&
                building.Fullname.find(person1Building) != string::npos) {
                foundP1 = true;
                P1Name = building.Fullname;
                P1Coords = building.Coords;
            }
            if (!foundP2 && building.Fullname.find(person2Building) != string::npos) {
                foundP2 = true;
                P2Name = building.Fullname;
                P2Coords = building.Coords;
            }
        }

        if (!foundP1) {
            cout << "Person 1's building not found" << endl;
        } else if (!foundP2) {
            cout << "Person 2's building not found" << endl;
        } else {
            cout << endl;
            cout << "Person 1's point:" << endl;
            cout << " " << P1Name << endl;
            cout << " (" << P1Coords.Lat << ", " << P1Coords.Lon << ")" << endl;
            cout << "Person 2's point:" << endl;
            cout << " " << P2Name << endl;
            cout << " (" << P2Coords.Lat << ", " << P2Coords.Lon << ")" << endl;

            string destName;
            Coordinates destCoords;

            Coordinates centerCoords = centerBetween2Points(
                P1Coords.Lat, P1Coords.Lon, P2Coords.Lat, P2Coords.Lon);

            double minDestDist = numeric_limits<double>::max();

            for (const BuildingInfo& building : Buildings) {
                double dist = distBetween2Points(
                    centerCoords.Lat, centerCoords.Lon,
                    building.Coords.Lat, building.Coords.Lon);
                if (dist < minDestDist) {
                    minDestDist = dist;
                    destCoords = building.Coords;
                    destName = building.Fullname;
                }
            }

            cout << "Destination Building:" << endl;
            cout << " " << destName << endl;
            cout << " (" << destCoords.Lat << ", " << destCoords.Lon << ")" << endl;

            vector<long long> P1Path = dijkstra(G, P1Coords.ID, destCoords.ID, buildingNodes);
            vector<long long> P2Path = dijkstra(G, P2Coords.ID, destCoords.ID, buildingNodes);

            // This should NEVER happen with how the graph is built
            if (P1Path.empty() || P2Path.empty()) {
                cout << endl;
                cout << "At least one person was unable to reach the destination building. Is an edge missing?" << endl;
                cout << endl;
            } else {
                cout << endl;
                cout << "Person 1's distance to dest: " << pathLength(G, P1Path);
                cout << " miles" << endl;
                cout << "Path: ";
                outputPath(P1Path);
                cout << endl;
                cout << "Person 2's distance to dest: " << pathLength(G, P2Path);
                cout << " miles" << endl;
                cout << "Path: ";
                outputPath(P2Path);
            }
        }

        //
        // another navigation?
        //
        cout << endl;
        cout << "Enter person 1's building (partial name or abbreviation), or #> ";
        getline(cin, person1Building);
    }
}