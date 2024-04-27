#pragma once

#include <map>
#include <set>
#include <vector>
#include <limits>
#include <algorithm>

using namespace std;

template <typename VertexT, typename WeightT>
class graph {
private:
    map<VertexT, map<VertexT, WeightT>> adjList;

public:
    graph() = default;

    bool addVertex(const VertexT& v) {
        if (adjList.find(v) != adjList.end()) {
            return false;
        }
        adjList[v];
        return true;
    }

    bool addEdge(const VertexT& from, const VertexT& to, const WeightT& weight) {
        if (adjList.find(from) == adjList.end() || adjList.find(to) == adjList.end()) {
            return false;
        }
        adjList[from][to] = weight;
        return true;
    }

    bool getWeight(const VertexT& from, const VertexT& to, WeightT& weight) const {
        auto it = adjList.find(from);
        if (it != adjList.end()) {
            auto it2 = it->second.find(to);
            if (it2 != it->second.end()) {
                weight = it2->second;
                return true;
            }
        }
        return false;
    }

    set<VertexT> neighbors(const VertexT& v) const {
        set<VertexT> result;
        auto it = adjList.find(v);
        if (it != adjList.end()) {
            for (auto& pair : it->second) {
                result.insert(pair.first);
            }
        }
        return result;
    }

 const map<VertexT, WeightT>& getAdjacent(const VertexT& v) const {
        static const map<VertexT, WeightT> empty; 
        auto it = adjList.find(v);
        if (it != adjList.end()) {
            return it->second;
        }
        return empty;
    }

    vector<VertexT> getVertices() const {
        vector<VertexT> vertices;
        for (auto& pair : adjList) {
            vertices.push_back(pair.first);
        }
        return vertices;
    }

    int NumVertices() const {
        return adjList.size();
    }

    int NumEdges() const {
        int edgeCount = 0;
        for (const auto& vertex : adjList) {
            edgeCount += vertex.second.size();
        }
        return edgeCount;
    }
};
