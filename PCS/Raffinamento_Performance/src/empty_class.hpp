#ifndef __EMPTY_H
#define __EMPTY_H
#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <queue>
#include <vector>


using namespace std;
using namespace Eigen;
namespace ProjectLibrary
{
// il vertice è formato da due punti nel piano: x e y
struct Vertex
{
    double x;
    double y;
    unsigned int id;

    Vertex(double x, double y, unsigned int id):
        x(x), y(y), id(id)
    {}
    Vertex():
        x(0), y(0), id(0)
    {}

    inline bool operator>(const Vertex& other) const
    {return id > other.id;}
    inline bool operator<(const Vertex& other) const
    {return id < other.id;}
    inline bool operator==(const Vertex& other) const
    {return id == other.id;}
    inline bool operator!=(const Vertex& other) const
    {return !(id == other.id);}

};

struct Edge
{
    bool active = true;
    static constexpr double geometricTol = 1e-12;
    Vertex start;
    Vertex finish ;
    vector<unsigned int> adjTriangles;
    unsigned int id;
    double length;
    Edge(Vertex start, Vertex finish, unsigned int id):
        start(start), finish(finish), id(id)
    {
        length = sqrt(pow(finish.y-start.y, 2) + pow(finish.x-start.x, 2));
    }
    Edge():
        id(0)
    {}
    void set_length()
    {
        length = sqrt(pow(finish.y-start.y, 2) + pow(finish.x-start.x, 2));
    }

    inline bool operator>(const Edge& other)
    {return length-geometricTol > other.length;}

    inline bool operator<(const Edge& other)
    {return length+geometricTol  < other.length;}

    inline bool operator==(const Edge& other)
    {return id == other.id;}

    inline bool operator!= (const Edge& other)
    {return !(id == other.id);}

};
// il triangolo è formato da area, id e lati
class Triangle
{
private:
public:
    static constexpr double geometricTol = 1e-12;
    bool active = true;
    double area;
    unsigned int id;
    vector<Edge> edges;
    Triangle (vector<Edge> edge, unsigned int id):
        id(id), edges(edge)
    {
        double s = (edges[0].length + edges[1].length + edges[2].length) / 2;
        area = sqrt(s * (s - edges[0].length) * (s - edges[1].length) * (s - edges[2].length));
    }

    Triangle():
        id(0)
    {}

    // quando confronto due triangoli, uno è maggiore dell'altro in base alle loro aree
    inline bool operator>(const Triangle& other)
    {return area-geometricTol > other.area;}

    inline bool operator<(const Triangle& other)
    {return area+ geometricTol < other.area;}

    inline bool operator==(const Triangle& other)
    {return id == other.id;}

    inline bool operator!=(const Triangle& other)
    {return !(id == other.id);}

    void set_area()
    {
        double s = (edges[0].length + edges[1].length + edges[2].length) / 2;
        area = sqrt(s * (s - edges[0].length) * (s - edges[1].length) * (s - edges[2].length));
    }
};

void Refine(vector<Triangle>& triangles, vector<Edge>& edges, vector<Vertex>& vertices, unsigned int &n);
//importa i vertici della mesh triangolare
bool ImportCell0Ds(vector<Vertex>& vertices, unsigned int n, string& test);
//importa i lati della mesh triangolare
bool ImportCell1Ds(vector<Edge>& edges, vector<Vertex>& vertices, unsigned int n, string& test);
//importa i triangoli della mesh
bool ImportCell2Ds(vector<Triangle>& triangles, vector<Edge>& edges, vector<Vertex>& vertices, unsigned int n, string &test);
//divido triangolo in 2
void split2(vector<Triangle> &triangles, vector<Edge> &edges, vector<Vertex> &vertices, unsigned int m,
            deque<unsigned int> &tempId, unsigned int &k, bool &permissible, deque<unsigned int> &tempId1);
//divido i triangoli successivi in 3 sottotriangoli
void split3(vector<Triangle> &triangles, vector<Edge> &edges, vector<Vertex> &vertices,
            deque<unsigned int> &tempId, unsigned int &k, bool &permissible, deque<unsigned int> &tempId1);
void split2again (vector<Triangle> &triangles, vector<Edge> &edges,vector<Vertex> &vertices,
                 unsigned int &k, bool &permissible, deque<unsigned int> &tempId, deque<unsigned int> &tempId1);

template <typename T>
inline void Pushback(vector<T> &edges, T &edge)
{
    edges.push_back(edge);
}


inline void insertionSort(vector<Edge> &edge)
{
    unsigned int n= edge.size();
    for (unsigned int i = 0; i < n; i++)
    {
        Edge key = edge[i];
        int j = i - 1;
        while (j >= 0 && edge[j] < key)
        {
            edge[j + 1] = edge[j];
            j = j - 1;
        }
        edge[j + 1] = key;
    }
}
template <typename T>
inline void Erase(vector<T> &triangles, unsigned int index)
{
    if (index < triangles.size())
    {
        triangles.erase(triangles.begin() + index);
    }
}

inline void massimoElementoAttivo(vector<Triangle> &vettore, unsigned int &m, unsigned int &p)
{
    Triangle massimo;
    massimo.area = -1;
    bool enteredActiveCondition = false;

    for (unsigned int i=p; i<vettore.size(); i++)
    {
        if (vettore[i].active)
        {
            enteredActiveCondition = true;
            if (vettore[i] > massimo)
                massimo = vettore[i];
        }
        if (!enteredActiveCondition)
        {
            p++;
        }
        if (vettore[m].area - massimo.area < 1.0e-12)
            break;
    }

    m = massimo.id;
}

inline Vertex getOppositeVertex(Triangle &triangle, Edge &edge)
{
    for (unsigned int i = 0; i < 3; i++)
    {
        if (triangle.edges[i].start != (edge.start) &&
            triangle.edges[i].start != (edge.finish))
            return triangle.edges[i].start;
        else if (triangle.edges[i].finish != (edge.start) &&
            triangle.edges[i].finish != (edge.finish))
            return triangle.edges[i].finish;
    }
    throw(runtime_error("Something went wrong getting the opposite vertex"));
}
// Prendo il punto medio di un lato
inline Vertex set_mid(Edge edge)
{
    Vertex mid;
    mid.x = (edge.finish.x + edge.start.x) / 2;
    mid.y = (edge.finish.y + edge.start.y) / 2;
    return mid; // Assegnare id e aggiornare vettore di vertici
}
}

#endif // __EMPTY_H

