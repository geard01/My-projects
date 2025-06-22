#ifndef __EMPTY_H
#define __EMPTY_H
#include "Eigen/Eigen"
#include <iostream>
#include <fstream>
#include <queue>
#include <set>
#include <vector>
#include <string>
#include <cmath>

using namespace std;
using namespace Eigen;
namespace ProjectLibrary
{
inline double evaluateExpressionf0(double x, double y, const string& function)
{
    if (function=="1")
        return abs(pow((x-0.5),(2)) + pow((y-0.5),(2)));
    else if (function=="2")
        return abs(exp(-(pow(x - 0.5, 2) + pow(y - 0.5, 2)) / (2.0 * pow(0.1, 2))));
    else if (function=="3")
        return abs(sqrt(pow(x - 0.5, 2) + pow(y - 0.5, 2)) - 0.25);
    else if (function=="4")
        return abs(pow((pow((x - 0.5),2) + pow((y - 0.5),2)),2) - 2*(pow((x - 0.5),2) - pow((y - 0.5),2)));
    return abs(sin(x) - cos(y));
}


inline double calculateGradientNorm(double (*evaluateExpressionf0)(double, double, const string&), double x, double y, const string& function)
{
    const double h = 0.0001;
    double dx = (evaluateExpressionf0(x + h, y, function) - evaluateExpressionf0(x - h, y, function)) / (2 * h);
    double dy = (evaluateExpressionf0(x, y + h, function) - evaluateExpressionf0(x, y - h, function)) / (2 * h);
    double gradientNorm = sqrt(pow(dx, 2) + pow(dy, 2));
    return gradientNorm;
}

inline double doubleIntegral(const string& function)
{
    if (function=="1")
        return 0.1666666666667;
    else if (function=="2")
        return 0.06285178;
    else if (function=="3")
        return 0.1653230000000;
    else if (function=="4")
        return 0.1726510000000;
    return 0.4003236626453;
}
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

    static constexpr double geometricTol = 1.0e-12;

    inline bool operator>(const Edge& other)
    {return length-geometricTol > other.length;}

    inline bool operator<(const Edge& other)
    {return length+geometricTol < other.length;}

    inline bool operator==(const Edge& other)
    {return id == other.id;}

    inline bool operator!= (const Edge& other)
    {return !(id == other.id);}

};
// il triangolo è formato da area, id, lati e vertici
class Triangle
{
private:
public:
    static string mode;
    static string function;
    bool active = true;
    double area;
    unsigned int id;
    vector<Edge> edges;
    vector<Vertex> vertices;
    Vertex center;
    double gradValue;
    Triangle (vector<Edge> edge, unsigned int id):
        id(id), edges(edge)
    {
        vertices.push_back(edges[0].start);
        vertices.push_back(edges[0].finish);
        if (edges[1].start != edges[0].start && edges[1].start != edges[0].finish)
            vertices.push_back(edges[1].start);
        else
            vertices.push_back(edges[1].finish);
        area = abs(vertices[0].x*(vertices[1].y-vertices[2].y) +
                   vertices[1].x*(vertices[2].y-vertices[0].y)  +
                   vertices[2].x*(vertices[0].y-vertices[1].y)) / 2;
        center.x = (vertices[0].x + vertices[1].x + vertices[2].x)/3;
        center.y = (vertices[0].y + vertices[1].y + vertices[2].y)/3;
        gradValue=calculateGradientNorm(evaluateExpressionf0, center.x, center.y, function);
    }

    Triangle():
        id(0)
    {}
    static constexpr double geometricTol = 1.0e-12;
    // INSERIRE NELLE FORMULE LA TOLLERANZA GEOMETRICA
    // quando confronto due triangoli, uno è maggiore dell'altro in base alle loro aree
    inline bool operator>(const Triangle& other)
    {
        if (mode=="1")
            return (area-geometricTol) > (other.area);
        else
            return (0.0001*gradValue+area-geometricTol) > (0.0001*(other.gradValue)+other.area);
    }

    inline bool operator<(const Triangle& other)
    {
        if (mode=="1")
            return (area+geometricTol) < (other.area);
        else
            return (0.0001*gradValue+area+geometricTol) < (0.0001*(other.gradValue)+other.area);
    }



    inline bool operator==(const Triangle& other)
    {return id == other.id;}

    inline bool operator!=(const Triangle& other)
    {return !(id == other.id);}
    void set_area()
    {
        area = abs(vertices[0].x*(vertices[1].y-vertices[2].y) +
                   vertices[1].x*(vertices[2].y-vertices[0].y)  +
                   vertices[2].x*(vertices[0].y-vertices[1].y)) / 2;
    }
    void set_grad()
    {
        center.x = (vertices[0].x + vertices[1].x + vertices[2].x)/3;
        center.y = (vertices[0].y + vertices[1].y + vertices[2].y)/3;
        gradValue=calculateGradientNorm(evaluateExpressionf0, center.x, center.y, function);
    }
};

//sfrutto l'algortimo di quicksort per ordinare i triangoli in base alla loro area
void Refine(vector<Triangle>& triangles, vector<Edge>& edges, vector<Vertex>& vertices, unsigned int &n, double &area, double& exactarea);
//importa i vertici della mesh triangolare
bool ImportCell0Ds(vector<Vertex>& vertices, unsigned int n, string& test);
//importa i lati della mesh triangolare
bool ImportCell1Ds(vector<Edge>& edges, vector<Vertex>& vertices, unsigned int n, string& test);
//importa i triangoli della mesh
bool ImportCell2Ds(vector<Triangle>& triangles, vector<Edge>& edges, vector<Vertex>& vertices, unsigned int n, string &test);
// sfrutto l'algoritmo di insertion sort per ordinare array di piccole dimensioni, oppure parzialmente ordinati
inline void insertionSort(vector<Edge> &edges)
{
    int n=edges.size();
    for (int i = 0; i < n; i++)
    {
        Edge key = edges[i];
        int j = i - 1;
        while (j >= 0 && edges[j] < key)
        {
            edges[j + 1] = edges[j];
            j = j - 1;
        }
        edges[j + 1] = key;
    }
}
// ricava il vertice opposto dato un lato e un triangolo
inline Vertex getOppositeVertex(Triangle &triangle, Edge &edge)
{
    for (unsigned int i = 0; i < 3; i++)
    {
        if (triangle.vertices[i] != (edge.start) &&
            triangle.vertices[i] != (edge.finish))
            return triangle.vertices[i];
    }
    throw(runtime_error("Something went wrong getting the opposite vertex"));
}
//divido triangolo in 2
void split2(vector<Triangle> &triangles, vector<Edge> &edges, vector<Vertex> &vertices, unsigned int m,
            deque<unsigned int> &tempId, unsigned int &k, bool &permissible, deque<unsigned int> &tempId1, double &area);
//divido i triangoli successivi in 3 sottotriangoli
void split3(vector<Triangle> &triangles, vector<Edge> &edges, vector<Vertex> &vertices,
            deque<unsigned int> &tempId, unsigned int &k, bool &permissible, deque<unsigned int> &tempId1, double &area);
void split2again (vector<Triangle> &triangles, vector<Edge> &edges,vector<Vertex> &vertices,
                 unsigned int &k, bool &permissible, deque<unsigned int> &tempId, deque<unsigned int> &tempId1, double &area);
template <typename T>
inline void Pushback(vector<T> &edges, T &edge)
{
    edges.push_back(edge);
}
template <typename T>
inline void Erase(vector<T> &triangles, unsigned int index)
{
    if (index < triangles.size())
    {
        triangles.erase(triangles.begin() + index);
    }
}
inline unsigned int massimoElementoAttivo(vector<Triangle> &vettore)
{
    Triangle massimo;
    massimo.area = 0;
    for (Triangle &elemento : vettore)
    {
        if (elemento.active) { // Verifica se l'elemento è attivo
            if (elemento > massimo)
                massimo = elemento;
        }
    }
    return massimo.id;
}
inline Vertex set_mid(Edge &edge)
{
    Vertex mid;
    mid.x = (edge.finish.x + edge.start.x) / 2;
    mid.y = (edge.finish.y + edge.start.y) / 2;
    return mid; // Assegnare id e aggiornare vettore di vertici
}

}

#endif // __EMPTY_H
