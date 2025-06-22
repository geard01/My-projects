#ifndef __TEST_EMPTY_H
#define __TEST_EMPTY_H

#include <gtest/gtest.h>
#include "iostream"
#include "empty_class.hpp"

using namespace testing;
using namespace ProjectLibrary;
vector<Vertex> Vertices;
vector<Edge> Edges;
vector<Triangle> Triangles;
unsigned int n=100;
string test= "1";
bool y=ImportCell0Ds(Vertices, n,test);
bool x=ImportCell1Ds(Edges,Vertices,n, test);
bool z=ImportCell2Ds(Triangles,Edges,Vertices,n, test);
bool permissible=true;
deque<unsigned int> tempId;
deque<unsigned int> tempId1;
double area=0;

TEST(TestData, TestImport1)
{
    int data = x;
    EXPECT_EQ(data, 1);
}
//il test seguente controlla che ogni file venga importato con successo
TEST(TestData, TestImport2)
{
    int data = y;
    EXPECT_EQ(data, 1);
}
//il test seguente controlla che ogni file venga importato con successo
TEST(TestData, TestImport3)
{
    int data = z;
    EXPECT_EQ(data, 1);
}
// il test seguente controlla che la somma delle aree totali venga uguale a 1
TEST(TestArea, TestSum)
{
    double sum=0;
    double tolerance = 1.0e-14;
    for(const Triangle &t : Triangles)
    {
        sum += t.area;
    }
    cout<<sum<<endl;
    ASSERT_TRUE((1 - sum) <= tolerance);
}

TEST(TestEdges, TestAllvertices)
{
    bool flag = true;
    for(unsigned int i=0; i<Triangles.size();i++)
    {
        vector<Vertex> tempVerticesT = Triangles[i].vertices;
        set<unsigned int> verticesT;
        for(unsigned int j=0; j<tempVerticesT.size(); j++)
            verticesT.insert(tempVerticesT[j].id);
        set<unsigned int> verticesE;
        for(unsigned int k=0; k<3; k++)
        {
            verticesE.insert(Triangles[i].edges[k].finish.id);
            verticesE.insert(Triangles[i].edges[k].start.id);
        }

        //controllo se il set dei vertici del triangolo è uguale al set formato dagli estremi dei lati
        if(verticesE != verticesT)
        {
            flag = false;
            break;
        }

    }
    ASSERT_TRUE(flag);
}

TEST(TestAdjacency, TestAdjacency)
{
    set<unsigned int> setT;
    set<unsigned int> setE;
    for(Triangle& triangle: Triangles)
    {
        setT.insert(triangle.id);
        for(unsigned int i=0; i< triangle.edges.size();i++)
        {
            if (triangle.edges[i].adjTriangles.size()>0)
            {
                Triangle pointedTriangle = Triangles[triangle.edges[i].adjTriangles[0]];
                setE.insert(pointedTriangle.id);
            }
        }
    }
    EXPECT_EQ(setT, setE);
}

TEST(Testinsertionsort, TestinsertionsortEdge)
{
    Edge edge1, edge2, edge3, edge4, edge5;
    edge1.length=1;
    edge1.id=1;
    edge2.length=10;
    edge2.id=2;
    edge3.length=0.5;
    edge3.id=3;
    edge4.length=15;
    edge4.id=4;
    edge5.length=2.2;
    edge5.id=5;
    vector<Edge> quickvector = {edge1, edge2, edge3, edge4, edge5};
    vector<Edge> expectedquick = {edge4, edge2, edge5, edge1, edge3};
    insertionSort(quickvector);
    unsigned int count=0;
    for (unsigned int i=0; i<5; i++)
    {
        cout<<quickvector[i].id<<endl;
        if(quickvector[i]==expectedquick[i])
            count++;
    }
    EXPECT_EQ(count, 5);
}

TEST(Testmassimo, TestmassimoTriangle)
{
    Triangle triangle1, triangle2, triangle3;
    triangle1.area=20;
    triangle1.id=30;
    triangle2.area=10;
    triangle2.id=10;
    triangle3.area=30;
    triangle3.id=20;
    triangle3.active=false;
    vector<Triangle> vector = {triangle1, triangle2, triangle3};
    unsigned int id = massimoElementoAttivo(vector);
    EXPECT_EQ(id, 30);
}

TEST(TestSetMid, TestSetMidEdge)
{
    Vertex vertex1;
    vertex1.x=4;
    vertex1.y=6;
    Vertex vertex2;
    vertex2.x=2;
    vertex2.y=7;
    Edge edge1;
    edge1.length=1;
    edge1.id=1;
    edge1.start= vertex1;
    edge1.finish=vertex2;
    Vertex mid = set_mid(edge1);
    Vertex expectedmid;
    expectedmid.x=3;
    expectedmid.y=6.5;
    EXPECT_EQ(mid.x, expectedmid.x);
    EXPECT_EQ(mid.y, expectedmid.y);
}

TEST(TestGetOppositeVertex, TestOppositeVertexTriangle)
{
    Edge edge1, edge2, edge3;
    Vertex vertex1(1.0,1.0,1);
    Vertex vertex2(2.0,1.0,2);
    Vertex vertex3(1.5,2.0,3);
    edge1.start=vertex1;
    edge1.finish=vertex2;
    edge1.id=1;
    edge2.start=vertex2;
    edge2.finish=vertex3;
    edge2.id=2;
    edge3.start=vertex3;
    edge3.finish=vertex1;
    edge3.id=3;
    vector<Edge> edges = {edge1, edge2, edge3};
    Triangle triangle(edges,1);
    Vertex opposite= getOppositeVertex(triangle, edge2);
    EXPECT_EQ(opposite.x, vertex1.x);
    EXPECT_EQ(opposite.y, vertex1.y);
}

TEST(TestSplit2, TestSplit2)
{
    Vertex vertex1(1.0,1.0,89);
    Vertex vertex2(2.0,1.0,90);
    Vertex vertex3(1.5,1.5,91);
    Edge edge1(vertex1, vertex2, 232);
    Edge edge2(vertex2, vertex3, 233);
    Edge edge3(vertex3, vertex1, 234);
    vector<Edge> edges = {edge1, edge2, edge3};
    Triangle triangle(edges,144);
    Pushback(Triangles,triangle);
    unsigned int k=0;
    split2(Triangles, Edges, Vertices, Triangles.size()-1,tempId, k, permissible, tempId1, area);
    EXPECT_EQ(Vertices[Vertices.size()-1].x, 1.5);
    EXPECT_EQ(Vertices[Vertices.size()-1].y, 1.0);
    EXPECT_EQ(Edges[Edges.size()-3].finish, Vertices[Vertices.size()-1]);
    EXPECT_EQ(Edges[Edges.size()-3].start, vertex3);
    EXPECT_EQ(Triangles[Triangles.size()-1].edges[0].adjTriangles[0], Triangles[Triangles.size()-2].id);
    EXPECT_EQ(Triangles[Triangles.size()-2].edges[0].adjTriangles[0], Triangles[Triangles.size()-1].id);
}

TEST(TestSplit2again, TestSplit2again)
{
    Vertex vertex1(1.0,1.0,93);
    Vertex vertex2(2.0,1.0,94);
    Vertex vertex3(1.5,1.5,95);
    Edge edge1(vertex1, vertex2, 236);
    Edge edge2(vertex2, vertex3, 237);
    Edge edge3(vertex3, vertex1, 238);
    vector<Edge> edges1 = {edge1, edge2, edge3};
    Triangle triangle1(edges1,147);
    Vertex vertex4(1.5,0.5,95);
    Edge edge4(vertex2, vertex4, 239);
    Edge edge5(vertex4, vertex1, 240);
    vector<Edge> edges2 = {edge1, edge4, edge5};
    Triangle triangle2(edges2,148);
    triangle1.edges[0].adjTriangles.push_back(triangle2.id);
    triangle2.edges[0].adjTriangles.push_back(triangle1.id);
    Pushback(Triangles,triangle1);
    Pushback(Triangles,triangle2);
    insertionSort(Triangles[Triangles.size()-1].edges);
    unsigned int k=0;
    split2(Triangles, Edges, Vertices, Triangles.size()-2, tempId, k, permissible, tempId1, area);
    split2again(Triangles, Edges, Vertices, k, permissible,tempId, tempId1, area);
    EXPECT_EQ(Vertices[Vertices.size()-1].x, 1.5);
    EXPECT_EQ(Vertices[Vertices.size()-1].y, 1.0);
    EXPECT_EQ(Edges[Edges.size()-1].finish, Vertices[Vertices.size()-1]);
    EXPECT_EQ(Edges[Edges.size()-1].start, vertex4);
    EXPECT_EQ(Triangles[Triangles.size()-1].edges[0].adjTriangles[0], Triangles[Triangles.size()-2].id);
    EXPECT_EQ(Triangles[Triangles.size()-2].edges[0].adjTriangles[0], Triangles[Triangles.size()-1].id);
}

TEST(TestSplit3, TestSplit3)
{
    Vertex vertex1(1.0,1.0,97);
    Vertex vertex2(2.0,1.0,98);
    Vertex vertex3(1.5,1.5,99);
    Edge edge1(vertex1, vertex2, 241);
    Edge edge2(vertex2, vertex3, 242);
    Edge edge3(vertex3, vertex1, 243);
    vector<Edge> edges1 = {edge1, edge2, edge3};
    Triangle triangle1(edges1,153);
    Vertex vertex4(1,0.5,100);
    Edge edge4(vertex2, vertex4, 244);
    Edge edge5(vertex4, vertex1, 245);
    vector<Edge> edges2 = {edge1, edge4, edge5};
    Triangle triangle2(edges2,154);
    triangle1.edges[0].adjTriangles.push_back(triangle2.id);
    triangle2.edges[0].adjTriangles.push_back(triangle1.id);
    Pushback(Triangles,triangle1);
    Pushback(Triangles,triangle2);
    insertionSort(Triangles[Triangles.size()-1].edges);
    unsigned int k=0;
    split2(Triangles, Edges, Vertices, Triangles.size()-2, tempId, k, permissible, tempId1, area);
    split3(Triangles, Edges, Vertices,tempId, k, permissible, tempId1, area);
    EXPECT_EQ(Vertices[Vertices.size()-1].x, 1.5);
    EXPECT_EQ(Vertices[Vertices.size()-1].y, 0.75);
    EXPECT_EQ(Edges[Edges.size()-1].start, Vertices[Vertices.size()-1]);
    EXPECT_EQ(Edges[Edges.size()-1].finish, Vertices[Vertices.size()-2]);
    EXPECT_EQ(Triangles[Triangles.size()-1].edges[0].adjTriangles[0], Triangles[Triangles.size()-2].id);
    EXPECT_EQ(Triangles[Triangles.size()-2].edges[0].adjTriangles[0], Triangles[Triangles.size()-1].id);
    EXPECT_EQ(Triangles[Triangles.size()-3].edges[0].adjTriangles[0], Triangles[Triangles.size()-4].id);
    EXPECT_EQ(Triangles[Triangles.size()-4].edges[0].adjTriangles[0], Triangles[Triangles.size()-3].id);
}

#endif // __TEST_EMPTY_H
