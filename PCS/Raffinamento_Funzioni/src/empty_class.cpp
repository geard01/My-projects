#include "empty_class.hpp"
#include <queue>

namespace ProjectLibrary
{

string ProjectLibrary::Triangle::mode;
string ProjectLibrary::Triangle::function;

void split2(vector<Triangle> &triangles, vector<Edge> &edges, vector<Vertex> &vertices, unsigned int m, deque<unsigned int> &tempId,
            unsigned int &k, bool &permissible, deque<unsigned int> &tempId1, double &area)
{
    insertionSort(triangles[m].edges); // Sistemo i lati in ordine decrescente
    tempId.push_back(triangles[m].edges[0].id); // Aggiungo nuovo lato lungo
    Vertex opposite =getOppositeVertex(triangles[m], triangles[m].edges[0]); // Lato opposto
    Vertex mid = set_mid(triangles[m].edges[0]);
    mid.id = vertices.size(); // Assegno nuovo id, +1 rispetto all'ultimo
    Pushback(vertices, mid);
    Edge bisection(opposite, mid,edges.size()); // Trovo il lato bisezione, gli id sono progressivi
    Edge newEdge1(triangles[m].edges[0].start, mid,edges.size() + 1); // Trovo primo lato nuovo, id progressivo
    Edge newEdge2(mid, triangles[m].edges[0].finish,edges.size() + 2); // Trovo secondo lato nuovo
    vector<Edge> newEdgesV1;
    vector<Edge> newEdgesV2;
    // nel modo seguente, creo i nuovi triangoli in modo corretto, assegnando il
    // lato giusto ai nuovi triangoli creati il flag risparmia tempo, per evitare
    // di ricontrollare
    bool flag = false;
    if (triangles[m].edges[1].start == triangles[m].edges[0].start ||
        triangles[m].edges[1].finish ==triangles[m].edges[0].start) // Lato sinistro con edge1 o con edge2
    {
        newEdgesV1 = {bisection, newEdge1,triangles[m].edges[1]}; // Creo nuovi triangoli con i lati giusti
        newEdgesV2 = {bisection, newEdge2, triangles[m].edges[2]};
        flag = true;
    }
    else
    {
        newEdgesV1 = {bisection, newEdge1,triangles[m].edges[2]}; // Creo nuovi triangoli con i lati giusti
        newEdgesV2 = {bisection, newEdge2, triangles[m].edges[1]};
    }
    Triangle newTriangle1(newEdgesV1, triangles.size());
    Triangle newTriangle2(newEdgesV2, triangles.size() + 1);
    // condizione di stop : se un triangolo ha il lato più lungo al bordo, raffino
    // una volta e mi fermo. La mesh è ammissibile
    if (triangles[m].edges[0].adjTriangles.size() == 0)
        permissible = true;
    else
        k = triangles[m].edges[0].adjTriangles[0];
    Pushback(edges, bisection);
    Pushback(edges, newEdge1);
    Pushback(edges, newEdge2);
    triangles[m].active =false; // Disattivo il triangolo dal vettore di triangolo
    // aggiorno l'adiacenza dei triangoli adiacenti (tranne a quelli adiacenti al
    // lato più lungo)
    if (triangles[m].edges[1].adjTriangles.size() > 0)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (triangles[triangles[m].edges[1].adjTriangles[0]].edges[j].adjTriangles.size() > 0 &&
                triangles[triangles[m].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] == triangles[m].id)
            {
                if (flag)
                {
                    triangles[triangles[m].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle1.id;
                    break;
                }
                else
                {
                    triangles[triangles[m].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle2.id;
                    break;
                }
            }
        }
    }
    if (triangles[m].edges[2].adjTriangles.size() > 0)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (triangles[triangles[m].edges[2].adjTriangles[0]].edges[j].adjTriangles.size() > 0 &&
                triangles[triangles[m].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] == triangles[m].id)
            {
                if (flag)
                {
                    triangles[triangles[m].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle2.id;
                    break;
                }
                else
                {
                    triangles[triangles[m].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle1.id;
                    break;
                }
            }
        }
    }
    // questi aggiornano le adiacenze del nuovo lato creato
    Pushback(newTriangle1.edges[0].adjTriangles, newTriangle2.id);
    Pushback(newTriangle2.edges[0].adjTriangles, newTriangle1.id);
    // inserisco i nuovi triangoli nel vettore
    Pushback(triangles, newTriangle1);
    Pushback(triangles, newTriangle2);
    // aggiorno l'id dei due nuovi triangoli creati
    tempId1.push_back(newTriangle1.id);
    tempId1.push_back(newTriangle2.id);
    area-=(triangles[m].area*evaluateExpressionf0(triangles[m].center.x,triangles[m].center.y, triangles[m].function));
    area+=(newTriangle1.area*evaluateExpressionf0(newTriangle1.center.x,newTriangle1.center.y, newTriangle1.function));
    area+=(newTriangle2.area*evaluateExpressionf0(newTriangle2.center.x,newTriangle2.center.y, newTriangle2.function));
}

void split2again(vector<Triangle> &triangles, vector<Edge> &edges,
                 vector<Vertex> &vertices, unsigned int &k, bool &permissible,
                 deque<unsigned int> &tempId, deque<unsigned int> &tempId1, double &area)
{
    Vertex newOpposite = getOppositeVertex(triangles[k], triangles[k].edges[0]);
    Edge newEdge(newOpposite, vertices[vertices.size() - 1],
                 edges.size()); // Creo nuovo lato
    vector<Edge> EdgesV1;
    vector<Edge> EdgesV2;
    bool flag = false;
    if (triangles[k].edges[1].start == triangles[k].edges[0].start ||triangles[k].edges[1].finish ==triangles[k].edges[0].start) // Lato sinistro con edge1 o con edge2
    {
        EdgesV1 = {newEdge, triangles[tempId1.back()].edges[1],triangles[k].edges[2]}; // Creo nuovi triangoli con i lati giusti
        EdgesV2 = {newEdge, triangles[tempId1.front()].edges[1],triangles[k].edges[1]};
        flag = true;
    }
    else
    {
        EdgesV1 = {newEdge, triangles[tempId1.back()].edges[1],triangles[k].edges[1]}; // Creo nuovi triangoli con i lati giusti
        EdgesV2 = {newEdge, triangles[tempId1.front()].edges[1],triangles[k].edges[2]};
    }
    Triangle newTriangle(EdgesV1, triangles.size());
    Triangle newTriangle0(EdgesV2, triangles.size() + 1);
    Pushback(edges, newEdge);
    triangles[k].active = false; // Cancello il vecchio triangolo grande
    // aggiorno l'adiacenza dei triangoli adiacenti (tranne a quelli adiacenti al
    // lato più lungo)
    if (triangles[k].edges[1].adjTriangles.size() > 0)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (triangles[triangles[k].edges[1].adjTriangles[0]].edges[j].adjTriangles.size() > 0 &&
                triangles[triangles[k].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] == triangles[k].id)
            {
                if (flag)
                {
                    triangles[triangles[k].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle0.id;
                    break;
                }
                else
                {
                    triangles[triangles[k].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle.id;
                    break;
                }
            }
        }
    }
    if (triangles[k].edges[2].adjTriangles.size() > 0)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (triangles[triangles[k].edges[2].adjTriangles[0]].edges[j].adjTriangles.size() > 0 &&
                triangles[triangles[k].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] == triangles[k].id)
            {
                if (flag)
                {
                    triangles[triangles[k].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle.id;
                    break;
                }
                else
                {
                    triangles[triangles[k].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle0.id;
                    break;
                }
            }
        }
    }

    // Sistemo le adiacenze del nuovo lato creato, come in split2
    Pushback(newTriangle.edges[0].adjTriangles, newTriangle0.id);
    Pushback(newTriangle0.edges[0].adjTriangles, newTriangle.id);
    // riaggiorno le adiacenze dei lati dei nuovi triangoli creati rispetto a lati
    // spezzati in split2 (penultimo triangolo spezzato)
    Pushback(newTriangle.edges[1].adjTriangles, triangles[tempId1.back()].id);
    Pushback(triangles[tempId1.back()].edges[1].adjTriangles, newTriangle.id);
    Pushback(newTriangle0.edges[1].adjTriangles, triangles[tempId1.front()].id);
    Pushback(triangles[tempId1.front()].edges[1].adjTriangles, newTriangle0.id);
    // Inserisco i triangoli nella lista
    Pushback(triangles, newTriangle);
    Pushback(triangles, newTriangle0);
    edges[tempId.front()].active = false;
    permissible = true;
    area-=(triangles[k].area*evaluateExpressionf0(triangles[k].center.x,triangles[k].center.y, triangles[k].function));
    area+=(newTriangle.area*evaluateExpressionf0(newTriangle.center.x,newTriangle.center.y, newTriangle.function));
    area+=(newTriangle0.area*evaluateExpressionf0(newTriangle0.center.x,newTriangle0.center.y, newTriangle0.function));

}
void split3(vector<Triangle> &triangles, vector<Edge> &edges, vector<Vertex> &vertices, deque<unsigned int> &tempId,
            unsigned int &k, bool &permissible, deque<unsigned int> &tempId1, double &area)
{
    split2(triangles, edges, vertices, k, tempId, k, permissible,tempId1,area); // Divido prima il triangolo adiacente in 2
    Edge newEdge3(vertices[vertices.size() - 1], vertices[vertices.size() - 2],edges.size()); // Creo lato tra i 2 punti medi
    Pushback(edges, newEdge3);
    int z = triangles.size() - 2;
    if (triangles[z + 1].edges[2].id == tempId.front()) //Trovo quale triangolo da dividere è quello tra i 2 appena creati
        z++;
    vector<Edge> EdgesV3;
    vector<Edge> EdgesV4;
    for (unsigned int j = 1; j < 3; j++) //Sistemo il lato più lungo in 0
    {
        if (triangles[z].edges[j].id == tempId.front())
        {
            Edge temp = triangles[z].edges[j];
            triangles[z].edges[j] = triangles[z].edges[0];
            triangles[z].edges[0] = temp;
            break;
        }
    }
    bool flag = false;
    if (triangles[z].edges[1].start == triangles[z].edges[0].start ||
        triangles[z].edges[1].finish ==triangles[z].edges[0].start) // solito controllo sui lati
    {
        EdgesV3 = {newEdge3, triangles[z].edges[1],triangles[tempId1.front()].edges[1]}; // Creo nuovi triangoli con i lati giusti
        EdgesV4 = {newEdge3, triangles[z].edges[2], triangles[tempId1[1]].edges[1]};
        flag = true;
    }
    else
    {
        EdgesV3 = {newEdge3, triangles[z].edges[2],triangles[tempId1.front()].edges[1]}; // Creo nuovi triangoli con i lati giusti
        EdgesV4 = {newEdge3, triangles[z].edges[1], triangles[tempId1[1]].edges[1]};
    }
    triangles[z].active = false;
    Triangle newTriangle3(EdgesV3, triangles.size());
    Triangle newTriangle4(EdgesV4, triangles.size() + 1);
    if (triangles[z].edges[1].adjTriangles.size() > 0) //Solito controllo adiacenze
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (triangles[triangles[z].edges[1].adjTriangles[0]].edges[j].adjTriangles.size() > 0 &&
                triangles[triangles[z].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] == triangles[z].id)
            {
                if (flag)
                {
                    triangles[triangles[z].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle3.id;
                    break;
                }
                else
                {
                    triangles[triangles[z].edges[1].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle4.id;
                    break;
                }
            }
        }
    }
    if (triangles[z].edges[2].adjTriangles.size() > 0)
    {
        for (unsigned int j = 0; j < 3; j++)
        {
            if (triangles[triangles[z].edges[2].adjTriangles[0]].edges[j].adjTriangles.size() > 0 &&
                triangles[triangles[z].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] == triangles[z].id)
            {
                if (flag)
                {
                    triangles[triangles[z].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle4.id;
                    break;
                }
                else
                {
                    triangles[triangles[z].edges[2].adjTriangles[0]].edges[j].adjTriangles[0] = newTriangle3.id;
                    break;
                }
            }
        }
    }
    Pushback(newTriangle3.edges[2].adjTriangles, triangles[tempId1.front()].id);
    Pushback(triangles[tempId1.front()].edges[1].adjTriangles, newTriangle3.id);
    tempId1.pop_front();
    Pushback(newTriangle4.edges[2].adjTriangles, triangles[tempId1.front()].id);
    Pushback(triangles[tempId1.front()].edges[1].adjTriangles, newTriangle4.id);
    Pushback(newTriangle3.edges[0].adjTriangles,newTriangle4.id); // Sistemo le adiacenze dei lati nuovi
    Pushback(newTriangle4.edges[0].adjTriangles, newTriangle3.id);
    tempId1.pop_front();
    bool flag1 = true;
    if (triangles[z].edges[1] == triangles[tempId1.front()].edges[1]) //Sistemo i nuovi triangoli creati dentro il mio tempId
        flag1 = false;
    if (newTriangle3.edges[1] == triangles[z].edges[1]) {
        if (flag1)
            tempId1.back() = newTriangle3.id;
        else
            tempId1.front() = newTriangle3.id;
    }
    else
    {
        if (flag1)
            tempId1.back() = newTriangle4.id;
        else
            tempId1.front() = newTriangle4.id;
    }
    Pushback(triangles, newTriangle3);
    Pushback(triangles, newTriangle4);
    edges[tempId.front()].active = false;
    tempId.pop_front(); // Mi dimentico il vecchio lato lungo e vado avanti
    area-=(triangles[z].area*evaluateExpressionf0(triangles[z].center.x,triangles[z].center.y, triangles[z].function));
    area+=(newTriangle3.area*evaluateExpressionf0(newTriangle3.center.x,newTriangle3.center.y, newTriangle3.function));
    area+=(newTriangle4.area*evaluateExpressionf0(newTriangle4.center.x,newTriangle4.center.y, newTriangle4.function));
}
// raffina tutto l'array di triangoli, partendo da un vettore ordinato in base
// all'area
void Refine(vector<Triangle> &triangles, vector<Edge> &edges, vector<Vertex> &vertices, unsigned int &n, double &area, double &exactarea)
{
    static constexpr double geometricTol = 1.0e-5;
    unsigned int h=0;
    while((area+geometricTol)<exactarea && h<n)
    {
        h++;
        bool permissible = false;
        unsigned int m = massimoElementoAttivo(triangles);
        deque<unsigned int> tempId;  // Salvo gli ultimi 2 lati più lunghi
        deque<unsigned int> tempId1; // Salvo gli id dei nuovi lati spezzati
        unsigned int k = 0;          // K-esimo triangolo da raffinare
        split2(triangles, edges, vertices, m, tempId, k, permissible, tempId1, area);
        while (!permissible)
        {
            // l'insertionSort è l'algoritmo più efficiente nel caso di vettori di
            // dimensioni molto piccole, come questo
            insertionSort(triangles[k].edges);
            if (triangles[k].edges[0].id == tempId.back())
            {
                split2again(triangles, edges, vertices, k, permissible, tempId,tempId1, area);
                permissible = true;
            }
            else
            {
                split3(triangles, edges, vertices, tempId, k, permissible, tempId1, area);
            }
        }
    }
    cout<<h<<endl;
}
// assegno a ogni triangolo le sue proprietà chiamando le seguenti tre funzioni
bool ImportCell0Ds(vector<Vertex> &vertices, unsigned int n, string &test)
{
    ifstream file;
    string inFile = "./Dataset/Test" + test + "/Cell0Ds.csv";
    file.open(inFile);

    if (file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);
    file.close();
    listLines.pop_front();
    vertices.reserve(listLines.size() * n);
    vertices.resize(listLines.size());
    if (vertices.size() == 0)
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }
    for (const string &line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        Vector2d coord;
        converter >> id >> marker >> coord(0) >> coord(1);
        vertices[id].id = id;
        vertices[id].x = coord(0);
        vertices[id].y = coord(1);
    }
    file.close();
    return true;
}
bool ImportCell1Ds(vector<Edge> &edges, vector<Vertex> &vertices, unsigned int n, string &test)
{
    ifstream file;
    string inFile = "./Dataset/Test" + test + "/Cell1Ds.csv";
    file.open(inFile);

    if (file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);
    listLines.pop_front();
    edges.reserve(2 * listLines.size() * n);
    edges.resize(listLines.size());
    if (edges.size() == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }
    for (const string &line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        unsigned int marker;
        unsigned int start;
        unsigned int finish;

        converter >> id >> marker >> start >> finish;

        edges[id].id = id;
        edges[id].start = vertices[start];
        edges[id].finish = vertices[finish];
        edges[id].set_length();
    }

    file.close();
    return true;
}
bool ImportCell2Ds(vector<Triangle> &triangles, vector<Edge> &edges, vector<Vertex> &vertices, unsigned int n, string &test)
{

    ifstream file;
    string inFile = "./Dataset/Test" + test + "/Cell2Ds.csv";
    file.open(inFile);

    if (file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);
    listLines.pop_front();
    triangles.reserve(3 * listLines.size() * n);
    triangles.resize(listLines.size());
    if (triangles.size() == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }
    vector<int> tempEdges;
    tempEdges.resize(3 * listLines.size());
    unsigned int z = 0;
    for (const string &line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        array<unsigned int, 3> tempVertices;

        converter >> id;
        for (unsigned int i = 0; i < 3; i++)
        {
            converter >> tempVertices[i];
        }
        for (unsigned int i = z; i < (z + 3); i++)
        {
            converter >> tempEdges[i];
        }

        triangles[id].id = id;
        Pushback(triangles[id].vertices, vertices[tempVertices[0]]);
        Pushback(triangles[id].vertices, vertices[tempVertices[1]]);
        Pushback(triangles[id].vertices, vertices[tempVertices[2]]);
        // assegno a ogni lato un puntatore al triangolo che lo costituisce
        for (unsigned int j = z; j < (z + 3); j++)
        {
            edges[tempEdges[j]].adjTriangles.reserve(3);
            Pushback(edges[tempEdges[j]].adjTriangles, triangles[id].id);
        }
        triangles[id].set_area(); // assegno l'area al triangolo
        triangles[id].set_grad(); // assegno il modulo del gradiente valutato nel baricentro del triangolo a ogni triangolo
        z = z + 3;
    }
    z = 0;
    for (unsigned int i = 0; i < listLines.size(); i++)
    {
        for (unsigned int j = z; j < (z + 3); j++)
        {
            Pushback(triangles[i].edges, edges[tempEdges[j]]);
            if (triangles[i].edges[j - z].adjTriangles[0] == triangles[i].id)
                Erase(triangles[i].edges[j - z].adjTriangles,0); // Se no levo l'altro// Levo il primo triangolo a cui punto
            else if (triangles[i].edges[j - z].adjTriangles[1] == triangles[i].id)
                Erase(triangles[i].edges[j - z].adjTriangles,triangles[i].edges[j - z].adjTriangles.size() -1); // Se no levo l'altro
        }
        z = z + 3;
    }
    file.close();
    return true;
}
}
