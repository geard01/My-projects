#include "empty_class.hpp"
#include <chrono>
using namespace std;
using namespace ProjectLibrary;
// Il seguente programma raffina una mesh triangolare, sono da inserire 2 oppure 3 input sotto command line argument
// separati da uno spazio. Il primo valore può essere 1 oppure 2 ed è il dataset che si utilizzerà. Il secondo valore
// è 1 oppure 2 ed è la modalità di raffinamento: in base all'area oppure in base al gradiente di una funzione a scelta.
// Nel caso si utilizzi la modalità 2, è necessario inserire un numero che sarà il numero della funzione che utilizzeremo
// per il calcolo dell'integrale sotteso.
int main(int argc, char **argv)
{
    if (argc < 4)
    {
        cerr << "Test number and mode shall be passed to the program" << endl;
        return -1;
    }
    string test = argv[1];
    Triangle::mode = argv[2];
    Triangle::function = argv[3];
    string function1=argv[3];
    auto start = chrono::high_resolution_clock::now();
    vector<Vertex> ciro;
    vector<Edge> marco;
    vector<Triangle> cosimo;
    unsigned int n = stoul(argv[4]);
    ImportCell0Ds(ciro, n, test);
    ImportCell1Ds(marco, ciro, n, test);
    ImportCell2Ds(cosimo, marco, ciro, n, test);
    double exactarea= doubleIntegral(function1);
    double area=0.0;
    for (unsigned int i=0; i<cosimo.size(); i++)
    {
        area+=(cosimo[i].area*evaluateExpressionf0(cosimo[i].center.x,cosimo[i].center.y, cosimo[i].function));
    }
    Refine(cosimo, marco, ciro, n, area, exactarea);
    ofstream outputFile("C:/Users/Gentian/Downloads/outputPunti.csv");
    for (unsigned int i = 0; i < ciro.size(); i++)
    {
        outputFile << ciro[i].x << ";" << ciro[i].y << "\n";
    }
    outputFile.close();

    ofstream outputFile1("C:/Users/Gentian/Downloads/outputLati.csv");
    for (unsigned int i = 0; i < marco.size(); i++)
    {
        if (marco[i].active)
            outputFile1 << marco[i].start.x << ";" << marco[i].start.y << ";"<< marco[i].finish.x << ";" << marco[i].finish.y << "\n";
    }

    auto end = chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Tempo di esecuzione: " << duration.count() << "s" << endl;
    return 0;
}
