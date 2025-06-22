#include "empty_class.hpp"
#include <chrono>
using namespace std;
using namespace ProjectLibrary;
int main(int argc, char** argv)
{
    if (argc < 2)
    {
        cerr<< "Test number shall be passed to the program"<< endl;
        return -1;
    }
    string test = argv[1];
    auto start = chrono::high_resolution_clock::now();
    vector<Vertex> ciro;
    vector<Edge> marco;
    vector<Triangle> cosimo;
    if (argc < 3)
    {
        cerr<< "Number of iterations shall be passed to the program"<< endl;
        return -1;
    }
    unsigned int n= stoul(argv[2]);
    ImportCell0Ds(ciro, n, test);
    ImportCell1Ds(marco, ciro, n, test);
    ImportCell2Ds(cosimo, marco, ciro, n, test);
    ProjectLibrary::Refine(cosimo, marco, ciro, n);
    ofstream outputFile("C:/Users/geard/Downloads/outputPunti.csv");
    for (unsigned int i=0; i<ciro.size(); i++)
    {
        outputFile<<ciro[i].x<<";"<<ciro[i].y<<"\n";
    }
    outputFile.close();

    ofstream outputFile1("C:/Users/geard/Downloads/outputLati.csv");
    for (unsigned int i=0; i<marco.size(); i++)
    {
        if (marco[i].active)
            outputFile1<<marco[i].start.x<<";"<<marco[i].start.y<<";"<<marco[i].finish.x<<";"<<marco[i].finish.y<<"\n";
    }

    auto end = std::chrono::high_resolution_clock::now();
    chrono::duration<double> duration = end - start;
    cout << "Tempo di esecuzione: " << duration.count() << "s" << endl;
    cout << "Numero di iterazioni: " << n << endl;
    return 0;
}
