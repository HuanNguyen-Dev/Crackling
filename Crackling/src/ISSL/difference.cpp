/*

Faster and better CRISPR guide RNA design with the Crackling method.
Jacob Bradford, Timothy Chappell, Dimitri Perrin
bioRxiv 2020.02.14.950261; doi: https://doi.org/10.1101/2020.02.14.950261


To compile:

g++ -o difference difference.cpp 

*/

#include <iostream>
#include <unordered_map>
#include <string>
#include <cstdint>
#include <vector>
#include <fstream>
#include <cmath>
using namespace std;


struct MapperResult {
    std::string querySignature;
    double mitScore;
    double cfdScore;
};
int seqLength;

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "Usage: %s [file-1]] [file-2] [seq-length]\n", argv[0]);
        exit(1);
    }
    
    seqLength = atoi(argv[3]);


    // store signature -> {mit, cfd}
    std::unordered_map<std::string, std::pair<double,double>> result;

    // load file 1
    ifstream in1(argv[1]);
    if (!in1) {
        cerr << "Failed to open " << argv[1] << endl;
        return 1;
    }
    std::string seq;
    double mit,cfd;
    while (in1 >> seq >> mit >> cfd) {
        result[seq] = {mit, cfd};
    }

    in1.close();

    ifstream in2(argv[2]);

    if (!in2) {
        cerr << "Failed to open " << argv[2] << endl;
        return 1;
    };

        while (in2 >> seq >> mit >> cfd) {
            auto tmp = result.find(seq);
            if (tmp != result.end()){
                double mitDif = tmp->second.first - mit;
                double cfdDif = tmp->second.second - cfd;
                const double epsilon = 1e-4;
                if (std::abs(mitDif) > epsilon || std::abs(cfdDif) > epsilon)
                {cout << seq << " " << mitDif << " " << cfdDif << "\n";}
                else{
                    cout << seq << " " << "difference < 1e-4 for both\n" ;
                }

            }else{
                cout << "===== Missed ======\n";
                cout << seq << "\n";
                cout << "===================\n";
            }
        };
        in2.close();


    return 0;
}