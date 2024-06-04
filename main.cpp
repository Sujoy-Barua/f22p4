// Project Identifier: 5949F553E20B650AB0FB2266D3C0822B13D248B0

#include <getopt.h>
#include <iostream>
#include <string>
#include <vector>
#include <limits>
#include <cmath>
#include <iomanip>
#include <algorithm>

using namespace std;

enum class Loc : uint8_t { Coast, Land, Sea };

struct PokemonMST {
    int xCoOrd;
    int yCoOrd;
    Loc loc = Loc::Land;
};

struct PokemonTSP {
    int xCoOrd;
    int yCoOrd;
    bool visited = false;
};

struct primDatum {
    bool k = false;
    double dSq = numeric_limits<double>::infinity(); // square of the distance
    uint32_t p = 0;
};

struct coOrd {
    int x;
    int y;
};


class PlanetAB {
private:
    vector<PokemonMST> pokeListMST;
    vector<primDatum> primData;
    vector<PokemonTSP> pokeListTSP;
    uint32_t totVertex;

    void inputMST() {
        pokeListMST.reserve(totVertex);
        primData.resize(totVertex);
        int inputX;
        int inputY;
        Loc inputLoc;
        PokemonMST pokemon;
        while (cin >> inputX >> inputY) {
            if (inputX > 0 || inputY > 0) {
                inputLoc = Loc::Land;
            }
            else if (inputX == 0 || inputY == 0) {
                inputLoc = Loc::Coast;
            }
            else {
                inputLoc = Loc::Sea;
            }

            pokemon = { inputX, inputY, inputLoc };
            pokeListMST.push_back(pokemon);

        }
    }

    void inputTSP() {
        pokeListTSP.reserve(totVertex);
        //primData.resize(totVertex);
        int inputX;
        int inputY;
        PokemonTSP pokemon;
        while (cin >> inputX >> inputY) {

            pokemon = { inputX, inputY, false };
            pokeListTSP.push_back(pokemon);

        }
    }

    double distanceSqMST(uint32_t indexFrom, uint32_t indexTo) {
        double xdiff = numeric_limits<double>::infinity();
        double ydiff = numeric_limits<double>::infinity();

        if ((pokeListMST[indexFrom].loc == Loc::Land && pokeListMST[indexTo].loc == Loc::Sea)
            || (pokeListMST[indexFrom].loc == Loc::Sea && pokeListMST[indexTo].loc == Loc::Land)) {

        }
        else {
            xdiff = pokeListMST[indexFrom].xCoOrd - pokeListMST[indexTo].xCoOrd;
            ydiff = pokeListMST[indexFrom].yCoOrd - pokeListMST[indexTo].yCoOrd;
        }

        return (xdiff * xdiff) + (ydiff * ydiff);
    }

    double distanceTSP(uint32_t indexFrom, uint32_t indexTo) {
        double xdiff = pokeListTSP[indexFrom].xCoOrd - pokeListTSP[indexTo].xCoOrd;
        double ydiff = pokeListTSP[indexFrom].yCoOrd - pokeListTSP[indexTo].yCoOrd;

        return sqrt((xdiff * xdiff) + (ydiff * ydiff));
    }

    double distanceSqTSP(uint32_t indexFrom, uint32_t indexTo) {
        double xdiff = pokeListTSP[indexFrom].xCoOrd - pokeListTSP[indexTo].xCoOrd;
        double ydiff = pokeListTSP[indexFrom].yCoOrd - pokeListTSP[indexTo].yCoOrd;

        return (xdiff * xdiff) + (ydiff * ydiff);
    }

    double primAlg() {
        primData[0].dSq = 0;
        double totWeight = 0;
        for (uint32_t i = 0; i < totVertex; ++i) {
            uint32_t curInd = 0;
            double curDis = numeric_limits<double>::infinity();
            for (uint32_t j = 0; j < totVertex; ++j) {
                if (primData[j].k == false) {
                    if (primData[j].dSq < curDis) {
                        curDis = primData[j].dSq;
                        curInd = j;
                    }
                }
            }
            primData[curInd].k = true;
            totWeight += sqrt(primData[curInd].dSq);

            for (uint32_t j = 0; j < totVertex; ++j) {
                if (primData[j].k == false) {
                    double dSqNew = distanceSqMST(curInd, j);
                    if (dSqNew < primData[j].dSq) {
                        primData[j].dSq = dSqNew;
                        primData[j].p = curInd;
                    }
                }
            }
        }

        return totWeight;
    }

    double nNeighborHeu(vector<uint32_t>& travelList) {

        uint32_t minInd = 0;
        uint32_t curInd = 0;
        double totWeight = 0;
        double curWeight = 0;
        pokeListTSP[curInd].visited = true;
        travelList.push_back(0);

        for (uint32_t i = 0; i < pokeListTSP.size() - 1; ++i) {
            double prevDis = numeric_limits<double>::infinity();
            double curDis = 0.0;
            for (uint32_t j = 0; j < pokeListTSP.size(); ++j) {
                if (pokeListTSP[j].visited == false) {
                    curDis = distanceSqTSP(curInd, j);
                    if (curDis < prevDis) {
                        prevDis = curDis;
                        minInd = j;
                        curWeight = curDis;
                    }
                }
            }
            pokeListTSP[minInd].visited = true;
            travelList.push_back(minInd);
            curInd = minInd;
            totWeight += sqrt(curWeight);
        }
        travelList.push_back(0);
        return totWeight += distanceTSP(curInd, 0);
    }

    double twoOptHeu(vector<uint32_t>& travelList) {
        double totWeight = nNeighborHeu(travelList);
        for (uint32_t i = 1; i < travelList.size() - 2; i++) {
            for (uint32_t j = i + 1; j < travelList.size() - 1; j++) {
                double d1 = distanceTSP(travelList[i - 1], travelList[i]) + distanceTSP(travelList[j], travelList[j + 1]);
                double d2 = distanceTSP(travelList[i - 1], travelList[j]) + distanceTSP(travelList[i], travelList[j + 1]);

                if (d1 > d2) {
                    reverse(travelList.begin() + i, travelList.begin() + j + 1);
                    totWeight -= (d1 - d2);
                    j = i + 1;
                }
            }
        }

        return totWeight;
    }


public:
    PlanetAB(uint32_t totVertex_In) : totVertex(totVertex_In) {}

    void mst() {
        inputMST();
        double totWeight = primAlg();
        cout << totWeight << "\n";
        for (uint32_t i = 1; i < primData.size(); ++i) {
            if (primData[i].p > i) {
                cout << i << " " << primData[i].p << "\n";
            }
            else {
                cout << primData[i].p << " " << i << "\n";
            }
        }
    }

    void fasttsp() {
        inputTSP();
        vector<uint32_t> travelList;
        double totWeight = 0.0;
        totWeight = twoOptHeu(travelList);

        cout << totWeight << "\n";
        for (uint32_t i = 0; i < travelList.size() - 1; ++i) {
            cout << travelList[i] << " ";
        }
        cout << "\n";
    }
};


class PlanetC {
private:
    vector<vector<double>> distanceMatrix;
    vector<primDatum> primData;
    vector<PokemonTSP> pokeListTSP;

    vector<uint32_t> path;
    vector<uint32_t> bestPath;
    double bestPathLength = 0.0;
    double curPathLength = 0.0;

    uint32_t totVertex;

    void inputTSP() {
        pokeListTSP.reserve(totVertex);
        //primData.resize(totVertex);
        int inputX;
        int inputY;
        PokemonTSP pokemon;
        while (cin >> inputX >> inputY) {
            pokemon = { inputX, inputY, false };
            pokeListTSP.push_back(pokemon);

        }

        initiateMatrix();
    }

    void initiateMatrix() {

        for (uint32_t i = 0; i < totVertex; i++) {
            for (uint32_t j = 0; j < totVertex; j++) {
                distanceMatrix[i][j] = distanceSqTSP(i, j);
            }
        }
    }

    double distanceTSP(uint32_t indexFrom, uint32_t indexTo) { // returns square root of distance
        double xdiff = pokeListTSP[indexFrom].xCoOrd - pokeListTSP[indexTo].xCoOrd;
        double ydiff = pokeListTSP[indexFrom].yCoOrd - pokeListTSP[indexTo].yCoOrd;

        return sqrt((xdiff * xdiff) + (ydiff * ydiff));
    }

    double distanceSqTSP(uint32_t indexFrom, uint32_t indexTo) { // returns squared distance
        double xdiff = pokeListTSP[indexFrom].xCoOrd - pokeListTSP[indexTo].xCoOrd;
        double ydiff = pokeListTSP[indexFrom].yCoOrd - pokeListTSP[indexTo].yCoOrd;

        return (xdiff * xdiff) + (ydiff * ydiff);
    }

    double distanceSqCoOrd(const coOrd& coOrd1, const coOrd& coOrd2) { // returns squared distance
        double xdiff = coOrd1.x - coOrd2.x;
        double ydiff = coOrd1.y - coOrd2.y;

        return (xdiff * xdiff) + (ydiff * ydiff);
    }

    void nNeighborHeu() {

        uint32_t minInd = 0;
        uint32_t curInd = 0;
        double curWeight = 0.0;
        pokeListTSP[curInd].visited = true;
        bestPath.push_back(0);

        for (uint32_t i = 0; i < pokeListTSP.size() - 1; ++i) {
            double prevDis = numeric_limits<double>::infinity();
            double curDis = 0.0;
            for (uint32_t j = 0; j < pokeListTSP.size(); ++j) {
                if (pokeListTSP[j].visited == false) {
                    curDis = distanceMatrix[curInd][j];
                    if (curDis < prevDis) {
                        prevDis = curDis;
                        minInd = j;
                        curWeight = curDis;
                    }
                }
            }
            pokeListTSP[minInd].visited = true;
            bestPath.push_back(minInd);
            curInd = minInd;
            bestPathLength += sqrt(curWeight);
        }
        bestPath.push_back(0);
        bestPathLength += sqrt(distanceMatrix[curInd][0]);
        twoOptHeu();
    }

    void twoOptHeu() {

        for (uint32_t i = 1; i < bestPath.size() - 2; i++) {
            for (uint32_t j = i + 1; j < bestPath.size() - 1; j++) {
                double d1 = sqrt(distanceMatrix[bestPath[i - 1]][bestPath[i]]) + sqrt(distanceMatrix[bestPath[j]][bestPath[j + 1]]);
                double d2 = sqrt(distanceMatrix[bestPath[i - 1]][bestPath[j]]) + sqrt(distanceMatrix[bestPath[i]][bestPath[j + 1]]);

                if (d1 > d2) {
                    reverse(bestPath.begin() + i, bestPath.begin() + j + 1);
                    bestPathLength -= (d1 - d2);
                    j = i + 1;
                }
            }
        }
    }

    void genPerms(uint32_t permLength) {

        if (permLength == path.size()) {
            // add update() code to run if better cycle
            curPathLength += sqrt(distanceMatrix[path[path.size() - 1]][path[0]]);
            if (curPathLength < bestPathLength) {
                bestPathLength = curPathLength;
                bestPath = path;
            }
            curPathLength -= sqrt(distanceMatrix[path[path.size() - 1]][path[0]]);
            return;
        }  // if ..complete path

        if (!promising(permLength)) {
            return;
        }  // if ..not promising

        for (uint32_t i = permLength; i < path.size(); ++i) {
            double curEdge = sqrt(distanceMatrix[path[permLength - 1]][path[i]]);
            swap(path[permLength], path[i]);
            curPathLength += curEdge;
            genPerms(permLength + 1);
            curPathLength -= curEdge;
            swap(path[permLength], path[i]);
        }  // for ..unpermuted elements
    }  // genPerms()

    bool promising(uint32_t permLength) {

        if ((path.size() - permLength) < 5) {
            return true;
        } // if

        // do the MST code for the unconnected portion of the graph

        vector<coOrd> coOrds;
        coOrds.reserve(path.size() - permLength);

        for (uint32_t i = permLength; i < path.size(); ++i) {
            coOrds.push_back({ pokeListTSP[path[i]].xCoOrd, pokeListTSP[path[i]].yCoOrd });
        }

        double lb = primAlgMod(coOrds);
        // look for the lowest two lengths to connect MST to the fixed graph
        uint32_t i = 0;
        while (i < permLength) {
            double minDist = numeric_limits<double>::infinity();
            for (uint32_t j = permLength; j < path.size(); ++j) {
                double dist = distanceMatrix[path[i]][path[j]];
                if (dist < minDist) {
                    minDist = dist;
                }
            }

            if (permLength > 1) {
                lb += sqrt(minDist);
                i += (permLength - 1);
            }
            else {
                lb += 2 * sqrt(minDist);
                i += 1;
            }
        }
        lb += curPathLength;

        // compute the lower bound for the remaining vertices
        if (lb < bestPathLength) {
            return true;
        }
        else {
            return false;
        }

    }

    double primAlgMod(const vector<coOrd>& coOrds) {
        primData.clear();
        primData.resize(coOrds.size());
        primData[0].dSq = 0.0;
        double totWeight = 0.0;
        for (uint32_t i = 0; i < coOrds.size(); ++i) {
            uint32_t curInd = 0;
            double curDis = numeric_limits<double>::infinity();
            for (uint32_t j = 0; j < coOrds.size(); ++j) {
                if (primData[j].k == false) {
                    if (primData[j].dSq < curDis) {
                        curDis = primData[j].dSq;
                        curInd = j;
                    }
                }
            }
            primData[curInd].k = true;
            totWeight += sqrt(primData[curInd].dSq);

            for (uint32_t j = 0; j < coOrds.size(); ++j) {
                if (primData[j].k == false) {
                    double dSqNew = distanceSqCoOrd(coOrds[curInd], coOrds[j]);
                    if (dSqNew < primData[j].dSq) {
                        primData[j].dSq = dSqNew;
                        primData[j].p = curInd;
                    }
                }
            }
        }

        return totWeight;
    }

public:
    PlanetC(uint32_t totVertex_In) : totVertex(totVertex_In) {
        distanceMatrix.resize(totVertex);
        for (uint32_t i = 0; i < totVertex; i++) {
            distanceMatrix[i].resize(totVertex);
        }
    }

    void opttsp() {
        inputTSP();
        bestPath.reserve(totVertex);
        path.reserve(totVertex);
        nNeighborHeu();
        bestPath.pop_back();
        path.clear();
        path = bestPath;
        genPerms(1);
        cout << bestPathLength << "\n";
        for (uint32_t i = 0; i < bestPath.size(); ++i) {
            cout << bestPath[i] << " ";
        }
        cout << "\n";
    }
};

void getMode(int argc, char* argv[]) {
    // These are used with getopt_long()
    opterr = false; // Let us handle all error output for command line options
    int choice;
    int index = 0;
    struct option long_options[] = {
        // Fill in two lines, for the "mode" ('m') and
        // the "help" ('h') options.
        // ./project0 -m nosize
        // ./project0 --help
        { "mode", required_argument, nullptr, 'm'  },
        { "help", no_argument,       nullptr, 'h'  },
        { nullptr,0,                 nullptr, '\0' },
    };  // long_options[]

    // Fill in the double quotes, to match the mode and help options.
    while ((choice = getopt_long(argc, argv, "hm:", long_options, &index)) != -1) {
        switch (choice) {

        case 'h':
            // help output
            break;

        case 'm':
            string arg{ optarg };
            // deal with the specified mode
            if (arg == "MST") {
                uint32_t totVertex;
                cin >> totVertex;
                PlanetAB planet(totVertex);
                planet.mst();
            }
            else if (arg == "FASTTSP") {
                uint32_t totVertex;
                cin >> totVertex;
                PlanetAB planet(totVertex);
                planet.fasttsp();
            }
            else if (arg == "OPTTSP") {
                uint32_t totVertex;
                cin >> totVertex;
                PlanetC planet(totVertex);
                planet.opttsp();
            }
            else {
                cerr << "Error \n";
                exit(1);
            }
            break;
        }  // switch ..choice
    }  // while
}  // getMode()

int main(int argc, char* argv[]) {
    ios::sync_with_stdio(false);
    cout << fixed << setprecision(2);
    getMode(argc, argv);
 }