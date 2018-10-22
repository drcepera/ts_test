#include <iostream>
#include <vector>
#include <unordered_map>
#include <sstream>
#include <iterator>
#include <algorithm>
#include <cassert>
#include <memory>
#include <list>
#include <cstdlib>
#include <fstream>
#include <chrono>

using namespace std;

#define MAX_N 300

const string filename = "../datasets/3.in";

//#define DEBUG_OUT
#ifdef DEBUG_OUT
    #define DEBUG(x) do { cout << x << endl; } while (0)
#else
    #define DEBUG(x) ;
#endif

int N_cities;
int N_zone;

typedef uint16_t cost_t;

struct city {
    int n;
    int zone;
    string name;
    
    bool operator <(const city& other) { return this->n < other.n; }
};

struct zone {
    int n;
    string name;
    int firstCity;
    int lastCity;
};

vector<zone> zones;
vector<city> cities;
unordered_map<string, int> citiesName2Num;
int startZone;
int startCity;

// flight table [day][departure][arriving] of costs
// dimension: Zones x Cities x Cities
typedef vector<vector<vector<cost_t >>> FlightTable;
unique_ptr<FlightTable> flightTable;

struct node {
    int day;  // from 0 to N_zones - 1
    int city;
    node* nextNode;
    int sum; // cost
    
    node(int d, int c) : day(d), city(c), nextNode(nullptr), sum(0) {}
    
    ~node() {
        if( nextNode )
            delete nextNode;
    }
    
    void dumpAnswer() {
        cout << sum << endl;
        
        node* n = this;
        while( n->nextNode ) {
            cout << cities[n->city].name << " "
                 << cities[n->nextNode->city].name << " "
                 << n->day + 1 << " "
                 << flightTable->at(n->day)[n->city][n->nextNode->city] << endl;
    
            n = n->nextNode;
        }
    }
};

// full search preserving zone arrange
bool tryFoundPath_greedy(node* start, vector<int>::iterator path) {
    // vector of next possible cities from path iterator
    struct todayFlight {
        int nextCity;
        cost_t cost;
    };
    
    vector<todayFlight> nextCities(zones[*path].lastCity - zones[*path].firstCity + 1);
    for( int i = 0; i < nextCities.size(); i++) {
        nextCities[i] = todayFlight{zones[*path].firstCity + i, flightTable->at(start->day)[start->city][zones[*path].firstCity + i]};
    }
    
    // thats why little bit greedy
    sort(nextCities.begin(), nextCities.end(), [](todayFlight& f1, todayFlight& f2){
       return (( f2.cost == 0 ) || (f1.cost < f2.cost));
    });
    
    for( auto &f : nextCities ) {
        if( !f.cost )
            break;
    
        start->nextNode = new node(start->day + 1, f.nextCity);
    
        if( *path == startZone ) { // start zone achieved
            start->sum += f.cost;
            return true;
        }
        
        vector<int>::iterator nextZone(path);
        if( tryFoundPath_greedy(start->nextNode, ++nextZone) ) {
            start->sum = f.cost + start->nextNode->sum;
            return true;    // success, path found!
        }
        else {
            delete start->nextNode;
            start->nextNode = nullptr;
            continue;
        }
    }
    
//    if( start->day > 5)
//        cout << "break on day " << start->day << "\n";
    
    return false;
}

void dumpFlightDayTable(const vector<vector<cost_t>> &costs) {
    // out table head
    cout << "\t";
    for( auto c : cities )
        cout << c.name << "\t";
    cout << endl;
    
    // rows
    for( auto c_dep : cities ) {
        cout << c_dep.name << "\t";
        for( auto cost : costs[c_dep.n] )
            cout << cost << "\t";
        cout << endl;
    }
}

void init()
{
    istream *input = filename.empty() ? &cin : new ifstream(filename);

    string startCityName;
    *input >> N_zone >> startCityName;
    
    DEBUG("num of zones: " << N_zone << ", start city: " << startCityName);
    
    // read cities in zones
    for( int i=0; i< N_zone; i++ ) {
        zone newZone{i, "", N_cities, N_cities-1};
        newZone.n = i;
        *input >> newZone.name;
        input->ignore();
        
        string citiesLine;
        getline(*input, citiesLine);
        istringstream iss(citiesLine);
        vector<string> cityNames{istream_iterator<string>{iss}, istream_iterator<string>{}};
        
        for( auto cityName : cityNames ) {
            city newCity({N_cities++, i, cityName});
            cities.push_back(newCity);
            DEBUG("city " << cityName << " number: " << cities[newCity.n].n << " zone: " << cities[newCity.n].zone);

            newZone.lastCity++;

            if( cityName == startCityName ) {
                startCity = newCity.n;
                startZone = newCity.zone;
                DEBUG("start city found: " << cityName << " number: " << cities[newCity.n].n << " zone: " << cities[newCity.n].zone);
            }
    
            citiesName2Num[newCity.name] = newCity.n;
        }

        zones.push_back(newZone);

        DEBUG("zone " << zones[i].name << " #" << zones[i].n <<
              " firstCity: " << zones[i].firstCity << " lastCity: " << zones[i].lastCity);
    }
    
    // create flight table
    flightTable.reset(new vector<vector<vector<cost_t >>>(N_zone, vector<vector<cost_t>>(N_cities, vector<cost_t>(N_cities, 0))));
    
    // read flights
    while( input->good() ) {
        string dep;
        string arr;
        int day;
        cost_t cost;
        
        *input >> dep >> arr >> day >> cost;
        if( !input->good() )
            break;
        
        DEBUG("departure: " << dep << "; arrives: " << arr << "; day: " << day << "; cost: " << cost);
        if( day )
            flightTable->at(day-1)[citiesName2Num[dep]][citiesName2Num[arr]] = cost;
        else
            for( int i=0; i<flightTable->size(); i++ )
                flightTable->at(i)[citiesName2Num[dep]][citiesName2Num[arr]] = cost;
    }
}

vector<int> randomPermutation(const vector<int> inp) {
    vector<int> out = inp;
    
    for( int i = out.size()-1; i>0; i-- ) {
        int pos = rand() % i;
        swap(out[i], out[pos]);
    }
    
//#ifdef DEBUG_OUT
//    cout << "random perm: \n";
//    for( auto c : out )
//        cout << c << " ";
//    cout << endl;
//#endif
    
    return out;
}

//unsigned int findPermutationCost(const vector<int> &perm) {
//    unsigned int sum = 0;
    
//    for(int i=0; i<perm.size(); i++) {
//        if( flightTable[])
//    }
//}

int main(int argc, char **argv)
{
    chrono::time_point<chrono::system_clock> startTime(chrono::system_clock::now());
    
    init();
    
//    dumpFlightDayTable(flightTable->at(0));
    
    vector<int> zonesWithoutStart;
    zonesWithoutStart.reserve(zones.size() - 1);
    for( int i=0; i<zones.size(); i++ ) {
        if( i != startZone )
            zonesWithoutStart.push_back(i);
    }
    
    node start(0, startCity);
    int counter = 0;
    // timer on random search
    while( chrono::system_clock::now() - startTime < std::chrono::milliseconds(14900) ) {
        counter++;
        
        auto perm = randomPermutation(zonesWithoutStart);
        perm.push_back(startZone);
    
        if( tryFoundPath_greedy(&start, perm.begin( )) )
            break;
    }
    
    start.dumpAnswer( );
    
    cout << "Num of permutations: " << counter;
    
    return 0;
}
