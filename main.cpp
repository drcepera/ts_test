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

//const string filename = "../datasets/3.in";
const string filename = "";

//#define DEBUG_OUT
#ifdef DEBUG_OUT
    #define DEBUG(x) do { cout << x << endl; } while (0)
#else
    #define DEBUG(x) ;
#endif

int N_cities;
int N_zone;

typedef uint16_t cost_t;

// flight table [day][departure][arriving] of costs
// dimension: Zones x Cities x Cities
typedef vector<vector<vector<cost_t >>> FlightTable;
unique_ptr<FlightTable> flightTable;
#define FLIGHT(day, from, to) (flightTable->at(day)[from][to])

struct city {
    int n;
    int zone;
    string name;
    
    bool operator <(const city& other) { return this->n < other.n; }
};

auto Timeout = chrono::milliseconds(3000-100);
chrono::time_point<chrono::system_clock> startTime;

struct zone {
    int n;
    string name;
    int firstCity;
    int lastCity;
    
    bool isConnected(int day, const zone &nextZone) {
        for( int i = firstCity; i <= lastCity; i++ ) {
            for( int j = nextZone.firstCity; j <= nextZone.lastCity; j++ ) {
                if( FLIGHT(day, i, j) )
                    return true;
            }
        }
        return false;
    }
};

vector<zone> zones;
vector<city> cities;
unordered_map<string, int> citiesName2Num;
int startZone;
int startCity;

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
                 << FLIGHT(n->day, n->city, n->nextNode->city) << endl;
    
            n = n->nextNode;
        }
    }
};

// vector of next possible cities from path iterator
struct todayFlight {
    int nextCity;
    cost_t cost;
    vector<int>::iterator zoneIt;
    
    todayFlight(int city, cost_t c) : nextCity(city), cost(c) {}
    todayFlight(int city, cost_t c, vector<int>::iterator zIt) : nextCity(city), cost(c), zoneIt(zIt) {}
    
    bool operator < (const todayFlight &f2) {
        if( !f2.cost ) return true;
        if( !cost ) return false;
        return cost < f2.cost;
    }
};

// just search in greedy manner through all cities but with respect to zones
// zonesLeftToVisit - iterator on vector part with not visited yet cities
bool greedyOnCities(node *start, vector<int> &zonesLeftToVisit) {
    // hey! end of search!
    if( start->day == N_zone )
        return true;
    
    // Timeout
    if( chrono::system_clock::now() - startTime >= Timeout ) {
//        cout << "Timeout!";
        return false;
    }
    
    // vector of cities
    vector<todayFlight> nextCities;
    nextCities.reserve(N_cities - start->day);  // upper limit
    for( auto zoneIt = zonesLeftToVisit.begin() + start->day; zoneIt != zonesLeftToVisit.end(); zoneIt++) {
        if(( *zoneIt != startZone ) || (  start->day == N_zone-1 )) {   // fly to start zone at the end only
            for ( int i = zones[*zoneIt].firstCity; i <= zones[*zoneIt].lastCity; i++ ) {
                if ( cost_t cost = FLIGHT(start->day, start->city, i) )
                    nextCities.push_back(todayFlight(i, cost, zoneIt));
            }
        }
    }
    
    while( nextCities.size() ) {
        auto minCostFlight = min_element( nextCities.begin(), nextCities.end() );
        
        start->nextNode = new node(start->day + 1, minCostFlight->nextCity );
        
        // change zones to leave in right part only not visited
        iter_swap( minCostFlight->zoneIt, zonesLeftToVisit.begin() + start->day );
        if( greedyOnCities(start->nextNode, zonesLeftToVisit) ) {
            start->sum = minCostFlight->cost + start->nextNode->sum;
            return true;
        }
        else {
            delete start->nextNode;
            start->nextNode = nullptr;
            // replace zones back
            iter_swap( minCostFlight->zoneIt, zonesLeftToVisit.begin() + start->day  );
            nextCities.erase(minCostFlight);
            continue;
        }
    }
    
    return false;
}

// random search through all cities but with respect to zones
// zonesLeftToVisit - iterator on vector part with not visited yet cities
bool randomOnCities(node *start, vector<int> &zonesLeftToVisit) {
    // hey! end of search!
    if( start->day == N_zone )
        return true;
    
    // Timeout
    if( chrono::system_clock::now() - startTime >= Timeout ) {
//        cout << "Timeout!";
        return false;
    }
    
    // vector of cities
    vector<todayFlight> nextCities;
    nextCities.reserve(N_cities - start->day);  // upper limit
    for( auto zoneIt = zonesLeftToVisit.begin() + start->day; zoneIt != zonesLeftToVisit.end(); zoneIt++) {
        if(( *zoneIt != startZone ) || (  start->day == N_zone-1 )) {   // fly to start zone at the end only
            for ( int i = zones[*zoneIt].firstCity; i <= zones[*zoneIt].lastCity; i++ ) {
                if ( cost_t cost = FLIGHT(start->day, start->city, i) )
                    nextCities.push_back(todayFlight(i, cost, zoneIt));
            }
        }
    }
    
    while( nextCities.size() ) {
        unsigned int random = rand() % nextCities.size();
        auto minCostFlight = nextCities.begin() + random;
        
        start->nextNode = new node(start->day + 1, minCostFlight->nextCity );
        
        // change zones to leave in right part only not visited
        iter_swap( minCostFlight->zoneIt, zonesLeftToVisit.begin() + start->day );
        if( greedyOnCities(start->nextNode, zonesLeftToVisit) ) {
            start->sum = minCostFlight->cost + start->nextNode->sum;
            return true;
        }
        else {
            delete start->nextNode;
            start->nextNode = nullptr;
            // replace zones back
            iter_swap( minCostFlight->zoneIt, zonesLeftToVisit.begin() + start->day  );
            nextCities.erase(minCostFlight);
            continue;
        }
    }
    
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
    
//    DEBUG("num of zones: " << N_zone << ", start city: " << startCityName);
    
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
//            DEBUG("city " << cityName << " number: " << cities[newCity.n].n << " zone: " << cities[newCity.n].zone);

            newZone.lastCity++;

            if( cityName == startCityName ) {
                startCity = newCity.n;
                startZone = newCity.zone;
//                DEBUG("start city found: " << cityName << " number: " << cities[newCity.n].n << " zone: " << cities[newCity.n].zone);
            }
    
            citiesName2Num[newCity.name] = newCity.n;
        }

        zones.push_back(newZone);

//        DEBUG("zone " << zones[i].name << " #" << zones[i].n <<
//              " firstCity: " << zones[i].firstCity << " lastCity: " << zones[i].lastCity);
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
        
//        DEBUG("departure: " << dep << "; arrives: " << arr << "; day: " << day << "; cost: " << cost);
        if( day ) {
            cost_t* c = &( FLIGHT(day - 1, citiesName2Num[dep], citiesName2Num[arr]) );
            *c = (*c) ? min(*c, cost) : cost;
        }
        else
            for( int i=0; i<flightTable->size(); i++ ) {
                cost_t* c = &( FLIGHT(i, citiesName2Num[dep], citiesName2Num[arr]) );
                *c = (*c) ? min(*c, cost) : cost;
            }
    }
    
    const int reserve_ms = 200;
    if( N_zone <= 20 && N_cities < 50 )
        Timeout = chrono::milliseconds(3000 - reserve_ms);
    else if( N_zone <= 100 && N_cities < 200 )
        Timeout = chrono::milliseconds(5000 - reserve_ms);
    else
        Timeout = chrono::milliseconds(15000 - reserve_ms);
}

vector<int> randomPermutation(const vector<int> inp) {
    vector<int> out = inp;
    
    for( int i = out.size()-1; i>0; i-- ) {
        int pos = rand() % i;
        swap(out[i], out[pos]);
    }
    
#ifdef DEBUG_OUT
    cout << "random perm: \n";
    for( auto c : out )
        cout << c << " ";
    cout << endl;
#endif
    
    return out;
}

int main(int argc, char **argv)
{
    startTime = chrono::system_clock::now();
    
    srand(time(NULL));
    
    init();
    
//    dumpFlightDayTable(flightTable->at(0));
//    dumpFlightDayTable(flightTable->at(1));
//    dumpFlightDayTable(flightTable->at(2));
    
    vector<int> zonesWithoutStart;
    zonesWithoutStart.reserve(zones.size() - 1);
    for( int i=0; i<zones.size(); i++ ) {
        if( i != startZone )
            zonesWithoutStart.push_back(i);
    }
    
    node start(0, startCity);
    
    vector<int> zonesWithStart = zonesWithoutStart;
    zonesWithStart.push_back(startZone);
    greedyOnCities(&start, zonesWithStart);
    
    node bestNode = start;
    start = node(0, startCity);
    
    int counter = 0;
    // timer on random search
    while( chrono::system_clock::now() - startTime < Timeout ) {
        counter++;
        
        if( randomOnCities(&start, zonesWithStart) ) {
           DEBUG("new sum : " << start.sum);
            if( start.sum < bestNode.sum ) {
                if( bestNode.nextNode ) {
                    delete bestNode.nextNode;
                    bestNode.nextNode = nullptr;
                }
                bestNode = start;
                start = node(0, startCity);
            }
        }
    }
    
    bestNode.dumpAnswer( );
    
    DEBUG("Num of permutations: " << counter);
    
    return 0;
}
