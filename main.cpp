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

//const string filename = "../datasets/4.in";
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

// full search preserving zone arrange
bool greedyOnPredefinedZoneSequence(node* start, vector<int>::iterator path) {
    vector<todayFlight> nextCities(zones[*path].lastCity - zones[*path].firstCity + 1, todayFlight(0, 0));
    for( int i = 0; i < nextCities.size(); i++) {
        nextCities[i] = todayFlight(zones[*path].firstCity + i, FLIGHT(start->day, start->city, zones[*path].firstCity + i));
    }
    
    // thats why little bit greedy
    // TODO: maybe better search min element than sort anytime ??
    // and not to paste empty (cost=0) flights ?
    sort(nextCities.begin(), nextCities.end());
    for( auto &f : nextCities ) {
        if( !f.cost )
            break;
    
        start->nextNode = new node(start->day + 1, f.nextCity);
    
        if( *path == startZone ) { // start zone achieved
            start->sum += f.cost;
            return true;
        }
        
        vector<int>::iterator nextZone(path);
        if( greedyOnPredefinedZoneSequence(start->nextNode, ++nextZone) ) {
            start->sum = f.cost + start->nextNode->sum;
            return true;    // success, path found!
        }
        else {
            delete start->nextNode;
            start->nextNode = nullptr;
            continue;
        }
    }
    
//    if( start->day > 1)
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

vector<int> randomPermutation_zonesConnected(vector<int> inp) {
    vector<int> out;
    int curZone = startZone;
    int day = 0;
    while( inp.size() ) {
        vector<int> temp = inp;
        while ( temp.size( ) ) {
            int randIndex = rand() % temp.size( );
            int nextZone = temp[randIndex];
            if ( zones[curZone].isConnected(day, zones[nextZone]) ) {
                out.push_back(nextZone);
                curZone = nextZone;
                break;
            }
            else
                temp.erase(temp.begin() + randIndex);
        }
        if( !temp.size() ) { // no next zone found
            break;
        }
        else {
            for( auto i = inp.begin(); i<inp.end(); i++ ) {
                if( *i == *(out.end() - 1) ) {
                    inp.erase(i);
                    break;
                }
            }
            continue;
        }
    }
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
    
    node bestNode(0, startCity);
    bestNode.sum = std::numeric_limits<int>::max();
    
    node start(0, startCity);
    
    zonesWithoutStart.push_back(startZone);
    if( greedyOnCities(&start, zonesWithoutStart) )
        start.dumpAnswer();
    return 0;
    
    int counter = 0;
    // timer on random search
    while( chrono::system_clock::now() - startTime < Timeout ) {
        counter++;
        
//        auto perm = randomPermutation(zonesWithoutStart);
//        perm.push_back(startZone);
    
        int findZoneCounter = 0;
        vector<int>perm;
        while( chrono::system_clock::now() - startTime < Timeout ) {
            findZoneCounter++;
            perm = randomPermutation_zonesConnected(zonesWithoutStart);
            if( (perm.size() == N_zone - 1 ) && zones[perm[perm.size()-1]].isConnected(N_zone - 1, zones[startZone]) ) {
                perm.push_back(startZone);
                break;
            }
        }
    
        DEBUG("Perm found counter: " << findZoneCounter << "\n");
    
        if( greedyOnPredefinedZoneSequence(&start, perm.begin( )) ) {
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
