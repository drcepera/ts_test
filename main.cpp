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

using namespace std;

#define MAX_N 300

#define DEBUG_OUT
#ifdef DEBUG_OUT
    #define DEBUG(x) do { std::cout << x << std::endl; } while (0)
#else
    #define DEBUG(x) ;
#endif

int N_cities;
int N_zone;

typedef uint16_t cost_t;

struct city {
    int n;
    int zone;
//    string name;
    
    bool operator <(const city& other) { return this->n < other.n; }
};

struct zone {
    int n;
    string name;
    int firstCity;
    int lastCity;
};

unordered_map<int, zone> N_zones;
unordered_map<string, city> cities;
int startZone;
int startCity;

// flight table [day][departure][arriving] of costs
// dimension: Zones x Cities x Cities
typedef vector<vector<vector<cost_t >>> FlightTable;
std::unique_ptr<FlightTable> flightTable;

void dumpFlightDayTable(const vector<vector<cost_t>> &costs) {
    // col-row names
    std::vector<std::pair<std::string, city>> citiesOrdered;
    citiesOrdered.resize(cities.size());
    std::copy(cities.begin(), cities.end(), citiesOrdered.begin());
    std::sort(citiesOrdered.begin(), citiesOrdered.end(), [](const std::pair<std::string, city>& c1, const std::pair<std::string, city>& c2){
        return c1.second.n < c2.second.n;
    });
    
    // out table head
    cout << "\t";
    for( auto c : citiesOrdered )
        cout << c.first << "\t";
    cout << endl;
    
    // rows
    for( auto c_dep : citiesOrdered ) {
        cout << c_dep.first << "\t";
        for( auto c_arr : citiesOrdered )
            cout << costs[c_dep.second.n][c_arr.second.n] << "\t";
        cout << endl;
    }
}

const std::string filename = "";

void init()
{
    std::istream *input = filename.empty() ? &cin : new ifstream(filename);

    std::string startCityName;
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
            city newCity({N_cities++, i});
            cities[cityName] = newCity;
            DEBUG("city " << cityName << " number: " << cities[cityName].n << " zone: " << cities[cityName].zone);

            newZone.lastCity++;

            if( cityName == startCityName ) {
                startCity = newCity.n;
                startZone = newCity.zone;
                DEBUG("start city found: " << cityName << " number: " << cities[cityName].n << " zone: " << cities[cityName].zone);
            }
        }

        N_zones[i] = newZone;

        DEBUG("zone " << N_zones[i].name << " #" << N_zones[i].n <<
              " firstCity: " << N_zones[i].firstCity << " lastCity: " << N_zones[i].lastCity);
    }
    
    // create flight table
    flightTable.reset(new vector<vector<vector<cost_t >>>(N_zone, vector<vector<cost_t>>(N_cities, vector<cost_t>(N_cities, 0))));
    
    // read flights
    while( input->good() ) {
        std::string dep;
        std::string arr;
        int day;
        cost_t cost;
        
        *input >> dep >> arr >> day >> cost;
        if( !input->good() )
            break;
        
        DEBUG("departure: " << dep << "; arrives: " << arr << "; day: " << day << "; cost: " << cost);
        if( day )
            flightTable->at(day-1)[cities[dep].n][cities[arr].n] = cost;
        else
            for( int i=0; i<flightTable->size(); i++ )
                flightTable->at(i)[cities[dep].n][cities[arr].n] = cost;
    }
    
//    dumpFlightDayTable(flightTable->at(0));
}

vector<int> randomPermutation(int start, int stop) {
    assert( stop >= start);
    
    vector<int> out(stop - start);
    for( int i = start; i<stop; i++ )
        out[i-start] = i;
    
    for( int i = out.size()-1; i>0; i-- ) {
        int pos = rand() % i;
        std::swap(out[i], out[pos]);
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
    init();
    
    auto perm = randomPermutation(1, N_zone);
    cout << "random perm: \n";
    for( auto c : perm )
        cout << c << " ";
    cout << endl;
    
    perm.push_back(0);
    
    return 0;
}
