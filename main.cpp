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
#include <iomanip>
#include <list>

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

#define TIMEOUT_ENABLE
#ifdef TIMEOUT_ENABLE
#define TIMEOUT (chrono::system_clock::now() - startTime >= Timeout)
#else
#define TIMEOUT (false)
#endif

int N_cities;
int N_zone;

typedef uint16_t cost_t;

// flight table [day][departure][arriving] of costs
// dimension: Zones x Cities x Cities
typedef vector<vector<vector<cost_t >>> FlightTable;
unique_ptr<FlightTable> flightTable;
unique_ptr<FlightTable> pseudoFlightTable;
#define FLIGHT(day, from, to) (flightTable->at(day)[from][to])
#define PSEUDO_FLIGHT(day, from, to) (pseudoFlightTable->at(day)[from][to])

auto Timeout = chrono::milliseconds(3000 - 100);
chrono::time_point<chrono::system_clock> startTime;

template<typename T>
void clearPtrList(list<T*> &l)
{
    for( auto i : l )
        delete i;
    l.clear();
}

template<typename T>
bool lessPtr(T a1, T a2) {
    return *a1 < *a2;
}

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

struct city {
    int n;
    int zone;
    string name;

    bool operator <(const city& other) { return this->n < other.n; }

    bool belongsToZone(int z);
};

vector<zone> zones;
vector<city> cities;
unordered_map<string, int> citiesName2Num;
int startZone;
int startCity;

bool city::belongsToZone(int z) {
    return( !( n < zones[z].firstCity ) && !( n > zones[z].lastCity ));
}

struct node {
    int day;  // from 0 to N_zones - 1
    int city;
    node* nextNode;
    node* prevNode;
    int sum;
    int costTillDay;    // cost of all path flights till day this->day

    node(int d = 0, int c = 0, int s = 0, int cost = 0) : day(d), city(c), sum(s), costTillDay(cost)
      , nextNode(nullptr), prevNode(nullptr) {}

    explicit node (const node &n) : day(n.day), city(n.city), sum(n.sum), costTillDay(n.costTillDay), nextNode(nullptr), prevNode(nullptr)
    {
        // forward
        node* from = n.nextNode;
        node** to = &nextNode;
        node *temp = this;
        while( from ) {
            node* n = new node(from->day, from->city, from->sum, from->costTillDay);
            *to = n;
            n->prevNode = temp;
            from = from->nextNode;
            to = &((*to)->nextNode);
            temp = n;
        }

        // backward
        from = n.prevNode;
        to = &prevNode;
        temp = this;
        while( from ) {
            node* n = new node(from->day, from->city, from->sum, from->costTillDay);
            *to = n;
            n->nextNode = temp;
            from = from->prevNode;
            to = &((*to)->prevNode);
            temp = n;
        }
    }

    ~node() {
        if( nextNode ) {
            nextNode->prevNode = nullptr;
            delete nextNode;
        }
        if( prevNode ) {
            prevNode->nextNode = nullptr;
            delete prevNode;
        }
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

    void dumpPath() {
#ifdef DEBUG_OUT
        cout << sum << ":\t" << endl;
        node* n = this;
        while( n ) {
            cout << cities[n->city].n << " ";
            n = n->nextNode;
        }
        cout << endl;
#endif
    }

    // with respect to zones visited
    bool checkFlightTo(int to)
    {
        node* n = this;
        while(n) {
            if( cities[n->city].zone == (cities[to].zone) )
                if( day != N_zone - 1 || cities[to].zone != startZone )
                    return false;
            n = n->prevNode;
        }
        return true;
    }

    void applyFlight(int to, cost_t cost)
    {
        // copy this and insert between this and prevNode
        if( prevNode ) {
            node* temp = prevNode;
            prevNode           = new node(day, city, sum, costTillDay);
            temp->nextNode = prevNode;
            prevNode->prevNode = temp;
            prevNode->nextNode = this;
        }
        else {
            prevNode = new node(day, city, sum, costTillDay);
            prevNode->nextNode = this;
        }

        // fill this fields as next
        this->day = prevNode->day+1;
        this->costTillDay = prevNode->costTillDay + cost;
        this->city = to;
    }

    bool operator <(const node &other)
    {
        return costTillDay < other.costTillDay ;
    }
};

namespace pathLimits
{
// -- Limits for random evaluation --
// worstDayCosts consists of biggest costs for today
vector<cost_t> worstDayCosts;
// worstCostsToTheEnd - upper limit for all cost from today to the end
vector<int>    worstCostsToTheEnd;

// bestDayCosts consists of chippest costs for today
vector<cost_t> bestDayCosts;
// worstCostsToTheEnd - bottom limit for all cost from today to the end
vector<int>    bestCostsToTheEnd;

// bestDayCosts consists of mean costs for today
vector<cost_t> meanDayCosts;
// worstCostsToTheEnd - mean summ of all cost from today to the end
vector<int>    meanCostsForDay;

// pathCostBestForDay consists of min found for this day cost
vector<int> bestPathCostForDay;

node *bestNode;

void init( )
{
    worstDayCosts.reserve(N_zone);
    bestDayCosts.reserve(N_zone);
    meanDayCosts.reserve(N_zone);
    for ( int i        = 0; i < N_zone; i++ ) {
        cost_t min_cost = numeric_limits<cost_t>::max();
        cost_t max_cost = 0;
        int flightsDayCostSumm = 0;
        int flightsDayCounter = 0;

        for ( int j = 0; j < N_cities; j++ ) {
            max_cost = max(max_cost, *max_element(flightTable->at(i)[j].begin( ), flightTable->at(i)[j].end( )));

            for(int k = 0; k < N_cities; k ++) {
                cost_t flightCost = FLIGHT(i, j, k);
                min_cost = !flightCost ? min_cost : min(flightCost, min_cost);

                if( flightCost ) {
                    flightsDayCostSumm += flightCost;
                    flightsDayCounter++;
                }
            }
        }
        worstDayCosts.push_back(max_cost);
        bestDayCosts.push_back(min_cost);
        meanDayCosts.push_back(flightsDayCostSumm / flightsDayCounter);
    }

    worstCostsToTheEnd.resize(N_zone);
    bestCostsToTheEnd.resize(N_zone);
    worstCostsToTheEnd[N_zone - 1] = worstDayCosts[N_zone - 1];
    bestCostsToTheEnd[N_zone - 1]  = bestDayCosts[N_zone - 1];

    for ( int i = N_zone - 2; i >= 0; i-- ) {
        worstCostsToTheEnd[i] = worstCostsToTheEnd[i + 1] + worstDayCosts[i];
        bestCostsToTheEnd[i]  = bestCostsToTheEnd[i + 1] + bestDayCosts[i];
        bestCostsToTheEnd[i]  = bestCostsToTheEnd[i + 1] + bestDayCosts[i];
    }

    meanCostsForDay.resize(N_zone);
    meanCostsForDay[0] = meanDayCosts[0];
    for( int i=1; i < N_zone; i++ ) {
        meanCostsForDay[i] = meanCostsForDay[i - 1] + meanDayCosts[i];
    }

    bestPathCostForDay.resize(N_zone+1);
    bestPathCostForDay[0] = 0;
    for( int i=0; i < N_zone; i++ ) {
        bestPathCostForDay[i+1] = worstDayCosts[i] + bestPathCostForDay[i];
    }
}

bool updateBestDayPaths(node* n) {
    bool updated = false;
    while( n ) {
        if( bestPathCostForDay[n->day] > n->costTillDay ) {
            bestPathCostForDay[n->day] = n->costTillDay;
            updated = true;
        }
        n = n->nextNode;
    }
#ifdef DEBUG_OUT
    if( updated )
        DEBUG("Best Updated!");
#endif
    return updated;
}

#define LIMITS_STATISTICS
static int checkGood = 0;
static int checkBad = 0;

inline bool checkForLimits( int day, int currentCost ) {
#ifdef LIMITS_STATISTICS
    if( currentCost + bestCostsToTheEnd[day] < bestNode->sum ) {
        checkGood++;
        return true;
    }
    else {
        checkBad++;
        return false;
    }
#else
    return currentCost + bestCostsToTheEnd[day] < bestNode.sum;
#endif
}

}   // namespace pathLimits

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
    if( start->day == N_zone ) {
        start->sum = start->costTillDay;
        return true;
    }

    // Timeout
    if( TIMEOUT ) {
        DEBUG("Timeout!");
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

    while( nextCities.size()  && !TIMEOUT) {
        auto minCostFlight = min_element( nextCities.begin(), nextCities.end() );

        start->nextNode = new node(start->day + 1, minCostFlight->nextCity );
        start->nextNode->costTillDay = start->costTillDay + minCostFlight->cost;

        // change zones to leave in right part only not visited
        iter_swap( minCostFlight->zoneIt, zonesLeftToVisit.begin() + start->day );
        if( greedyOnCities(start->nextNode, zonesLeftToVisit) ) {
            start->sum = start->nextNode->sum;
            start->nextNode->prevNode = start;
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

// search backward greedy with respect to zones
node* greedyBackward(node *last, vector<int> &zonesLeftToVisit) {
    // hey! end of search!
    if( last->day == 1 ) {
        if( cost_t cost = FLIGHT(0, startCity, last->city) ) {
            node* start = new node(0, startCity);
            start->sum = last->sum + cost;
            start->costTillDay = 0;
            start->nextNode = last;
            last->costTillDay = cost;
            return start;
        }
        else
            return nullptr;
    }

    // Timeout
    if( TIMEOUT ) {
        DEBUG("Timeout!");
        return nullptr;
    }

    int numOfIter = N_zone - last->day;

    // vector of cities
    vector<todayFlight> prevCities;
    prevCities.reserve(N_cities - numOfIter);  // upper limit
    for( auto zoneIt = zonesLeftToVisit.begin() + numOfIter; zoneIt != zonesLeftToVisit.end(); zoneIt++) {
        for ( int i = zones[*zoneIt].firstCity; i <= zones[*zoneIt].lastCity; i++ ) {
            if ( cost_t cost = FLIGHT(last->day-1, i, last->city) )
                prevCities.push_back(todayFlight(i, cost, zoneIt));
        }
    }

    while( prevCities.size() ) {
        auto minCostFlight = min_element( prevCities.begin(), prevCities.end() );

        node* prevNode = new node(last->day-1, minCostFlight->nextCity);
        prevNode->sum = last->sum + minCostFlight->cost;
        prevNode->nextNode = last;

        // change zones to leave in right part only not visited
        iter_swap( minCostFlight->zoneIt, zonesLeftToVisit.begin() + numOfIter );
        if( node* start = greedyBackward(prevNode, zonesLeftToVisit) ) {
            last->costTillDay = prevNode->costTillDay + minCostFlight->cost;
            last->prevNode = prevNode;
            return start;
        }
        else {
            prevNode->nextNode = nullptr;
            delete prevNode;
            // replace zones back
            iter_swap( minCostFlight->zoneIt, zonesLeftToVisit.begin() + numOfIter  );
            prevCities.erase(minCostFlight);
            continue;
        }
    }

    return nullptr;
}

// random search through all cities but with respect to zones
#define RANDOM_TO_GREEDY // if on some day cost for this day is better than best known - go on with greedy
bool randomOnCities(node *start, vector<int> &zonesLeftToVisit) {
    // hey! end of search!
    if( start->day == N_zone )
        return true;

    // Timeout
    if( TIMEOUT ) {
        return false;
    }

    // vector of cities
    vector<todayFlight> nextCities;
    nextCities.reserve(N_cities - start->day);  // upper limit
    for( auto zoneIt = zonesLeftToVisit.begin() + start->day; zoneIt != zonesLeftToVisit.end(); zoneIt++) {
        if(( *zoneIt != startZone ) || (  start->day == N_zone-1 )) {   // fly to start zone at the end only
            for ( int i = zones[*zoneIt].firstCity; i <= zones[*zoneIt].lastCity; i++ ) {
                if ( cost_t cost = FLIGHT(start->day, start->city, i) ) {
                    if ( pathLimits::checkForLimits(start->day + 1, start->sum + cost) )
                        nextCities.push_back(todayFlight(i, cost, zoneIt));
                }
            }
        }
    }

    while( nextCities.size() && !TIMEOUT) {

        unsigned int random = rand() % nextCities.size();
        auto minCostFlight = nextCities.begin() + random;

        start->nextNode = new node(start->day + 1, minCostFlight->nextCity );
        start->nextNode->costTillDay = start->costTillDay + minCostFlight->cost;
        start->nextNode->sum = start->nextNode->costTillDay;

        // change zones to leave in right part only not visited
        iter_swap( minCostFlight->zoneIt, zonesLeftToVisit.begin() + start->day );
#ifdef RANDOM_TO_GREEDY
        if( start->day > 0 &&  // greedy filtered out flights less than found path for day 0
                start->nextNode->costTillDay < pathLimits::bestPathCostForDay[start->nextNode->day] ) {

            DEBUG("Random - found better cost for day. Go greedy! Day: " << start->nextNode->day << " -> " <<  start->nextNode->costTillDay << " vs "
                                                                         << pathLimits::bestPathCostForDay[start->nextNode->day]);
            pathLimits::bestPathCostForDay[start->nextNode->day] = start->nextNode->costTillDay;

            if( greedyOnCities(start->nextNode, zonesLeftToVisit) ) {
                start->sum = start->nextNode->sum;
                return true;
            }
        }
        else
#endif
        if( randomOnCities(start->nextNode, zonesLeftToVisit) ) {
            start->sum = start->nextNode->sum;
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

    pathLimits::init();
    
    pseudoFlightTable.reset(new FlightTable);
    *pseudoFlightTable.get() = *flightTable.get();
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

vector<int> flightsPerDay;
vector<double> growthCoefficients; // coefficients of growth of possible path if everything is equally distributed

void updateGrowthCoefficients() {
    flightsPerDay = vector<int>(N_zone, 0);
    for( int i = 0; i < N_zone; i++ ) {
        for( auto from : flightTable->at(i) )
            flightsPerDay[i] += count_if( from.begin(), from.end(), [](cost_t cost){
                return cost > 0;
            });
    }

    growthCoefficients = vector<double>(N_zone, 0);
    for( int i = 0; i < N_zone; i++ ) {
        double k = (float) flightsPerDay[i] / (N_cities * (N_cities - 1) ); // k - density of flights per day coefficient
        double l = N_cities / N_zone; // l - cities per zone coefficient
        if( i == N_zone - 1)
            growthCoefficients[i] = k * l * 1;
        else
            growthCoefficients[i] = k * l * ( N_zone - i - 1 );
    }

#ifdef DEBUG_OUT
    DEBUG("\nGrowth coef-s:");
    for( int i = 0; i < N_zone; i++ )
        DEBUG(i << " -> \t" << setprecision(6) << growthCoefficients[i]);
    DEBUG("");
#endif
}

void filterListProbabilistic(list<node*>& paths, int maxLimit ) {
    if( paths.size() <= maxLimit )
        return;

#ifdef DEBUG_OUT
    cout << endl << "Input costs: ";
    for( auto p : paths )
        cout << p->costTillDay << ", ";
    cout << endl;
#endif

    // mean cost at input
#ifdef DEBUG_OUT
    int mean = 0;
    for( auto p : paths)
        mean += p->costTillDay;
    cout << "mean on input: " << (float) mean / paths.size() << endl;
#endif

    auto worstPath = max_element(paths.begin(), paths.end(), lessPtr<node*>);
    // find coefficient knowing biggest cost such that ( 1 / cost_biggest * coef = 100 )
    float coef = 100 * (*worstPath)->costTillDay;

    struct pathNumberAndProbability {
        int path;
        int probability;
    };

    list<pathNumberAndProbability> pr(paths.size());   // probabilities
    auto it = paths.begin();
    int i=0;
    for( auto &prIt : pr ) {
        prIt = {i, (int) ceil( 1.f / (*it)->costTillDay * coef)};
        it++;
        i++;
    }

    // calculate probabilities sum
    long allPathsCostSum = 0;
    for( auto &p : pr )
        allPathsCostSum += p.probability;
    DEBUG("all points: " << allPathsCostSum);

    // sort probabilities in descending
    pr.sort([](pathNumberAndProbability &pr1, pathNumberAndProbability &pr2){
        return pr1.probability > pr2.probability;
    });

#ifdef DEBUG_OUT
    cout << endl << "Probabilities: ";
    for( auto p : pr )
        cout << p.path << " - " << p.probability << ", ";
    cout << endl;
#endif

    long randMod = allPathsCostSum;
    for( int i=0; i < maxLimit; i++ ) {
        long randPoints = (rand() * (long) RAND_MAX + rand()) % randMod;
        DEBUG("random points: " << randPoints);
        long partSum = 0;
        for( auto p = pr.begin(); p != pr.end(); p++ ) {
            partSum += p->probability;
            if( partSum >= randPoints ) {
                DEBUG("path: " << p->path);
                randMod -= p->probability;
                pr.erase(p);
                break;
            }
            assert( p != pr.end() );
        }
    }

    pr.sort([&](pathNumberAndProbability &pr1, pathNumberAndProbability &pr2){
        return pr1.path < pr2.path;
    });

    // delete extra paths
    auto pathIt = paths.begin();
    auto erasePathIt = pr.begin();
    int inputPathsSize = paths.size();
    for( int i=0; i < inputPathsSize; i++ ) {
        if( i == (*erasePathIt).path ) {
            delete *pathIt;
            paths.erase(pathIt++);
            erasePathIt++;
        }
        else
            pathIt++;
    }

    // mean cost at input
#ifdef DEBUG_OUT
    mean = 0;
    for( auto p : paths)
        mean += p->costTillDay;
    cout << "mean on output: " << (float) mean / paths.size() << endl;
#endif

#ifdef DEBUG_OUT
    cout << endl << "Output costs: ";
    for( auto p : paths )
        cout << p->costTillDay << ", ";
    cout << endl;
#endif
}

void filterListRandom(list<node*>& paths, int maxLimit ) {
    if( paths.size( ) <= maxLimit )
        return;

    int n_inp = paths.size();
    for( auto it = paths.begin(); it != paths.end();  ) {
        if( rand() % n_inp < maxLimit )
            it++;
        else {
            delete *it;
            paths.erase(it++);
        }
    }
}

void recalcPseudoFlightTable(list<node*> &paths) {  // paths - pointers on last nodes !
    const float mult = 0.95;
    for( node* p : paths ) {
        node* n = p->prevNode;
        while( n ) {
            cost_t* cost = &PSEUDO_FLIGHT(n->day, n->city, n->nextNode->city);
            int newCost = round((*cost) * mult);
            if( newCost == *cost && *cost > 1 )
                *cost -= 1;
            else
                *cost = newCost;
            
            n = n->prevNode;
        }
    }
}

void recalcPathsFromPseudo(list<node*> &paths) {  // paths - pointers on last nodes !
    const float mult = 0.95;
    for( node* p : paths ) {
        node* n = p->prevNode;
        while( n->prevNode )
            n = n->prevNode;
        
        while( n->nextNode ) {
            n->nextNode->costTillDay = n->costTillDay + FLIGHT(n->day, n->city, n->nextNode->city);
            n = n->nextNode;
        }
    }
}

void probabilisticDynamic() {
    const int M = 100;  // desired num of paths on exit
    updateGrowthCoefficients();

    vector<int> desiredPathNumForDay = vector<int>(N_zone, M);
    desiredPathNumForDay.back() = max((int) ceil(desiredPathNumForDay.back() /  growthCoefficients.back()), M);
    for( int i = N_zone - 2; i > 0; i-- ) {
        desiredPathNumForDay[i] = max((int) ceil(desiredPathNumForDay[i+1] /  growthCoefficients[i]), M);
    }
#ifdef DEBUG_OUT
    DEBUG("\nDesired paths for day:");
    for( int i = 0; i < N_zone; i++ )
        DEBUG(i << " -> \t" << desiredPathNumForDay[i]);
    DEBUG("");
#endif

    list<node*> inputList;
    inputList.push_back(new node(0, startCity));
    list<node*> outputList;
    for(int i=0; i<N_zone; i++) {

        if( TIMEOUT ) {
            clearPtrList(inputList);
            clearPtrList(outputList);
            return;
        }
    
        int nextDay = i + 1;
        int n_inp = growthCoefficients[i] * inputList.size();
        int n_out = desiredPathNumForDay[nextDay];
        bool mustFilter = (( n_inp > n_out ) && ( i < N_zone - 1 ));
        int worstPath, bestPath;
        int meanPath = 0;
        if( mustFilter ) {
            worstPath = (*max_element(inputList.begin( ), inputList.end( ), lessPtr<node*>))->costTillDay;
            worstPath += pathLimits::worstDayCosts[i];
            bestPath = (*min_element(inputList.begin( ), inputList.end( ), lessPtr<node*>))->costTillDay;
            bestPath += pathLimits::bestDayCosts[i];
            for( auto p : inputList )
                meanPath += p->costTillDay;
            meanPath /= inputList.size();
            meanPath += pathLimits::meanDayCosts[i];
        }
        
        // put all paths after this day flights into outputList with some probability
        auto &dayFlights = pseudoFlightTable->at(i);
        for( auto path : inputList ) {
            for( int to = 0; to < N_cities; to++ ) {
                if( path->city == to )
                    continue;
                if( cost_t cost = dayFlights[path->city][to] ) {
                    if( path->checkFlightTo(to) ) {
                        if( mustFilter ) {
                            int costNext = path->costTillDay + cost;
                            double factor = 1;
                            if( costNext > meanPath ) {
                                factor = (float) (worstPath - costNext) / ( worstPath - meanPath );
                            }
                            else if( costNext < meanPath ) {
                                factor = (costNext - bestPath + (float) n_inp / n_out * ( meanPath - costNext ))
                                         / ( meanPath - bestPath );
                            }
                            long randMod = (factor > 0) ? (long) (n_inp / factor) : numeric_limits<long>::max();
                            long num = (rand() * (long) RAND_MAX + rand()) % randMod;
                            if( num > n_out )
                                continue;
                        }
                        outputList.push_back(new node(*path));
                        outputList.back()->applyFlight(to, cost);
                    }
                }
            }
        }
        DEBUG("day " << i << " out paths not filtered: " << outputList.size());

        if( TIMEOUT ) {
            clearPtrList(inputList);
            clearPtrList(outputList);
            return;
        }

        if( i != N_zone - 1 ) {
//            filterListProbabilistic(outputList, n_out);
            filterListRandom(outputList, n_out);
            DEBUG("out paths filtered: " << outputList.size());
            swap(inputList, outputList);
            clearPtrList(outputList);
        }
    }

    DEBUG("paths on output: " << outputList.size());
    
    recalcPseudoFlightTable(outputList);
    recalcPathsFromPseudo(outputList);
    
    node* n = *(min_element(outputList.begin(), outputList.begin(), lessPtr<node*>));
    n->sum = n->costTillDay;
    while( n->prevNode ) {
        n->prevNode->sum = n->sum;
        n = n->prevNode;
    }
    DEBUG("best cost: " << n->sum);
//    n->dumpAnswer();

    if( pathLimits::bestNode ) {
        if( pathLimits::bestNode->sum > n->sum ) {
            delete pathLimits::bestNode;
            pathLimits::bestNode = new node(*n);
        }
    }
    else
        pathLimits::bestNode = new node(*n);
    
    clearPtrList(inputList);
    clearPtrList(outputList);
}

int main(int argc, char **argv)
{
    startTime = chrono::system_clock::now();

    srand(time(NULL));

    init();

    vector<int> zonesWithoutStart;
    zonesWithoutStart.reserve(zones.size() - 1);
    for( int i=0; i<zones.size(); i++ ) {
        if( i != startZone )
            zonesWithoutStart.push_back(i);
    }

    node start(0, startCity);

    // greedy forward
    vector<int> zonesWithStart = zonesWithoutStart;
    zonesWithStart.push_back(startZone);
    greedyOnCities(&start, zonesWithStart);
    pathLimits::bestNode = new node(start);
    pathLimits::updateBestDayPaths(pathLimits::bestNode);
    pathLimits::bestNode->dumpPath();

    // greedy backward
    for( int i = zones[startZone].firstCity; i<= zones[startZone].lastCity; i++ ) {

        if( TIMEOUT )
            break;

        node* lastNode = new node(N_zone, i);
        node* startNode = greedyBackward(lastNode, zonesWithoutStart);
        if( startNode ) {
            DEBUG("backward:");
            startNode->dumpPath();

            pathLimits::updateBestDayPaths(startNode);

            if( startNode->sum < pathLimits::bestNode->sum ) {
                DEBUG("new sum (backward): " << startNode->sum);
                delete pathLimits::bestNode;
                node *n = new node(*startNode);
                pathLimits::bestNode = new node(*startNode);
            }
        }
        delete lastNode;
    }
    
    int counter = 0;
    while( !TIMEOUT ) {
        probabilisticDynamic( );
        counter++;
    }
    
    pathLimits::bestNode->dumpAnswer( );
    delete pathLimits::bestNode;
    
    DEBUG("Num of loops: " << counter);
    return 0;
}
