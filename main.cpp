#include <iostream>
#include <vector>
#include <unordered_map>

#define MAX_N 300

int N_cities;
int N_zone;

struct city {
    int n;
    int zone;
};

std::unordered_map<int, std::string> zones2Numbers;
std::unordered_map<std::string, city> cities;

struct SingleFlightDay {
    struct tableColumn {
        tableColumn() {
            column.fill(0);
        }
        std::array<int, MAX_N> column;
    };
    std::array<tableColumn, MAX_N> table;
};

void init()
{
    std::cin >> N_zone;
    for( int i=0; i< N_zone; i++ ) {
        std::string zoneName;
        std::getline(std::cin, zoneName);
        zones2Numbers[i] = zoneName;
        
        std::string cieties;
        std::getline(std::cin, cities)
    }
    
}

int main(int argc, char **argv)
{
    std::cout << "Hello, World!" << std::endl;
    return 0;
}