#include "ecef/ecef.hpp"
#include <iostream>
#include <vector>


using namespace std;

int main() {
    ECEF ecef;
    vector<double> bih{38.13579617, 140.91581617, 41.940};
    vector<double> xyz = ecef.bih2ecef(bih);
    cout << fixed << xyz[0] << " " << fixed << xyz[1] << " " << fixed << xyz[2] << " " << endl;
    return 0;
}