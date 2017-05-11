#include "ecef/ecef.hpp"
#include <iostream>
#include <vector>
#include <cmath>


using namespace std;

int main() {
    ECEF ecef;
    vector<double> origin{38.13877338, 140.89872429, 44.512};
    vector<double> dest{38.14227288, 140.93265738, 45.664};

    //vector<double> ecef_origin = ecef.bih2ecef(origin);
    //vector<double> ecef_dest = ecef.bih2ecef(dest);

    //cout << fixed << ecef_origin[0] << " " << fixed << ecef_origin[1] << " " << fixed << ecef_origin[2] << endl;
    //cout << fixed << ecef_dest[0] << " " << fixed << ecef_dest[1] << " " << fixed << ecef_dest[2] << endl;

    //vector<double> xyz = ecef.bih2enu(dest, origin);
    //cout << fixed << xyz[0] << " " << fixed << xyz[1] << " " << fixed << xyz[2] << " " << endl;
    //cout << fixed << sqrt(pow(xyz[0], 2) + pow(xyz[1], 2)) << endl;

    double length = ecef.bih2length(dest, origin);
    double angle = ecef.bih2angle(dest, origin);
    double direction = ecef.bih2direction(dest, origin);

    cout << length << " " << angle << " " << direction <<  endl;
    return 0;
}