// Required: C++11 (for OSX: g++ -std=c++11 ecef.cpp)
// http://www.enri.go.jp/~fks442/K_MUSEN/1st/1st021118.pdf

#include "ecef.hpp"

// WGS84 to ECEF
// input :: phi = latitude, lambda = longitude, height = height
// output :: vector (x, y, z)
vector<double> ECEF::bih2ecef(double phi, double lambda, double height) {
    double x = (ECEF::n(phi) + height) * cos(phi * M_PI / 180.0) * cos(lambda * M_PI / 180.0);
    double y = (ECEF::n(phi) + height) * cos(phi * M_PI / 180.0) * sin(lambda * M_PI / 180.0);
    double z = (ECEF::n(phi) * (1 - e2) + height) * sin(phi * M_PI / 180.0);
    vector<double> xyz{x, y, z};
    return xyz;
    }

// WGS84 to ECEF
// input :: vector(latitude, longitude, height)
// output :: vector (x, y, z)
vector<double> ECEF::bih2ecef(vector<double> bih) {
    double phi = bih[0];
    double lambda = bih[1];
    double height = bih[2];
    double x = (n(phi) + height) * cos(phi * M_PI / 180.0) * cos(lambda * M_PI / 180.0);
    double y = (n(phi) + height) * cos(phi * M_PI / 180.0) * sin(lambda * M_PI / 180.0);
    double z = (n(phi) * (1 - e2) + height) * sin(phi * M_PI / 180.0);
    vector<double> xyz{x, y, z};
    return xyz;
}

// ECEF to WGS84
// input :: x = x, y = y, z = z
// output :: vector (latitude, longitude, height)
vector<double> ECEF::ecef2bih(double x, double y, double z) {
    double p = sqrt(x * x + y * y);
    double sita = (180.0 / M_PI) * atan2(z * a, p * b);
    double phi = (180.0 / M_PI) * atan2(z + ed2 * b * pow(sin(sita * M_PI / 180.0), 3), (p - e2 * a * pow(cos(sita * M_PI / 180.0), 3)));
    double lambda = (180.0 / M_PI) * atan2(y, x);
    double height = (p / cos(phi * M_PI / 180.0)) - n(phi);
    vector<double> bih{phi, lambda, height};
    return bih;
}

// ECEF to WGS84
// input :: vector (x, y, z)
// output :: vector (latitude, longitude, height)
vector<double> ECEF::ecef2bih(vector<double> xyz) {
    double x = xyz[0];
    double y = xyz[1];
    double z = xyz[2];
    double p = sqrt(x * x + y * y);
    double sita = (180.0 / M_PI) * atan2(z * ECEF::a, p * b);
    double phi = (180.0 / M_PI) * atan2(z + ed2 * b * pow(sin(sita * M_PI / 180.0), 3), (p - e2 * a * pow(cos(sita * M_PI / 180.0), 3)));
    double lambda = (180.0 / M_PI) * atan2(y, x);
    double height = (p / cos(phi * M_PI / 180.0)) - n(phi);
    vector<double> bih{phi, lambda, height};
    return bih;
}