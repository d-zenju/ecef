// Required: C++11 (for OSX: g++ -std=c++11 ecef.cpp)
// http://www.enri.go.jp/~fks442/K_MUSEN/1st/1st021118.pdf

#ifndef _ECEF_H_
#define _ECEF_H_

#include <iostream>
#include <vector>
#include <cmath>

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// ECEF(Earth Centered, Earth Fixed) <-> WGS84(Latitude, Longitude, Height) convert
class ECEF {

    private:
        const double a = 6378137.0;
        const double one_f = 298.257223563;
        const double b = a * (1.0 - 1.0 / one_f);
        const double e2 = (1.0 / one_f) * (2 - (1.0 / one_f));
        const double ed2 = e2 * a * a / (b * b);
        double n(double p) {
            return a / sqrt(1.0 - e2 * pow(sin(p * M_PI / 180.0), 2));
        }
        void rx(double theta, double (*mat)[3]);
        void ry(double theta, double (*mat)[3]);
        void rz(double theta, double (*mat)[3]);
        void mm(double (*mat_a)[3], double (*mat_b)[3], double (*mat_answer)[3]);
        void mv(double (*mat_a)[3], double *mat_b, double *mat_answer);

    public:
        // WGS84 to ECEF
        // input :: phi = latitude, lambda = longitude, height = height
        // output :: vector (x, y, z)
        vector<double> bih2ecef(double phi, double lambda, double height);

        // WGS84 to ECEF
        // input :: vector(latitude, longitude, height)
        // output :: vector (x, y, z)
        vector<double> bih2ecef(vector<double> bih);

        // ECEF to WGS84
        // input :: x = x, y = y, z = z
        // output :: vector (latitude, longitude, height)
        vector<double> ecef2bih(double x, double y, double z);

        // ECEF to WGS84
        // input :: vector (x, y, z)
        // output :: vector (latitude, longitude, height)
        vector<double> ecef2bih(vector<double> xyz);

        // ENU
        vector<double> ecef2enu(vector<double> ecef_dest, vector<double> ecef_origin);
        vector<double> bih2enu(vector<double> bih_dest, vector<double> bih_origin);

        // ENU to length, angle
        double enu2length(vector<double> enu);
        double enu2angle(vector<double> enu);
        double enu2direction(vector<double> enu);

        // BIH to length, angle
        double bih2length(vector<double> bih_dest, vector<double> bih_origin);
        double bih2angle(vector<double> bih_dest, vector<double> bih_origin);
        double bih2direction(vector<double> bih_dest, vector<double> bih_origin);
};

#endif