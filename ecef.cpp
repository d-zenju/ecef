// Required: C++11 (for OSX: g++ -std=c++11 ecef.cpp)

#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// ECEF(Earth Centered, Earth Fixed) <-> WGS84(Latitude, Longitude, Height) convert
class ECEF {
    private:
        double a;
        double one_f;
        double b;
        double e2;
        double ed2;
        double n(double p) {
            return a / sqrt(1.0 - e2 * pow(sin(p * M_PI / 180.0), 2));
        }

    public:
        // construct
        ECEF(void) {
            a = 6378137.0;
            one_f = 298.257223563;
            b = a * (1.0 - 1.0 / one_f);
            e2 = (1.0 / one_f) * (2 - (1.0 / one_f));
            ed2 = e2 * a * a / (b * b);
        }

        // WGS84 to ECEF
        // input :: phi = latitude, lambda = longitude, height = height
        // output :: vector (x, y, z)
        vector<double> bih2ecef(double phi, double lambda, double height) {
            double x = (n(phi) + height) * cos(phi * M_PI / 180.0) * cos(lambda * M_PI / 180.0);
            double y = (n(phi) + height) * cos(phi * M_PI / 180.0) * sin(lambda * M_PI / 180.0);
            double z = (n(phi) * (1 - e2) + height) * sin(phi * M_PI / 180.0);
            vector<double> xyz{x, y, z};
            return xyz;
        }

        // WGS84 to ECEF
        // input :: vector(latitude, longitude, height)
        // output :: vector (x, y, z)
        vector<double> bih2ecef(vector<double> bih) {
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
        vector<double> ecef2bih(double x, double y, double z) {
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
        vector<double> ecef2bih(vector<double> xyz) {
            double x = xyz[0];
            double y = xyz[1];
            double z = xyz[2];
            double p = sqrt(x * x + y * y);
            double sita = (180.0 / M_PI) * atan2(z * a, p * b);
            double phi = (180.0 / M_PI) * atan2(z + ed2 * b * pow(sin(sita * M_PI / 180.0), 3), (p - e2 * a * pow(cos(sita * M_PI / 180.0), 3)));
            double lambda = (180.0 / M_PI) * atan2(y, x);
            double height = (p / cos(phi * M_PI / 180.0)) - n(phi);
            vector<double> bih{phi, lambda, height};
            return bih;
        }
};


int main() {
    ECEF ecef;
    vector<double> bih{38.13579617, 140.91581617, 41.940};
    vector<double> xyz = ecef.bih2ecef(bih);
    vector<double> bih2 = ecef.ecef2bih(xyz);
    cout << fixed << xyz[0] << " " << fixed << xyz[1] << " " << fixed << xyz[2] << " " << endl;
    cout << fixed << bih2[0] << " " << fixed << bih2[1] << " " << fixed << bih2[2] << " " << endl;
    printf("%.8f, %.8f, %.3f\n", bih2[0], bih2[1], bih2[2]);
    return 0;
}