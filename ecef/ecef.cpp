// Required: C++11 (for OSX: g++ -std=c++11 ecef.cpp)
// http://www.enri.go.jp/~fks442/K_MUSEN/1st/1st021118.pdf

#include "ecef.hpp"

void ECEF::rx(double theta, double (*mat)[3]) {
    mat[0][0] = 1.0;
    mat[0][1] = 0.0;
    mat[0][2] = 0.0;
    mat[1][0] = 0.0;
    mat[1][1] = cos(theta * M_PI / 180.0);
    mat[1][2] = sin(theta * M_PI / 180.0);
    mat[2][0] = 0.0;
    mat[2][1] = -sin(theta * M_PI / 180.0);
    mat[2][2] = cos(theta * M_PI / 180.0);
}

void ECEF::ry(double theta, double (*mat)[3]) {
    mat[0][0] = cos(theta * M_PI / 180.0);
    mat[0][1] = 0.0;
    mat[0][2] = -sin(theta * M_PI / 180.0);
    mat[1][0] = 0.0;
    mat[1][1] = 1.0;
    mat[1][2] = 0.0;
    mat[2][0] = sin(theta * M_PI / 180.0);
    mat[2][1] = 0.0;
    mat[2][2] = cos(theta * M_PI / 180.0);
}

void ECEF::rz(double theta, double (*mat)[3]) {
    mat[0][0] = cos(theta * M_PI / 180.0);
    mat[0][1] = sin(theta * M_PI / 180.0);
    mat[0][2] = 0.0;
    mat[1][0] = -sin(theta * M_PI / 180.0);
    mat[1][1] = cos(theta * M_PI / 180.0);
    mat[1][2] = 0.0;
    mat[2][0] = 0.0;
    mat[2][1] = 0.0;
    mat[2][2] = 1.0;
}

void ECEF::mm(double (*mat_a)[3], double (*mat_b)[3], double (*mat_answer)[3]) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            for (int k = 0; k < 3; ++k) {
                mat_answer[i][j] += mat_a[i][k] * mat_b[k][j];
            }
        }
    }
}

void ECEF::mv(double (*mat_a)[3], double *mat_b, double *mat_answer) {
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            mat_answer[i] += mat_a[i][j] * mat_b[j];
        }
    }
}

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

vector<double> ECEF::ecef2enu(vector<double> ecef_dest, vector<double> ecef_origin) {
    double rot0[3][3];
    double rot1[3][3];
    double rot2[3][3];
    double m[3][3];
    double r[3][3];
    double p[3] = {ecef_dest[0] - ecef_origin[0], ecef_dest[1] - ecef_origin[1], ecef_dest[2] - ecef_origin[2]};
    double enu_a[3];

    vector<double> bih = ecef2bih(ecef_origin);

    rz(90.0 * M_PI / 180.0, rot0);
    ry((90.0 - bih[0]) * M_PI / 180.0, rot1);
    rz(bih[2] * M_PI / 180.0, rot2);

    mm(rot0, rot1, m);
    mm(m, rot2, r);
    mv(r, p, enu_a);

    vector<double> enu = {enu_a[0], enu_a[1], enu_a[2]};
    return enu;
}


vector<double> ECEF::bih2enu(vector<double> bih_dest, vector<double> bih_origin) {
    vector<double> ecef_dest = bih2ecef(bih_dest);
    vector<double> ecef_origin = bih2ecef(bih_origin);
    
    double rot0[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double rot1[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double rot2[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double m[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double r[3][3] = {{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}};
    double p[3];
    for (int i = 0; i < 3; i++) {
        p[i] = ecef_dest[i] - ecef_origin[i];
    }
    double enu_a[3] = {0.0, 0.0, 0.0};

    vector<double> bih = ecef2bih(ecef_origin);

    rz(90.0, rot0);
    ry(90.0 - bih[0], rot1);
    rz(bih[1], rot2);

    mm(rot0, rot1, m);
    mm(m, rot2, r);
    mv(r, p, enu_a);

    vector<double> enu = {enu_a[0], enu_a[1], enu_a[2]};
    return enu;
}

double ECEF::enu2length(vector<double> enu) {
    return sqrt(pow(enu[0], 2) + pow(enu[1], 2));
}

double ECEF::enu2angle(vector<double> enu) {
    return atan(enu[1] / enu[0]) * 180.0 / M_PI;
}

double ECEF::enu2direction(vector<double> enu) {
    double r = atan(enu[1] / enu[0]) * 180.0 / M_PI;
    double direction;
    if (enu[0] < 0 && enu[1] >= 0) direction = 270.0 + (180.0 - r);
    else direction = 90.0 - r;
    return direction;
}

double ECEF::bih2length(vector<double> bih_dest, vector<double> bih_origin) {
    vector<double> enu = bih2enu(bih_dest, bih_origin);
    return sqrt(pow(enu[0], 2) + pow(enu[1], 2));
}

double ECEF::bih2angle(vector<double> bih_dest, vector<double> bih_origin) {
    vector<double> enu = bih2enu(bih_dest, bih_origin);
    return atan(enu[1] / enu[0]) * 180.0 / M_PI;
}

double ECEF::bih2direction(vector<double> bih_dest, vector<double> bih_origin) {
    vector<double> enu = bih2enu(bih_dest, bih_origin);
    double r = atan(enu[1] / enu[0]) * 180.0 / M_PI;
    double direction;
    if (enu[0] < 0 && enu[1] >= 0) direction = 270.0 + (180.0 - r);
    else direction = 90.0 - r;
    return direction;
}