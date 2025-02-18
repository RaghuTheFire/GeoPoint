#include <math.h>
#include <stdio.h>
#define PI 3.14159265358979323846
// WGS-84 geodetic constants
const double a = 6378137.0;  // WGS-84 Earth semimajor axis (m)
const double b = 6356752.314245;  // Derived Earth semiminor axis (m)
const double f = (a - b) / a;  // Ellipsoid Flatness
const double f_inv = 1.0 / f;  // Inverse flattening
const double a_sq = a * a;
const double b_sq = b * b;
const double e_sq = f * (2 - f);  // Square of Eccentricity

void deg2rad(double degrees, double *radians) 
{
    *radians = degrees * (PI / 180.0);
}

void rad2deg(double radians, double *degrees) 
{
    *degrees = radians * (180.0 / PI);
}

int sign(float x) 
{
    return (x > 0) - (x < 0);
}

typedef struct 
{
    float xSouthing;
    float yWesting;
} Coordinate;

/******************************************************************************************************************
    Converts WGS-84 Geodetic point (lat, lon, h) to the Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z).
    :param lat: latitude in degrees
    :param lon: longitude in degrees
    :param height: meters
    :output: (x, y, z) in [meters, meters, meters]
*******************************************************************************************************************/
void geodetic_to_ecef(double lat, double lon, double height, double *x, double *y, double *z) 
{
    double lambd, phi, s, N;
    deg2rad(lat, &lambd);
    deg2rad(lon, &phi);
    s = sin(lambd);
    N = a / sqrt(1 - e_sq * s * s);

    double sin_lambda = sin(lambd);
    double cos_lambda = cos(lambd);
    double cos_phi = cos(phi);
    double sin_phi = sin(phi);

    *x = (height + N) * cos_lambda * cos_phi;
    *y = (height + N) * cos_lambda * sin_phi;
    *z = (height + (1 - e_sq) * N) * sin_lambda;
}
/******************************************************************************************************************
 Converts Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z) to the WGS-84 Geodetic point (lat, lon, height)
    
    :param x: meters
    :param y: meters
    :param z: meters
    :output: (lat, lon, height) in [degrees, degrees, meters]
*******************************************************************************************************************/
void ecef_to_geodetic(double x, double y, double z, double *lat, double *lon, double *height) 
{
    double eps = e_sq / (1.0 - e_sq);
    double p = sqrt(x * x + y * y);
    double q = atan2((z * a), (p * b));
    double sin_q = sin(q);
    double cos_q = cos(q);
    double sin_q_3 = sin_q * sin_q * sin_q;
    double cos_q_3 = cos_q * cos_q * cos_q;
    double phi = atan2((z + eps * b * sin_q_3), (p - e_sq * a * cos_q_3));
    double lambd = atan2(y, x);
    double v = a / sqrt(1.0 - e_sq * sin(phi) * sin(phi));
    *height = (p / cos(phi)) - v;

    rad2deg(phi, lat);
    rad2deg(lambd, lon);
}

/******************************************************************************************************************
    Converts the Earth-Centered Earth-Fixed (ECEF) coordinates (x, y, z) to 
    East-North-Up coordinates in a Local Tangent Plane that is centered at the 
    (WGS-84) Geodetic point (lat0, lon0, h0).
    
    :param x: meters
    :param y: meters
    :param z: meters
    :param lat0: latitude in degrees
    :param lon0: longitude in degrees
    :param height0: meters
    :output: (xEast, yNorth, zUp) in [meters, meters, meters]
*******************************************************************************************************************/
void ecef_to_enu(double x, double y, double z, double lat0, double lon0, double h0, double *xEast, double *yNorth, double *zUp) 
{
    double lamb, phi, s, N;
    deg2rad(lat0, &lamb);
    deg2rad(lon0, &phi);
    s = sin(lamb);
    N = a / sqrt(1 - e_sq * s * s);

    double sin_lambda = sin(lamb);
    double cos_lambda = cos(lamb);
    double sin_phi = sin(phi);
    double cos_phi = cos(phi);

    double x0 = (h0 + N) * cos_lambda * cos_phi;
    double y0 = (h0 + N) * cos_lambda * sin_phi;
    double z0 = (h0 + (1 - e_sq) * N) * sin_lambda;

    double xd = x - x0;
    double yd = y - y0;
    double zd = z - z0;

    *xEast = -sin_phi * xd + cos_phi * yd;
    *yNorth = -cos_phi * sin_lambda * xd - sin_lambda * sin_phi * yd + cos_lambda * zd;
    *zUp = cos_lambda * cos_phi * xd + cos_lambda * sin_phi * yd + sin_lambda * zd;
}
/******************************************************************************************************************
    Converts the geodetic WGS-84 coordinated (lat, lon, h) to 
    East-North-Up coordinates in a Local Tangent Plane that is centered at the 
    (WGS-84) Geodetic point (lat0, lon0, h0).
    
    :param lat: latitude of target in degrees
    :param lon: longitude of target in degrees
    :param h: height of target in meters
    :param lat0: latitude of origin in degrees
    :param lon0: longitude of origin in degrees
    :param h0: height of origin meters
    :output: (xEast, yNorth, zUp) in [meters, meters, meters]
*******************************************************************************************************************/
void geodetic_to_enu(float lat, float lon, float h, float lat0, float lon0, float h0, float *xEast, float *yNorth, float *zUp) 
{
    float x, y, z;
    geodetic_to_ecef(lat, lon, h, &x, &y, &z);
    ecef_to_enu(x, y, z, lat0, lon0, h0, xEast, yNorth, zUp);
}
/******************************************************************************************************************
    South African Coordinate Reference System (Hartebeesthoek94) to Geodetic
    From "CDNGI Coordinate Conversion Utility v1 Sep 2009.xls".
    CM = central meridian

    :param loMeridian: central meridian (Lo.) degrees
    :param yWesting: coordinates measured in meters from the CM of the respective zone, increasing from the CM (where Y=0) in a westerly direction. Y is +ve west of the CM and -ve east of the CM.
    :param xSouthing: coordinates measured in meters southwards from the equator, increasing from the equator (where X = 0m) towards the South Pole.
    :output: (lat, lon) in [degree, degree]
*******************************************************************************************************************/
void sacrs_to_geodetic(int loMeridian, int yWesting, int xSouthing, double *lat, double *lon) 
{
    double loMeridianRadians;
    deg2rad(loMeridian, &loMeridianRadians);
    
    double ec_sq = (a_sq - b_sq) / a_sq;   // e^2
    double ep_sq = (a_sq - b_sq) / b_sq;   // e'^2
    double n = (a - b) / (a + b);
    double n_2 = n * n;
    double n_3 = n_2 * n;
    double n_4 = n_2 * n_2;
    double n_5 = n_2 * n_3;
    double p2 = 3.0 / 2.0 * n - 27.0 / 32.0 * n_3 + 269.0 / 512.0 * n_5;
    double p4 = 21.0 / 16.0 * n_2 - 55.0 / 32.0 * n_4;
    double p6 = 151.0 / 96.0 * n_3 - 417.0 / 128.0 * n_5;
    double p8 = 1097.0 / 512.0 * n_4;
    double p10 = 8011.0 / 2560.0 * n_5;
    double a0 = 1.0 / (n + 1.0) * (1.0 + 1.0 / 4.0 * n_2 + 1.0 / 64.0 * n_4);
    double footBar = xSouthing / (a * a0);
    double foot = footBar + p2 * sin(2.0 * footBar) + p4 * sin(4.0 * footBar) + p6 * sin(6.0 * footBar) + p8 * sin(8.0 * footBar) + p10 * sin(10.0 * footBar);
    double Nf = a / sqrt(1.0 - ec_sq * sin(foot) * sin(foot));

    double b1 = 1.0 / (Nf * cos(foot));
    double b2 = tan(foot) / (2.0 * Nf * Nf * cos(foot));
    double b3 = (1.0 + 2.0 * tan(foot) * tan(foot) + ep_sq * cos(foot) * cos(foot)) / (6.0 * Nf * Nf * Nf * cos(foot));
    double b4 = (tan(foot) * (5.0 + 6.0 * tan(foot) * tan(foot) + ep_sq * cos(foot) * cos(foot))) / (24.0 * Nf * Nf * Nf * Nf * cos(foot));
    double b5 = (5.0 + 28.0 * tan(foot) * tan(foot) + 24.0 * tan(foot) * tan(foot) * tan(foot) * tan(foot)) / (120.0 * Nf * Nf * Nf * Nf * Nf * cos(foot));
    double d1 = cos(foot) * (1.0 + ep_sq * cos(foot) * cos(foot));
    double d2 = -1.0 / 2.0 * cos(foot) * cos(foot) * tan(foot) * (1.0 + 4.0 * ep_sq * cos(foot) * cos(foot));

    double latRadians = -(foot - b2 * d1 * yWesting * yWesting + (b4 * d1 + b2 * b2 * d2) * yWesting * yWesting * yWesting * yWesting);
    double lonRadians = loMeridianRadians - (b1 * yWesting - b3 * yWesting * yWesting * yWesting + b5 * yWesting * yWesting * yWesting * yWesting * yWesting);
    
    rad2deg(latRadians, lat);
    rad2deg(lonRadians, lon);
}
/******************************************************************************************************************
    Geodetic to South African Coordinate Reference System 
   
    :param lat: latitude in degrees
    :param lon: longitude in degrees
    :output: (loMeridian, yWesting, xSouthing) in [degree, meters, meters]
*******************************************************************************************************************/
Coordinate geodetic_to_sacrs(float lat, float lon) 
{
    float a = 6378137.0; // semi-major axis
    float b = 6356752.314245; // semi-minor axis
    float a_sq = a * a;
    float b_sq = b * b;

    float loMeridian = 2 * (int)(lon / 2) + sign(lon);
    float loMeridianRadians = deg2rad(loMeridian);
    float ec_sq = (a_sq - b_sq) / a_sq;
    float ep_sq = (a_sq - b_sq) / b_sq;
    float n = (a - b) / (a + b);
    float n_2 = n * n;
    float n_3 = n_2 * n;
    float n_4 = n_2 * n_2;
    float n_5 = n_2 * n_3;
    float A0 = 1.0 / (n + 1.0) * (1.0 + 1.0 / 4.0 * n_2 + 1.0 / 64.0 * n_4);
    float A2 = -1.0 / (n + 1.0) * (3.0 / 2.0 * n - 3.0 / 16.0 * n_3 - 3.0 / 128.0 * n_5);
    float A4 = 1.0 / (n + 1.0) * (15.0 / 16.0 * n_2 - 15.0 / 64.0 * n_4);
    float A6 = -1.0 / (n + 1.0) * (35.0 / 48.0 * n_3 - 175.0 / 768.0 * n_5);
    float A8 = 1.0 / (n + 1.0) * (315.0 / 512.0 * n_4);
    float A10 = 1.0 / (n + 1.0) * (693.0 / 1280.0 * n_5);

    float latRadians = deg2rad(-lat);
    float lonRadians = deg2rad(lon);
    float G = a * (A0 * latRadians + A2 * sin(2 * latRadians) + A4 * sin(4 * latRadians) + A6 * sin(6 * latRadians) + A8 * sin(8 * latRadians) + A10 * sin(10 * latRadians));
    float N = a / sqrt(1.0 - ec_sq * sin(latRadians) * sin(latRadians));

    float latCos = cos(latRadians);
    float latCos_2 = latCos * latCos;
    float latCos_3 = latCos_2 * latCos;
    float latCos_4 = latCos_2 * latCos_2;
    float latCos_5 = latCos_4 * latCos;
    float latTan = tan(latRadians);
    float latTan_2 = latTan * latTan;
    float latTan_4 = latTan_2 * latTan_2;
    float a1 = N * latCos;
    float a2 = -1.0 / 2.0 * N * latCos_2 * latTan;
    float a3 = -1.0 / 6.0 * N * latCos_3 * (1.0 - latTan_2 + ep_sq * latCos_2);
    float a4 = 1.0 / 24.0 * N * latCos_4 * latTan * (5 - latTan_2 + 9.0 * ep_sq * latCos_2);
    float a5 = 1.0 / 120.0 * N * latCos_5 * (5.0 - 18.0 * latTan_2 + latTan_4);
    float l = lonRadians - loMeridianRadians;
    float l_2 = l * l;
    float l_3 = l * l_2;
    float l_4 = l_2 * l_2;
    float l_5 = l_2 * l_3;

    Coordinate result;
    result.xSouthing = G - a2 * l_2 + a4 * l_4;
    result.yWesting = -(a1 * l - a3 * l_3 + a5 * l_5);

    return result;
}
