#include <stdio.h>
#include <math.h>

// Rotational matrices for rotating the coordinate system
void Rx(double deg, double result[3][3]) 
{
    double sin_x = sin(deg * M_PI / 180.0);
    double cos_x = cos(deg * M_PI / 180.0);
    result[0][0] = 1;     result[0][1] = 0;     result[0][2] = 0;
    result[1][0] = 0;     result[1][1] = cos_x; result[1][2] = sin_x;
    result[2][0] = 0;     result[2][1] = -sin_x; result[2][2] = cos_x;
}

void Ry(double deg, double result[3][3])
{
    double sin_x = sin(deg * M_PI / 180.0);
    double cos_x = cos(deg * M_PI / 180.0);
    result[0][0] = cos_x; result[0][1] = 0;     result[0][2] = -sin_x;
    result[1][0] = 0;     result[1][1] = 1;     result[1][2] = 0;
    result[2][0] = sin_x; result[2][1] = 0;     result[2][2] = cos_x;
}

void Rz(double deg, double result[3][3]) 
{
    double sin_x = sin(deg * M_PI / 180.0);
    double cos_x = cos(deg * M_PI / 180.0);
    result[0][0] = cos_x; result[0][1] = sin_x; result[0][2] = 0;
    result[1][0] = -sin_x; result[1][1] = cos_x; result[1][2] = 0;
    result[2][0] = 0;     result[2][1] = 0;     result[2][2] = 1;
}

/*******************************************************************************************
    Rotates the coordinate system along the Z-, Y-, and X-axis,
    and returns the new coordinates of (x,y,z) expressed in this
    new, shifted coordiante system
********************************************************************************************/
void rotate_frame_zyx(double x, double y, double z, double z_offset, double y_offset, double x_offset, double *x_new, double *y_new, double *z_new) 
{
    double vector[3][1] = {{x}, {y}, {z}};
    double temp[3][1];
    double rot_x[3][3], rot_y[3][3], rot_z[3][3];

    Rx(x_offset, rot_x);
    Ry(y_offset, rot_y);
    Rz(z_offset + 90, rot_z);

    // First apply Rz
    for (int i = 0; i < 3; i++) 
    {
        temp[i][0] = rot_z[i][0] * vector[0][0] + rot_z[i][1] * vector[1][0] + rot_z[i][2] * vector[2][0];
    }
    // Then apply Ry
    for (int i = 0; i < 3; i++) {
        vector[i][0] = temp[i][0];
        temp[i][0] = rot_y[i][0] * vector[0][0] + rot_y[i][1] * vector[1][0] + rot_y[i][2] * vector[2][0];
    }
    // Finally apply Rx
    for (int i = 0; i < 3; i++) 
    {
        vector[i][0] = temp[i][0];
        temp[i][0] = rot_x[i][0] * vector[0][0] + rot_x[i][1] * vector[1][0] + rot_x[i][2] * vector[2][0];
    }

    *x_new = temp[0][0];
    *y_new = temp[1][0];
    *z_new = temp[2][0];
}

void offset_north_pitch_roll(double x, double y, double z, double north_offset, double pitch_offset, double roll_offset, double *x_new, double *y_new, double *z_new) 
{
    rotate_frame_zyx(x, y, z, north_offset, pitch_offset, roll_offset, x_new, y_new, z_new);
}

