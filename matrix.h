//From :Eplankton Date: 2021/11/15(Start date)
#ifndef _MATRIX
#define _MATRIX

#include <iostream>
#include <vector>
#include <cmath>
#include <string>

#define PI 3.1415926
#define MAXO 50
#define Set Init

using namespace std;

class matrix
{
public:
    string name;
    int row;
    int column;
    int rank;
    double module;       //Only for vectors.
    double det;          //Only for cubes.
    vector<double> body; //Contain the actual elements of matrix.

    vector<double> &Init(const string init_name, const int init_row, const int init_column); //Initialize a matrix by name, row and column.
    vector<double> &Init(const int init_row, const int init_column);                         //Non-named initialization.
    void Head(void);                                                                         //Display main information.
    void Out(char key);                                                                      //Print by format.
    void Out(void);                                                                          //Print the matrix on screen.
    void Clear(void);                                                                        //Empty the matrix into void.
    double &operator()(int row, int column);                                                 //Access to elements.
    double Module(void);                                                                     //Calculate the module of a vector.
    int Rank(void);
    double Det(void);

    matrix Power(unsigned int pow);         //Calculate the power of matrix.
    matrix Tran(void);                      //Calculate the transpose of matrix.
    matrix Approx(void);                    //Calculate the  approx of matrix.
    matrix operator+(const matrix &latter); //Calculate the sum of matrices.
    matrix operator-(const matrix &latter); //Calculate the subtract of matrices.
    matrix operator*(const matrix &second); //Calculate the multiplication of matrices.
    matrix operator^(const int pow);        //Short for Power().
    matrix operator[](const int column);    //Short for getColumn().
};

class R_3d : public matrix
{
public:
    double theta_x, theta_y, theta_z;
    void Boost(double alpha, double beta, double gramma)
    {
        row = 3, column = 3;
        theta_x = alpha, theta_y = beta, theta_z = gramma;
        matrix Rx, Ry, Rz;
        Rx.Init(3, 3) = {1, 0, 0, 0, cos(theta_x), -sin(theta_x), 0, sin(theta_x), cos(theta_x)};
        Ry.Init(3, 3) = {cos(theta_y), 0, -sin(theta_y), 0, 1, 0, sin(theta_y), 0, cos(theta_y)};
        Rz.Init(3, 3) = {cos(theta_z), -sin(theta_z), 0, sin(theta_z), cos(theta_z), 0, 0, 0, 1};
        this->body = (Rx * Ry * Rz).body;
    }
};

class R_2d : public matrix
{
public:
    double theta_x;
    void Boost(double alpha)
    {
        row = 2, column = 2;
        theta_x = alpha;
        this->body = {cos(theta_x), -sin(theta_x), sin(theta_x), cos(theta_x)};
    }
};

typedef matrix mat;
typedef matrix vct;

vct getColumn(matrix Mt, int c);
vct getRow(matrix Mt, int r);

matrix Tran(matrix Mt);
matrix Power(const matrix Mt, unsigned int pow);
double Dot(vct &va, vct &vb);
double Angle(vct &va, vct &vb);
double standard_echelon(double matrix[MAXO][MAXO], int r, int c, int x, int y);
matrix Echelon(matrix Mt);
int Rank(matrix Mt);
matrix Inv(matrix Mt);

double determinant(double matrix[MAXO][MAXO], int order);                     //From: https://zhuanlan.zhihu.com/p/305328519
double laplace_expansion(double matrix[MAXO][MAXO], int r, int c, int order); //From: https://zhuanlan.zhihu.com/p/305328519
double Det(matrix Mt);

#endif