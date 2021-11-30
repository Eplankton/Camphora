#include <iostream>
#include <vector>
#include <cmath>
#include <string>
#include "matrix.h"

using namespace std;

vector<double> &matrix::Init(const int init_row, const int init_column)
{
    row = init_row;
    column = init_column;
    module = EOF;
    det = EOF;
    rank = EOF;
    body.clear();
    for (int i = 0; i < row * column; i++)
    {
        body.push_back(0);
    }
    /* body.resize(init_row);
    for (int i = 0; i < init_row; i++)
    {
        body[i].resize(init_column);
    } */
    return body;
}

vector<double> &matrix::Init(const string init_name, const int init_row, const int init_column)
{
    name = init_name;
    row = init_row;
    column = init_column;
    module = EOF;
    rank = EOF;
    det = EOF;
    body.clear();
    for (int i = 0; i < row * column; i++)
    {
        body.push_back(0);
    }
    return body;
}

void matrix::Head(void)
{
    //string tpid = typeid(*this).name();
    //string::iterator itr = tpid.begin();
    //tpid.erase(itr);
    //tpid = "[" + tpid + "]: ";
    //cout << tpid;
    cout << name << ' ' << '<' << row << ',' << column << '>' << endl;
}

void matrix::Out(char key)
{
    if (!body.empty())
    {
        Head();
        if (key == 'f') //Print in format.
        {
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < column; j++)
                {
                    if (j != 0)
                        printf("%g\t", (*this)(i, j));
                    else
                        printf("[ %g\t", (*this)(i, j));
                }
                printf("]\n");
            }
        }
        else
        {
            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < column; j++)
                {
                    if (j != 0)
                        printf(" %g ", (*this)(i, j));
                    else
                        printf("[ %g ", (*this)(i, j));
                }
                printf("]\n");
            }
        }
        cout << endl;
    }
    else
    {
        cout << "\nWrong with error:-1 No value exists!\n";
    }
}

void matrix::Out(void)
{
    if (!body.empty())
    {
        Head();
        for (int i = 0; i < row; i++)
        {
            for (int j = 0; j < column; j++)
            {
                if (j != 0)
                    printf(" %g ", (*this)(i, j));
                else
                    printf("[ %g ", (*this)(i, j));
            }
            printf("]\n");
        }

        cout << endl;
    }
    else
    {
        cout << "\nWrong with error:-1  No value exists!\n";
    }
}

void matrix::Clear(void)
{
    this->body.clear();
    this->name = "";
    this->row = 0;
    this->column = 0;
    this->rank = 0;
    this->module = 0;
    this->det = 0;
}

matrix matrix::Approx(void)
{
    for (int i = 0; i < this->row * this->column; i++)
    {
        if (fabs(this->body[i]) <= LIMIT)
        {
            this->body[i] = 0;
        }
    }
    return *this;
}

double &matrix::operator()(int i, int j)
{
    double &p = body[i * column + j];
    return p;
}

matrix matrix::operator+(const matrix &latter)
{
    if (this->row == latter.row && this->column == latter.column && (this->row * this->column != 0) && (latter.row * latter.column != 0))
    {
        matrix result;
        result.name = "(" + (this->name) + "+" + latter.name + ")";
        result.row = this->row;
        result.column = this->column;
        for (int i = 0; i < latter.row * latter.column; i++)
        {
            result.body.push_back((this->body)[i] + latter.body[i]);
        }
        return result;
    }
    else
    {
        matrix fail;
        fail.Clear();
        cout << "Wrong with error:-1\n";
        return fail;
    }
}

matrix matrix::operator-(const matrix &latter)
{
    if (this->row == latter.row && this->column == latter.column && (this->row * this->column != 0) && (latter.row * latter.column != 0))
    {
        matrix result;
        result.name = "(" + (this->name) + "-" + latter.name + ")";
        result.row = this->row;
        result.column = this->column;
        for (int i = 0; i < latter.row * latter.column; i++)
        {
            result.body.push_back((this->body)[i] - latter.body[i]);
        }
        return result;
    }
    else
    {
        matrix fail;
        fail.Clear();
        cout << "Wrong with error:-1\n";
        return fail;
    }
}

vct getColumn(matrix Mt, const int c)
{
    if (c < Mt.column)
    {
        vct a;
        a.Init("getColumn{}", Mt.row, 1);
        a.name = Mt.name + "[" + to_string(c) + "]";
        a.body.clear();
        for (int m = 0; m < Mt.row; m++)
        {
            a.body.push_back(Mt(m, c));
        }
        return a;
    }
    else
    {
        vct fail;
        fail.Clear();
        cout << "Wrong with error:-1\n";
        return fail;
    }
}

matrix matrix::operator[](const int c)
{
    return getColumn(*this, c);
}

vct getRow(matrix Mt, const int r)
{
    if (r < Mt.row)
    {
        vct a;
        a.Init("getRow{}", 1, Mt.column);
        a.name = Mt.name + ".getRow{" + to_string(r) + "}";
        a.body.clear();

        for (int m = 0; m < Mt.column; m++)
        {
            a.body.push_back(Mt(r, m));
        }
        return a;
    }
    else
    {
        vct fail;
        fail.Clear();
        cout << "Wrong with error:-1\n";
        return fail;
    }
}

matrix matrix::operator*(const matrix &second)
{
    if (column == second.row)
    {
        matrix result;
        result.name = "(" + this->name + "*" + second.name + ")";
        result.row = this->row;
        result.column = second.column;
        result.body.clear();

        for (int i = 0; i < result.row; i++)
        {
            for (int j = 0; j < result.column; j++)
            {
                result.body.push_back(0);
            }
        }

        for (int i = 0; i < this->row; i++)
        {
            for (int j = 0; j < second.column; j++)
            {
                for (int c = 0; c < this->column; c++)
                {
                    result.body[(i)*result.column + j] += (this->body[(i) * (this->column) + c] * second.body[(c)*second.column + j]);
                }
            }
        }
        return result;
    }
    else
    {
        matrix fail;
        fail.Clear();
        cout << "Wrong with error:-1\n";
        return fail;
    }
}

matrix matrix::Power(unsigned int pow)
{
    matrix result;
    string temp_name = name;
    if (row == column)
    {
        if (pow > 0)
        {
            result.row = row;
            result.column = column;
            result.body = body;
            for (int i = 1; i < pow; i++)
            {
                result = result * (*this);
            }
            result.name = "(" + temp_name + '^' + to_string(pow) + ")";
            return result;
        }
        else
        {
            result.name = "(" + temp_name + '^' + to_string(pow) + ")";
            result.row = row;
            result.column = column;
            result.body = body;

            for (int i = 0; i < row; i++)
            {
                for (int j = 0; j < column; j++)
                {
                    if (i == j)
                    {
                        result.body[(i)*column + j] = 1;
                    }
                    else
                    {
                        result.body[(i)*column + j] = 0;
                    }
                }
            }

            return result;
        }
    }
    else
    {
        cout << "\nWrong with error:-1\n";
        return *this;
    }
}

matrix Power(const matrix Mt, unsigned int pow)
{
    matrix result;
    string temp_name = Mt.name;
    if (Mt.row == Mt.column)
    {
        if (pow > 0)
        {
            result.row = Mt.row;
            result.column = Mt.column;
            result.body = Mt.body;
            for (int i = 1; i < pow; i++)
            {
                result = result * (Mt);
            }
            result.name = "(" + temp_name + '^' + to_string(pow) + ")";
            return result;
        }
        else
        {
            result.name = "(" + temp_name + '^' + to_string(pow) + ")";
            result.row = Mt.row;
            result.column = Mt.column;
            result.body = Mt.body;

            for (int i = 0; i < Mt.row; i++)
            {
                for (int j = 0; j < Mt.column; j++)
                {
                    if (i == j)
                    {
                        result.body[(i)*Mt.column + j] = 1;
                    }
                    else
                    {
                        result.body[(i)*Mt.column + j] = 0;
                    }
                }
            }
            return result;
        }
    }
    else
    {
        cout << "\nWrong with error:-1\n";
        matrix fail;
        fail.Clear();
        cout << "Wrong with error:-1\n";
        return fail;
    }
}

matrix matrix::operator^(const int pow)
{
    if (row == column)
    {
        switch (pow)
        {
        case -1:
            return Inv(*this);
            break;
        default:
            return this->Power(pow);
            break;
        }
    }
    else
    {
        matrix fail;
        fail.Clear();
        cout << "Wrong with error:-1\n";
        return fail;
    }
}

matrix matrix::Tran(void)
{
    matrix result;
    result.name = "(" + name + "^T)";
    result.row = column;
    result.column = row;
    result.body = body;
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < column; j++)
        {
            result.body[j + (i)*column] = body[i + (j)*column];
        }
    }
    return result;
}

matrix Tran(matrix Mt)
{
    matrix result;
    result.name = "(" + Mt.name + "^T)";
    result.row = Mt.column;
    result.column = Mt.row;
    result.body = Mt.body;
    for (int i = 0; i < result.row; i++)
    {
        for (int j = 0; j < result.column; j++)
        {
            result(i, j) = Mt(j, i);
        }
    }
    return result;
}

double vct::Module(void)
{
    if (this->row == 1 || this->column == 1)
    {
        double sum = 0;
        for (int i = 0; i < row * column; i++)
        {
            sum += body[i] * body[i];
        }
        module = sqrt(sum);
        sum = 0;
        return module;
    }
    else
    {
        cout << "Wrong with error:-1\n";
        return EOF;
    }
}

double Dot(vct &va, vct &vb)
{
    if (va.row * va.column == vb.row * vb.column && (va.row == 1 || va.column == 1) && (vb.row == 1 || vb.column == 1))
    {
        double ans;
        for (int i = 0; i < (va.row * va.column); i++)
        {
            ans += (va.body[i]) * (vb.body[i]);
        }
        return ans;
    }
    else
    {
        cout << "Wrong with error:-1\n";
        return EOF;
    }
}

double Angle(vct &va, vct &vb)
{
    if (va.row * va.column == vb.row * vb.column && (va.row == 1 || va.column == 1) && (vb.row == 1 || vb.column == 1))
    {
        double theta, dot, mod;
        mod = va.Module() * vb.Module();
        dot = Dot(va, vb);
        theta = acos(dot / mod);
        return theta;
    }
    else
    {
        cout << "Wrong with error:-1\n";
        return EOF;
    }
}

double standard_echelon(double matrix[MAXO][MAXO], int r, int c, int x, int y)
{
    int i, j, k, l, total[MAXO] = {0};
    double times, temp, result = 0, original_matrix[MAXO][MAXO];

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            original_matrix[i][j] = matrix[i][j];
    for (i = 0; i < r - 1; i++)
        for (k = i + 1; k < r; k++)
        {
            j = 0;
            while (matrix[i][j] == 0)
                j++;
            if (matrix[i][j] != 0)
            {
                times = matrix[k][j] / matrix[i][j];
                for (j = 0; j < c; j++)
                    matrix[k][j] -= matrix[i][j] * times;
            }
        }
    for (i = 0; i < r; i++)
    {
        j = 0;
        while (matrix[i][j] == 0)
            j++;
        if (matrix[i][j] != 0)
        {
            times = matrix[i][j];
            for (j = 0; j < c; j++)
                matrix[i][j] /= times;
        }
    }
    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            if (matrix[i][j] == 0)
                total[i]++;
            else
                break;
    for (l = r - 1; l > 0; l--)
        for (i = 0; i < l; i++)
            if (total[l] < total[i])
                for (j = 0; j < c; j++)
                {
                    temp = matrix[l][j];
                    matrix[l][j] = matrix[i][j];
                    matrix[i][j] = temp;
                }
    for (i = 0; i < r; i++)
    {
        j = 0;
        while (matrix[i][j] == 0)
            j++;
        if (matrix[i][j] != 0)
            for (k = 0; k < i; k++)
            {
                times = matrix[k][j] / matrix[i][j];
                for (l = 0; l < c; l++)
                    matrix[k][l] -= times * matrix[i][l];
            }
    }
    result = matrix[x][y];
    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            matrix[i][j] = original_matrix[i][j];
    if (fabs(result) <= LIMIT)
        result = 0;

    return result;
}

matrix Echelon(matrix Mt)
{
    matrix result = Mt;
    result.name = "e#" + Mt.name;

    double TEMP_matrix[MAXO][MAXO] = {0};

    for (int i = 0; i < result.row; i++)
    {
        for (int j = 0; j < result.column; j++)
        {
            TEMP_matrix[i][j] = Mt(i, j);
        }
    }

    for (int i = 0; i < result.row; i++)
    {
        for (int j = 0; j < result.column; j++)
        {
            result(i, j) = standard_echelon(TEMP_matrix, result.row, result.column, i, j);
        }
    }
    return result;
}

int Rank(matrix Mt)
{
    int none_zero = 0;
    Mt.rank = 0;
    matrix TEMP = Echelon(Mt);
    for (int i = 0; i < Mt.row; i++)
    {
        for (int j = 0; j < Mt.column; j++)
        {
            if (TEMP(i, j) != 0)
            {
                none_zero = 1;
                break;
            }
        }
        if (none_zero == 1)
        {
            Mt.rank++;
        }
        none_zero = 0;
    }

    return Mt.rank;
}

matrix Inv(matrix Mt)
{
    if (Mt.row == Mt.column && Mt.row == Rank(Mt))
    {
        matrix TEMP;
        TEMP.row = Mt.row;
        TEMP.column = Mt.column * 2;
        TEMP.body.clear();

        for (int i = 0; i < Mt.row; i++)
        {
            for (int j = 0; j < Mt.column * 2; j++)
            {
                if ((j - i) == Mt.row)
                {
                    TEMP.body.push_back(1);
                }
                else
                {
                    TEMP.body.push_back(0);
                }
            }
        }

        for (int i = 0; i < Mt.row; i++)
        {
            for (int j = 0; j < Mt.column; j++)
            {
                TEMP(i, j) = Mt(i, j);
            }
        }

        matrix result = Mt;
        matrix combi = Echelon(TEMP);

        for (int i = 0; i < Mt.row; i++)
        {
            for (int j = Mt.column; j < Mt.column * 2; j++)
            {
                result(i, j - Mt.column) = (combi)(i, j);
            }
        }
        result.name = "(" + Mt.name + "^-1" + ")";
        return result;
    }
    else
    {
        matrix fail;
        fail.Clear();
        cout << "Wrong with error:-1\n";
        return fail;
    }
}

double determinant(double matrix[MAXO][MAXO], int order)
{
    int sign = 1, i;
    double result = 0, cofactor;
    if (order == 1)
        result = matrix[0][0];
    else
        for (i = 0; i < order; i++)
        {
            cofactor = laplace_expansion(matrix, i, 0, order);
            result += sign * matrix[i][0] * cofactor;
            sign *= -1;
        }

    return result;
}

double laplace_expansion(double matrix[MAXO][MAXO], int r, int c, int order)
{
    double result = 0, cofactor[MAXO][MAXO];
    int original_i, original_j, i, j;

    for (i = 0; i < order; i++)
        for (j = 0; j < order; j++)
        {
            original_i = i;
            original_j = j;
            if (i == r || j == c)
                ;
            else
            {
                if (i > r)
                    i--;
                if (j > c)
                    j--;
                cofactor[i][j] = matrix[original_i][original_j];
                i = original_i;
                j = original_j;
            }
        }
    if (order >= 2)
        result = determinant(cofactor, order - 1);

    return result;
}

double Det(matrix Mt)
{
    if (Mt.row == Mt.column)
    {
        double TEMP_matrix[MAXO][MAXO];
        for (int i = 0; i < Mt.row; i++)
        {
            for (int j = 0; j < Mt.column; j++)
            {
                TEMP_matrix[i][j] = Mt(i, j);
            }
        }
        Mt.det = determinant(TEMP_matrix, Mt.row);
        return Mt.det;
    }
    else
    {
        cout << "Wrong with error:-1\n";
        return EOF;
    }
}