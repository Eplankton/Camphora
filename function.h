
#include <iostream>
#include <cmath>

#define LIMIT 1e-005

using namespace std;

double derivative(double (*)(double), double);
double derivative(double (*fun)(double), double x0)
{
    return ((*fun)(x0 + LIMIT) - (*fun)(x0 - LIMIT)) / (2 * LIMIT);
}

double integral(double (*)(double), double, double);
double integral(double (*fun)(double), double lb, double ub)
{
    double ptr = lb;
    double result = 0;
    while (ub - ptr >= LIMIT * 0.1)
    {
        result += 0.5 * (LIMIT) * ((*fun)(ptr) + (*fun)(ptr + LIMIT));
        ptr += LIMIT;
    }
    return result;
}
