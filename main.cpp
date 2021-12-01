#include <iostream>
#include "matrix.h"
#include "function.h"

int main()
{
    matrix A;
    A.Init("A", 2, 2) = {1, 0, 0, 1};
    A.Out();
    return 0;
}
