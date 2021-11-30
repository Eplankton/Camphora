#include <iostream>
#include "matrix.h"
#include "function.h"

int main()
{
    matrix A;
    A.Init("hello_world!", 2, 2) = {1, 2, 3, 4};
    A.Out();
    return 0;
}
