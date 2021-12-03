#include <iostream>
#include "matrix.h"
#include "function.h"

int main()
{
    matrix A;
    A.Init("A", 3, 3) = {0, 0, 0, 0, 0, 0, 0, 0, 0};
    Echelon(A ^ -1).Out('f');
    cout << (A ^ -1).Rank();
    return 0;
}
