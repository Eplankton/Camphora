# Camphora



## #Download

You can download the header files and a test example here:
https://gitee.com/Eplankton/Camphora.git

    $ git clone https://gitee.com/Eplankton/Camphora.git



## #How to use 

You can create a matrix by:

```C++
#include <iostream>
#include "matrix.h"

int main()
{
	matrix A;
	return 0;
}

```

And initialize it by:

```C++
#include <iostream>
#include "matrix.h"

int main()
{
	matrix A;
    A.Init("A", 2, 2) = {1, 0, 0, 1};
	return 0;
}
```

Print it on your screen:

```C++
#include <iostream>
#include "matrix.h"

int main()
{
	matrix A;
    A.Init("A", 2, 2) = {1, 0, 0, 1};
    A.Out('f'); //Print in format
	return 0;
}
```

```
A <2,2>
[ 1  0 ]
[ 0  1 ]
```



You may calculate matrices by:

```C++
#include <iostream>
#include "matrix.h"

int main()
{
	matrix A;
    A.Init("A", 2, 2) = {1, 0, 0, 1};
    matrix B = A; //To initialize by another existed matrix
    B.name = "B";
    
    (A * (A + B)).Out(); //Do some calculation
    cout << Det(A) << endl;
    cout << Rank(A) << endl;  
	return 0;
}

```

```
(A*(A+B)) <2,2>
[ 2  0 ]
[ 0  2 ]

3.30237e-317   //Approximately: 0
2
```





