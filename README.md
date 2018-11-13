# Minimal Matrix

A header-only minimal matrix implementation which accepts user-defined numerical types. Includes unary +/-, addition, subtraction, multiplication assignment/binary operators, and comparison operators.

# Usage
Here's a basic example:
```
#include "matrix_t.h"
#include <iostream>
#include <iomanip>

using namespace math_nerd::matrix_t;

template<typename T, std::size_t R, std::size_t C>
std::ostream &operator<<(std::ostream &os, matrix_t<T, R, C> matrix);

int main()
{
    matrix_t<double, 3, 7> A;

    for( auto i = 0u; i < 3; ++i )
    {
        for( auto j = 0u; j < 7; ++j )
        {
            A[i][j] = 2.71828*i + 3.14159*j;
        }
    }

    std::cout << A << '\n';

    return 0;
}

template<typename T, std::size_t R, std::size_t C>
std::ostream &operator<<(std::ostream &os, matrix_t<T, R, C> matrix)
{
    for( auto i = 0u; i < R; ++i )
    {
        for( auto j = 0u; j < C; ++j )
        {
            os << std::setw(10) << matrix[i][j] << ' ';
        }
        os << '\n';
    }

    return os;
}
```

Output:
```
         0    3.14159    6.28318    9.42477    12.5664     15.708    18.8495
   2.71828    5.85987    9.00146     12.143    15.2846    18.4262    21.5678
   5.43656    8.57815    11.7197    14.8613    18.0029    21.1445    24.2861
```