# Minimal Matrix

A header-only minimal matrix implementation which accepts user-defined numerical types. Includes unary +/-, addition, subtraction, multiplication assignment/binary operators, and comparison operators.

# Usage
Here's a basic example:
```
#include "matrix_t.h"
#include <iostream>
#include <iomanip>

using namespace math_nerd::matrix_t; // For demonstration purposes

template<typename T, std::size_t R, std::size_t C>
std::ostream &operator<<(std::ostream &os, matrix_t<T, R, C> matrix);

int main()
{
    matrix_t<double, 4, 7> A;
    matrix_t<double, 7, 6> B;

    for( auto i = 0u; i < 4; ++i )
    {
        for( auto j = 0u; j < 7; ++j )
        {
            A[i][j] = 2.71828*i + j;
        }
    }

    for( auto i = 0u; i < 7; ++i )
    {
        for( auto j = 0u; j < 6; ++j )
        {
            B[i][j] = i + 3.14159*j;
        }
    }

    try
    {
        std::cout << A * B << '\n';
    }
    catch( std::invalid_argument const &e )
    {
        std::cout << e.what();
        return EXIT_FAILURE;
    }
    

    return EXIT_SUCCESS;
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
        91    156.973    222.947     288.92    354.894    420.867
   148.084    273.835    399.587    525.338     651.09    776.841
   205.168    390.697    576.227    761.756    947.286    1132.82
   262.252    507.559    752.867    998.174    1243.48    1488.79
```