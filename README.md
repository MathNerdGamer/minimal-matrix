# Minimal Matrix

A header-only minimal matrix implementation which accepts user-defined numerical types. Includes unary +/-, addition, subtraction, multiplication assignment/binary operators, and comparison operators. This is minimal in the sense that only the most basic operations on a matrix are included.

# Usage
Here's a basic example:
```
#include "matrix_t.h"
#include <iostream>
#include <iomanip>

using namespace math_nerd::matrix_t; // For demonstration purposes

template<typename T>
std::ostream &operator<<(std::ostream &os, matrix_t<T> matrix);

int main()
{
    std::int64_t A_rows{ 4 }, A_columns{ 7 };
    std::int64_t B_rows{ 7 }, B_columns{ 6 };

    matrix_t<double> A{ A_rows, A_columns };
    matrix_t<double> B{ B_rows, B_columns };

    for( auto i = 0u; i < A_rows; ++i )
    {
        for( auto j = 0u; j < A_columns; ++j )
        {
            A[i][j] = 2.71828*i + j;
        }
    }

    for( auto i = 0u; i < B_rows; ++i )
    {
        for( auto j = 0u; j < B_columns; ++j )
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

template<typename T>
std::ostream &operator<<(std::ostream &os, matrix_t<T> matrix)
{
    for( auto i = 0u; i < matrix.row_count(); ++i )
    {
        for( auto j = 0u; j < matrix.column_count(); ++j )
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