# Minimal Matrix

A header-only minimal matrix implementation which accepts user-defined numerical types. Includes unary +/-, addition, subtraction, multiplication assignment/binary operators, and comparison operators. This is minimal in the sense that only the most basic operations on a matrix are included.

# Usage
Here's a basic example:
```
#include "matrix_t.h"
#include <iostream>
#include <iomanip>

using namespace math_nerd::matrix_t;

template<typename T>
std::ostream &operator<<(std::ostream &os, matrix_t<T> matrix);

int main(int argc, char **argv)
{
    constexpr std::int64_t A_rows{ 7 }, A_cols{ 6 }, B_rows{ 6 }, B_cols{ 9 };

    matrix_t<double> A{ A_rows, A_cols };
    matrix_t<double> B{ B_rows, B_cols };

    for( auto i = 0u; i < A_rows; ++i )
    {
        for( auto j = 0u; j < A_cols; ++j )
        {
            A[i][j] = 2.71828*i + j;
        }
    }

    for( auto i = 0u; i < B_rows; ++i )
    {
        for( auto j = 0u; j < B_cols; ++j )
        {
            B[i][j] = i + 3.14159*j;
        }
    }

    std::array<double, B_cols> first_only;

    first_only[0] = 1.0;

    for( auto i = 1u; i < B_cols; ++i )
    {
        first_only[i] = 0.0;
    }

    try
    {
        auto C = A * B;
        std::cout << C << '\n';

        std::cout << '\n' << '\n';

        auto D = C * first_only; // Vector with only first element of each row.
        std::cout << D << '\n';
    }
    catch( std::invalid_argument &e )
    {
        std::cout << e.what() << '\n';
    }

    return 0;
}

template<typename T>
std::ostream &operator<<(std::ostream &os, matrix_t<T> const matrix)
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