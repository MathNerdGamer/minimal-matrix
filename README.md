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
        55    102.124    149.248    196.372    243.495    290.619    337.743    384.867    431.991
   95.7742    194.136    292.499    390.861    489.223    587.585    685.947    784.309    882.672
   136.548    286.149    435.749     585.35     734.95    884.551    1034.15    1183.75    1333.35
   177.323    378.161        579    779.839    980.678    1181.52    1382.36    1583.19    1784.03
   218.097    470.174    722.251    974.328    1226.41    1478.48    1730.56    1982.64    2234.71
   258.871    562.186    865.502    1168.82    1472.13    1775.45    2078.76    2382.08    2685.39
   299.645    654.199    1008.75    1363.31    1717.86    2072.41    2426.97    2781.52    3136.08



        55
   95.7742
   136.548
   177.323
   218.097
   258.871
   299.645
```