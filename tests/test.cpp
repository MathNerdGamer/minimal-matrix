#include <sstream>

#include <math_nerd/matrix_t.h>


#define CATCH_DEFINE_MAIN
#include "catch.hpp"

namespace mt = math_nerd::matrix_t;

TEST_CASE("Testing Constructor")
{
    mt::matrix_t<int>   test_subject1;
    REQUIRE(test_subject1.row_count() == 2);
    REQUIRE(test_subject1.column_count() == 2);
    
    for( auto i{ 0 }; i < test_subject1.row_count(); ++i )
    {
        for( auto j{ 0 }; j < test_subject1.column_count(); ++j )
        {
            REQUIRE(test_subject1[i][j] == 0);
        }
    }

    mt::matrix_t<float> test_subject2{ 1, 5 };
    REQUIRE(test_subject2.row_count() == 1);
    REQUIRE(test_subject2.column_count() == 5);

    for( auto i{ 0u }; i < test_subject2.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject2.column_count(); ++j )
        {
            REQUIRE(test_subject2[i][j] == 0);
        }
    }

    mt::matrix_t<double> test_subject3{ 3, 1 };
    REQUIRE(test_subject3.row_count() == 3);
    REQUIRE(test_subject3.column_count() == 1);

    for( auto i{ 0u }; i < test_subject3.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject3.column_count(); ++j )
        {
            REQUIRE(test_subject3[i][j] == 0);
        }
    }
}

TEST_CASE("Testing Comparison Operators")
{
    mt::matrix_t<double> test_subject1{ 7, 6 };
    mt::matrix_t<double> test_subject2{ 7, 6 };
    mt::matrix_t<double> test_subject3{ 7, 6 };
    mt::matrix_t<double> test_subject4{ 6, 9 };

    for( auto i{ 0u }; i < test_subject1.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject1.column_count(); ++j )
        {
            test_subject1[i][j] = 2.71828 * i + j;
        }
    }


    for( auto i{ 0u }; i < test_subject2.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject2.column_count(); ++j )
        {
            test_subject2[i][j] = 2.71828 * i + j;
        }
    }

    REQUIRE(test_subject1 == test_subject2);

    for( auto i{ 0u }; i < test_subject3.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject3.column_count(); ++j )
        {
            test_subject3[i][j] = i + 3.14159 * j;
        }
    }

    REQUIRE(test_subject1 != test_subject3);

    for( auto i{ 0u }; i < test_subject4.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject4.column_count(); ++j )
        {
            test_subject4[i][j] = i + 3.14159 * j;
        }
    }

    REQUIRE(test_subject1 != test_subject4);
}

TEST_CASE("Testing Unary Operators")
{
    mt::matrix_t<double> test_subject{ 7, 6 };

    for( auto i{ 0u }; i < test_subject.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject.column_count(); ++j )
        {
            test_subject[i][j] = 2.71828 * i + j;
        }
    }

    SECTION("Positive")
    {
        REQUIRE(test_subject == +test_subject);
    }

    SECTION("Negative")
    {
        REQUIRE(test_subject != -test_subject);

        for( auto i{ 0u }; i < test_subject.row_count(); ++i )
        {
            for( auto j{ 0u }; j < test_subject.column_count(); ++j )
            {
                test_subject[i][j] *= -1;
            }
        }

        REQUIRE(test_subject == test_subject);
    }
}


TEST_CASE("Testing Assignment Operators")
{
    mt::matrix_t<double> test_subject1{ 7, 6 };
    mt::matrix_t<double> test_subject2{ 7, 6 };
    mt::matrix_t<double> test_subject3{ 6, 9 };
    mt::matrix_t<double> temp_subject { 7, 6 };

    for( auto i{ 0u }; i < test_subject1.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject1.column_count(); ++j )
        {
            test_subject1[i][j] = 2.71828 * i + j;
        }
    }

    for( auto i{ 0u }; i < test_subject3.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject3.column_count(); ++j )
        {
            test_subject3[i][j] = i + 3.14159 * j;
        }
    }
    
    SECTION("Assignment")
    {
        test_subject2 = test_subject1;
        REQUIRE(test_subject1 == test_subject2);

        try
        {
            test_subject1 = test_subject3;
        }
        catch( std::invalid_argument const &e )
        {
            REQUIRE(std::string(e.what()) == "Cannot assign matrices which are not the same dimension.\nLeft-hand matrix dimensions: 7 x 6\nRight-hand matrix dimensions: 6 x 9\n");
        }
    }

    for( auto i{ 0u }; i < test_subject2.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject2.column_count(); ++j )
        {
            test_subject2[i][j] = i + 3.14159 * j;
        }
    }

    SECTION("Addition Assignment")
    {
        for( auto i{ 0u }; i < temp_subject.row_count(); ++i )
        {
            for( auto j{ 0u }; j < temp_subject.column_count(); ++j )
            {
                temp_subject[i][j] = test_subject1[i][j] + test_subject2[i][j];
            }
        }

        test_subject1 += test_subject2;
        
        REQUIRE(test_subject1 == temp_subject);

        try
        {
            test_subject1 += test_subject3;
        }
        catch( std::invalid_argument const &e )
        {
            REQUIRE(std::string(e.what()) == "Cannot add matrices which are not the same dimension.\nLeft-hand matrix dimensions: 7 x 6\nRight-hand matrix dimensions: 6 x 9\n");
        }
    }

    SECTION("Subtraction Assignment")
    {
        for( auto i{ 0u }; i < temp_subject.row_count(); ++i )
        {
            for( auto j{ 0u }; j < temp_subject.column_count(); ++j )
            {
                temp_subject[i][j] = test_subject1[i][j] - test_subject2[i][j];
            }
        }

        test_subject1 -= test_subject2;

        REQUIRE(test_subject1 == temp_subject);

        try
        {
            test_subject1 -= test_subject3;
        }
        catch( std::invalid_argument const &e )
        {
            REQUIRE(std::string(e.what()) == "Cannot subtract matrices which are not the same dimension.\nLeft-hand matrix dimensions: 7 x 6\nRight-hand matrix dimensions: 6 x 9\n");
        }
    }
}

TEST_CASE("Testing Binary Operators")
{
    mt::matrix_t<double>  test_subject1{ 7, 6 };
    mt::matrix_t<double>  test_subject2{ 7, 6 };
    mt::matrix_t<double>  test_subject3{ 6, 9 };
    std::vector<double>   test_subject4(6, 0);
    std::array<double, 6> test_subject5;

    mt::matrix_t<double>  temp_subject1{ 7, 6 };
    mt::matrix_t<double>  temp_subject2{ 7, 9 };
    mt::matrix_t<double>  temp_subject3{ 7, 1 };

    for( auto i{ 0u }; i < test_subject1.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject1.column_count(); ++j )
        {
            test_subject1[i][j] = 2.71828 * i + j;
        }
    }

    for( auto i{ 0u }; i < test_subject2.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject2.column_count(); ++j )
        {
            test_subject2[i][j] = i + 3.14159 * j;
        }
    }

    SECTION("Addition")
    {
        for( auto i{ 0u }; i < temp_subject1.row_count(); ++i )
        {
            for( auto j{ 0u }; j < temp_subject1.column_count(); ++j )
            {
                temp_subject1[i][j] = test_subject1[i][j] + test_subject2[i][j];
            }
        }

        REQUIRE(test_subject1 + test_subject2 == temp_subject1);

        try
        {
            test_subject1 + test_subject3;
        }
        catch( std::invalid_argument const &e )
        {
            REQUIRE(std::string(e.what()) == "Cannot add matrices which are not the same dimension.\nLeft-hand matrix dimensions: 7 x 6\nRight-hand matrix dimensions: 6 x 9\n");
        }
    }

    SECTION("Subtraction")
    {
        for( auto i{ 0u }; i < temp_subject1.row_count(); ++i )
        {
            for( auto j{ 0u }; j < temp_subject1.column_count(); ++j )
            {
                temp_subject1[i][j] = test_subject1[i][j] - test_subject2[i][j];
            }
        }

        REQUIRE(test_subject1 - test_subject2 == temp_subject1);

        try
        {
            test_subject1 - test_subject3;
        }
        catch( std::invalid_argument const &e )
        {
            REQUIRE(std::string(e.what()) == "Cannot subtract matrices which are not the same dimension.\nLeft-hand matrix dimensions: 7 x 6\nRight-hand matrix dimensions: 6 x 9\n");
        }
    }

    SECTION("Multiplication")
    {
        try
        {
            test_subject1 * test_subject2;
        }
        catch( std::invalid_argument const &e )
        {
            REQUIRE(std::string(e.what()) == "Cannot multiply matrices because the left-hand matrix has 6 columns, which does not equal the number of rows of the right-hand matrix, 7.\n");
        }

        for( auto i{ 0u }; i < temp_subject2.row_count(); ++i )
        {
            for( auto j{ 0u }; j < temp_subject2.column_count(); ++j )
            {
                for( auto k{ 0u }; k < test_subject1.column_count(); ++k )
                {
                    temp_subject2[i][j] += test_subject1[i][k] * test_subject3[k][j];
                }
            }
        }

        REQUIRE(test_subject1 * test_subject3 == temp_subject2);

        for( auto i{ 0u }; i < 6; ++i )
        {
            test_subject4[i] = 1.0;
            test_subject5[i] = 1.0;

            for( auto j{ 0u }; j < test_subject1.row_count(); ++j )
            {
                temp_subject3[j][0] += 1.0 * test_subject1[j][i];
            }
        }

        REQUIRE(test_subject1 * test_subject4 == temp_subject3);
        REQUIRE(test_subject1 * test_subject5 == temp_subject3);
    }
}

TEST_CASE("Testing Scalar Multiplication Operators")
{
    mt::matrix_t<double> test_subject{ 7, 6 };

    for( auto i{ 0u }; i < test_subject.row_count(); ++i )
    {
        for( auto j{ 0u }; j < test_subject.column_count(); ++j )
        {
            test_subject[i][j] = 2.71828 * i + j;
        }
    }

    auto temp_subject{ test_subject };

    for( auto i{ 0u }; i < temp_subject.row_count(); ++i )
    {
        for( auto j{ 0u }; j < temp_subject.column_count(); ++j )
        {
            temp_subject[i][j] *= 1.7;
        }
    }

    SECTION("Left Multiplication by Scalar")
    {
        REQUIRE(1.7 * test_subject == temp_subject);
    }

    SECTION("Right Multiplication by Scalar")
    {
        REQUIRE(test_subject * 1.7 == temp_subject);
    }
}

