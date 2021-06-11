/*
MIT License

Copyright (c) 2019 Math Nerd

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
*/
#pragma once
#ifndef MATH_NERD_MATRIX_T_H
#define MATH_NERD_MATRIX_T_H
#include <array>
#include <cstdint>
#include <iostream>
#include <sstream>
#include <type_traits>
#include <vector>

/** \file matrix_t.h
    \brief Minimal matrix implemetation which takes user-defined numerical types.
 */

 /** \namespace math_nerd
      \brief Namespace for all of my projects.
   */
namespace math_nerd
{
    /** \namespace math_nerd::matrix_t
        \brief Namespace for matrix_t implementation.
     */
    namespace matrix_t
    {
        /** \namespace math_nerd::matrix_t::impl_details
            \brief Contains implementation details.
         */
        namespace impl_details
        {
            /** \name Type Traits
             */
            // Modified from T.C.'s answer on Stack Overflow: https://stackoverflow.com/a/26434209
            template<class...> struct voidify { using type = void; };
            template<class... Ts> using void_t = typename voidify<Ts...>::type;

            template<class T, class = void>
            struct supports_arithmetic : std::false_type {};

            template<class T>
            struct supports_arithmetic<T,
                void_t<decltype(std::declval<T>() + std::declval<T>()),
                decltype(std::declval<T>() - std::declval<T>()),
                decltype(std::declval<T>() *std::declval<T>()),
                decltype(std::declval<T>() / std::declval<T>())>>
                : std::true_type {};

            // Added to emulate std::is_arithmetic_v<T> syntax.
            template<class T>
            inline constexpr bool supports_arithmetic_v = supports_arithmetic<T>::value;
        }
        

        /** \name Signed 64-bit integer
         */
        using s64 = std::int64_t;

        /** \class matrix_t<T>
            \brief A very minimal matrix implementation.
         */
        template<typename T>
        class matrix_t
        {
            // Check for arithmetic support.
            static_assert(impl_details::supports_arithmetic_v<T>, "Matrix requires a type which supports arithmetic operations (+, -, *, /).");

            private:
                std::vector<std::vector<T>> mat;
                s64 rows_{ 2 }, columns_{ 2 };

            public:
                constexpr matrix_t()
                {
                    std::vector<T> zero(2, 0);

                    for( auto i = 0; i < 2; ++i )
                    {
                        mat.push_back(zero);
                    }
                }

                constexpr matrix_t(s64 rows = 2, s64 columns = 2)
                {
                    if( rows < 1 )
                    {
                        rows = 1;
                    }

                    if( columns < 1 )
                    {
                        columns = rows;
                    }

                    rows_ = rows;
                    columns_ = columns;

                    std::vector<T> zero(static_cast<std::uint32_t>(columns_), 0);

                    for( auto i = 0; i < rows_; ++i )
                    {
                        mat.push_back(zero);
                    }
                }

                /** \name Array operators */
                /** \fn auto operator[](std::size_t pos) -> std::vector<T> &
                    \brief Returns a reference to the column array at the specified position.
                 */
                constexpr auto operator[](std::size_t pos) -> std::vector<T> &
                {
                    return mat[pos];
                }

                /** \fn auto operator[](std::size_t pos) const -> std::vector<T>
                    \brief Returns a copy to the column array at the specified position.
                 */
                constexpr auto operator[](std::size_t pos) const -> std::vector<T>
                {
                    return mat[pos];
                }

                /** \name Row and Column Count */
                /** \fn auto row_count() const noexcept -> s64
                    \brief Returns the number of rows.
                 */
                constexpr auto row_count() const noexcept -> s64
                {
                    return rows_;
                }

                /** \fn auto column_count() const noexcept -> s64
                    \brief Returns the number of columns.
                 */
                constexpr auto column_count() const noexcept -> s64
                {
                    return columns_;
                }

                /** \name Matrix Inverse */
                /** \fn auto inverse() const -> matrix_t<T>
                    \brief Returns the inverse matrix. Throws if not invertible. Defined by the user.
                 */
                auto inverse() const -> matrix_t<T>;

                /** \name Unary operators */
                /** \fn constexpr auto operator+() const noexcept -> matrix_t<T>
                    \brief Returns the matrix as-is.
                 */
                constexpr auto operator+() const noexcept -> matrix_t<T>;

                /** \fn constexpr auto operator-() const noexcept -> matrix_t<T>
                    \brief Returns the matrix with the signs of all elements flipped.
                 */
                constexpr auto operator-() const noexcept -> matrix_t<T>;

                /** \name Assignment operators */
                /** \fn constexpr auto operator+=(matrix_t<T> const &rhs) -> matrix_t<T> &
                    \brief Adds rhs to the matrix.
                 */
                constexpr auto operator+=(matrix_t<T> const rhs) -> matrix_t<T> &;

                /** \fn constexpr auto operator-=(matrix_t<T> const &rhs) -> matrix_t<T> &
                    \brief Subtracts rhs from the matrix.
                 */
                constexpr auto operator-=(matrix_t<T> const rhs) -> matrix_t<T> &;

                /** \name Comparison operators */
                /** \fn constexpr auto operator==(matrix_t<T> const rhs) const noexcept -> bool
                    \brief Compares the values and returns true if they are all equal.
                 */
                constexpr auto operator==(matrix_t<T> const &rhs) const noexcept -> bool;

                /** \fn constexpr auto operator!=(matrix_t<T> const rhs) const noexcept -> bool
                    \brief Compares the values and returns false if they are all equal.
                 */
                constexpr auto operator!=(matrix_t<T> const &rhs) const noexcept -> bool;
        };

        // Unary operators
        template<typename T>
        constexpr auto matrix_t<T>::operator+() const noexcept -> matrix_t<T>
        {
            return *this;
        }

        template<typename T>
        constexpr auto matrix_t<T>::operator-() const noexcept -> matrix_t<T>
        {
            matrix_t<T> new_mat{ rows_, columns_  };

            for( auto i = 0u; i < rows_; ++i )
            {
                for( auto j = 0u; j < columns_; ++j )
                {
                    new_mat[i][j] = -new_mat[i][j];
                }
            }

            return new_mat;
        }

        // Assignment operators
        template<typename T>
        constexpr auto matrix_t<T>::operator+=(matrix_t<T> const rhs) -> matrix_t<T> &
        {
            if( rows_ != rhs.row_count() || columns_ != rhs.column_count() )
            {
                std::stringstream error_stream;

                error_stream << "Cannot add matrices which are not the same dimension.\n"
                             << "Left-hand matrix dimensions: " << rows_ << " x " << columns_ << '\n'
                             << "Right-hand matrix dimensions: "
                             << rhs.row_count() << " x " << rhs.column_count() << '\n';

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            for( auto i = 0u; i < rows_; ++i )
            {
                for( auto j = 0u; j < columns_; ++j )
                {
                    mat[i][j] += rhs[i][j];
                }
            }

            return *this;
        }

        template<typename T>
        constexpr auto matrix_t<T>::operator-=(matrix_t<T> const rhs) -> matrix_t<T> &
        {
            if( rows_ != rhs.row_count() || columns_ != rhs.column_count() )
            {
                std::stringstream error_stream;

                error_stream << "Cannot subtract matrices which are not the same dimension.\n"
                             << "Left-hand matrix dimensions: " << rows_ << " x " << columns_ << '\n'
                             << "Right-hand matrix dimensions: "
                             << rhs.row_count() << " x " << rhs.column_count() << '\n';

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            for( auto i = 0u; i < rows_; ++i )
            {
                for( auto j = 0u; j < columns_; ++j )
                {
                    mat[i][j] -= rhs[i][j];
                }
            }

            return *this;
        }

        // Comparison operators
        template<typename T>
        constexpr auto matrix_t<T>::operator==(matrix_t<T> const &rhs) const noexcept -> bool
        {
            bool same = true;

            if( rows_ != rhs.row_count() || columns_ != rhs.column_count() )
            {
                same = false;
            }
            else
            {
                for( auto i = 0u; i < rows_; ++i )
                {
                    for( auto j = 0u; j < columns_; ++j )
                    {
                        if( mat[i][j] != rhs[i][j] )
                        {
                            same = false;
                        }
                    }
                }
            }

            return same;
        }

        template<typename T>
        constexpr auto matrix_t<T>::operator!=(matrix_t<T> const &rhs) const noexcept -> bool
        {
            return !(*this == rhs);
        }

        /** \name Binary operators */
        /** \fn auto operator+(matrix_t<T> lhs, matrix_t<T> const &rhs) noexcept -> matrix_t<T> &
            \brief Adds two matrices of the same dimension together.
         */
        template<typename T>
        constexpr auto operator+(matrix_t<T> lhs, matrix_t<T> const &rhs) -> matrix_t<T>
        {
            try
            {
                lhs += rhs;
            }
            catch( std::invalid_argument &e )
            {
                throw;
            }

            return lhs;
        }

        /** \fn auto operator-(matrix_t<T> lhs, matrix_t<T> const &rhs) noexcept -> matrix_t<T> &
            \brief Subtracts one matrix from another of the same dimension.
         */
        template<typename T>
        constexpr auto operator-(matrix_t<T> lhs, matrix_t<T> const &rhs) noexcept -> matrix_t<T> &
        {
            try
            {
                lhs -= rhs;
            }
            catch( std::invalid_argument &e )
            {
                throw;
            }

            return lhs;
        }

        /** \fn constexpr auto operator*(matrix_t<T> const &lhs, matrix_t<T> const &rhs) -> matrix_t<T>
            \brief Multiples rhs to matrix, throws if dimensions aren't compatible. Returns a matrix of possibly different dimensions.
         */
        template<typename T>
        constexpr auto operator*(matrix_t<T> const &lhs, matrix_t<T> const &rhs) -> matrix_t<T>
        {
            s64 left_rows = lhs.row_count(), left_columns = lhs.column_count();
            s64 right_rows = rhs.row_count(), right_columns = rhs.column_count();

            if( right_rows != left_columns )
            {
                std::stringstream error_stream;

                error_stream << "Cannot multiply matrices because the left-hand matrix has "
                             << left_columns << " columns, which does not equal the number of rows of the right-hand matrix, "
                             << right_rows << ".\n";

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            matrix_t<T> new_mat{ left_rows, right_columns };

            for( auto i = 0u; i < left_rows; ++i )
            {
                for( auto j = 0u; j < right_columns; ++j )
                {
                    new_mat[i][j] = 0;
                }
            }

            for( auto i = 0u; i < left_rows; ++i )
            {
                for( auto j = 0u; j < right_columns; ++j )
                {
                    for( auto k = 0u; k < left_columns; ++k )
                    {
                        new_mat[i][j] += lhs[i][k] * rhs[k][j];
                    }
                }
            }

            return new_mat;
        }
        
        /** \fn constexpr auto operator*(matrix_t<T> const &lhs, std::array<T, len> const &rhs) -> matrix_t<T>
            \brief Matrix-Vector multiplication with std::array.
         */
        template<typename T, std::size_t len>
        constexpr auto operator*(matrix_t<T> const &lhs, std::array<T, len> const &rhs) -> matrix_t<T>
        {
            s64 left_rows = lhs.row_count(), left_columns = lhs.column_count();
            s64 right_rows = len, right_columns = 1;

            if( right_rows != left_columns )
            {
                std::stringstream error_stream;

                error_stream << "Cannot multiply matrices because the left-hand matrix has "
                             << left_columns << " columns, which does not equal the number of rows of the right-hand matrix, "
                             << right_rows << ".\n";

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            matrix_t<T> new_mat{ left_rows, right_columns };

            for( auto i = 0u; i < left_rows; ++i )
            {
                new_mat[i][0] = 0;
            }

            for( auto i = 0u; i < left_rows; ++i )
            {
                for( auto j = 0u; j < left_columns; ++j )
                {
                    new_mat[i][0] += lhs[i][j] * rhs[j];
                }
            }

            return new_mat;
        }
        
        /** \fn constexpr auto operator*(matrix_t<T> const &lhs, std::vector<T> const &rhs) -> matrix_t<T>
            \brief Matrix-Vector multiplication with std::vector.
         */
        template<typename T>
        constexpr auto operator*(matrix_t<T> const &lhs, std::vector<T> const &rhs) -> matrix_t<T>
        {
            s64 left_rows = lhs.row_count(), left_columns = lhs.column_count();
            s64 right_rows = rhs.size(), right_columns = 1;

            if( right_rows != left_columns )
            {
                std::stringstream error_stream;

                error_stream << "Cannot multiply matrices because the left-hand matrix has "
                             << left_columns << " columns, which does not equal the number of rows of the right-hand matrix, "
                             << right_rows << ".\n";

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            matrix_t<T> new_mat{ left_rows, right_columns };

            for( auto i = 0u; i < left_rows; ++i )
            {
                new_mat[i][0] = 0;
            }

            for( auto i = 0u; i < left_rows; ++i )
            {
                for( auto j = 0u; j < left_columns; ++j )
                {
                    new_mat[i][0] += lhs[i][j] * rhs[j];
                }
            }

            return new_mat;
        }

        /** \name Scalar multiplication operators */
        /** \fn constexpr auto operator*(T const &lhs, matrix_t<T> const &rhs) -> matrix_t<T>
            \brief Scales matrix rhs by lhs.
         */
        template<typename T>
        constexpr auto operator*(T const &lhs, matrix_t<T> &rhs) -> matrix_t<T>
        {
            for( auto i = 0u; i < rhs.row_count(); ++i )
            {
                for( auto j = 0u; j < rhs.column_count(); ++j )
                {
                    rhs[i][j] *= lhs;
                }
            }

            return rhs;
        }

        /** \fn constexpr auto operator*(matrix_t<T> const &lhs, T const &rhs) -> matrix_t<T>
            \brief Scales matrix lhs by rhs.
         */
        template<typename T>
        constexpr auto operator*(matrix_t<T> const &lhs, T const &rhs) -> matrix_t<T>
        {
            return rhs * lhs;
        }
        
    } // namespace matrix_t
} // namespace math_nerd

/** \mainpage Minimal Matrix Implemntation
    \section gitlab_link GitLab Link
    View the source code at <a href="https://gitlab.com/mathnerd/minimal-matrix">GitLab</a>.
 */
#endif

