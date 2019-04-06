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
#include <vector>
#include <cstdint>
#include <sstream>
#include <iostream>

/** \file matrix_t.h
    \brief Minimal matrix implemetation which takes user-defined numerical types.
 */

namespace math_nerd
{
    namespace matrix_t
    {
        using s64 = std::int64_t;

        template<typename T>
        class matrix_t
        {
            private:
                std::vector<std::vector<T>> mat;
                s64 r = 2, c = 2;

            public:
                matrix_t(s64 rows = 2, s64 columns = 0)
                {
                    if( rows < 1 )
                    {
                        rows = 1;
                    }

                    if( columns < 1 )
                    {
                        columns = rows;
                    }

                    r = rows;
                    c = columns;

                    std::vector<T> zero;

                    for( auto i = 0; i < c; ++i )
                    {
                        zero.push_back(0);
                    }

                    for( auto i = 0; i < r; ++i )
                    {
                        mat.push_back(zero);
                    }
                }

                /** \name Array operators */
                /** \fn std::vector<T> &operator[](std::size_t pos)
                    \brief Returns a reference to the column array at the specified position.
                 */
                std::vector<T> &operator[](std::size_t pos)
                {
                    return mat[pos];
                }

                /** \fn std::vector<T> operator[](std::size_t pos) const
                    \brief Returns a copy to the column array at the specified position.
                 */
                std::vector<T> operator[](std::size_t pos) const
                {
                    return mat[pos];
                }

                /** \fn s64 row_count() const noexcept
                    \brief Returns the number of rows.
                 */
                s64 row_count() const noexcept
                {
                    return r;
                }

                /** \fn s64 column_count() const noexcept
                    \brief Returns the number of columns.
                 */
                s64 column_count() const noexcept
                {
                    return c;
                }

                /** \fn matrix_t<T> inverse() const;
                    \brief Returns the inverse matrix. Throws if not invertible. Defined by the user.
                 */
                matrix_t<T> inverse() const;

                /** \name Unary operators */
                /** \fn matrix_t<T> operator+() const noexcept
                    \brief Returns the matrix as-is.
                 */
                matrix_t<T> operator+() const noexcept;

                /** \fn matrix_t<T> operator-() const noexcept
                    \brief Returns the matrix with the signs of all elements flipped.
                 */
                matrix_t<T> operator-() const noexcept;

                /** \name Assignment operators */
                /** \fn matrix_t<T> &operator+=(matrix_t<T> const &rhs)
                    \brief Adds rhs to the matrix.
                 */
                matrix_t<T> &operator+=(matrix_t<T> const rhs);

                /** \fn matrix_t<T> &operator-=(matrix_t<T> const &rhs)
                    \brief Subtracts rhs from the matrix.
                 */
                matrix_t<T> &operator-=(matrix_t<T> const rhs);

                /** \name Comparison operators */
                /** \fn constexpr bool operator==(matrix_t<T> const rhs) const noexcept
                    \brief Compares the values and returns true if they are all equal.
                 */
                constexpr bool operator==(matrix_t<T> const &rhs) const noexcept;

                /** \fn constexpr bool operator!=(matrix_t<T> const rhs) const noexcept
                    \brief Compares the values and returns false if they are all equal.
                 */
                constexpr bool operator!=(matrix_t<T> const &rhs) const noexcept;
        };

        // Unary operators
        template<typename T>
        matrix_t<T> matrix_t<T>::operator+() const noexcept
        {
            return *this;
        }

        template<typename T>
        matrix_t<T> matrix_t<T>::operator-() const noexcept
        {
            s64 R = this->row_count, C = this->column_count();
            matrix_t<T> new_mat{ R, C };

            for( auto i = 0u; i < R; ++i )
            {
                for( auto j = 0u; j < C; ++j )
                {
                    new_mat[i][j] = -new_mat[i][j];
                }
            }

            return new_mat;
        }

        // Assignment operators
        template<typename T>
        matrix_t<T> &matrix_t<T>::operator+=(matrix_t<T> const rhs)
        {
            if( r != rhs.row_count() || c != rhs.column_count() )
            {
                std::stringstream error_stream;

                error_stream << "Cannot add matrices which are not the same dimension.\n"
                             << "Left-hand matrix dimensions: " << r << " x " << c << '\n'
                             << "Right-hand matrix dimensions: "
                             << rhs.row_count() << " x " << rhs.column_count() << '\n';

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            for( auto i = 0u; i < this->row_count(); ++i )
            {
                for( auto j = 0u; j < this->column_count(); ++j )
                {
                    mat[i][j] += rhs[i][j];
                }
            }

            return *this;
        }

        template<typename T>
        matrix_t<T> &matrix_t<T>::operator-=(matrix_t<T> const rhs)
        {
            if( r != rhs.row_count() || c != rhs.column_count() )
            {
                std::stringstream error_stream;

                error_stream << "Cannot subtract matrices which are not the same dimension.\n"
                             << "Left-hand matrix dimensions: " << r << " x " << c << '\n'
                             << "Right-hand matrix dimensions: "
                             << rhs.row_count() << " x " << rhs.column_count() << '\n';

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            for( auto i = 0u; i < this->row_count(); ++i )
            {
                for( auto j = 0u; j < this->column_count(); ++j )
                {
                    mat[i][j] -= rhs[i][j];
                }
            }

            return *this;
        }

        // Comparison operators
        template<typename T>
        constexpr bool matrix_t<T>::operator==(matrix_t<T> const &rhs) const noexcept
        {
            bool same = true;

            if( r != rhs.row_count() || c != rhs.column_count() )
            {
                same = false;
            }
            else
            {
                for( auto i = 0u; i < this->row_count(); ++i )
                {
                    for( auto j = 0u; j < this->column_count(); ++j )
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
        constexpr bool matrix_t<T>::operator!=(matrix_t<T> const &rhs) const noexcept
        {
            return !(*this == rhs);
        }

        // Binary operators
        /** \name Binary operators */
        /** \fn matrix_t<T> &operator+(matrix_t<T> lhs, matrix_t<T> const &rhs) noexcept
            \brief Adds two matrices of the same dimension together.
         */
        template<typename T>
        matrix_t<T> operator+(matrix_t<T> lhs, matrix_t<T> const &rhs)
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

        /** \fn matrix_t<T> &operator-(matrix_t<T> lhs, matrix_t<T> const &rhs) noexcept
            \brief Subtracts one matrix from another of the same dimension.
         */
        template<typename T>
        matrix_t<T> &operator-(matrix_t<T> lhs, matrix_t<T> const &rhs) noexcept
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

        /** \fn constexpr matrix_t<T> operator*(matrix_t<T> const &lhs, matrix_t<T> const &rhs)
            \brief Multiples rhs to matrix, throws if dimensions aren't compatible. Returns a matrix of possibly different dimensions.
         */
        template<typename T>
        constexpr matrix_t<T> operator*(matrix_t<T> const &lhs, matrix_t<T> const &rhs)
        {
            s64 R = lhs.row_count(), C = lhs.column_count();
            s64 R2 = rhs.row_count(), C2 = rhs.column_count();

            if( R2 != C )
            {
                std::stringstream error_stream;

                error_stream << "Cannot multiply matrices because the left-hand matrix has "
                             << C << " columns, which does not equal the number of rows of the right-hand matrix, "
                             << R2 << ".\n";

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            matrix_t<T> new_mat{ R, C2 };

            for( auto i = 0u; i < R; ++i )
            {
                for( auto j = 0u; j < C2; ++j )
                {
                    new_mat[i][j] = 0;
                }
            }

            for( auto i = 0u; i < R; ++i )
            {
                for( auto j = 0u; j < C2; ++j )
                {
                    for( auto k = 0u; k < C; ++k )
                    {
                        new_mat[i][j] += lhs[i][k] * rhs[k][j];
                    }
                }
            }

            return new_mat;
        }

        // Scalar multiplication operator
        /** \name Scalar multiplication operators */
        /** \fn constexpr matrix_t<T> operator*(T const &lhs, matrix_t<T> const &rhs)
            \brief Scales matrix rhs by lhs.
         */
        template<typename T>
        constexpr matrix_t<T> operator*(T const &lhs, matrix_t<T> const &rhs)
        {
            auto new_mat = rhs;

            for( auto i = 0u; i < rhs.row_count(); ++i )
            {
                for( auto j = 0u; j < rhs.column_count(); ++j )
                {
                    new_mat[i][j] *= lhs;
                }
            }

            return new_mat;
        }

        /** \fn constexpr matrix_t<T> operator*(matrix_t<T, R2, C2> const &lhs, T const &rhs)
            \brief Scales matrix lhs by rhs.
         */
        template<typename T>
        constexpr matrix_t<T> operator*(matrix_t<T> const &lhs, T const &rhs)
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
