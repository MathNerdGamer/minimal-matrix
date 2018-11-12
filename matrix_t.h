/*
MIT License

Copyright (c) 2018 Math Nerd

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
#ifndef MATRIX_T_H
#define MATRIX_T_H
#include <array>
#include <cstdint>
#include <sstream>

/** \file matrix_t.h
    \brief Minimal matrix implemetation which takes user-defined numerical types.
 */

namespace math_nerd
{
    namespace matrix_t
    {
        using std::size_t;

        template<typename T, size_t R, size_t C>
        class matrix_t
        {
            private:
                static_assert(R > 0 && C > 0, "The row and column count must be positive.\n");
                std::array<std::array<T, C>, R> mat{ 0 };

            public:
                matrix_t() = default;

                /** \name Array operators */
                /** \fn std::array<T, C> &operator[](std::size_t pos)
                    \brief Returns a reference to the column array at the specified position.
                 */
                std::array<T, C> &operator[](std::size_t pos)
                {
                    return mat[pos];
                }

                /** \fn std::array<T, C> operator[](std::size_t pos) const
                    \brief Returns a copy to the column array at the specified position.
                 */
                std::array<T, C> operator[](std::size_t pos) const
                {
                    return mat[pos];
                }

                /** \name Unary operators */
                /** \fn matrix_t<T, R, C> operator+() const noexcept
                    \brief Returns the matrix as-is.
                 */
                matrix_t<T, R, C> operator+() const noexcept;

                /** \fn matrix_t<T, R, C> operator-() const noexcept
                    \brief Returns the matrix with the signs of all elements flipped.
                 */
                matrix_t<T, R, C> operator-() const noexcept;

                /** \name Assignment operators */
                /** \fn matrix_t<T, R, C> &operator=(matrix_t<T, R, C> const rhs) noexcept
                    \brief Assigns rhs to the matrix.
                 */
                matrix_t<T, R, C> &operator=(matrix_t<T, R, C> const rhs) noexcept;

                /** \fn matrix_t<T,R,C> &operator+=(matrix_t<T,R,C> const &rhs) noexcept
                    \brief Adds rhs to the matrix.
                 */
                matrix_t<T, R, C> &operator+=(matrix_t<T, R, C> const rhs) noexcept;

                /** \fn matrix_t<T, R, C> &operator-=(matrix_t<T, R, C> const &rhs) noexcept
                    \brief Subtracts rhs from the matrix.
                 */
                matrix_t<T, R, C> &operator-=(matrix_t<T, R, C> const rhs) noexcept;

                /** \name Comparison operators */
                /** \fn constexpr bool operator==(matrix_t<T,R,C> const rhs) const noexcept
                    \brief Compares the values and returns true if they are all equal.
                 */
                constexpr bool operator==(matrix_t<T, R, C> const rhs) const noexcept;

                /** \fn constexpr bool operator!=(matrix_t<T,R,C> const rhs) const noexcept
                    \brief Compares the values and returns false if they are all equal.
                 */
                constexpr bool operator!=(matrix_t<T, R, C> const rhs) const noexcept;
        };

        // Unary operators
        template<typename T, size_t R, size_t C>
        matrix_t<T, R, C> matrix_t<T, R, C>::operator+() const noexcept
        {
            return *this;
        }

        template<typename T, size_t R, size_t C>
        matrix_t<T, R, C> matrix_t<T, R, C>::operator-() const noexcept
        {
            matrix_t<T, R, C> new_mat;

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
        template<typename T, size_t R, size_t C>
        matrix_t<T, R, C> &matrix_t<T, R, C>::operator=(matrix_t<T, R, C> const rhs) noexcept
        {
            for( auto i = 0u; i < R; ++i )
            {
                for( auto j = 0u; j < C; ++i )
                {
                    mat[i][j] = rhs[i][j];
                }
            }

            return *this;
        }

        template<typename T, size_t R, size_t C>
        matrix_t<T, R, C> &matrix_t<T, R, C>::operator+=(matrix_t<T, R, C> const rhs) noexcept
        {
            for( auto i = 0u; i < R; ++i )
            {
                for( auto j = 0u; j < C; ++i )
                {
                    mat[i][j] += rhs[i][j];
                }
            }

            return *this;
        }

        template<typename T, size_t R, size_t C>
        matrix_t<T, R, C> &matrix_t<T, R, C>::operator-=(matrix_t<T, R, C> const rhs) noexcept
        {
            for( auto i = 0u; i < R; ++i )
            {
                for( auto j = 0u; j < C; ++i )
                {
                    mat[i][j] -= rhs[i][j];
                }
            }

            return *this;
        }

        // Comparison operators
        template<typename T, size_t R, size_t C>
        constexpr bool matrix_t<T, R, C>::operator==(matrix_t<T, R, C> const rhs) const noexcept
        {
            bool same = true;

            for( auto i = 0u; i < R; ++i )
            {
                for( auto j = 0u; j < C; ++j )
                {
                    if( mat[i][j] != rhs[i][j] )
                    {
                        same = false;
                    }
                }
            }

            return same;
        }

        template<typename T, size_t R, size_t C>
        constexpr bool matrix_t<T, R, C>::operator!=(matrix_t<T, R, C> const rhs) const noexcept
        {
            return !(*this == rhs);
        }

        // Binary operators
        /** \name Binary operators */
        /** \fn matrix_t<T, R, C> &operator+(matrix_t<T, R, C> lhs, matrix_t<T, R, C> const &rhs) noexcept
            \brief Adds two matrices of the same dimension together.
         */
        template<typename T, size_t R, size_t C>
        matrix_t<T, R, C> &operator+(matrix_t<T, R, C> lhs, matrix_t<T, R, C> const &rhs) noexcept
        {
            lhs += rhs;

            return lhs;
        }

        /** \fn matrix_t<T, R, C> &operator-(matrix_t<T, R, C> lhs, matrix_t<T, R, C> const &rhs) noexcept
            \brief Subtracts one matrix from another of the same dimension.
         */
        template<typename T, size_t R, size_t C>
        matrix_t<T, R, C> &operator-(matrix_t<T, R, C> lhs, matrix_t<T, R, C> const &rhs) noexcept
        {
            lhs -= rhs;

            return lhs;
        }

        /** \fn constexpr matrix_t<T, R, C2> operator*(matrix_t<T, R2, C2> &lhs, matrix_t<T, R2, C2> const &rhs)
            \brief Multiples rhs to matrix, throws if dimensions aren't compatible. Returns a matrix of possibly different dimensions.
         */
        template<typename T, size_t R, size_t C, size_t R2, size_t C2>
        constexpr matrix_t<T, R, C2> operator*(matrix_t<T, R, C> &lhs, matrix_t<T, R2, C2> const &rhs)
        {
            if( C != R2 )
            {
                std::stringstream error_stream;

                error_stream << "Cannot multiply matrices because the left-hand matrix has "
                             << C << " columns, which does not equal the number of rows of the right-hand matrix, "
                             << R << ".\n";

                std::string error_text = error_stream.str();
                throw std::invalid_argument(error_text);
            }

            matrix_t<T, R, C2> new_mat;

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
        /** \fn constexpr matrix_t<T, R, C> operator*(T const &lhs, matrix_t<T, R, C> const &rhs)
            \brief Scales matrix rhs by lhs.
         */
        template<typename T, size_t R, size_t C>
        constexpr matrix_t<T, R, C> operator*(T const &lhs, matrix_t<T, R, C> const &rhs)
        {
            auto new_mat = rhs;

            for( auto i = 0u; i < R; ++i )
            {
                for( auto j = 0u; j < C; ++j )
                {
                    new_mat[i][j] *= lhs;
                }
            }

            return new_mat;
        }

        /** \fn constexpr matrix_t<T, R, C> operator*(matrix_t<T, R2, C2> const &lhs, T const &rhs)
            \brief Scales matrix lhs by rhs.
         */
        template<typename T, size_t R, size_t C>
        constexpr matrix_t<T, R, C> operator*(matrix_t<T, R, C> const &lhs, T const &rhs)
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

