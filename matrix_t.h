#pragma once
#ifndef MATH_NERD_MATRIX_T_H
#define MATH_NERD_MATRIX_T_H

#include <array>
#include <cstdint>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>

#if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
    #include <omp.h>
    #define maybe_const
    // maybe_const means: const when OMP not enabled; not const when OMP enabled.
#else
    #define maybe_const const
#endif

#ifdef MATH_NERD_OMP_THREADCOUNT
    // Check if the value is empty. If so, set to 0.
    #if (0 - MATH_NERD_OMP_THREADCOUNT - 1) == 1 && (MATH_NERD_OMP_THREADCOUNT + 0) != -2
    // If macro is empty, this evaluates to (0 - - 1) == 1 && ( + 0) != -2.
        #define MATH_NERD_OMP_THREADCOUNT 0
    #endif
#endif


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
            constexpr matrix_t(s64 rows = 2, s64 columns = 0)
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

                for( auto i{ 0 }; i < rows_; ++i )
                {
                    mat.push_back(zero);
                }
            }

            constexpr matrix_t(matrix_t const &rhs) = default;

            constexpr matrix_t(matrix_t &&rhs) = default;

            constexpr ~matrix_t() = default;

            /** \name Array operators */
            /** \fn auto operator[](std::size_t pos) -> std::vector<T> &
                \brief Returns a reference to the row array at the specified position.
             */
            constexpr auto operator[](std::size_t pos) -> std::vector<T> &
            {
                return mat[pos];
            }

            /** \fn auto operator[](std::size_t pos) const -> std::vector<T>
                \brief Returns a copy to the row array at the specified position.
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
            /** \fn constexpr auto operator=(matrix_t<T> maybe_const &rhs) -> matrix_t<T> &
                \brief Sets matrix to rhs.
             */
            constexpr auto operator=(matrix_t<T> maybe_const  &rhs) -> matrix_t<T> &;

            /** \fn constexpr auto operator=(matrix_t<T> &&rhs) -> matrix_t<T> &
                \brief Move assignment.
             */
            constexpr auto operator=(matrix_t<T> &&rhs) -> matrix_t<T> &;

            /** \fn constexpr auto operator+=(matrix_t<T> maybe_const &rhs) -> matrix_t<T> &
                \brief Adds rhs to the matrix.
             */
            constexpr auto operator+=(matrix_t<T> maybe_const &rhs) -> matrix_t<T> &;

            /** \fn constexpr auto operator-=(matrix_t<T> maybe_const &rhs) -> matrix_t<T> &
                \brief Subtracts rhs from the matrix.
             */
            constexpr auto operator-=(matrix_t<T> maybe_const &rhs) -> matrix_t<T> &;

            /** \name Comparison operators */
            /** \fn constexpr auto operator==(matrix_t<T> maybe_const rhs) const noexcept -> bool
                \brief Compares the values and returns true if they are all equal.
             */
            constexpr auto operator==(matrix_t<T> maybe_const  &rhs) const noexcept -> bool;

            /** \fn constexpr auto operator!=(matrix_t<T> maybe_const rhs) const noexcept -> bool
                \brief Compares the values and returns false if they are all equal.
             */
            constexpr auto operator!=(matrix_t<T> maybe_const  &rhs) const noexcept -> bool;
        };

        /** \namespace math_nerd::matrix_t::vectorized
            \brief Contains multithreaded and vectorized versions of each operation for when not in a constexpr context.
        */
        namespace vectorized
        {
            /** \fn auto negation(matrix_t<T> maybe_const &old) -> matrix_t<T>
                \brief Returns the matrix whose entries are the negation of old.
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T>
            auto negation(matrix_t<T> maybe_const &old) -> matrix_t<T>;

            /** \fn auto assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void
                \brief Assigns rhs to lhs.
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T>
            auto assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void;

            /** \fn auto add_assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void
                \brief Assigns lhs to the sum of lhs and rhs.
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T>
            auto add_assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void;

            /** \fn auto subtract_assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void
                \brief Assigns lhs to the difference of lhs and rhs.
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T>
            auto subtract_assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void;

            /** \fn auto is_equal(matrix_t<T> maybe_const &lhs, matrix_t<T> maybe_const &rhs) -> bool
                \brief Returns whether the matrices are equal or not.
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T>
            auto is_equal(matrix_t<T> maybe_const &lhs, matrix_t<T> maybe_const &rhs) -> bool;

            /** \fn auto multiply_matrix(matrix_t<T> &C, matrix_t<T> maybe_const &A, matrix_t<T> maybe_const &B) -> void
                \brief Multiplies matrices A and B and assigns the result to .
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T>
            auto multiply_matrix(matrix_t<T> &C, matrix_t<T> maybe_const &A, matrix_t<T> maybe_const &B) -> void;

            /** \fn auto multiply_array(matrix_t<T> &C, matrix_t<T> maybe_const &A, std::array<T, len> maybe_const &B) -> void
                \brief Multiplies matrix A by array B and assigns the result to C.
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T, std::size_t len>
            auto multiply_array(matrix_t<T> &C, matrix_t<T> maybe_const &A, std::array<T, len> maybe_const &B) -> void;

            /** \fn auto multiply_vector(matrix_t<T> &C, matrix_t<T> maybe_const &A, std::vector<T> maybe_const &B) -> void
                \brief Multiplies matrix A by vector B and assigns the result to C.
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T>
            auto multiply_vector(matrix_t<T> &C, matrix_t<T> maybe_const &A, std::vector<T> maybe_const &B) -> void;

            /** \fn auto scalar_multiply(T maybe_const &lhs, matrix_t<T> &rhs) -> void
                \brief Returns the matrix whose entries are those of rhs scaled by lhs.
                       When OpenMP is enabled, it uses vector instructions.
            */
            template<typename T>
            auto scalar_multiply(T maybe_const &lhs, matrix_t<T> &rhs) -> void;
        } // namespace vectorized

        // Unary operators
        template<typename T>
        constexpr auto matrix_t<T>::operator+() const noexcept -> matrix_t<T>
        {
            return *this;
        }

        template<typename T>
        constexpr auto matrix_t<T>::operator-() const noexcept -> matrix_t<T>
        {
            if( std::is_constant_evaluated() )
            {
                matrix_t<T> new_mat{ rows_, columns_ };

                for( auto i = 0; i < rows_; ++i )
                {
                    for( auto j = 0; j < columns_; ++j )
                    {
                        new_mat[i][j] = -new_mat[i][j];
                    }
                }

                return new_mat;
            }
            else
            {
                return vectorized::negation<T>(*this);
            }
        }

        // Assignment operators
        template<typename T>
        constexpr auto matrix_t<T>::operator=(matrix_t<T> maybe_const &rhs) -> matrix_t<T> &
        {
            if( rows_ != rhs.row_count() || columns_ != rhs.column_count() )
            {
                throw std::invalid_argument(std::string("Cannot assign matrices which are not the same dimension.\n")
                    + "Left-hand matrix dimensions: " + std::to_string(rows_) + " x " + std::to_string(columns_) + "\n"
                    + "Right-hand matrix dimensions: " + std::to_string(rhs.row_count()) + " x "
                    + std::to_string(rhs.column_count()) + "\n");
            }

            if( std::is_constant_evaluated() )
            {
                for( auto i = 0; i < rows_; ++i )
                {
                    for( auto j = 0; j < columns_; ++j )
                    {
                        mat[i][j] = rhs[i][j];
                    }
                }
            }
            else
            {
                vectorized::assign(*this, rhs);
            }

            return *this;
        }

        template<typename T>
        constexpr auto matrix_t<T>::operator=(matrix_t<T> &&rhs) -> matrix_t<T> &
        {
            if( rows_ != rhs.row_count() || columns_ != rhs.column_count() )
            {
                throw std::invalid_argument(std::string("Cannot assign matrices which are not the same dimension.\n")
                    + "Left-hand matrix dimensions: " + std::to_string(rows_) + " x " + std::to_string(columns_) + "\n"
                    + "Right-hand matrix dimensions: " + std::to_string(rhs.row_count()) + " x "
                    + std::to_string(rhs.column_count()) + "\n");
            }

            mat = std::move(rhs);

            return *this;
        }

        template<typename T>
        constexpr auto matrix_t<T>::operator+=(matrix_t<T> maybe_const &rhs) -> matrix_t<T> &
        {
            if( rows_ != rhs.row_count() || columns_ != rhs.column_count() )
            {
                throw std::invalid_argument(std::string("Cannot add matrices which are not the same dimension.\n")
                    + "Left-hand matrix dimensions: " + std::to_string(rows_) + " x " + std::to_string(columns_) + "\n"
                    + "Right-hand matrix dimensions: " + std::to_string(rhs.row_count()) + " x "
                    + std::to_string(rhs.column_count()) + "\n");
            }

            if( std::is_constant_evaluated() )
            {
                for( auto i = 0; i < rows_; ++i )
                {
                    for( auto j = 0; j < columns_; ++j )
                    {
                        mat[i][j] += rhs[i][j];
                    }
                }
            }
            else
            {
                vectorized::add_assign(*this, rhs);
            }

            return *this;
        }

        template<typename T>
        constexpr auto matrix_t<T>::operator-=(matrix_t<T> maybe_const &rhs) -> matrix_t<T> &
        {
            if( rows_ != rhs.row_count() || columns_ != rhs.column_count() )
            {
                throw std::invalid_argument(std::string("Cannot subtract matrices which are not the same dimension.\n")
                    + "Left-hand matrix dimensions: " + std::to_string(rows_) + " x " + std::to_string(columns_) + "\n"
                    + "Right-hand matrix dimensions: " + std::to_string(rhs.row_count()) + " x "
                    + std::to_string(rhs.column_count()) + "\n");
            }

            if( std::is_constant_evaluated() )
            {
                for( auto i = 0; i < rows_; ++i )
                {
                    for( auto j = 0; j < columns_; ++j )
                    {
                        mat[i][j] -= rhs[i][j];
                    }
                }
            }
            else
            {
                vectorized::subtract_assign(*this, rhs);
            }

            return *this;
        }

        // Comparison operators
        template<typename T>
        constexpr auto matrix_t<T>::operator==(matrix_t<T> maybe_const &rhs) const noexcept -> bool
        {
            if( rows_ != rhs.row_count() || columns_ != rhs.column_count() )
            {
                return false;
            }
            else
            {
                if( std::is_constant_evaluated() )
                {
                    for( auto i = 0; i < rows_; ++i )
                    {
                        for( auto j = 0; j < columns_; ++j )
                        {
                            if( mat[i][j] != rhs[i][j] )
                            {
                                return false;
                            }
                        }
                    }
                }
                else
                {
                    return vectorized::is_equal(*this, rhs);
                }
            }

            return true;
        }

        template<typename T>
        constexpr auto matrix_t<T>::operator!=(matrix_t<T> maybe_const &rhs) const noexcept -> bool
        {
            return !(*this == rhs);
        }

        /** \name Binary operators */
        /** \fn auto operator+(matrix_t<T> lhs, matrix_t<T> maybe_const  &rhs) -> matrix_t<T> &
            \brief Adds two matrices of the same dimension together.
         */
        template<typename T>
        constexpr auto operator+(matrix_t<T> lhs, matrix_t<T> maybe_const &rhs) -> matrix_t<T>
        {
            try
            {
                lhs += rhs;
            }
            catch( std::invalid_argument & )
            {
                throw;
            }

            return lhs;
        }

        /** \fn auto operator-(matrix_t<T> lhs, matrix_t<T> maybe_const  &rhs) -> matrix_t<T> &
            \brief Subtracts one matrix from another of the same dimension.
         */
        template<typename T>
        constexpr auto operator-(matrix_t<T> lhs, matrix_t<T> maybe_const &rhs) -> matrix_t<T>
        {
            try
            {
                lhs -= rhs;
            }
            catch( std::invalid_argument & )
            {
                throw;
            }

            return lhs;
        }

        /** \fn constexpr auto operator*(matrix_t<T> maybe_const  &lhs, matrix_t<T> maybe_const  &rhs) -> matrix_t<T>
            \brief Multiples rhs to matrix, throws if dimensions aren't compatible. Returns a matrix of possibly different dimensions.
                   Non-const when OMP enabled for SIMD.
         */
        template<typename T>
        constexpr auto operator*(matrix_t<T> maybe_const &lhs, matrix_t<T> maybe_const &rhs) -> matrix_t<T>
        {
            s64 left_rows = lhs.row_count(), left_columns = lhs.column_count();
            s64 right_rows = rhs.row_count(), right_columns = rhs.column_count();

            if( right_rows != left_columns )
            {
                throw std::invalid_argument("Cannot multiply matrices because the left-hand matrix has "
                    + std::to_string(left_columns)
                    + " columns, which does not equal the number of rows of the right-hand matrix, "
                    + std::to_string(right_rows) + ".\n");
            }

            matrix_t<T> new_mat{ left_rows, right_columns };

            if( std::is_constant_evaluated() )
            {
                for( auto i = 0; i < left_rows; ++i )
                {
                    for( auto k = 0; k < left_columns; ++k )
                    {
                        for( auto j = 0; j < right_columns; ++j )
                        {
                            new_mat[i][j] += lhs[i][k] * rhs[k][j];
                        }
                    }
                }
            }
            else
            {
                vectorized::multiply_matrix(new_mat, lhs, rhs);
            }

            return new_mat;
        }

        /** \fn constexpr auto operator*(matrix_t<T> maybe_const  &lhs, std::array<T, len> maybe_const  &rhs) -> matrix_t<T>
            \brief Matrix-Vector multiplication with std::array.
         */
        template<typename T, std::size_t len>
        constexpr auto operator*(matrix_t<T> maybe_const &lhs, std::array<T, len> maybe_const &rhs) -> matrix_t<T>
        {
            s64 left_rows = lhs.row_count(), left_columns = lhs.column_count();
            s64 right_rows = len;

            if( right_rows != left_columns )
            {
                throw std::invalid_argument("Cannot multiply matrices because the left-hand matrix has "
                    + std::to_string(left_columns)
                    + " columns, which does not equal the number of rows of the right-hand matrix, "
                    + std::to_string(right_rows) + ".\n");
            }

            matrix_t<T> new_mat{ left_rows, 1 };

            if( std::is_constant_evaluated() )
            {
                for( auto i = 0; i < left_rows; ++i )
                {
                    for( auto k = 0; k < left_columns; ++k )
                    {
                        new_mat[i][0] += lhs[i][k] * rhs[k];
                    }
                }
            }
            else
            {
                vectorized::multiply_array(new_mat, lhs, rhs);
            }

            return new_mat;
        }

        /** \fn constexpr auto operator*(matrix_t<T> maybe_const &lhs, std::vector<T> maybe_const  &rhs) -> matrix_t<T>
            \brief Matrix-Vector multiplication with std::vector.
         */
        template<typename T>
        constexpr auto operator*(matrix_t<T> maybe_const &lhs, std::vector<T> maybe_const &rhs) -> matrix_t<T>
        {
            s64 left_rows = lhs.row_count(), left_columns = lhs.column_count();
            s64 right_rows = rhs.size();

            if( right_rows != left_columns )
            {
                throw std::invalid_argument("Cannot multiply matrices because the left-hand matrix has "
                    + std::to_string(left_columns)
                    + " columns, which does not equal the number of rows of the right-hand matrix, "
                    + std::to_string(right_rows) + ".\n");
            }

            matrix_t<T> new_mat{ left_rows, 1 };

            if( std::is_constant_evaluated() )
            {
                for( auto i = 0; i < left_rows; ++i )
                {
                    for( auto k = 0; k < left_columns; ++k )
                    {
                        new_mat[i][0] += lhs[i][k] * rhs[k];
                    }
                }
            }
            else
            {
                vectorized::multiply_vector(new_mat, lhs, rhs);
            }

            return new_mat;
        }

        /** \name Scalar multiplication operators */
        /** \fn constexpr auto operator*(T maybe_const &lhs, matrix_t<T> rhs) -> matrix_t<T>
            \brief Scales matrix rhs by lhs.
         */
        template<typename T>
        constexpr auto operator*(T maybe_const &lhs, matrix_t<T> rhs) -> matrix_t<T>
        {
            if( std::is_constant_evaluated() )
            {
                for( auto i = 0; i < rhs.row_count(); ++i )
                {
                    for( auto j = 0; j < rhs.column_count(); ++j )
                    {
                        rhs[i][j] *= lhs;
                    }
                }
            }
            else
            {
                vectorized::scalar_multiply(lhs, rhs);
            }

            return rhs;
        }

        /** \fn constexpr auto operator*(matrix_t<T> &lhs, T maybe_const &rhs) -> matrix_t<T>
            \brief Scales matrix lhs by rhs.
         */
        template<typename T>
        constexpr auto operator*(matrix_t<T> &lhs, T maybe_const &rhs) -> matrix_t<T>
        {
            return rhs * lhs;
        }


        // Vectorized implementations.
        namespace vectorized
        {
            template<typename T>
            auto negation(matrix_t<T> maybe_const &old) -> matrix_t<T>
            {
                auto rows_{ old.row_count() };
                auto columns_{ old.column_count() };

                matrix_t<T> new_mat{ rows_, columns_ };

                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel for
                    #endif
                #endif
                for( auto i = 0; i < rows_; ++i )
                {
                    #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                        #pragma omp simd
                    #endif
                    for( auto j = 0; j < columns_; ++j )
                    {
                        new_mat[i][j] = -old[i][j];
                    }
                }

                return new_mat;
            }

            template<typename T>
            auto assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void
            {
                auto rows_{ lhs.row_count() }, columns_{ lhs.column_count() };

                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel for
                    #endif
                #endif
                for( auto i = 0; i < rows_; ++i )
                {
                    #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                    #pragma omp simd
                    #endif
                    for( auto j = 0; j < columns_; ++j )
                    {
                        lhs[i][j] = rhs[i][j];
                    }
                }
            }

            template<typename T>
            auto add_assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void
            {
                auto rows_{ lhs.row_count() }, columns_{ lhs.column_count() };

                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel for
                    #endif
                #endif
                for( auto i = 0; i < rows_; ++i )
                {
                    #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                        #pragma omp simd
                    #endif
                    for( auto j = 0; j < columns_; ++j )
                    {
                        lhs[i][j] += rhs[i][j];
                    }
                }
            }

            template<typename T>
            auto subtract_assign(matrix_t<T> &lhs, matrix_t<T> maybe_const &rhs) -> void
            {
                auto rows_{ lhs.row_count() }, columns_{ lhs.column_count() };

                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel for
                    #endif
                #endif
                for( auto i = 0; i < rows_; ++i )
                {
                    #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                        #pragma omp simd
                    #endif
                    for( auto j = 0; j < columns_; ++j )
                    {
                        lhs[i][j] -= rhs[i][j];
                    }
                }
            }

            template<typename T>
            auto is_equal(matrix_t<T> maybe_const &lhs, matrix_t<T> maybe_const &rhs) -> bool
            {
                auto rows_{ lhs.row_count() }, columns_{ lhs.column_count() };
                bool res = true;

                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel for
                    #endif
                #endif
                for( auto i = 0; i < rows_ && res; ++i )
                {
                    #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                        #pragma omp simd
                    #endif
                    for( auto j = 0; j < columns_; ++j )
                    {
                        if( lhs[i][j] != rhs[i][j] )
                        {
                            res = false;
                        }
                    }
                }

                return res;
            }

            template<typename T>
            auto multiply_matrix(matrix_t<T> &C, matrix_t<T> maybe_const &A, matrix_t<T> maybe_const &B) -> void
            {
                auto left_rows{ A.row_count() }, left_columns{ A.column_count() };
                auto right_columns{ B.column_count() };

                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel
                    #endif
                #endif
                for( auto i = 0; i < left_rows; ++i )
                {
                    for( auto k = 0; k < left_columns; ++k )
                    {
                        #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                            #pragma omp simd
                        #endif
                        for( auto j = 0; j < right_columns; ++j )
                        {
                            C[i][j] += A[i][k] * B[k][j];
                        }
                    }
                }
            }

            template<typename T, std::size_t len>
            auto multiply_array(matrix_t<T> &C, matrix_t<T> maybe_const &A, std::array<T, len> maybe_const &B) -> void
            {
                auto left_rows{ A.row_count() }, left_columns{ A.column_count() };

                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel for
                    #endif
                #endif
                for( auto i = 0; i < left_rows; ++i )
                {
                    #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                        #pragma omp simd
                    #endif
                    for( auto k = 0; k < left_columns; ++k )
                    {
                        C[i][0] += A[i][k] * B[k];
                    }
                }
            }


            template<typename T>
            auto multiply_vector(matrix_t<T> &C, matrix_t<T> maybe_const &A, std::vector<T> maybe_const &B) -> void
            {
                auto left_rows{ A.row_count() }, left_columns{ A.column_count() };

                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel for
                    #endif
                #endif
                for( auto i = 0; i < left_rows; ++i )
                {
                    #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                        #pragma omp simd
                    #endif
                    for( auto k = 0; k < left_columns; ++k )
                    {
                        C[i][0] += A[i][k] * B[k];
                    }
                }
            }

            template<typename T>
            auto scalar_multiply(T maybe_const &lhs, matrix_t<T> &rhs) -> void
            {
                #ifdef MATH_NERD_OMP_THREADCOUNT
                    #if MATH_NERD_OMP_THREADCOUNT > 0
                        #pragma omp parallel for num_threads(MATH_NERD_OMP_THREADCOUNT)
                    #else
                        #pragma omp parallel for
                    #endif
                #endif
                for( auto i = 0; i < rhs.row_count(); ++i )
                {
                    #if defined(MATH_NERD_OMP_THREADCOUNT) || defined(MATH_NERD_OMP)
                        #pragma omp simd
                    #endif
                    for( auto j = 0; j < rhs.column_count(); ++j )
                    {
                        rhs[i][j] *= lhs;
                    }
                }
            }
        } // namespace vectorized

    } // namespace matrix_t
} // namespace math_nerd

/** \mainpage Minimal Matrix Implementation
    \section minmat Minimal Matrix
    This is a very minimal implementation of a matrix class. It takes arithmetic-like types (int, float, or custom types with +, -, *, /)
    and allows for adding/subtracting matrices, as well as multiplying. The inverse function is not implemented here.

    The main use of this was simply to implement the <a href="../hill_cipher/index.html">Hill Cipher (modulo 97)</a>, which uses another type I made,
    the <a href="../int_mod/index.html">Integers Modulo N</a>.
    \section gitlab_link GitLab Link
    View the source code at <a href="https://gitlab.com/mathnerd/minimal-matrix">GitLab</a>.
 */
#endif
