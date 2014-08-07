/*
 * Copyright 2013-2014 Milian Wolff <mail@milianw.de>
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Library General Public License as
 * published by the Free Software Foundation; either version 2 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public
 * License along with this program; if not, write to the
 * Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#ifndef VECTOR3D_H
#define VECTOR3D_H

#include <iostream>
#include <cmath>
#include <cassert>

#if defined(__INTEL_COMPILER)
#define FUNC_CONSTEXPR
#else
#define FUNC_CONSTEXPR constexpr
#endif

template<typename S, typename E1, typename E2>
class Vector3D_CrossExpr;

/**
 * Expression template pattern implementation using the CRTP pattern.
 *
 * see also: http://en.wikipedia.org/wiki/Expression_templates
 */
template <typename S, typename E>
class Vector3D_Expr
{
public:
    /**
     * \return the x-component of this vector expression.
     */
    inline FUNC_CONSTEXPR S x() const noexcept
    {
        return static_cast<const E&>(*this).x();
    }

    /**
     * \return the y-component of this vector expression.
     */
    inline FUNC_CONSTEXPR S y() const noexcept
    {
        return static_cast<const E&>(*this).y();
    }

    /**
     * \return the z-component of this vector expression.
     */
    inline FUNC_CONSTEXPR S z() const noexcept
    {
        return static_cast<const E&>(*this).z();
    }

    /**
     * \return the dot product of this vector expression with the other \p vector expression.
     */
    template<typename E2>
    inline FUNC_CONSTEXPR S dot(const Vector3D_Expr<S, E2>& vector) const noexcept
    {
        return x() * vector.x()
             + y() * vector.y()
             + z() * vector.z();
    }

    /**
     * \return An expression template equivalent to the cross product of this vector expression with the \p other expression.
     */
    template<typename E2>
    inline FUNC_CONSTEXPR Vector3D_CrossExpr<S, E, E2> cross(const Vector3D_Expr<S, E2>& b) const noexcept;

    /**
     * \return the squared norm of this vector expression, i.e. the value of the dot-product with itself.
     */
    inline FUNC_CONSTEXPR S squaredNorm() const noexcept
    {
        return dot(*this);
    }

    /**
     * \return the norm or length of this vector, i.e. the square root of the result of the dot-product with itself.
     */
    inline FUNC_CONSTEXPR S norm() const noexcept
    {
        return std::sqrt(squaredNorm());
    }
};

/**
 * A 3D vector to be used in the calculations.
 */
template<typename S>
class Vector3D : public Vector3D_Expr<S, Vector3D<S> >
{
public:
    /**
     * Create a vector with all elements initialized to zero.
     */
    inline FUNC_CONSTEXPR Vector3D() noexcept
    : m_x(0), m_y(0), m_z(0)
    {
    }

    /**
     * Create a vector initialized with the element values \p x, \p y and \p z.
     */
    inline FUNC_CONSTEXPR Vector3D(S x, S y, S z) noexcept
    : m_x(x), m_y(y), m_z(z)
    {
    }

    /**
     * Create a vector with all elements initialized with \p value.
     */
    inline FUNC_CONSTEXPR Vector3D(S value) noexcept
    : m_x(value), m_y(value), m_z(value)
    {
    }

    /**
     * Create a vector from another \p vector or a temporary vector expression.
     */
    template<typename S2, typename E>
    inline FUNC_CONSTEXPR Vector3D(const Vector3D_Expr<S2, E>& vector) noexcept
    : m_x(vector.x()), m_y(vector.y()), m_z(vector.z())
    {
    }

    /**
     * Set all elements in this vector to \p scalar.
     *
     * \return reference to this vector
     */
    inline Vector3D<S>& setConstant(S scalar) noexcept
    {
        m_x = m_y = m_z = scalar;
        return *this;
    }

    /**
     * \return true if the \p other vector equals this vector, false otherwise.
     */
    inline FUNC_CONSTEXPR bool operator==(const Vector3D<S>& vector) const noexcept
    {
        return m_x == vector.x && m_y == vector.y && m_z == vector.z;
    }

    /**
     * Add \p vector element-wise to this vector.
     *
     * \return reference to this vector.
     */
    template<typename E>
    inline Vector3D<S>& operator+=(const Vector3D_Expr<S, E>& vector) noexcept
    {
        m_x += vector.x();
        m_y += vector.y();
        m_z += vector.z();
        return *this;
    }

    /**
     * Subtract \p vector element-wise from this vector.
     *
     * \return reference to this vector.
     */
    template<typename E>
    inline Vector3D<S>& operator-=(const Vector3D_Expr<S, E>& vector) noexcept
    {
        m_x -= vector.x();
        m_y -= vector.y();
        m_z -= vector.z();
        return *this;
    }

    /**
     * Devide every element in this vector by \p scalar.
     *
     * \return reference to this vector.
     */
    inline Vector3D& operator/=(S scalar) noexcept
    {
        m_x /= scalar;
        m_y /= scalar;
        m_z /= scalar;
        return *this;
    }

    /**
     * Multiply every element in this vector with \p scalar.
     *
     * \return reference to this vector.
     */
    inline Vector3D& operator*=(S scalar) noexcept
    {
        m_x *= scalar;
        m_y *= scalar;
        m_z *= scalar;
        return *this;
    }

    inline FUNC_CONSTEXPR S x() const noexcept
    {
        return m_x;
    }
    inline FUNC_CONSTEXPR S y() const noexcept
    {
        return m_y;
    }
    inline FUNC_CONSTEXPR S z() const noexcept
    {
        return m_z;
    }

    inline S& x() noexcept
    {
        return m_x;
    }
    inline S& y() noexcept
    {
        return m_y;
    }
    inline S& z() noexcept
    {
        return m_z;
    }

    /**
     * \return a copy of this vector with all elements divided by the norm of this vector.
     */
    inline Vector3D<S> normalized() const noexcept
    {
        const S n = this->norm();
        return Vector3D<S>(m_x / n, m_y / n, m_z / n);
    }

private:
    S m_x;
    S m_y;
    S m_z;
};

/**
 * Output space-separated elements of \p vector to \p stream.
 */
template<typename S, typename E>
inline std::ostream& operator<<(std::ostream& stream, const Vector3D_Expr<S, E>& vector)
{
    stream << vector.x() << ' ' << vector.y() << ' ' << vector.z();
    return stream;
}

/**
 * Read three scalars from \p stream into \p vector.
 */
template<typename S>
inline std::istream& operator>>(std::istream& stream, Vector3D<S>& vector)
{
    stream >> vector.x() >> vector.y() >> vector.z();
    return stream;
}

/**
 * Expression template for the sum of two vectors.
 */
template<typename S, typename E1, typename E2>
class Vector3D_SumExpr : public Vector3D_Expr<S, Vector3D_SumExpr<S, E1, E2> >
{
public:
    FUNC_CONSTEXPR Vector3D_SumExpr(const Vector3D_Expr<S, E1>& l, const Vector3D_Expr<S, E2>& r) noexcept
    : m_l(l)
    , m_r(r)
    {}
    inline FUNC_CONSTEXPR S x() const noexcept
    {
        return m_l.x() + m_r.x();
    }
    inline FUNC_CONSTEXPR S y() const noexcept
    {
        return m_l.y() + m_r.y();
    }
    inline FUNC_CONSTEXPR S z() const noexcept
    {
        return m_l.z() + m_r.z();
    }
private:
    const Vector3D_Expr<S, E1>& m_l;
    const Vector3D_Expr<S, E2>& m_r;
};

/**
 * \return An expression template equivalent to the elementwise addition of \p augend and \p addend.
 */
template<typename S, typename E1, typename E2>
FUNC_CONSTEXPR Vector3D_SumExpr<S, E1, E2> operator+(const Vector3D_Expr<S, E1>& augend, const Vector3D_Expr<S, E2>& addend) noexcept
{
    return Vector3D_SumExpr<S, E1, E2>(augend, addend);
}

/**
 * Expression template for the element-wise difference of two vectors.
 */
template<typename S, typename E1, typename E2>
class Vector3D_DiffExpr : public Vector3D_Expr<S, Vector3D_DiffExpr<S, E1, E2> >
{
public:
    FUNC_CONSTEXPR Vector3D_DiffExpr(const Vector3D_Expr<S, E1>& l, const Vector3D_Expr<S, E2>& r) noexcept
    : m_l(l)
    , m_r(r)
    {}
    inline FUNC_CONSTEXPR S x() const noexcept
    {
        return m_l.x() - m_r.x();
    }
    inline FUNC_CONSTEXPR S y() const noexcept
    {
        return m_l.y() - m_r.y();
    }
    inline FUNC_CONSTEXPR S z() const noexcept
    {
        return m_l.z() - m_r.z();
    }
private:
    const Vector3D_Expr<S, E1>& m_l;
    const Vector3D_Expr<S, E2>& m_r;
};

/**
 * \return An expression template equivalent to the elementwise subtraction of  \p minuend and \p subtrahend.
 */
template<typename S, typename E1, typename E2>
FUNC_CONSTEXPR const Vector3D_DiffExpr<S, E1, E2>
inline operator-(const Vector3D_Expr<S, E1>& minuend, const Vector3D_Expr<S, E2>& subtrahend) noexcept
{
    return Vector3D_DiffExpr<S, E1, E2>(minuend, subtrahend);
}

/**
 * Expression template for the negation of a vector.
 */
template<typename S, typename E>
class Vector3D_NegateExpr : public Vector3D_Expr<S, Vector3D_NegateExpr<S, E> >
{
public:
    FUNC_CONSTEXPR Vector3D_NegateExpr(const Vector3D_Expr<S, E>& vec) noexcept
    : m_vec(vec)
    {}
    FUNC_CONSTEXPR inline S x() const noexcept
    {
        return -m_vec.x();
    }
    FUNC_CONSTEXPR inline S y() const noexcept
    {
        return -m_vec.y();
    }
    FUNC_CONSTEXPR inline S z() const noexcept
    {
        return -m_vec.z();
    }
private:
    const E& m_vec;
};

/**
 * \return An expression template equivalent to a copy of \p vector with every element negated.
 */
template<typename S, typename E>
FUNC_CONSTEXPR inline const Vector3D_NegateExpr<S, E> operator-(const Vector3D_Expr<S, E>& vector) noexcept
{
    return Vector3D_NegateExpr<S, E>(v);
}

/**
 * Expression template for the element-wise multiplication of a vector with a scalar.
 */
template<typename S, typename E, typename S2>
class Vector3D_MultExpr : public Vector3D_Expr<S, Vector3D_MultExpr<S, E, S2> >
{
public:
    FUNC_CONSTEXPR Vector3D_MultExpr(const Vector3D_Expr<S, E>& vec, S2 s) noexcept
    : m_vec(vec)
    , m_s(s)
    {}
    FUNC_CONSTEXPR inline S x() const noexcept
    {
        return m_vec.x() * m_s;
    }
    FUNC_CONSTEXPR inline S y() const noexcept
    {
        return m_vec.y() * m_s;
    }
    FUNC_CONSTEXPR inline S z() const noexcept
    {
        return m_vec.z() * m_s;
    }
private:
    const Vector3D_Expr<S, E>& m_vec;
    S2 m_s;
};

/**
 * \return An expression template equivalent to a copy of \p vector with every element multiplied by \p scalar.
 */
template<typename S, typename E, typename S2>
FUNC_CONSTEXPR inline const Vector3D_MultExpr<S, E, S2> operator*(const Vector3D_Expr<S, E>&vector, S2 scalar) noexcept
{
    return Vector3D_MultExpr<S, E, S2>(vector, scalar);
}

/**
 * \return An expression template equivalent to a copy of \p vector with every element multiplied by \p scalar.
 */
template<typename S, typename E, typename S2>
FUNC_CONSTEXPR inline const Vector3D_MultExpr<S, E, S2> operator*(S2 scalar, const Vector3D_Expr<S, E>&vector) noexcept
{
    return Vector3D_MultExpr<S, E, S2>(vector, scalar);
}

/**
 * Expression template for the element-wise division of a vector with a scalar.
 */
template<typename S, typename E, typename S2>
class Vector3D_DivExpr : public Vector3D_Expr<S, Vector3D_DivExpr<S, E, S2> >
{
public:
    FUNC_CONSTEXPR Vector3D_DivExpr(const Vector3D_Expr<S, E>& vec, S2 s) noexcept
    : m_vec(vec)
    , m_s(s)
    {}
    FUNC_CONSTEXPR inline S x() const noexcept
    {
        return m_vec.x() / m_s;
    }
    FUNC_CONSTEXPR inline S y() const noexcept
    {
        return m_vec.y() / m_s;
    }
    FUNC_CONSTEXPR inline S z() const noexcept
    {
        return m_vec.z() / m_s;
    }
private:
    const Vector3D_Expr<S, E>& m_vec;
    S2 m_s;
};

/**
 * \return An expression template equivalent to a copy of \p vector with every element devided by \p scalar.
 */
template<typename S, typename E, typename S2>
FUNC_CONSTEXPR inline const Vector3D_DivExpr<S, E, S2> operator/(const Vector3D_Expr<S, E>& vector, S2 scalar) noexcept
{
    return Vector3D_DivExpr<S, E, S2>(vector, scalar);
}

/**
 * Expression template for the cross product between two vectors.
 */
template<typename S, typename E1, typename E2>
class Vector3D_CrossExpr : public Vector3D_Expr<S, Vector3D_CrossExpr<S, E1, E2> >
{
public:
    FUNC_CONSTEXPR Vector3D_CrossExpr(const Vector3D_Expr<S, E1>& a, const Vector3D_Expr<S, E2>& b) noexcept
    : a(a)
    , b(b)
    {}
    FUNC_CONSTEXPR inline S x() const noexcept
    {
        return a.y() * b.z() - b.y() * a.z();
    }
    FUNC_CONSTEXPR inline S y() const noexcept
    {
        return a.z() * b.x() - b.z() * a.x();
    }
    FUNC_CONSTEXPR inline S z() const noexcept
    {
        return a.x() * b.y() - b.x() * a.y();
    }
private:
    const Vector3D_Expr<S, E1>& a;
    const Vector3D_Expr<S, E2>& b;
};

/**
 * \return An expression template equivalent to the cross product of this vector expression with the \p other expression.
 */
template<typename S, typename E>
template<class E2>
inline FUNC_CONSTEXPR Vector3D_CrossExpr<S, E, E2> Vector3D_Expr<S, E>::cross(const Vector3D_Expr<S, E2>& other) const noexcept
{
    return Vector3D_CrossExpr<S, E, E2>(*this, other);
}

#undef FUNC_CONSTEXPR

typedef Vector3D<double> dvec;

#endif // VECTOR3D_H
