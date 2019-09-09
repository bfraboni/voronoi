/****************************************************************************
Copyright (C) 2010-2020 Alexandre Meyer

This file is part of  library.

gkit2light is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

gkit2light is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with gkit2light.  If not, see <http://www.gnu.org/licenses/>.
*****************************************************************************/

#ifndef _QUATERNION_H
#define _QUATERNION_H


#ifndef _USE_MATH_DEFINES
#define _USE_MATH_DEFINES
#endif
#include <cmath>
#include <iostream>
#include <cstdlib>


	/*! \brief A Quaternion class
	 *
	 */
	template<typename Real, typename Vec3Real>
	class TQuaternion
	{
	public:
		/*! @name Defining a Quaternion */
		//@{
		/*! Default constructor, builds an identity rotation. */
		TQuaternion()
		{ q[0]=q[1]=q[2]=0.0;  q[3]=1.0; }

		/*! Constructor from rotation axis (non null) and angle (in radians). See also setAxisAngle(). */
		TQuaternion(const Vec3Real& axis, const Real angle)
		{
			setAxisAngle(axis, angle);
		}

		TQuaternion(const Vec3Real& from, const Vec3Real& to)
		{
			const Real epsilon = 1E-10f;

			const Real fromSqNorm = from.squaredNorm();
			const Real toSqNorm   = to.squaredNorm();
			// Identity TQuaternion when one vector is null
			if ((fromSqNorm < epsilon) || (toSqNorm < epsilon))
			{
				q[0]=q[1]=q[2]=0.0;
				q[3]=1.0;
			}
			else
			{
				Vec3Real axis = cross(from, to);
				const Real axisSqNorm = axis.squaredNorm();

				// Aligned vectors, pick any axis, not aligned with from or to
				if (axisSqNorm < epsilon)
				{
					axis = Vec3Real(1.0, 0.0, 0.0);
					if (axis*from < 0.1*sqrt(fromSqNorm))
						axis = Vec3Real(0.0, 1.0, 0.0);
				}

				Real angle = asin(sqrt(axisSqNorm / (fromSqNorm * toSqNorm)));

				if (from*to < 0.0)
					angle = M_PI-angle;

				setAxisAngle(axis, angle);
			}
		}


		/*! Constructor from the four values of a TQuaternion. First three values are axis*sin(angle/2) and
		  last one is cos(angle/2).
		\attention The identity TQuaternion is TQuaternion(0,0,0,1) and \e not TQuaternion(0,0,0,0) (which is
		not unitary). The default TQuaternion() creates such identity TQuaternion. */
		TQuaternion(Real q0, Real q1, Real q2, Real q3)
		{ q[0]=q0;    q[1]=q1;    q[2]=q2;    q[3]=q3; }

		/*! Copy constructor. */
		TQuaternion(const TQuaternion& Q)
		{ for (int i=0; i<4; ++i) q[i] = Q.q[i]; }

		/*! Equal operator. */
		TQuaternion& operator=(const TQuaternion& Q)
		{
			for (int i=0; i<4; ++i)
				q[i] = Q.q[i];
			return (*this);
		}

		TQuaternion& operator+=(const TQuaternion& Q)
		{
			for (int i = 0; i < 4; ++i)
				q[i] += Q.q[i];
			return (*this);
		}


		/*! Sets the TQuaternion as a rotation of axis \p axis and angle \p angle (in radians).

		\p axis does not need to be normalized. A null \p axis will result in an identity TQuaternion. */
		void setAxisAngle(const Vec3Real& axis, const Real angle)
		{
			const Real norm = length(axis); // axis.norm();
			if (norm < 1E-8)
			{
				// Null rotation
				q[0] = 0.0;
				q[1] = 0.0;
				q[2] = 0.0;
				q[3] = 1.0;
			}
			else
			{
				const Real sin_half_angle = sin(angle / 2.0);
				q[0] = sin_half_angle*axis.x/norm;
				q[1] = sin_half_angle*axis.y/norm;
				q[2] = sin_half_angle*axis.z/norm;
				q[3] = cos(angle / 2.0);
			}
		}

		/*! Sets the TQuaternion value. See the TQuaternion(Real, Real, Real, Real) constructor documentation. */
		void setValue(Real q0, Real q1, Real q2, Real q3)
		{ q[0]=q0;    q[1]=q1;    q[2]=q2;    q[3]=q3; }


		template<typename MAT>
		void setFromRotationMatrix(const MAT& m)
        {
			// First, find largest diagonal in matrix:
			int i = 2;
			if (m[0][0] > m[1][1])
			{
				if (m[0][0] > m[2][2])
				{
					i = 0;
				}
			}
			else
			{
				if (m[1][1] > m[2][2])
				{
					i = 1;
				}
			}

			if (m[0][0]+m[1][1]+m[2][2] > m[i][i])
			{
				// Compute w first:
				q[3] = sqrt(m[0][0]+m[1][1]+m[2][2]+1.0)/2.0;
				// And compute other values:
				q[0] = (m[2][1]-m[1][2])/(4.0*q[3]);
				q[1] = (m[0][2]-m[2][0])/(4.0*q[3]);
				q[2] = (m[1][0]-m[0][1])/(4.0*q[3]);
			}
			else
			{
				// Compute x, y, or z first:
				int j = (i+1)%3;
				int k = (i+2)%3;

				// Compute first value:
				q[i] = sqrt(m[i][i]-m[j][j]-m[k][k]+1.0)/2.0;

				// And the others:
				q[j] = (m[i][j]+m[j][i])/(4.0*q[i]);
				q[k] = (m[i][k]+m[k][i])/(4.0*q[i]);

				q[3] = (m[k][j]-m[j][k])/(4.0*q[i]);
			}
		}

		void setFromRotatedBase(const Vec3Real& X,  const Vec3Real& Y,    const Vec3Real& Z)
		{
			Real m[3][3];
			Real normX = X.norm();
			Real normY = Y.norm();
			Real normZ = Z.norm();

			for (int i=0; i<3; ++i)
			{
				m[i][0] = X[i] / normX;
				m[i][1] = Y[i] / normY;
				m[i][2] = Z[i] / normZ;
			}
			setFromRotationMatrix(m);
		}

		//@}


		/*! @name Accessing values */
		//@{
		Vec3Real axis() const
		{
			Vec3Real res = Vec3Real(q[0], q[1], q[2]);
			const Real sinus = res.norm();
			if (sinus > 1E-8)
				res /= sinus;
			return (acos(q[3]) <= M_PI/2.0) ? res : -res;
		}


		Real angle() const
		{
			const Real angle = 2.0 * acos(q[3]);
			return (angle <= M_PI) ? angle : 2.0*M_PI - angle;
		}


		void getAxisAngle(Vec3Real& axis, Real& angle) const
		{
			angle = 2.0*acos(q[3]);
			axis = Vec3Real(q[0], q[1], q[2]);
			const Real sinus = axis.norm();
			if (sinus > 1E-8)
				axis /= sinus;

			if (angle > M_PI)
			{
				angle = 2.0*M_PI - angle;
				axis = -axis;
			}
		}


		/*! Bracket operator, with a constant return value. \p i must range in [0..3]. See the TQuaternion(Real, Real, Real, Real) documentation. */
		Real operator[](int i) const { return q[i]; }

		/*! Bracket operator returning an l-value. \p i must range in [0..3]. See the TQuaternion(Real, Real, Real, Real) documentation. */
		Real& operator[](int i) { return q[i]; }
		//@}



        /*! multiplication  */
        friend TQuaternion operator*(const Real a, const TQuaternion& b)
        {
            return TQuaternion<Real,Vec3Real>(a*b[0],a*b[1],a*b[2],a*b[3]);
        }


		/*! @name Rotation computations */
		//@{
		/*! Returns the composition of the \p a and \p b rotations.
		The order is important. When applied to a Vec \c v (see operator*(const TQuaternion&, const Vec&)
		and rotate()) the resulting TQuaternion acts as if \p b was applied first and then \p a was
		applied. This is obvious since the image \c v' of \p v by the composited rotation satisfies: \code
		v'= (a*b) * v = a * (b*v) \endcode
		Note that a*b usually differs from b*a.
		\attention For efficiency reasons, the resulting TQuaternion is not normalized. Use normalize() in
		case of numerical drift with small rotation composition. */
		friend TQuaternion operator*(const TQuaternion& a, const TQuaternion& b)
		{
			return TQuaternion(a.q[3]*b.q[0] + b.q[3]*a.q[0] + a.q[1]*b.q[2] - a.q[2]*b.q[1],
								a.q[3]*b.q[1] + b.q[3]*a.q[1] + a.q[2]*b.q[0] - a.q[0]*b.q[2],
								a.q[3]*b.q[2] + b.q[3]*a.q[2] + a.q[0]*b.q[1] - a.q[1]*b.q[0],
								a.q[3]*b.q[3] - b.q[0]*a.q[0] - a.q[1]*b.q[1] - a.q[2]*b.q[2]);
		}

		/*! TQuaternion rotation is composed with \p q.

		See operator*(), since this is equivalent to \c this = \c this * \p q.

		\note For efficiency reasons, the resulting TQuaternion is not normalized.
		You may normalize() it after each application in case of numerical drift. */
		TQuaternion& operator*=(const TQuaternion &q)
		{
			*this = (*this)*q;
			return *this;
		}

		/*! Returns the image of \p v by the rotation \p q.

		Same as q.rotate(v). See rotate() and inverseRotate(). */
		friend Vec3Real operator*(const TQuaternion& q, const Vec3Real& v)
		{
			return q.rotate(v);
		}



		Vec3Real rotate(const Vec3Real& v) const
		{
			const Real q00 = 2.0l * q[0] * q[0];
			const Real q11 = 2.0l * q[1] * q[1];
			const Real q22 = 2.0l * q[2] * q[2];

			const Real q01 = 2.0l * q[0] * q[1];
			const Real q02 = 2.0l * q[0] * q[2];
			const Real q03 = 2.0l * q[0] * q[3];

			const Real q12 = 2.0l * q[1] * q[2];
			const Real q13 = 2.0l * q[1] * q[3];

			const Real q23 = 2.0l * q[2] * q[3];

			return Vec3Real((1.0 - q11 - q22)*v.x + (      q01 - q23)*v.y + (      q02 + q13)*v.z,
							(      q01 + q23)*v.x + (1.0 - q22 - q00)*v.y + (      q12 - q03)*v.z,
							(      q02 - q13)*v.x + (      q12 + q03)*v.y + (1.0 - q11 - q00)*v.z );
		}


		Vec3Real inverseRotate(const Vec3Real& v) const
		{
			return inverse().rotate(v);
		}



		//@}


		/*! @name Inversion */
		//@{
		/*! Returns the inverse TQuaternion (inverse rotation).

		Result has a negated axis() direction and the same angle(). A composition (see operator*()) of a
		TQuaternion and its inverse() results in an identity function.

		Use invert() to actually modify the TQuaternion. */
		TQuaternion inverse() const { return TQuaternion(-q[0], -q[1], -q[2], q[3]); }

		/*! Inverses the TQuaternion (same rotation angle(), but negated axis()).

		See also inverse(). */
		void invert() { q[0] = -q[0]; q[1] = -q[1]; q[2] = -q[2]; }

		/*! Negates all the coefficients of the TQuaternion.

		This results in an other representation of the \e same rotation (opposite rotation angle, but with
		a negated axis direction: the two cancel out). However, note that the results of axis() and
		angle() are unchanged after a call to this method since angle() always returns a value in [0,pi].

		This method is mainly useful for TQuaternion interpolation, so that the spherical
		interpolation takes the shortest path on the unit sphere. See slerp() for details. */
		void negate() { invert(); q[3] = -q[3]; }

		/*! Normalizes the TQuaternion coefficients.

		This method should not need to be called since we only deal with unit TQuaternions. This is however
		useful to prevent numerical drifts, especially with small rotational increments. */
		Real normalize()
		{
			const Real norm = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
			for (int i=0; i<4; ++i)
				q[i] /= norm;
			return norm;
		}
		//@}


		/*! @name Associated matrix */
		//@{
		/*! Returns the TQuaternion associated 4x4 OpenGL rotation matrix.

		 Use \c glMultMatrixd(q.matrix()) or \c glLoadMatrixd(q.matrix()) to apply the TQuaternion rotation to the current OpenGL matrix.

		 See also getMatrix(), getRotationMatrix() and inverseMatrix().

		 \attention The result is only valid until the next call to matrix(). Use it immediately (as shown
		 above) or consider using getMatrix() instead.

		 \attention The matrix is given in OpenGL format (row-major order) and is the transpose of the
		 actual mathematical European representation. Consider using getRotationMatrix() instead. */
		const Real* matrix() const
		{
			static Real m[4][4];
			getMatrix44(m);
			return (const Real*)(m);
		}


		/*! Fills \p m with the OpenGL representation of the TQuaternion rotation.

		Use matrix() if you do not need to store this matrix and simply want to alter the current OpenGL
		matrix. See also getInverseMatrix() and Frame::getMatrix(). */
        //void getMatrix(Real m[4][4]) const
        template<typename MAT>
        void getMatrix44(MAT& m) const
		{
			const Real q00 = 2.0l * q[0] * q[0];
			const Real q11 = 2.0l * q[1] * q[1];
			const Real q22 = 2.0l * q[2] * q[2];

			const Real q01 = 2.0l * q[0] * q[1];
			const Real q02 = 2.0l * q[0] * q[2];
			const Real q03 = 2.0l * q[0] * q[3];

			const Real q12 = 2.0l * q[1] * q[2];
			const Real q13 = 2.0l * q[1] * q[3];

			const Real q23 = 2.0l * q[2] * q[3];

			m[0][0] = 1.0l - q11 - q22;
			m[1][0] =        q01 - q23;
			m[2][0] =        q02 + q13;

			m[0][1] =        q01 + q23;
			m[1][1] = 1.0l - q22 - q00;
			m[2][1] =        q12 - q03;

			m[0][2] =        q02 - q13;
			m[1][2] =        q12 + q03;
			m[2][2] = 1.0l - q11 - q00;

            m[0][3] = 0.0l;
            m[1][3] = 0.0l;
            m[2][3] = 0.0l;

            m[3][0] = 0.0l;
            m[3][1] = 0.0l;
            m[3][2] = 0.0l;
            m[3][3] = 1.0l;
        }

        /*! Fills \p m with the OpenGL representation of the TQuaternion rotation.

        Use matrix() if you do not need to store this matrix and simply want to alter the current OpenGL
        matrix. See also getInverseMatrix() and Frame::getMatrix(). */
        //void getMatrix(Real m[4][4]) const
        template<typename MAT>
        void getMatrix33(MAT& m) const
        {
            const Real q00 = 2.0l * q[0] * q[0];
            const Real q11 = 2.0l * q[1] * q[1];
            const Real q22 = 2.0l * q[2] * q[2];

            const Real q01 = 2.0l * q[0] * q[1];
            const Real q02 = 2.0l * q[0] * q[2];
            const Real q03 = 2.0l * q[0] * q[3];

            const Real q12 = 2.0l * q[1] * q[2];
            const Real q13 = 2.0l * q[1] * q[3];

            const Real q23 = 2.0l * q[2] * q[3];

            m[0][0] = 1.0l - q11 - q22;
            m[1][0] =        q01 - q23;
            m[2][0] =        q02 + q13;

            m[0][1] =        q01 + q23;
            m[1][1] = 1.0l - q22 - q00;
            m[2][1] =        q12 - q03;

            m[0][2] =        q02 - q13;
            m[1][2] =        q12 + q03;
            m[2][2] = 1.0l - q11 - q00;
        }

        void getMatrix16(Real m[16]) const
        {
            static Real mat[4][4];
            getMatrix44(mat);
            int count = 0;
            for (int i=0; i<4; ++i)
                for (int j=0; j<4; ++j)
                    m[count++] = mat[i][j];
        }


		/*! Fills \p m with the 3x3 rotation matrix associated with the TQuaternion.
		  See also getInverseRotationMatrix().
		  \attention \p m uses the European mathematical representation of the rotation matrix. Use matrix()
		  and getMatrix() to retrieve the OpenGL transposed version. */
		void getRotationMatrix(Real m[3][3]) const
		{
			static Real mat[4][4];
			getMatrix(mat);
			for (int i=0; i<3; ++i)
				for (int j=0; j<3; ++j)
					// Beware of transposition
					m[i][j] = mat[j][i];
		}

		/*! Returns the associated 4x4 OpenGL \e inverse rotation matrix. This is simply the matrix() of the
		  inverse().
		  \attention The result is only valid until the next call to inverseMatrix(). Use it immediately (as
		  in \c glMultMatrixd(q.inverseMatrix())) or use getInverseMatrix() instead.
		  \attention The matrix is given in OpenGL format (row-major order) and is the transpose of the
		  actual mathematical European representation. Consider using getInverseRotationMatrix() instead. */
		const Real* inverseMatrix() const
		{
			static Real m[4][4];
			getInverseMatrix(m);
			return (const Real*)(m);
		}

		/*! Fills \p m with the OpenGL matrix corresponding to the inverse() rotation.

		Use inverseMatrix() if you do not need to store this matrix and simply want to alter the current
		OpenGL matrix. See also getMatrix(). */
		void getInverseMatrix(Real m[4][4]) const
		{
			inverse().getMatrix(m);
		}


		void getInverseMatrix(Real m[16]) const
		{
			inverse().getMatrix(m);
		}


		/*! \p m is set to the 3x3 \e inverse rotation matrix associated with the TQuaternion.

		 \attention This is the classical mathematical rotation matrix. The OpenGL format uses its
		 transposed version. See inverseMatrix() and getInverseMatrix(). */
		void getInverseRotationMatrix(Real m[3][3]) const
		{
			static Real mat[4][4];
			getInverseMatrix(mat);
			for (int i=0; i<3; ++i)
				for (int j=0; j<3; ++j)
					// Beware of transposition
					m[i][j] = mat[j][i];
		}
		//@}


		/*! @name Slerp interpolation */
		//@{
		static TQuaternion slerp(const TQuaternion& a, const TQuaternion& b, Real t, bool allowFlip=true)
		{
			Real cosAngle = TQuaternion::dot(a, b);

			Real c1, c2;
			// Linear interpolation for close orientations
			if ((1.0 - fabs(cosAngle)) < 0.01)
			{
				c1 = 1.0 - t;
				c2 = t;
			}
			else
			{
				// Spherical interpolation
				Real angle    = acos(fabs(cosAngle));
				Real sinAngle = sin(angle);
				c1 = sin(angle * (1.0 - t)) / sinAngle;
				c2 = sin(angle * t) / sinAngle;
			}

			// Use the shortest path
			if (allowFlip && (cosAngle < 0.0))
				c1 = -c1;

			return TQuaternion(c1*a[0] + c2*b[0], c1*a[1] + c2*b[1], c1*a[2] + c2*b[2], c1*a[3] + c2*b[3]);
		}


		static TQuaternion squad(const TQuaternion& a, const TQuaternion& tgA, const TQuaternion& tgB, const TQuaternion& b, Real t)
		{
			TQuaternion ab = TQuaternion::slerp(a, b, t);
			TQuaternion tg = TQuaternion::slerp(tgA, tgB, t, false);
			return TQuaternion::slerp(ab, tg, 2.0*t*(1.0-t), false);
		}



		/*! Returns the "dot" product of \p a and \p b: a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]. */
		static Real dot(const TQuaternion& a, const TQuaternion& b) { return a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3]; }

		TQuaternion log()
		{
			Real len = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

			if (len < 1E-6)
				return TQuaternion(q[0], q[1], q[2], 0.0);
			else
			{
				Real coef = acos(q[3]) / len;
				return TQuaternion(q[0]*coef, q[1]*coef, q[2]*coef, 0.0);
			}
		}


		TQuaternion exp()
		{
			Real theta = sqrt(q[0]*q[0] + q[1]*q[1] + q[2]*q[2]);

			if (theta < 1E-6)
				return TQuaternion(q[0], q[1], q[2], cos(theta));
			else
			{
				Real coef = sin(theta) / theta;
				return TQuaternion(q[0]*coef, q[1]*coef, q[2]*coef, cos(theta));
			}
		}


		static TQuaternion lnDif(const TQuaternion& a, const TQuaternion& b)
		{
			TQuaternion dif = a.inverse()*b;
			dif.normalize();
			return dif.log();
		}


		static TQuaternion squadTangent(const TQuaternion& before, const TQuaternion& center, const TQuaternion& after)
		{
			TQuaternion l1 = TQuaternion::lnDif(center,before);
			TQuaternion l2 = TQuaternion::lnDif(center,after);
			TQuaternion e;
			for (int i=0; i<4; ++i)
				e.q[i] = -0.25 * (l1.q[i] + l2.q[i]);
			e = center*(e.exp());

			// if (TQuaternion::dot(e,b) < 0.0)
			// e.negate();

			return e;
		}
		//@}

		/*! @name Random TQuaternion */
		//@{
		static TQuaternion randomQuaternion()
		{
			// The rand() function is not very portable and may not be available on your system.
			// Add the appropriate include or replace by an other random function in case of problem.
			Real seed = rand()/(Real)RAND_MAX;
			Real r1 = sqrt(1.0 - seed);
			Real r2 = sqrt(seed);
			Real t1 = 2.0 * M_PI * (rand()/(Real)RAND_MAX);
			Real t2 = 2.0 * M_PI * (rand()/(Real)RAND_MAX);
			return TQuaternion(sin(t1)*r1, cos(t1)*r1, sin(t2)*r2, cos(t2)*r2);
		}
		//@}

		/*! @name XML representation */
		//@{
		//~ explicit TQuaternion(const QDomElement& element);
		//~ QDomElement domElement(const QString& name, QDomDocument& document) const;
		//~ void initFromDOMElement(const QDomElement& element);
		//@}

	#ifdef DOXYN
		/*! @name Output stream */
		//@{
		/*! Output stream operator. Enables debugging code like:
		\code
		TQuaternion rot(...);
		cout << "Rotation=" << rot << endl;
		\endcode */
		inline std::ostream& operator<<(std::ostream& o, const qglviewer::Vec&);
		//@}
	#endif

		//template<typename C>
		friend inline std::ostream& operator<<(std::ostream& o, const TQuaternion& Q)
		{
			o << "(" << Q[0] << ',' << Q[1] << ',' << Q[2] << ',' << Q[3] << ")";
			return o;
		}


	private:
		/*! The internal data representation is private, use operator[] to access values. */
		Real q[4];
	};

	typedef TQuaternion<float, Vector> Quaternion;
	//typedef TQuaternion<double, Vector> Quaterniond;




#endif // _QUATERNION_H

