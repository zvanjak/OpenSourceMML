#if !defined  MML_VECTORN_H
#define MML_VECTORN_H

#include "MMLBase.h"
#include "base/Geometry.h"

namespace MML
{
	template<class Type, int N>
	class VectorN
	{
	protected:
		Type  _val[N] = { 0 };

	public:
		///////////////////////          Constructors and destructor       //////////////////////
		VectorN() {}
		explicit VectorN(const Type& init_val) {
			for (int i = 0; i < N; ++i)
				_val[i] = init_val;
		}
		explicit VectorN(std::initializer_list<Type> list)
		{
			int count = 0;
			for (auto element : list)
			{
				_val[count] = element;
				++count;

				if (count >= N)
					break;
			}
		}
		explicit VectorN(std::vector<Type> list)
		{
			int count{ 0 };
			for (auto element : list)
			{
				_val[count] = element;
				++count;

				if (count >= N)
					break;
			}
		}
		explicit VectorN(Type* vals)
		{
			for (int i = 0; i < N; ++i)
				_val[i] = vals[i];
		}
		static VectorN GetUnitVector(int indUnit)
		{
			VectorN ret;
			ret[indUnit] = 1.0;
			return ret;
		}

		////////////////////////            Standard stuff             ////////////////////////
		int size() const { return N; }
		void clear() {
			for (int i = 0; i < N; ++i)
				_val[i] = Type{ 0 };
		}

		VectorN GetAsUnitVector() const
		{
			return VectorN{ (*this) / NormL2() };
		}
		void MakeUnitVector()
		{
			(*this) = GetAsUnitVector();
		}

		bool IsEqualTo(const VectorN& b, Type eps = Defaults::VectorEqualityPrecision) const
		{
			for (int i = 0; i < N; i++)
			{
				if (Abs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}
		static bool AreEqual(const VectorN& a, const VectorN& b, Type eps = Defaults::VectorEqualityPrecision)
		{
			return a.IsEqualTo(b, eps);
		}
		bool IsNullVec() const
		{
			for (int i = 0; i < N; i++)
				if (_val[i] != 0.0)
					return false;

			return true;
		}
		///////////////////////////            Operators             ///////////////////////////
		Type& operator[](int n) { return _val[n]; }
		Type  operator[](int n) const { return _val[n]; }

		// checked access
		Type& at(int n) {
			if (n < 0 || n >= N)
				throw VectorDimensionError("VectorN::at - index out of bounds", N, n);
			else
				return _val[n];
		}
		Type  at(int n) const {
			if (n < 0 || n >= N)
				throw VectorDimensionError("VectorN::at - index out of bounds", N, n);
			else
				return _val[n];
		}

		VectorN operator-() const        // unary minus
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = Type{ -1 } *_val[i];
			return ret;
		}
		VectorN operator+(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] + b._val[i];
			return ret;
		}
		VectorN operator-(const VectorN& b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] - b._val[i];
			return ret;
		}
		bool operator==(const VectorN& b) const
		{
			for (int i = 0; i < size(); i++)
				if ((*this)[i] != b[i])
					return false;

			return true;
		}

		VectorN operator*(const Type &b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] * b;
			return ret;
		}
		VectorN operator/(const Type &b) const
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = _val[i] / b;
			return ret;
		}		
		friend VectorN operator*(Type a, const VectorN<Type, N>& b)
		{
			VectorN ret;
			for (int i = 0; i < N; i++)
				ret._val[i] = a * b[i];
			return ret;
		}

		//////////////////////                 Operations                 ///////////////////////
		Type ScalarProductCartesian(const VectorN& b) const
		{
			Type product{ 0.0 };
			for (int i = 0; i < N; i++)
				product += (*this)[i] * b[i];
			return product;
		}
		Type NormL2() const
		{
			Type norm{ 0.0 };
			for (int i = 0; i < N; i++)
				norm += (*this)[i] * (*this)[i];
			return std::sqrt(norm);
		}
		Type AngleToVector(const VectorN& b) const
		{
			Type cosAngle = ScalarProductCartesian(b) / (NormL2() * b.NormL2());
			return std::acos(cosAngle);
		}

		///////////////////////////               I/O                 ///////////////////////////
		std::string to_string(int width, int precision) const
		{
			std::stringstream str;

			Print(str, width, precision);

			return str.str();
		}
		std::ostream& Print(std::ostream& stream, int width, int precision) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _val)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				stream << std::setw(width) << std::setprecision(precision) << x;
			}
			stream << "]";

			return stream;
		}
		std::ostream& PrintLine(std::ostream& stream, std::string msg, int width, int precision) const
		{
			stream << msg;
			Print(stream, width, precision);
			stream << std::endl;

			return stream;
		}

		std::ostream& Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _val)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				if (Abs(x) > zeroThreshold)
					stream << std::setw(width) << std::setprecision(precision) << x;
				else
					stream << std::setw(width) << std::setprecision(precision) << 0.0;
			}
			stream << "]";

			return stream;
		}
		friend std::ostream& operator<<(std::ostream& stream, const VectorN<Type, N>& a)
		{
			a.Print(stream, 15, 10);

			return stream;
		}
	};

	class Vector2Cartesian : public VectorN<Real, 2>
	{
	public:
		Vector2Cartesian() {}
		Vector2Cartesian(Real x, Real y)
		{
			_val[0] = x;
			_val[1] = y;
		}
		Vector2Cartesian(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}
		Vector2Cartesian(const Point2Cartesian& a, const Point2Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
		}

		Real  X() const { return _val[0]; }
		Real& X() { return _val[0]; }
		Real  Y() const { return _val[1]; }
		Real& Y() { return _val[1]; }
			
		// For Cartesian vector, we will enable operator* to represent standard scalar product
		Real operator*(const Vector2Cartesian& b)
		{
			return this->ScalarProductCartesian(b);
		}

		friend Vector2Cartesian operator*(const Vector2Cartesian& a, Real b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a[i] * b;
			return ret;
		}
		friend Vector2Cartesian operator*(Real a, const Vector2Cartesian& b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a * b[i];
			return ret;
		}
		friend Vector2Cartesian operator/(const Vector2Cartesian& a, Real b)
		{
			Vector2Cartesian ret;
			for (int i = 0; i < 2; i++)
				ret._val[i] = a[i] / b;
			return ret;
		}

		Vector2Cartesian GetUnitVector() const
		{
			VectorN<Real, 2> res = (*this) / NormL2();

			return Vector2Cartesian(res[0], res[1]);
		}
		Vector2Cartesian GetAsUnitVectorAtPos(const Vector2Cartesian& pos) const
		{
			return Vector2Cartesian{ (*this) / NormL2() };
		}

		friend Point2Cartesian operator+(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() + b[0], a.Y() + b[1]); }
		friend Point2Cartesian operator-(const Point2Cartesian& a, const Vector2Cartesian& b) { return Point2Cartesian(a.X() - b[0], a.Y() - b[1]); }
	};

	class Vector2Polar : public VectorN<Real, 2>
	{
	public:
		Real  R() const { return _val[0]; }
		Real& R() { return _val[0]; }
		Real  Phi() const { return _val[1]; }
		Real& Phi() { return _val[1]; }

		Vector2Polar() {}
		Vector2Polar(Real r, Real phi)
		{
			_val[0] = r;
			_val[1] = phi;
		}
		Vector2Polar(const VectorN<Real, 2>& b) : VectorN<Real, 2>{ b[0], b[1] } {}

		Vector2Polar GetAsUnitVectorAtPos(const Vector2Polar& pos) const
		{
			// TODO 1.1 - BIG!!!
			return Vector2Polar{ (*this) / NormL2() };
		}
	};

	class Vector3Cartesian : public VectorN<Real, 3>
	{
	public:
		Real  X() const { return _val[0]; }
		Real& X() { return _val[0]; }
		Real  Y() const { return _val[1]; }
		Real& Y() { return _val[1]; }
		Real  Z() const { return _val[2]; }
		Real& Z() { return _val[2]; }

		Vector3Cartesian() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Cartesian(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b } {}
		Vector3Cartesian(Real x, Real y, Real z) : VectorN<Real, 3>{ x, y, z } {}
		Vector3Cartesian(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }
		Vector3Cartesian(const Point3Cartesian& a, const Point3Cartesian& b)
		{
			_val[0] = b.X() - a.X();
			_val[1] = b.Y() - a.Y();
			_val[2] = b.Z() - a.Z();
		}

		// For Cartesian vector, we will enable operator* to represent standard scalar product
		Real operator*(const Vector3Cartesian& b)
		{
			return this->ScalarProductCartesian(b);
		}

		friend Vector3Cartesian operator*(const Vector3Cartesian& a, Real b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a[i] * b;
			return ret;
		}
		friend Vector3Cartesian operator*(Real a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a * b[i];
			return ret;
		}
		friend Vector3Cartesian operator/(const Vector3Cartesian& a, Real b)
		{
			Vector3Cartesian ret;
			for (int i = 0; i < 3; i++)
				ret._val[i] = a[i] / b;
			return ret;
		}

		friend Point3Cartesian operator+(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() + b[0], a.Y() + b[1], a.Z() + b[2]); }
		friend Point3Cartesian operator-(const Point3Cartesian& a, const Vector3Cartesian& b) { return Point3Cartesian(a.X() - b[0], a.Y() - b[1], a.Z() - b[2]); }

		Point3Cartesian getAsPoint()
		{
			return Point3Cartesian(_val[0], _val[1], _val[2]);
		}

		bool IsParallelTo(const Vector3Cartesian& b, Real eps = 1e-15) const
		{
			Real norm1 = NormL2();
			Real norm2 = b.NormL2();

			return std::abs(X() / norm1 - b.X() / norm2) < eps &&
				std::abs(Y() / norm1 - b.Y() / norm2) < eps &&
				std::abs(Z() / norm1 - b.Z() / norm2) < eps;
		}
		bool IsPerpendicularTo(const Vector3Cartesian& b, Real eps = 1e-15) const
		{
			if (std::abs(ScalarProd(*this, b)) < eps)
				return true;
			else
				return false;
		}
		Real AngleToVector(const Vector3Cartesian& b)
		{
			Real cos_phi = ScalarProd(*this, b) / (NormL2() * b.NormL2());

			return acos(cos_phi);
		}

		Vector3Cartesian GetAsUnitVector() const
		{
			return Vector3Cartesian{ (*this) / NormL2() };
		}
		Vector3Cartesian GetAsUnitVectorAtPos(const Vector3Cartesian& pos) const
		{
			return Vector3Cartesian{ (*this) / NormL2() };
		}

		friend Real ScalarProd(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			return a.ScalarProductCartesian(b);
		}
		friend Vector3Cartesian VectorProd(const Vector3Cartesian& a, const Vector3Cartesian& b)
		{
			Vector3Cartesian ret;

			ret.X() = a.Y() * b.Z() - a.Z() * b.Y();
			ret.Y() = a.Z() * b.X() - a.X() * b.Z();
			ret.Z() = a.X() * b.Y() - a.Y() * b.X();

			return ret;
		}
	};

	class Vector3Spherical : public VectorN<Real, 3>
	{
	public:
		Real  R()     const { return _val[0]; }
		Real& R() { return _val[0]; }
		Real  Theta() const { return _val[1]; }
		Real& Theta() { return _val[1]; }
		Real  Phi()   const { return _val[2]; }
		Real& Phi() { return _val[2]; }

		Vector3Spherical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Spherical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Spherical(Real r, Real theta, Real phi) : VectorN<Real, 3>{ r, theta, phi } {}
		Vector3Spherical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		// TODO - HIGH, HARD, osnovne operacije +, -
		// sve operaciju pretpostavlju da se odvijaju NA ISTOJ TOCKI U PROSTORU
		Vector3Spherical GetAsUnitVectorAtPos(const Vector3Spherical& pos) const
		{
			// TODO 1.1 - VERIFY this!!!
			return Vector3Spherical{ R(), Theta() / pos.R(), Phi() / (pos.R() * sin(pos.Theta())) };
		}

		std::ostream& PrintDeg(std::ostream& stream, int width, int precision) const
		{
			stream << "[ ";
			stream << std::fixed << std::setw(width) << std::setprecision(precision);
			stream << R();
			stream << ", " << Theta() * 180.0 / Constants::PI;
			stream << ", " << Phi() * 180.0 / Constants::PI << " ]" << std::endl;

			return stream;
		}
	};

	class Vector3Cylindrical : public VectorN<Real, 3>
	{
	public:
		Real    R()   const { return _val[0]; }
		Real& R() { return _val[0]; }
		Real    Phi() const { return _val[1]; }
		Real& Phi() { return _val[1]; }
		Real    Z()   const { return _val[2]; }
		Real& Z() { return _val[2]; }

		Vector3Cylindrical() : VectorN<Real, 3>{ 0.0, 0.0, 0.0 } {}
		Vector3Cylindrical(const VectorN<Real, 3>& b) : VectorN<Real, 3>{ b[0], b[1], b[2] } {}
		Vector3Cylindrical(Real r, Real phi, Real z) : VectorN<Real, 3>{ r, phi, z } {}
		Vector3Cylindrical(std::initializer_list<Real> list) : VectorN<Real, 3>(list) { }

		// TODO - MED, implement Vector3Cylindrical GetAsUniTVector()
		// TODO - MED, razmisliti generalno vektor i vektor at point
		Vector3Cylindrical GetAsUnitVectorAtPos(const Vector3Cylindrical& pos) const
		{
			return Vector3Cylindrical{ R(), Phi() / pos.R(), Z() };
		}
	};

	typedef Vector2Cartesian    Vec2Cart;
	typedef Vector3Cartesian    Vec3Cart;
	typedef Vector3Spherical    Vec3Sph;
	typedef Vector3Cylindrical  Vec3Cyl;

	typedef VectorN<float, 2> Vec2Flt;
	typedef VectorN<float, 3> Vec3Flt;
	typedef VectorN<float, 4> Vec4Flt;

	typedef VectorN<double, 2> Vec2Dbl;
	typedef VectorN<double, 3> Vec3Dbl;
	typedef VectorN<double, 4> Vec4Dbl;

	typedef VectorN<Complex, 2> Vec2Complex;
	typedef VectorN<Complex, 3> Vec3Complex;
	typedef VectorN<Complex, 4> Vec4Complex;

	typedef VectorN<float, 2> Vec2F;
	typedef VectorN<float, 3> Vec3F;
	typedef VectorN<float, 4> Vec4F;

	typedef VectorN<double, 2> Vec2D;
	typedef VectorN<double, 3> Vec3D;
	typedef VectorN<double, 4> Vec4D;

	typedef VectorN<Complex, 2> Vec2C;
	typedef VectorN<Complex, 3> Vec3C;
	typedef VectorN<Complex, 4> Vec4C;
}

#endif