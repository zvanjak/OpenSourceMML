#if !defined MML_GEOMETRY_H
#define MML_GEOMETRY_H

#include "MMLBase.h"


namespace MML
{
	class Point2Cartesian
	{
	private:
		Real _x, _y;

	public:
		Real  X() const { return _x; }
		Real& X() { return _x; }
		Real  Y() const { return _y; }
		Real& Y() { return _y; }

		Point2Cartesian() : _x(0), _y(0) {}
		Point2Cartesian(Real x, Real y) : _x(x), _y(y) {}

		Real Dist(const Point2Cartesian& b) const { return sqrt(POW2(b.X() - X()) + POW2(b.Y() - Y())); }

		Point2Cartesian operator+(const Point2Cartesian& b) const { return Point2Cartesian(X() + b.X(), Y() + b.Y()); }
		Point2Cartesian operator-(const Point2Cartesian& b) const { return Point2Cartesian(X() - b.X(), Y() - b.Y()); }

		friend Point2Cartesian operator*(const Point2Cartesian& a, Real b) { return Point2Cartesian(a.X() * b, a.Y() * b); }
		friend Point2Cartesian operator*(Real a, const Point2Cartesian& b) { return Point2Cartesian(a * b.X(), a * b.Y()); }
		friend Point2Cartesian operator/(const Point2Cartesian& a, Real b) { return Point2Cartesian(a.X() / b, a.Y() / b); }
	};	

	class Point2Polar
	{
	private:
		Real _r, _phi;

	public:
		Real  R() const { return _r; }
		Real& R() { return _r; }
		Real  Phi() const { return _phi; }
		Real& Phi() { return _phi; }

		Point2Polar() : _r(0), _phi(0) {}
		Point2Polar(Real r, Real phi) : _r(r), _phi(phi) {}

		static Point2Polar GetFromCartesian(const Point2Cartesian &pnt)
		{
			Real r = sqrt(POW2(pnt.X()) + POW2(pnt.Y()));
			Real phi = atan2(pnt.Y(), pnt.X());
			return Point2Polar(r, phi);
		}
		Point2Cartesian GetCartesian() const {  return Point2Cartesian(R() * cos(Phi()), R() * sin(Phi())); }

		Real Dist(const Point2Polar& b) const { return sqrt(R() * R() + b.R() * b.R() - 2 * R() * b.R() * cos(b.Phi() - Phi())); }
	};

	class Point3Cartesian
	{
	private:
		Real _x, _y, _z;

	public:
		Real  X() const { return _x; }
		Real& X() { return _x; }
		Real  Y() const { return _y; }
		Real& Y() { return _y; }
		Real  Z() const { return _z; }
		Real& Z() { return _z; }

		Point3Cartesian() : _x(0), _y(0), _z(0) {}
		Point3Cartesian(Real x, Real y, Real z) : _x(x), _y(y), _z(z) {}

		Real Dist(const Point3Cartesian& b) const { return sqrt(POW2(b.X() - X()) + POW2(b.Y() - Y()) + POW2(b.Z() - Z())); }

		bool	operator==(const Point3Cartesian& b) const {return (X() == b.X()) && (Y() == b.Y()) && (Z() == b.Z()); }
		bool	operator!=(const Point3Cartesian& b) const {return (X() != b.X()) || (Y() != b.Y()) || (Z() != b.Z()); }

		Point3Cartesian operator+(const Point3Cartesian& b) const { return Point3Cartesian(X() + b.X(), Y() + b.Y(), Z() + b.Z()); }
		Point3Cartesian operator-(const Point3Cartesian& b) const { return Point3Cartesian(X() - b.X(), Y() - b.Y(), Z() - b.Z()); }

		friend Point3Cartesian operator*(const Point3Cartesian& a, Real b) { return Point3Cartesian(a.X() * b, a.Y() * b, a.Z() * b); }
		friend Point3Cartesian operator*(Real a, const Point3Cartesian& b) { return Point3Cartesian(a * b.X(), a * b.Y(), a * b.Z()); }
		friend Point3Cartesian operator/(const Point3Cartesian& a, Real b) { return Point3Cartesian(a.X() / b, a.Y() / b, a.Z() / b); }
	};

	class Triangle
	{
	private:
		Real _a, _b, _c;

	public:
		Real  A() const { return _a; }
		Real& A() { return _a; }
		Real  B() const { return _b; }
		Real& B() { return _b; }
		Real  C() const { return _c; }
		Real& C() { return _c; }

		Triangle() : _a(0.0), _b(0.0), _c(0.0){}
		Triangle(Real a, Real b, Real c) : _a(a), _b(b), _c(c) {}

		Real Area() const
		{
			Real s = (_a + _b + _c) / 2.0;
			return sqrt(s * (s - _a) * (s - _b) * (s - _c));
		}
		bool IsRight() const
		{
			return ( POW2(_a) + POW2(_b) == POW2(_c)) || 
							(POW2(_a) + POW2(_c) == POW2(_b)) || 
							(POW2(_b) + POW2(_c) == POW2(_a) );
		}
		bool IsIsosceles() const
		{
			return (_a == _b) || (_a == _c) || (_b == _c);
		}
		bool IsEquilateral() const
		{
			return (_a == _b) && (_a == _c);
		}
	};

	typedef Point2Cartesian Pnt2Cart;
	typedef Point2Polar			Pnt2Pol;
	typedef Point3Cartesian Pnt3Cart;
}
#endif // MML_GEOMETRY_H
