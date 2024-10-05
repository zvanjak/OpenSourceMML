#if !defined  MML_Vector_H
#define MML_Vector_H

#include "MMLBase.h"

namespace MML
{
	template<class Type>
	class	Vector
	{
	private:
		std::vector<Type> _elems;

	public:
		///////////////////////          Constructors and destructor       //////////////////////
		Vector() {}
		explicit Vector(int n) {
			if(n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n);
		}
		explicit Vector(int n, const Type &val) {
			if (n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n, val);
		}
		explicit Vector(int n, Type* vals) 
		{
			if (n < 0)
				throw VectorInitializationError("Vector::Vector - negative size", n);

			_elems.resize(n);
			for (int i = 0; i < n; ++i)
				_elems[i] = vals[i];
		}
		explicit Vector(std::vector<Type> values) : _elems(values) {}
		explicit Vector(std::initializer_list<Type> list) : _elems(list) {}

		// not really needed, but let's be explicit
		Vector(const Vector& b) = default; 
		Vector(Vector&& b) = default;
		Vector& operator=(const Vector& b) = default; 
		Vector& operator=(Vector&& b) = default;

		////////////////            Standard std::vector stuff             ////////////////////
		int  size() const { return (int)_elems.size(); }
		bool empty() const { return _elems.empty(); }

		void clear() { _elems.clear(); }
		void resize(int newLen)		{ _elems.resize(newLen); }

		////////////////////////            Standard stuff             ////////////////////////
		static Vector GetUnitVector(int dimVec, int indUnit)
		{
			if (indUnit < 0 || indUnit >= dimVec)
				throw VectorDimensionError("Vector::GetUnitVector - wrong unit index", dimVec, indUnit);

			Vector ret(dimVec);
			ret[indUnit] = Type{ 1.0 };
			return ret;
		}
		
		static bool AreEqual(const Vector& a, const Vector& b, Type eps = Defaults::VectorEqualityPrecision)
		{
			return a.IsEqualTo(b, eps);
		}
		bool IsEqualTo(const Vector& b, Real eps = Defaults::VectorEqualityPrecision) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::IsEqual - vectors must be equal size", size(), b.size());

			for (int i = 0; i < size(); i++)
			{
				if (Abs((*this)[i] - b[i]) > eps)
					return false;
			}
			return true;
		}

		///////////////////////////            Operators             ///////////////////////////
		Type&       operator[](int n)       { return _elems[n]; }
		const Type& operator[](int n) const { return _elems[n]; }

		// checked access
		Type& at(int n)	{
			if(n < 0 || n >= size())
				throw VectorDimensionError("Vector::at - index out of bounds", size(), n);
			else
				return _elems[n];
		}
		Type  at(int n) const { 
			if(n < 0 || n >= size())
				throw VectorDimensionError("Vector::at - index out of bounds", size(), n);
			else
				return _elems[n];
		}

		Vector operator-() const         // unary minus
		{
			Vector ret(size());
			for (int i = 0; i < size(); i++)
				ret._elems[i] = Type{ -1 } * (*this)[i];
			return ret;
		}
		Vector operator+(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator+() - vectors must be equal size", size(), b.size());

			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = (*this)[i] + b._elems[i];
			return ret;
		}
		Vector operator-(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator-() - vectors must be equal size", size(), b.size());

			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = (*this)[i] - b._elems[i];
			return ret;
		}
		bool   operator==(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::operator==() - vectors must be equal size", size(), b.size());

			for (int i = 0; i < size(); i++)
				if ((*this)[i] != b[i])
					return false;

			return true;
		}

		Vector operator*(Type b)
		{
			Vector ret(size());;
			for (int i = 0; i < size(); i++)
				ret._elems[i] = b * _elems[i];
			return ret;
		}
		Vector operator/(Type b)
		{
			Vector ret(size());
			for (int i = 0; i < size(); i++)
				ret._elems[i] = _elems[i] / b;
			return ret;
		}
		
		friend Vector operator*(Type a, const Vector& b)
		{
			Vector ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret._elems[i] = a * b._elems[i];
			return ret;
		}
		
		//////////////////////                 Operations                 ///////////////////////
		Type ScalarProductCartesian(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::ScalarProductCartesian - vectors must be equal size", size(), b.size());

			Type product{ 0.0 };
			for (int i = 0; i < size(); i++)
				product += (*this)[i] * b[i];
			return product;
		}
		Type NormL2() const
		{
			Type norm{ 0.0 };
			for (int i = 0; i < size(); i++)
				norm += (*this)[i] * (*this)[i];
			return std::sqrt(norm);
		}
		Type AngleToVector(const Vector& b) const
		{
			if (size() != b.size())
				throw VectorDimensionError("Vector::AngleToVector - vectors must be equal size", size(), b.size());

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
		void Print(std::ostream& stream, int width, int precision) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _elems)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

				stream << std::setw(width) << std::setprecision(precision) << x;
			}
			stream << "]";
		}
		std::ostream& Print(std::ostream& stream, int width, int precision, Real zeroThreshold) const
		{
			stream << "[";
			bool first = true;
			for (const Type& x : _elems)
			{
				if (!first)
					stream << ", ";
				else
					first = false;

        if( Abs(x) > zeroThreshold )
				  stream << std::setw(width) << std::setprecision(precision) << x;
        else
          stream << std::setw(width) << std::setprecision(precision) << 0.0;
			}
			stream << "]";

			return stream;
		}      
		friend std::ostream& operator<<(std::ostream& stream, const Vector& a)
		{
			a.Print(stream, 15, 10);

			return stream;
		}
	};

	typedef Vector<int>     VectorInt;
	typedef Vector<float>   VectorFlt;
	typedef Vector<double>  VectorDbl;
	typedef Vector<Complex> VectorComplex;

	typedef Vector<int>     VecI;
	typedef Vector<float>   VecF;
	typedef Vector<double>  VecD;
	typedef Vector<Complex> VecC;
}

#endif // MML_Vector_H