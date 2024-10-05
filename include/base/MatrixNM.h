#if !defined  MML_MATRIXNM_H
#define MML_MATRIXNM_H

#include "MMLBase.h"

#include "VectorN.h"

namespace MML
{
	template <class Type, int N, int M>
	class MatrixNM
	{
	public:
		Type _vals[N][M] = { {0} };

	public:
		//////////////////////////             Constructors           /////////////////////////
		MatrixNM() {}
		MatrixNM(std::initializer_list<Type> values)
		{
			auto val = values.begin();
			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					if (val != values.end())
						_vals[i][j] = *val++;
					else
						_vals[i][j] = 0.0;
		}
		MatrixNM(const MatrixNM& m)
		{
			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m._vals[i][j];
		}
		MatrixNM(const Type& m)        // initialize as diagonal matrix
		{
			for (int i = 0; i < N; i++)
				_vals[i][i] = Type{ m };
		}

		typedef Type value_type;      // make T available externally

		////////////////////////            Standard stuff             ////////////////////////
		int RowNum() const { return N; }
		int ColNum() const { return M; }

		static MatrixNM GetUnitMatrix()
		{
			MatrixNM unitMat;

			for (int i = 0; i < N; i++)
				unitMat._vals[i][i] = 1.0;

			return unitMat;
		}
		void   MakeUnitMatrix(void)
		{
			if (RowNum() == ColNum())
			{
				for (int i = 0; i < RowNum(); i++)
					for (int j = 0; j < ColNum(); j++)
						if (i == j)
							_vals[i][j] = 1;
						else
							_vals[i][j] = 0;
			}
			else
				throw MatrixDimensionError("MatrixNM::MakeUnitMatrix - must be square matrix", N, M, -1, -1);
		}

		MatrixNM GetLower(bool includeDiagonal = true) const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetLower - must be square matrix", N, M, -1, -1);

			MatrixNM ret;
			for (int i = 0; i < RowNum(); i++)
			{
				if (includeDiagonal)
					for (int j = 0; j < i + 1; j++)
						ret[i][j] = _vals[i][j];
				else
					for (int j = 0; j < i; j++)
						ret[i][j] = _vals[i][j];
			}

			return ret;
		}
		MatrixNM GetUpper(bool includeDiagonal = true) const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("Matrix::GetUpper - must be square matrix", N, M, -1, -1);

			MatrixNM ret;
			for (int i = 0; i < RowNum(); i++)
			{
				if (includeDiagonal)
					for (int j = i; j < ColNum(); j++)
						ret[i][j] = _vals[i][j];
				else
					for (int j = i + 1; j < ColNum(); j++)
						ret[i][j] = _vals[i][j];
			}

			return ret;
		}

		/////////////////////          Vector-Matrix conversion           /////////////////////
		static MatrixNM<Type, 1, M> RowMatrixFromVector(const VectorN<Type, M>& b)
		{
			MatrixNM<Type, 1, M>  ret;
			for (int j = 0; j < M; j++)
				ret._vals[0][j] = b[j];

			return ret;
		}
		static MatrixNM<Type, N, 1> ColumnMatrixFromVector(const VectorN<Type, N>& b)
		{
			MatrixNM<Type, N, 1>  ret;
			for (int i = 0; i < M; i++)
				ret._vals[i][0] = b[i];

			return ret;
		}
		static MatrixNM<Type, N, N> DiagonalMatrixFromVector(const VectorN<Type, N>& b)
		{
			MatrixNM<Type, N, N> ret;
			for (int i = 0; i < N; i++)
				ret[i][i] = b[i];

			return ret;
		}
		static VectorN<Type, M> VectorFromRow(const MatrixNM<Type, N, M>& a, int rowInd)
		{
			VectorN<Type, M> ret;
			for (int j = 0; j < M; j++)
				ret[j] = a._vals[rowInd][j];

			return ret;
		}
		static VectorN<Type, N> VectorFromColumn(const MatrixNM<Type, N, M>& a, int colInd)
		{
			VectorN<Type, N> ret;
			for (int i = 0; i < N; i++)
				ret[i] = a._vals[i][colInd];

			return ret;
		}
		static VectorN<Type, N> VectorFromDiagonal(const MatrixNM<Type, N, N>& a)
		{
			VectorN<Type, N> ret;
			for (int i = 0; i < N; i++)
				ret[i] = a._vals[i][i];

			return ret;
		}

		/////////////////////            Assignment operators             ////////////////////
		MatrixNM& operator=(const MatrixNM& m)
		{
			if (this == &m)
				return *this;

			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m._vals[i][j];

			return *this;
		}
		MatrixNM& operator=(const Type& m)
		{
			if (this == &m)
				return *this;

			for (size_t i = 0; i < RowNum(); ++i)
				for (size_t j = 0; j < ColNum(); ++j)
					_vals[i][j] = m;

			return *this;
		}

		////////////////////            Access operators             ///////////////////////
		Type* operator[](int i) { return _vals[i]; }
		const Type* operator[](const int i) const { return _vals[i]; }

		Type  operator()(int i, int j) const { return _vals[i][j]; }
		Type& operator()(int i, int j) { return _vals[i][j]; }

		// version with checking bounds
		Type  ElemAt(int i, int j) const
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixNM::ElemAt", i, j, RowNum(), ColNum());

			return _vals[i][j];
		}
		Type& ElemAt(int i, int j)
		{
			if (i < 0 || i >= RowNum() || j < 0 || j >= ColNum())
				throw MatrixAccessBoundsError("MatrixNM::ElemAt", i, j, RowNum(), ColNum());

			return _vals[i][j];
		}

		////////////////////            Arithmetic operators             ////////////////////
		MatrixNM operator-()            // unary minus
		{
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = -_vals[i][j];
			return temp;
		}
		MatrixNM operator+(const MatrixNM& b) const
		{
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = b._vals[i][j] + _vals[i][j];
			return temp;
		}
		MatrixNM operator-(const MatrixNM& b) const
		{
			MatrixNM temp;
			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = 0; j < ColNum(); j++)
					temp._vals[i][j] = _vals[i][j] - b._vals[i][j];
			return temp;
		}
		template<int K>
		MatrixNM<Type, N, K>  operator*(const MatrixNM<Type, M, K>& b) const
		{
			MatrixNM<Type, N, K>	ret;

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++)
				{
					ret._vals[i][j] = 0;
					for (int k = 0; k < ColNum(); k++)
						ret._vals[i][j] += _vals[i][k] * b._vals[k][j];
				}

			return	ret;
		}

		MatrixNM operator*(Type b)
		{
			int	i, j;
			MatrixNM	ret(*this);

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret._vals[i][j] *= b;

			return ret;
		}
		MatrixNM& operator*=(Type b)
		{
			int	i, j;

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					_vals[i][j] *= b;

			return *this;
		}

		MatrixNM operator/(const Type &b) const
		{
			int	i, j;
			MatrixNM	ret(*this);

			for (i = 0; i < RowNum(); i++)
				for (j = 0; j < ColNum(); j++)
					ret._vals[i][j] /= b;

			return ret;
		}
		VectorN<Type, N> operator*(const VectorN<Type, M>& b) const
		{
			int	i, j;
			VectorN<Type, N>	ret;

			for (i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (j = 0; j < M; j++)
					ret[i] += _vals[i][j] * b[j];
			}

			return ret;
		}
		
		friend MatrixNM operator*(Type a, const MatrixNM& b)
		{
			int	i, j;
			MatrixNM	ret;

			for (i = 0; i < b.RowNum(); i++)
				for (j = 0; j < b.ColNum(); j++)
					ret._vals[i][j] = a * b._vals[i][j];

			return ret;
		}
		friend VectorN<Type, M> operator*(const VectorN<Type, N>& a, const MatrixNM<Type, N, M>& b)
		{
			int	i, j;
			VectorN<Type, M>	ret;

			for (i = 0; i < M; i++)
			{
				ret[i] = 0;
				for (j = 0; j < N; j++)
					ret[i] += a[j] * b._vals[j][i];
			}

			return ret;
		}

		///////////////////////            Equality operations             //////////////////////
		bool operator==(const MatrixNM& b) const
		{
			for (int i = 0; i < N; i++)
				for (int j = 0; j < M; j++)
					if (_vals[i][j] != b._vals[i][j])
						return false;

			return true;
		}
		bool operator!=(const MatrixNM& b) const
		{
			return !(*this == b);
		}

		bool IsEqual(const MatrixNM& b, Type eps = Defaults::MatrixEqualityPrecision) const
		{
			for (int i = 0; i < RowNum(); i++)
				for (int j = 0; j < ColNum(); j++)
					if (std::abs(_vals[i][j] - b._vals[i][j]) > eps)
						return false;

			return true;
		}
		bool AreEqual(const MatrixNM& a, const MatrixNM& b, Type eps = Defaults::MatrixEqualityPrecision) const
		{
			return a.IsEqual(b, eps);
		}

		///////////////////            Trace, Inverse & Transpose             ///////////////////
		Type   Trace() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Trace - must be square matrix", N, M, -1, -1);

			Type sum = 0;
			for (int i = 0; i < RowNum(); i++)
				sum += _vals[i][i];

			return sum;
		}

		void Invert()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Invert - must be square matrix", N, M, -1, -1);

			MatrixNM& a = *this;
			MatrixNM<Type, N, 1>  b;      // dummy rhs

			b(0, 0) = 1.0;

			int i, icol, irow, j, k, l, ll;
			Type big, dum, pivinv;

			int n = RowNum();
			int m = b.ColNum();
			std::vector<int> indxc(n), indxr(n), ipiv(n);
			for (j = 0; j < n; j++) ipiv[j] = 0;
			for (i = 0; i < n; i++) {
				big = 0.0;
				for (j = 0; j < n; j++)
					if (ipiv[j] != 1)
						for (k = 0; k < n; k++) {
							if (ipiv[k] == 0) {
								if (std::abs(a._vals[j][k]) >= big) {
									big = std::abs(a._vals[j][k]);
									irow = j;
									icol = k;
								}
							}
						}
				++(ipiv[icol]);
				if (irow != icol) {
					for (l = 0; l < n; l++) std::swap(a._vals[irow][l], a._vals[icol][l]);
					for (l = 0; l < m; l++) std::swap(b._vals[irow][l], b._vals[icol][l]);
				}
				indxr[i] = irow;
				indxc[i] = icol;

				if (a._vals[icol][icol] == 0.0)
					throw SingularMatrixError("MatrixNM::Invert, gaussj: Singular Matrix");

				pivinv = 1.0 / a._vals[icol][icol];
				a._vals[icol][icol] = 1.0;
				for (l = 0; l < n; l++) a._vals[icol][l] *= pivinv;
				for (l = 0; l < m; l++) b._vals[icol][l] *= pivinv;
				for (ll = 0; ll < n; ll++)
					if (ll != icol) {
						dum = a._vals[ll][icol];
						a._vals[ll][icol] = 0.0;
						for (l = 0; l < n; l++) a._vals[ll][l] -= a._vals[icol][l] * dum;
						for (l = 0; l < m; l++) b._vals[ll][l] -= b._vals[icol][l] * dum;
					}
			}
			for (l = n - 1; l >= 0; l--) {
				if (indxr[l] != indxc[l])
					for (k = 0; k < n; k++)
						std::swap(a._vals[k][indxr[l]], a._vals[k][indxc[l]]);
			}
		}
		MatrixNM GetInverse() const
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::GetInverse - must be square matrix", N, M, -1, -1);

			MatrixNM a(*this);              // making a copy, where inverse will be stored at the end

			a.Invert();

			return a;
		}

		void Transpose()
		{
			if (RowNum() != ColNum())
				throw MatrixDimensionError("MatrixNM::Transpose - inplace Transpose possible only for square  matrix", N, M, -1, -1);

			for (size_t i = 0; i < RowNum(); i++)
				for (size_t j = i + 1; j < ColNum(); j++)
					std::swap(_vals[i][j], _vals[j][i]);
		}
		MatrixNM<Type, M, N> GetTranspose() const
		{
			MatrixNM<Type, M, N> ret;

			for (size_t i = 0; i < ColNum(); i++)
				for (size_t j = 0; j < RowNum(); j++)
					ret._vals[i][j] = _vals[j][i];

			return ret;
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
			stream << "Rows: " << RowNum() << "  Cols: " << ColNum() << std::endl;

			for (size_t i = 0; i < RowNum(); i++)
			{
				stream << "[ ";
				for (size_t j = 0; j < ColNum(); j++)
				{
					stream << std::setw(width) << std::setprecision(precision) << _vals[i][j] << ", ";
				}
				stream << " ]" << std::endl;
			}
		}
		friend std::ostream& operator<<(std::ostream& stream, const MatrixNM& a)
		{
			a.Print(stream, 10, 3);

			return stream;
		}
	};

	///////////////////               Default MatrixNM typdefs                 ///////////////////
	typedef MatrixNM<float, 2, 2> Matrix22Flt;
	typedef MatrixNM<float, 3, 3> Matrix33Flt;
	typedef MatrixNM<float, 4, 4> Matrix44Flt;

	typedef MatrixNM<Real, 2, 2> Matrix22Dbl;
	typedef MatrixNM<Real, 3, 3> Matrix33Dbl;
	typedef MatrixNM<Real, 4, 4> Matrix44Dbl;

	typedef MatrixNM<Complex, 2, 2> Matrix22Complex;
	typedef MatrixNM<Complex, 3, 3> Matrix33Complex;
	typedef MatrixNM<Complex, 4, 4> Matrix44Complex;

	typedef MatrixNM<float, 2, 2> Mat22F;
	typedef MatrixNM<float, 3, 3> Mat33F;
	typedef MatrixNM<float, 4, 4> Mat44F;

	typedef MatrixNM<Real, 2, 2> Mat22D;
	typedef MatrixNM<Real, 3, 3> Mat33D;
	typedef MatrixNM<Real, 4, 4> Mat44D;

	typedef MatrixNM<Complex, 2, 2> Mat22C;
	typedef MatrixNM<Complex, 3, 3> Mat33C;
	typedef MatrixNM<Complex, 4, 4> Mat44C;
}

#endif