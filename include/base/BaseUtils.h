#if !defined MML_COREUTILS_H
#define MML_COREUTILS_H

#include "MMLBase.h"

#include "Vector.h"
#include "VectorN.h"
#include "Matrix.h"
#include "MatrixNM.h"

namespace MML
{
	namespace Utils
	{
		static int LeviCivita(int i, int j, int k)
		{
			if (i == j || j == k || i == k)
				return 0;

			if (i == 1 && j == 2 && k == 3)
				return 1;
			if (i == 2 && j == 3 && k == 1)
				return 1;
			if (i == 3 && j == 1 && k == 2)
				return 1;
			if (i == 3 && j == 2 && k == 1)
				return -1;
			if (i == 2 && j == 1 && k == 3)
				return -1;
			if (i == 1 && j == 3 && k == 2)
				return -1;

			return 0;
		}
		static int LeviCivita(int i, int j, int k, int l)
		{
			int a[4] = { i, j, k, l };
			int ret = 1;

			// TODO check if 1, 2, 3, 4 are present
			for (int i = 0; i < 4; i++)
				for (int j = i + 1; j < 4; j++)
					if (a[i] > a[j])
					{
						int tmp = a[i];
						a[i] = a[j];
						a[j] = tmp;
						ret *= -1;
					}
			return ret;
		}

		static Real DegToRad(Real angleDeg) { return angleDeg * Constants::PI / 180.0; }
		static Real RadToDeg(Real angleRad) { return angleRad * 180.0 / Constants::PI; }

		static void AngleDegToExplicit(Real angle, Real& deg, Real& min, Real& sec)
		{
			deg = floor(angle);
			min = floor((angle - deg) * 60.0);
			sec = (angle - deg - min / 60.0) * 3600.0;
		}
		static void AngleRadToExplicit(Real angleRad, Real& deg, Real& min, Real& sec) { AngleDegToExplicit(angleRad * 180.0 / Constants::PI, deg, min, sec); }
		static Real ExplicitToAngleDeg(Real deg, Real min, Real sec) { return deg + min / 60.0 + sec / 3600.0; }
		static Real ExplicitToAngleRad(Real deg, Real min, Real sec) { return ExplicitToAngleDeg(deg, min, sec) * Constants::PI / 180.0; }

		///////////////////                     Complex helpers                   ///////////////////
		static bool AreEqual(const Complex& a, const Complex& b, double eps = Defaults::ComplexEqualityPrecision)
		{
			if (std::abs(a.real() - b.real()) > eps || std::abs(a.imag() - b.imag()) > eps)
				return false;
			return true;
		}
		static bool AreEqualAbs(const Complex& a, const Complex& b, double eps = Defaults::ComplexAbsEqualityPrecision)
		{
			if (Abs(a - b) > eps)
				return false;
			return true;
		}
		static bool AreEqual(const Vector<Complex>& a, const Vector<Complex>& b, double eps = Defaults::ComplexEqualityPrecision)
		{
			if (a.size() != b.size())
				return false;

			for (int i = 0; i < a.size(); i++)
				if( !AreEqual(a[i], b[i], eps) )
					return false;

			return true;
		}
		static bool AreEqualAbs(const Vector<Complex>& a, const Vector<Complex>& b, double eps = Defaults::ComplexAbsEqualityPrecision)
		{
			if (a.size() != b.size())
				return false;

			for (int i = 0; i < a.size(); i++)
				if ( !AreEqualAbs(a[i], b[i], eps) )
					return false;

			return true;
		}

		///////////////////                     Vector helpers                     ///////////////////
		static Vector<Real> VectorProjectionParallelTo(const Vector<Real>& orig, const Vector<Real>& b)
		{
			return orig.ScalarProductCartesian(b) / b.NormL2() * b;
		}
		static Vector<Real> VectorProjectionPerpendicularTo(const Vector<Real>& orig, const Vector<Real>& b)
		{
			return orig - VectorProjectionParallelTo(orig, b);
		}

		static Real ScalarProduct(const Vector<Real>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("ScalarProduct - must be same dim", a.size(), b.size());

			Real ret = 0;
			for (int i = 0; i < a.size(); i++)
				ret += a[i] * b[i];

			return ret;
		}
		static Complex ScalarProduct(const Vector<Complex>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("ScalarProduct - must be same dim", a.size(), b.size());

			Complex ret = 0;
			for (int i = 0; i < a.size(); i++)
				ret += a[i] * std::conj(b[i]);

			return ret;
		}

		template<class Type> static Matrix<Type> OuterProduct(const Vector<Type>& a, const Vector<Type>& b)
		{
			Matrix<Type> ret(a.size(), b.size());

			for (int i = 0; i < a.size(); i++)
				for (int j = 0; j < b.size(); j++)
					ret[i][j] = a[i] * b[j];

			return ret;
		}

		///////////////////       Vector<Complex> - Vector<Real> operations       ///////////////////
		static Vector<Complex> AddVec(const Vector<Complex>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("AddVec(Complex, Real) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] + b[i];
			return ret;
		}
		static Vector<Complex> AddVec(const Vector<Real>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("AddVec(Real, Complex) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] + b[i];
			return ret;
		}

		static Vector<Complex> SubVec(const Vector<Complex>& a, const Vector<Real>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("SubVec(Complex, Real) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] - b[i];
			return ret;
		}
		static Vector<Complex> SubVec(const Vector<Real>& a, const Vector<Complex>& b)
		{
			if (a.size() != b.size())
				throw VectorDimensionError("SubVec(Real, Complex) - must be same dim", a.size(), b.size());

			Vector<Complex> ret(b.size());;
			for (int i = 0; i < b.size(); i++)
				ret[i] = a[i] - b[i];
			return ret;
		}
	};

	namespace MatrixUtils
	{
		const static inline MatrixNM<Complex, 2, 2> Pauli[] = {
				MatrixNM<Complex, 2, 2>{ 0, 1,
																 1, 0 },
				MatrixNM<Complex, 2, 2>{ 0,              Complex(0, -1),
																 Complex(0, 1),  0},
				MatrixNM<Complex, 2, 2>{ 1,  0,
																 0, -1 }
		};
		const static inline MatrixNM<Complex, 4, 4> DiracGamma[] = {
				MatrixNM<Complex, 4, 4>{ 1,  0,  0,  0,
																 0,  1,  0,  0,
																 0,  0, -1,  0,
																 0,  0,  0, -1 },
				MatrixNM<Complex, 4, 4>{ 0,  0,  0,  1,
																 0,  0,  1,  0,
																 0, -1,  0,  0,
																-1,  0,  0,  0 },
				MatrixNM<Complex, 4, 4>{              0,             0,             0, Complex(0, -1),
																              0,             0, Complex(0, 1),              0,
																              0, Complex(0, 1),             0,              0,
																 Complex(0, -1),             0,             0,              0 },
				MatrixNM<Complex, 4, 4>{ 0,  0,  1,  0,
																 0,  0,  0, -1,
																-1,  0,  0,  0,
																 0,  1,  0,  0 }
		};
		const static inline MatrixNM<Complex, 4, 4> DiracGamma5{ 0,  0,  1,  0,
																														 0,  0,  0,  1,
																														 1,  0,  0,  0,
																														 0,  1,  0,  0 };

		///////////////////             Creating Matrix from Vector              ///////////////////
		template<class Type> static Matrix<Type> RowMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret(1, (int)b.size());
			for (int i = 0; i < b.size(); i++)
				ret[0][i] = b[i];

			return ret;
		}
		template<class Type> static Matrix<Type> ColumnMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret((int)b.size(), 1);
			for (int i = 0; i < b.size(); i++)
				ret[i][0] = b[i];

			return ret;
		}
		template<class Type> static Matrix<Type> DiagonalMatrixFromVector(const Vector<Type>& b)
		{
			Matrix<Type> ret((int)b.size(), (int)b.size());
			for (int i = 0; i < b.size(); i++)
				ret[i][i] = b[i];

			return ret;
		}

		///////////////////                   Matrix helpers                     ///////////////////
		template<class Type> static Matrix<Type> Commutator(const Matrix<Type>& a, const Matrix<Type>& b)
		{
			return a * b - b * a;
		}
		template<class Type> static Matrix<Type> AntiCommutator(const Matrix<Type>& a, const Matrix<Type>& b)
		{
			return a * b + b * a;
		}

		template<class Type> static void MatrixDecomposeToSymAntisym(const Matrix<Type>& orig, Matrix<Type>& outSym, Matrix<Type>& outAntiSym)
		{
			if (orig.RowNum() != orig.ColNum())
				throw MatrixDimensionError("MatrixDecompose - matrix must be square", orig.RowNum(), orig.ColNum(), -1, -1);

			outSym = (orig + orig.GetTranspose()) * 0.5;
			outAntiSym = (orig - orig.GetTranspose()) * 0.5;
		}

		///////////////////                  Matrix functions                    ///////////////////
		template<class Type> static Matrix<Type> Exp(const Matrix<Type>& a, int n = 10)
		{
			Matrix<Type> ret(a.RowNum(), a.ColNum());
			Matrix<Type> a_pow_n(a);

			double fact = 1.0;
			ret.MakeUnitMatrix();

			for (int i = 1; i <= n; i++)
			{
				ret += a_pow_n / fact;

				a_pow_n = a_pow_n * a;
				fact *= (double)i;
			}

			return ret;
		}

		///////////////////                Real matrix helpers                   ///////////////////
		static bool IsOrthogonal(const Matrix<Real>& mat, double eps = Defaults::IsMatrixOrthogonalPrecision)
		{
			if (mat.RowNum() != mat.ColNum())
				throw MatrixDimensionError("IsOrthogonal - matrix must be square", mat.RowNum(), mat.ColNum(), -1, -1);

			Matrix<Real> matProd = mat * mat.GetTranspose();

			return matProd.IsUnit(eps);
		}

		///////////////////               Complex matrix helpers                 ///////////////////
		static Matrix<Real> GetRealPart(const Matrix<Complex>& a)
		{
			Matrix<Real> ret(a.RowNum(), a.ColNum());

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = a[i][j].real();
				}

			return	ret;
		}
		static Matrix<Real> GetImagPart(const Matrix<Complex>& a)
		{
			Matrix<Real> ret(a.RowNum(), a.ColNum());

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = a[i][j].imag();
				}

			return	ret;
		}

		static Matrix<Complex> GetConjugateTranspose(const Matrix<Complex>& mat)
		{
			Matrix<Complex> ret(mat.ColNum(), mat.RowNum());

			for (int i = 0; i < mat.RowNum(); i++)
				for (int j = 0; j < mat.ColNum(); j++)
					ret[j][i] = std::conj(mat[i][j]);

			return ret;
		}
		static Matrix<Complex> CmplxMatFromRealMat(const Matrix<Real>& mat)
		{
			Matrix<Complex> mat_cmplx(mat.RowNum(), mat.ColNum());

			for (int i = 0; i < mat.RowNum(); i++)
				for (int j = 0; j < mat.ColNum(); j++)
					mat_cmplx[i][j] = Complex(mat(i, j), 0.0);

			return mat_cmplx;
		}

		static bool IsComplexMatReal(const Matrix<Complex>& mat)
		{
			for (int i = 0; i < mat.RowNum(); i++)
				for (int j = 0; j < mat.ColNum(); j++)
					if (mat[i][j].imag() != 0.0)
						return false;

			return true;
		}
		static bool IsHermitian(const Matrix<Complex>& mat)
		{
			if (mat.RowNum() != mat.ColNum())
				throw MatrixDimensionError("IsHermitian - matrix must be square", mat.RowNum(), mat.ColNum(), -1, -1);

			for (int i = 0; i < mat.RowNum(); i++)
				for (int j = i + 1; j < mat.ColNum(); j++)
					if (mat[i][j] != std::conj(mat[j][i]))
						return false;
			return true;
		}
		static bool IsUnitary(const Matrix<Complex>& mat)
		{
			// IsUnitary - complex square matrix U is unitary if its conjugate transpose U* is also its inverse
			if (mat.RowNum() != mat.ColNum())
				throw MatrixDimensionError("IsUnitary - matrix must be square", mat.RowNum(), mat.ColNum(), -1, -1);

			Matrix<Complex> matProd = mat * GetConjugateTranspose(mat);

			return matProd.IsUnit();
		}

		///////////////////       Matrix<Complex> - Matrix<Real>  operations     ///////////////////
		static Matrix<Complex> AddMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
				throw MatrixDimensionError("AddMat(Complex, Real) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex> ret(a);
			for (int i = 0; i < b.RowNum(); i++)
				for (int j = 0; j < b.ColNum(); j++)
					ret[i][j] += b[i][j];
			return ret;
		}
		static Matrix<Complex> AddMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
				throw MatrixDimensionError("AddMat(Real, Complex) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex> ret(b);
			for (int i = 0; i < a.RowNum(); i++)
				for (int j = 0; j < a.ColNum(); j++)
					ret[i][j] += a[i][j];
			return ret;
		}

		static Matrix<Complex> SubMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
				throw MatrixDimensionError("AddMat(Complex, Real) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex> ret(a);
			for (int i = 0; i < b.RowNum(); i++)
				for (int j = 0; j < b.ColNum(); j++)
					ret[i][j] -= b[i][j];
			return ret;
		}
		static Matrix<Complex> SubMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.RowNum() != b.RowNum() || a.ColNum() != b.ColNum())
				throw MatrixDimensionError("AddMat(Real, Complex) - must be same dim", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex> ret(b);
			for (int i = 0; i < a.RowNum(); i++)
				for (int j = 0; j < a.ColNum(); j++)
					ret[i][j] = a[i][j] - b[i][j];
			return ret;
		}

		static Matrix<Complex> MulMat(const Complex& a, const Matrix<Real>& b)
		{
			Matrix<Complex>	ret(b.RowNum(), b.ColNum());

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = a * b[i][j];
				}

			return	ret;
		}
		static Matrix<Complex> MulMat(const Matrix<Real>& a, const Complex& b)
		{
			Matrix<Complex>	ret(a.RowNum(), a.ColNum());

			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = b * a[i][j];
				}

			return	ret;
		}
		
		static Matrix<Complex> MulMat(const Matrix<Complex>& a, const Matrix<Real>& b)
		{
			if (a.ColNum() != b.RowNum())
				throw MatrixDimensionError("Matrix::operator*(Complex, Real)) - a.colNum must be equal to b.rowNum", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex>	ret(a.RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = 0.0;
					for (int k = 0; k < a.ColNum(); k++)
						ret[i][j] += a[i][k] * b[k][j];
				}

			return	ret;
		}
		static Matrix<Complex> MulMat(const Matrix<Real>& a, const Matrix<Complex>& b)
		{
			if (a.ColNum() != b.RowNum())
				throw MatrixDimensionError("Matrix::operator*(Real, Complex)) - a.colNum must be equal to b.rowNum", a.RowNum(), a.ColNum(), b.RowNum(), b.ColNum());

			Matrix<Complex>	ret(a.RowNum(), b.ColNum());
			for (int i = 0; i < ret.RowNum(); i++)
				for (int j = 0; j < ret.ColNum(); j++) {
					ret[i][j] = 0.0;
					for (int k = 0; k < a.ColNum(); k++)
						ret[i][j] += a[i][k] * b[k][j];
				}

			return	ret;
		}

		static Vector<Complex> MulMatVec(const Matrix<Real>& a, const Vector<Complex>& b)
		{
			if (a.ColNum() != b.size())
				throw MatrixDimensionError("operator*(Mat a, Vec b) - a.colNum must be equal to vector size", a.RowNum(), a.ColNum(), (int)b.size(), -1);

			Vector<Complex>	ret(a.RowNum());
			for (int i = 0; i < a.RowNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.ColNum(); j++)
					ret[i] += a[i][j] * b[j];
			}
			return ret;
		}
		static Vector<Complex> MulMatVec(const Matrix<Complex>& a, const Vector<Real>& b)
		{
			if (a.ColNum() != b.size())
				throw MatrixDimensionError("operator*(Mat a, Vec b) - a.colNum must be equal to vector size", a.RowNum(), a.ColNum(), (int)b.size(), -1);

			Vector<Complex>	ret(a.RowNum());
			for (int i = 0; i < a.RowNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < a.ColNum(); j++)
					ret[i] += a[i][j] * b[j];
			}
			return ret;
		}

		static Vector<Complex> MulVecMat(const Vector<Complex>& a, const Matrix<Real>& b)
		{
			if (a.size() != b.RowNum())
				throw MatrixDimensionError("operator*(Vec a, Mat b) - vector size must be equal to b.rowNum", (int)a.size(), -1, b.RowNum(), b.ColNum());

			Vector<Complex>	ret(b.ColNum());
			for (int i = 0; i < b.ColNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.RowNum(); j++)
					ret[i] += a[j] * b[j][i];
			}
			return ret;
		}
		static Vector<Complex> MulVecMat(const Vector<Real>& a, const Matrix<Complex>& b)
		{
			if (a.size() != b.RowNum())
				throw MatrixDimensionError("operator*(Vec a, Mat b) - vector size must be equal to b.rowNum", (int)a.size(), -1, b.RowNum(), b.ColNum());

			Vector<Complex>	ret(b.ColNum());
			for (int i = 0; i < b.ColNum(); i++)
			{
				ret[i] = 0;
				for (int j = 0; j < b.RowNum(); j++)
					ret[i] += a[j] * b[j][i];
			}
			return ret;
		}
	}
}
#endif