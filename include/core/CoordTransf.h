#if !defined MML_COORD_TRANSF_H
#define MML_COORD_TRANSF_H

#include "MMLBase.h"

#include "interfaces/ICoordTransf.h"

#include "base/Vector.h"
#include "base/VectorN.h"
#include "base/Matrix.h"
#include "base/MatrixNM.h"
#include "base/Tensor.h"
#include "base/Geometry.h"

#include "core/Function.h"
#include "core/Derivation.h"

namespace MML
{
	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransf : public virtual ICoordTransf<VectorFrom, VectorTo, N>
	{
	public:
		virtual VectorTo getCovariantBasisVec(int ind, const VectorFrom& pos)
		{
			VectorTo ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->coordTransfFunc(i), ind, pos);

			return ret;
		}

		virtual VectorTo getUnitBasisVector(int ind, const VectorFrom& pos)
		{
			return getCovariantBasisVec(ind, pos).GetAsUnitVectorAtPos(pos);
		}

		MatrixNM<Real, N, N> jacobian(const VectorN<Real, N>& pos)
		{
			MatrixNM<Real, N, N> jac;

			for (int i = 0; i < N; ++i)
				for (int j = 0; j < N; ++j)
				{
					jac(i, j) = Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos);
				}
            
      return jac;
		}

		VectorTo   transfVecContravariant(const VectorFrom& vec, const VectorFrom& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++) {
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(i), j, pos) * vec[j];
			}
			return ret;
		}

		VectorFrom transfInverseVecCovariant(const VectorTo& vec, const VectorFrom& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->coordTransfFunc(j), i, pos, 1e-8) * vec[j];
			}
			return ret;
		}
	};

	template<typename VectorFrom, typename VectorTo, int N>
	class CoordTransfWithInverse : public virtual CoordTransf<VectorFrom, VectorTo, N>,
																 public virtual ICoordTransfWithInverse<VectorFrom, VectorTo, N>
	{
	public:
		virtual VectorFrom getContravariantBasisVec(int ind, const VectorTo& pos)
		{
			VectorFrom ret;

			for (int i = 0; i < N; i++)
				ret[i] = Derivation::NDer4Partial(this->inverseCoordTransfFunc(ind), i, pos);

			return ret;
		}

		virtual VectorFrom getUnitBasisVectorInverse(int ind, const VectorTo& pos)
		{
			return getContravariantBasisVec(ind, pos).GetAsUnitVectorAtPos(pos);
		}

		VectorTo transfVecCovariant(const VectorFrom& vec, const VectorTo& pos)
		{
			VectorTo ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->inverseCoordTransfFunc(j), i, pos) * vec[j];
			}

			return ret;
		}

		VectorFrom transfInverseVecContravariant(const VectorTo& vec, const VectorTo& pos)
		{
			VectorFrom ret;
			for (int i = 0; i < N; i++)
			{
				ret[i] = 0;
				for (int j = 0; j < N; j++)
					ret[i] += Derivation::NDer4Partial(this->inverseCoordTransfFunc(i), j, pos, 1e-8) * vec[j];
			}

			return ret;
		}

		Tensor2<N> transfTensor2(const Tensor2<N>& tensor, const VectorFrom& pos)
		{
			Tensor2<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
				{
					ret.Component(i, j) = 0;
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
						{
							double coef1, coef2;
							if (tensor._isContravar[0])
								coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), k, pos);
							else
								coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(k), i, pos);

							if (tensor._isContravar[1])
								coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), l, pos);
							else
								coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(l), j, pos);

							ret.Component(i, j) += coef1 * coef2 * tensor.Component(k, l);
						}
				}

			return ret;
		}

		Tensor3<N> transfTensor3(const Tensor3<N>& tensor, const VectorFrom& pos)
		{
			Tensor3<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
					{
						ret.Component(i, j, k) = 0;
						for (int l = 0; l < N; l++)
							for (int m = 0; m < N; m++)
								for (int n = 0; n < N; n++)
								{
									double coef1, coef2, coef3;
									if (tensor._isContravar[0])
										coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), l, pos);
									else
										coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(l), i, pos);

									if (tensor._isContravar[1])
										coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), m, pos);
									else
										coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(m), j, pos);

									if (tensor._isContravar[2])
										coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), n, pos);
									else
										coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), k, pos);

									ret.Component(i, j, k) += coef1 * coef2 * coef3 * tensor.Component(l, m, n);
								}
					}

			return ret;
		}

		Tensor4<N> transfTensor4(const Tensor4<N>& tensor, const VectorFrom& pos)
		{
			Tensor4<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
						{
							ret[i][j][k][l] = 0;
							for (int m = 0; m < N; m++)
								for (int n = 0; n < N; n++)
									for (int o = 0; o < N; o++)
										for (int p = 0; p < N; p++)
										{
											double coef1, coef2, coef3, coef4;
											if (tensor._isContravar[0])
												coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), m, pos);
											else
												coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(m), i, pos);

											if (tensor._isContravar[1])
												coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), n, pos);
											else
												coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), j, pos);

											if (tensor._isContravar[2])
												coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), o, pos);
											else
												coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(o), k, pos);

											if (tensor._isContravar[3])
												coef4 = Derivation::NDer1Partial(this->coordTransfFunc(l), p, pos);
											else
												coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(p), l, pos);

											ret[i][j][k][l] += coef1 * coef2 * coef3 * coef4 * tensor[m][n][o][p];
										}
						}

			return ret;
		}

		Tensor5<N> transfTensor5(const Tensor5<N>& tensor, const VectorFrom& pos)
		{
			Tensor5<N> ret(tensor.NumContravar(), tensor.NumCovar());

			for (int i = 0; i < N; i++)
				for (int j = 0; j < N; j++)
					for (int k = 0; k < N; k++)
						for (int l = 0; l < N; l++)
							for (int m = 0; m < N; m++)
							{
								ret[i][j][k][l][m] = 0;
								for (int n = 0; n < N; n++)
									for (int o = 0; o < N; o++)
										for (int p = 0; p < N; p++)
											for (int q = 0; q < N; q++)
												for (int r = 0; r < N; r++)
												{
													double coef1, coef2, coef3, coef4, coef5;
													if (tensor._isContravar[0])
														coef1 = Derivation::NDer1Partial(this->coordTransfFunc(i), n, pos);
													else
														coef1 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(n), i, pos);

													if (tensor._isContravar[1])
														coef2 = Derivation::NDer1Partial(this->coordTransfFunc(j), o, pos);
													else
														coef2 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(o), j, pos);

													if (tensor._isContravar[2])
														coef3 = Derivation::NDer1Partial(this->coordTransfFunc(k), p, pos);
													else
														coef3 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(p), k, pos);

													if (tensor._isContravar[3])
														coef4 = Derivation::NDer1Partial(this->coordTransfFunc(l), q, pos);
													else
														coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(q), l, pos);

													if (tensor._isContravar[4])
														coef4 = Derivation::NDer1Partial(this->coordTransfFunc(m), r, pos);
													else
														coef4 = Derivation::NDer1Partial(this->inverseCoordTransfFunc(r), m, pos);

													ret[i][j][k][l][m] += coef1 * coef2 * coef3 * coef4 * coef5 * tensor[n][o][p][q][r];
												}
							}

			return ret;
		}
	};

	class CoordTransfPolarToCartesian2D : public CoordTransfWithInverse<Vector2Polar, Vector2Cartesian, 2>
	{
		// q[0] = r     - radial distance
		// q[1] = phi   - polar angle
	public:
		static Real func1(const VectorN<Real, 2>& q) { return q[0] * cos(q[1]); }
		static Real func2(const VectorN<Real, 2>& q) { return q[0] * sin(q[1]); }

		// q[0] = x
		// q[1] = y
		static Real funcInverse1(const VectorN<Real, 2>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real funcInverse2(const VectorN<Real, 2>& q) { return atan2(q[1], q[0]); }

		inline static ScalarFunction<2> _func[2] = { ScalarFunction<2>{func1},
																								 ScalarFunction<2>{func2}
		};

		inline static ScalarFunction<2> _funcInverse[2] = { ScalarFunction<2>{funcInverse1},
																												ScalarFunction<2>{funcInverse2}
		};

		Vector2Cartesian     transf(const Vector2Polar& q) const { return Vector2Cartesian{ func1(q), func2(q) }; }
		Vector2Polar         transfInverse(const Vector2Cartesian& q) const { return Vector2Polar{ funcInverse1(q), funcInverse2(q) }; }

		const IScalarFunction<2>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	class CoordTransfCart2DRotation : public CoordTransfWithInverse<Vector2Cartesian, Vector2Cartesian, 2>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 2, 2>  _transf;
		MatrixNM<Real, 2, 2>  _inverse;

		const ScalarFunctionFromStdFunc<2> _f1;
		const ScalarFunctionFromStdFunc<2> _f2;

		const ScalarFunctionFromStdFunc<2> _fInverse1;
		const ScalarFunctionFromStdFunc<2> _fInverse2;

	public:

		CoordTransfCart2DRotation(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::func2, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 2>&)> { std::bind(&CoordTransfCart2DRotation::funcInverse2, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][1] = -sin(_angle);
			_transf[1][0] = sin(_angle);
			_transf[1][1] = cos(_angle);

			_inverse[0][0] = cos(_angle);
			_inverse[0][1] = sin(_angle);
			_inverse[1][0] = -sin(_angle);
			_inverse[1][1] = cos(_angle);
		}

		Real func1(const VectorN<Real, 2>& q) const { return _transf[0][0] * q[0] + _transf[0][1] * q[1]; }
		Real func2(const VectorN<Real, 2>& q) const { return (_transf * q)[1]; }

		Real funcInverse1(const VectorN<Real, 2>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 2>& q) const { return (_inverse * q)[1]; }

		Vector2Cartesian    transf(const Vector2Cartesian& q) const { return Vector2Cartesian{ func1(q), func2(q) }; }
		Vector2Cartesian    transfInverse(const Vector2Cartesian& q) const { return Vector2Cartesian{ funcInverse1(q), funcInverse2(q) }; }

		const IScalarFunction<2>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else return _f2;
		}
		const IScalarFunction<2>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else return _fInverse2;
		}
	};

	class CoordTransfCart3DRotationXAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationXAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationXAxis::funcInverse3, this, std::placeholders::_1) })
		{
			_transf[0][0] = 1.0;
			_transf[1][1] = cos(_angle);
			_transf[1][2] = -sin(_angle);
			_transf[2][1] = sin(_angle);
			_transf[2][2] = cos(_angle);

			_inverse[0][0] = 1.0;
			_inverse[1][1] = cos(_angle);
			_inverse[1][2] = sin(_angle);
			_inverse[2][1] = -sin(_angle);
			_inverse[2][2] = cos(_angle);
		}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i)const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};
	class CoordTransfCart3DRotationYAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationYAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationYAxis::funcInverse3, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][2] = sin(_angle);
			_transf[1][1] = 1.0;
			_transf[2][0] = -sin(_angle);
			_transf[2][2] = cos(_angle);

			_inverse[0][0] = cos(_angle);
			_inverse[0][2] = -sin(_angle);
			_inverse[1][1] = 1.0;
			_inverse[2][0] = sin(_angle);
			_inverse[2][2] = cos(_angle);
		}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};
	class CoordTransfCart3DRotationZAxis : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Real    _angle;
		MatrixNM<Real, 3, 3>  _transf;
		MatrixNM<Real, 3, 3>  _inverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

	public:
		CoordTransfCart3DRotationZAxis(Real inAngle) : _angle(inAngle),
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransfCart3DRotationZAxis::funcInverse3, this, std::placeholders::_1) })
		{
			_transf[0][0] = cos(_angle);
			_transf[0][1] = -sin(_angle);
			_transf[1][0] = sin(_angle);
			_transf[1][1] = cos(_angle);
			_transf[2][2] = 1.0;

			_inverse[0][0] = cos(_angle);
			_inverse[0][1] = sin(_angle);
			_inverse[1][0] = -sin(_angle);
			_inverse[1][1] = cos(_angle);
			_inverse[2][2] = 1.0;
		}

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_inverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_inverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_inverse * q)[2]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}
	};

	// TODO 0.9 - LOW, TESKO, Euler angles, finish CoordTransfCart3DRotationGeneralAxis
	// class CoordTransfCart3DRotationGeneralAxis  : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	// {

	// };

	// Performs tranformation from original (Cartesian) system to orthogonal system defined by 
	// its base vectors expressed in original system.
	class CoordTransf3DCartOrthogonal : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _transfInverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[2]; }

	public:
		CoordTransf3DCartOrthogonal(const VectorN<Real, 3> &b1, const VectorN<Real, 3> &b2, const VectorN<Real, 3> &b3) : 
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOrthogonal::funcInverse3, this, std::placeholders::_1) })
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			for (int i = 0; i < 3; i++)	{
				for (int j = 0; j < 3; j++)
				{
					_transf(i, j)  = _base[i][j];
					_transfInverse(i, j) = _base[j][i]; 
				}
			}
		}

		MatrixNM<Real, 3, 3> getTransfMatrix() { return _transf; }
		MatrixNM<Real, 3, 3> getInvTransfMatrix() { return _transfInverse; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProd(_base[0], _base[1]);
			if (ScalarProd(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};

	// General 3D Cartesian transformation, given by matrix
class CoordTransf3DCartGeneral : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];

		MatrixNM<Real, 3, 3> _transf;
		MatrixNM<Real, 3, 3> _transfInverse;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		Real func1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real func2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real func3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transfInverse * q)[2]; }

	public:
		CoordTransf3DCartGeneral(const MatrixNM<Real, 3, 3> &transfMat) : 
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartGeneral::funcInverse3, this, std::placeholders::_1) })
		{
			_transf = transfMat;

			for (int i = 0; i < 3; i++)	{
				for (int j = 0; j < 3; j++)
				{
					_base[i][j]  = _transf(i, j); 
					_transfInverse(i, j) = _transf[j][i];
				}
			}
		}

		MatrixNM<Real, 3, 3> getTransfMatrix() { return _transf; }
		MatrixNM<Real, 3, 3> getInvTransfMatrix() { return _transfInverse; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProd(_base[0], _base[1]);
			if (ScalarProd(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};
	class CoordTransf3DCartOblique : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cartesian, 3>
	{
	private:
		Vector3Cartesian _base[3];
		Vector3Cartesian _dual[3];

		MatrixNM<Real, 3, 3> _alpha;
		MatrixNM<Real, 3, 3> _transf;

		const ScalarFunctionFromStdFunc<3> _f1, _f2, _f3;
		const ScalarFunctionFromStdFunc<3> _fInverse1, _fInverse2, _fInverse3;

		// TODO - ovo popraviti i napraviti kako spada
		Real func1(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[0])); }
		Real func2(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[1])); }
		Real func3(const VectorN<Real, 3>& q) const { return ScalarProd(q, MML::Vector3Cartesian(_dual[2])); }

		Real funcInverse1(const VectorN<Real, 3>& q) const { return (_transf * q)[0]; }
		Real funcInverse2(const VectorN<Real, 3>& q) const { return (_transf * q)[1]; }
		Real funcInverse3(const VectorN<Real, 3>& q) const { return (_transf * q)[2]; }

	public:
		CoordTransf3DCartOblique(VectorN<Real, 3> b1, VectorN<Real, 3> b2, VectorN<Real, 3> b3) : 
			_f1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func1, this, std::placeholders::_1) }),
			_f2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func2, this, std::placeholders::_1) }),
			_f3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::func3, this, std::placeholders::_1) }),
			_fInverse1(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse1, this, std::placeholders::_1) }),
			_fInverse2(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse2, this, std::placeholders::_1) }),
			_fInverse3(std::function<Real(const VectorN<Real, 3>&)> { std::bind(&CoordTransf3DCartOblique::funcInverse3, this, std::placeholders::_1) })
		{
			_base[0] = b1;
			_base[1] = b2;
			_base[2] = b3;

			Real V = ScalarProd(_base[0], VectorProd(_base[1], _base[2]));

			_dual[0] = (1 / V) * VectorProd(_base[1], _base[2]);
			_dual[1] = (1 / V) * VectorProd(_base[2], _base[0]);
			_dual[2] = (1 / V) * VectorProd(_base[0], _base[1]);

			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					_alpha(i, j) = _base[i][j];
					_transf(i, j) = _base[j][i];     // transponirano
				}
			}
		}

		MatrixNM<Real, 3, 3> getAlpha() { return _alpha; }
		MatrixNM<Real, 3, 3> getTransf() { return _alpha; }

		Vector3Cartesian    Base(int i) { return _base[i]; }
		Vector3Cartesian    Dual(int i) { return _dual[i]; }

		Vector3Cartesian    transf(const Vector3Cartesian& q) const { return Vector3Cartesian{ func1(q), func2(q), func3(q) }; }
		Vector3Cartesian    transfInverse(const Vector3Cartesian& q) const { return Vector3Cartesian{ funcInverse1(q), funcInverse2(q), funcInverse3(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const
		{
			if (i == 0) return _f1;
			else if (i == 1) return _f2;
			else return _f3;
		}
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const
		{
			if (i == 0) return _fInverse1;
			else if (i == 1) return _fInverse2;
			else return _fInverse3;
		}

		bool IsRightHanded()
		{
			Vector3Cartesian cross = VectorProd(_base[0], _base[1]);
			if (ScalarProd(cross, _base[2]) > 0.0)
				return true;
			else
				return false;
		}
	};

	// TODO 0.9 - VIDJETI STO S MATH vs PHY konvencijama o redoslijedu koordinata
	class CoordTransfSphericalToCartesian : public CoordTransfWithInverse<Vector3Spherical, Vector3Cartesian, 3>
	{
	private:
		// q[0] = r     - radial distance
		// q[1] = theta - inclination
		// q[2] = phi   - azimuthal angle
		static Real x(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * cos(q[2]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * sin(q[2]); }
		static Real z(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }

		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]); }
		static Real theta(const VectorN<Real, 3>& q) { return acos(q[2] / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{x},
																								 ScalarFunction<3>{y},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{r},
																												ScalarFunction<3>{theta},
																												ScalarFunction<3>{phi}
		};
	public:
		Vector3Cartesian     transf(const Vector3Spherical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }
		Vector3Spherical     transfInverse(const Vector3Cartesian& q) const { return Vector3Spherical{ r(q), theta(q), phi(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }

		// TODO 0.9 - overload covar i contravar transf. funkcije s analiticki izracunatim jakobijanom derivacija
	};

	class CoordTransfCartesianToSpherical : public CoordTransfWithInverse<Vector3Cartesian, Vector3Spherical, 3>
	{
	private:
		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2]); }
		static Real theta(const VectorN<Real, 3>& q) { return acos(q[2] / sqrt(q[0] * q[0] + q[1] * q[1] + q[2] * q[2])); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }

		// q[0] = r     - radial distance
		// q[1] = theta - inclination
		// q[2] = phi   - azimuthal angle
		static Real x(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * cos(q[2]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]) * sin(q[2]); }
		static Real z(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{r},
																								 ScalarFunction<3>{theta},
																								 ScalarFunction<3>{phi}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{x},
																												ScalarFunction<3>{y},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Spherical     transf(const Vector3Cartesian& q) const { return Vector3Spherical{ r(q), theta(q), phi(q) }; }
		Vector3Cartesian     transfInverse(const Vector3Spherical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	class CoordTransfCylindricalToCartesian : public CoordTransfWithInverse<Vector3Cylindrical, Vector3Cartesian, 3>
	{
	private:
		// q1 = r   - distance from symmetry axis
		// q2 = phi - angle to symmetry axis
		// q3 = z   - z
		static Real x(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]); }
		static Real z(const VectorN<Real, 3>& q) { return q[2]; }

		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }
		//static Real funcInverse3(const VectorN<Real, 3> &q) { return q[2]; }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{x},
																								 ScalarFunction<3>{y},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{r},
																												ScalarFunction<3>{phi},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Cartesian     transf(const Vector3Cylindrical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }
		Vector3Cylindrical   transfInverse(const Vector3Cartesian& q) const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	class CoordTransfCartesianToCylindrical : public CoordTransfWithInverse<Vector3Cartesian, Vector3Cylindrical, 3>
	{
	private:
		// q[0] = x
		// q[1] = y
		// q[2] = z
		static Real r(const VectorN<Real, 3>& q) { return sqrt(q[0] * q[0] + q[1] * q[1]); }
		static Real phi(const VectorN<Real, 3>& q) { return atan2(q[1], q[0]); }
		static Real z(const VectorN<Real, 3>& q) { return q[2]; }

		// q1 = r   - distance from symmetry axis
		// q2 = phi - angle to symmetry axis
		// q3 = z   - z
		static Real x(const VectorN<Real, 3>& q) { return q[0] * cos(q[1]); }
		static Real y(const VectorN<Real, 3>& q) { return q[0] * sin(q[1]); }

		inline static ScalarFunction<3> _func[3] = { ScalarFunction<3>{r},
																								 ScalarFunction<3>{phi},
																								 ScalarFunction<3>{z}
		};

		inline static ScalarFunction<3> _funcInverse[3] = { ScalarFunction<3>{x},
																												ScalarFunction<3>{y},
																												ScalarFunction<3>{z}
		};
	public:
		Vector3Cylindrical transf(const Vector3Cartesian& q) const { return Vector3Cylindrical{ r(q), phi(q), z(q) }; }
		Vector3Cartesian   transfInverse(const Vector3Cylindrical& q) const { return Vector3Cartesian{ x(q), y(q), z(q) }; }

		const IScalarFunction<3>& coordTransfFunc(int i) const { return _func[i]; }
		const IScalarFunction<3>& inverseCoordTransfFunc(int i) const { return _funcInverse[i]; }
	};

	static CoordTransfSphericalToCartesian      CoordTransfSpherToCart;
	static CoordTransfCylindricalToCartesian    CoordTransfCylToCart;
	static CoordTransfCartesianToSpherical      CoordTransfCartToSpher;
	static CoordTransfCartesianToCylindrical    CoordTransfCartToCyl;
}
#endif
