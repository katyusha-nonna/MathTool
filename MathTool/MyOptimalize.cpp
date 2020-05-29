/*! \file MyOptimalize.cpp */

#include "MyOptimalize.h"
//矩阵运算包
#include <Eigen/Eigen>
#define DEIGEN_MPL2_ONLY

using namespace Utility;

constexpr double zeroAct = 1e-10;

double LPBarrierInteriorPointKernel(const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, double t, double mu, double error);
double LPPrimalDualInteriorPointKernel(const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, double alpha, double sigma, double error);
double QPConjugateGradientKernel(const Eigen::MatrixXd& G, const Eigen::VectorXd& p, Eigen::VectorXd& x);
double QPASMKernel(const Eigen::MatrixXd& G, const Eigen::MatrixXd& A, const Eigen::VectorXd& p, const Eigen::VectorXd& b, Eigen::VectorXd& x);
double QPLagrangeKernel(const Eigen::MatrixXd& G, const Eigen::MatrixXd& A, const Eigen::VectorXd& p, const Eigen::VectorXd& b, Eigen::VectorXd& x);

/*!
	障碍内点法的核心函数
	\li 用于求解不等式约束线性优化问题

	\f[
	\begin{aligned}
	\min_{x,y}\quad&c^{T}x\\
	s.t.\quad&A^{T}x\geq b
	\end{aligned}
	\f]

	\param c 目标函数系数向量(Eigen向量)
	\param A 约束条件系数矩阵(Eigen矩阵)
	\param b 约束条件常数向量(Eigen向量)
	\param x 决策变量优化结果向量(Eigen向量)

	\return 目标函数优化结果

	\todo 暂无
*/
double LPBarrierInteriorPointKernel(
	const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, double t, double mu, double error)
{
	// Min: c'*x
	// s.t A'*x >= b
	assert(A.cols() == b.rows());
	assert(A.rows() == x.rows());
	// 两个主要的维度(变量数目/约束数目)
	const auto n = A.rows();
	const auto m = A.cols();
	// 临时变量-梯度向量
	Eigen::VectorXd gradV;
	// 临时变量-海森矩阵
	Eigen::MatrixXd hessianM;
	// 临时变量-决策变量增量
	Eigen::VectorXd delta;
	// 当前最大误差
	double curGlobalError = 1.0;
	// lambda对象，用途：初始化
	auto init = [&]() {
		// 这里不对x进行检查，直接认为传入函数的x长度等于n且具有合适的初值
	};
	// lambda对象，用途：牛顿迭代求解非线性方程组
	auto newton = [&]() {
		int maxCount = 8, curCount = 0;
		double maxError = 1e-4, curError = 1.0;
		while (curError > maxError && curCount < maxCount)
		{
			// 更新梯度向量
			gradV = c * t - A * (A.transpose()*x - b).cwiseInverse();
			// 更新海森矩阵
			hessianM = A * (A.transpose()*x - b).cwiseAbs2().cwiseInverse().asDiagonal()*A.transpose();
			// 更新当前解
			delta = hessianM.colPivHouseholderQr().solve(gradV);
			x -= delta;
			// 更新计数器与最大误差
			curError = delta.norm() / x.norm();
			curCount++;
		}
	};

	init();
	// 开始迭代
	while (curGlobalError > error)
	{
		newton();
		curGlobalError = m / t;
		t *= mu;
	}
	// 返回最优解
	auto f = c.transpose()*x;
	return f.value();
}

/*!
	原始对偶内点法的核心函数
	\li 用于求解等式约束线性优化问题

	\f[
	\begin{aligned}
	\min_{x,y}\quad&c^{T}x\\
	s.t.\quad&A^{T}x=b
	\end{aligned}
	\f]

	\param c 目标函数系数向量(Eigen向量)
	\param A 约束条件系数矩阵(Eigen矩阵)
	\param b 约束条件常数向量(Eigen向量)
	\param x 决策变量优化结果向量(Eigen向量)
	\param alpha 牛顿法步长(0<alpha<1，若alpha<=0将自动选择步长)
	\param sigma 牛顿法下降因子(0<=sigma<1，若sigma<0将自动选择下降因子)
	\param error 牛顿法最大迭代误差(原问题-对偶问题的对偶间隙)

	\return 目标函数优化结果(当无解时返回INFITY)

	\todo 暂无
*/
double LPPrimalDualInteriorPointKernel(
	const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, double alpha, double sigma, double error)
{
	// Min: c'*x
	// s.t A'*x = b
	assert(A.cols() == b.rows());
	assert(A.rows() == x.rows());
	// 两个主要的维度(变量数目/约束数目)
	const auto n = A.rows();
	const auto m = A.cols();
	// 一些标志位
	const bool autoAlpha = alpha <= 0;
	const bool autoSigma = sigma < 0;
	bool isUnbounded = false;
	// 临时变量-决策变量 var=[x; lambda; s]
	Eigen::VectorXd var;
	// 临时变量-雅克比矩阵(稀疏)
	Eigen::SparseMatrix<double, Eigen::ColMajor> jacobiM(2 * n + m, 2 * n + m);
	// 临时变量-KKT方程组右端向量
	Eigen::VectorXd right;
	// 临时变量迭代增量求解结果(负值) delta=[xDelta; lambdaDelta; sDelta]
	Eigen::VectorXd delta;
	// 当前的对偶间隙
	double mu;
	// lambda对象，用途：初始化
	auto init = [&]() {
		// 初始化决策变量var
		var.setOnes(2 * n + m);
		var.head(n) = x;
		// 初始化雅可比矩阵
		auto xx = Eigen::VectorXd::Map(var.data(), n);
		auto s = Eigen::VectorXd::Map(var.data() + n + m, n);
		std::list<Eigen::Triplet<double> > tripletList;
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				tripletList.emplace_back(i, n + j, A.coeff(i, j));
		for (int i = 0; i < n; i++)
			for (int j = 0; j < m; j++)
				tripletList.emplace_back(n + j, i, A.coeff(i, j));
		for (int k = 0; k < n; k++)
			tripletList.emplace_back(n + m + k, k, s.coeff(k));
		for (int k = 0; k < n; k++)
			tripletList.emplace_back(n + m + k, n + m + k, xx.coeff(k));
		for (int k = 0; k < n; k++)
			tripletList.emplace_back(k, n + m + k, 1.0);
		tripletList.sort([](const auto& p1, const auto&p2) {return p1.col() == p2.col() ? p1.row() < p2.row() : p1.col() < p2.col(); });
		jacobiM.setFromTriplets(tripletList.begin(), tripletList.end());
		jacobiM.makeCompressed();
		// 初始化右端向量
		right.setZero(2 * n + m);
		// 初始化增量向量
		delta.setZero(2 * n + m);
		// 初始化对偶间隙
		mu = xx.dot(s) / n;
	};
	// lambda对象，用途：更新决策变量
	auto updataVar = [&]() {
		// 判断alpha决定更新步长
		if (autoAlpha)
		{
			// 自动决定步长更新决策变量
			auto xDelta = Eigen::VectorXd::Map(delta.data(), n);
			auto sDelta = Eigen::VectorXd::Map(delta.data() + n + m, n);
			auto xx = Eigen::VectorXd::Map(var.data(), n);
			auto s = Eigen::VectorXd::Map(var.data() + n + m, n);
			Eigen::VectorXd tempV;
			tempV.resize(2 * n);
			tempV.head(n) = -xDelta.cwiseQuotient(xx);
			tempV.tail(n) = -sDelta.cwiseQuotient(s);
			alpha = std::min(0.9999 / tempV.maxCoeff(), 1.0);
		}
		var -= alpha * delta;
	};
	// lambda对象，用途：计算下降因子
	auto updateSigma = [&]() {
		if (autoSigma)
		{
			// 暂时固定为1.0
			sigma = 1.0;
		}
	};
	// lambda对象，用途：更新对偶间隙
	auto checkDualMeasure = [&]() {
		// mu=x^{T}s/n
		auto xx = Eigen::VectorXd::Map(var.data(), n);
		auto s = Eigen::VectorXd::Map(var.data() + n + m, n);
		mu = xx.dot(s) / n;
	};
	// lambda对象，用途：判断是否无界
	auto checkUnbounded = [&]() {
		auto xNorm = Eigen::VectorXd::Map(var.data(), n).norm();
		auto lNorm = Eigen::VectorXd::Map(var.data() + n, m).norm();
		// 判断x和lambda的二范数是否是NAN或者INFITY(判断无界)
		return isinf(xNorm) || isinf(lNorm) || isnan(xNorm) || isnan(lNorm);
	};
	// lambda对象，用途：牛顿法迭代KKT条件方程求解原始-对偶问题
	auto newton = [&]() {
		// 雅克比矩阵为一个稀疏矩阵，每次迭代时更新左下角和右下角区块
		// [      0        A        I
		//      A^{T}   0        0
		//        S        0        X    ]
		int maxCount = 100, curCount = 0;
		double maxError = error, curError = 1.0;
		auto xx = Eigen::VectorXd::Map(var.data(), n);
		auto lambda = Eigen::VectorXd::Map(var.data() + n, m);
		auto s = Eigen::VectorXd::Map(var.data() + n + m, n);
		auto rx = Eigen::VectorXd::Map(right.data(), n);
		auto rl = Eigen::VectorXd::Map(right.data() + n, m);
		auto rs = Eigen::VectorXd::Map(right.data() + n + m, n);
		double tempProd = mu * sigma;
		Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
		while (curError > maxError && curCount < maxCount)
		{
			// 更新下降因子
			updateSigma();
			// 更新雅克比矩阵
			for (int k = 0; k < n; k++)
				jacobiM.coeffRef(n + m + k, k) = s.coeff(k);
			for (int k = 0; k < n; k++)
				jacobiM.coeffRef(n + m + k, n + m + k) = xx.coeff(k);
			// 更新右端向量
			rx = A * lambda + s - c;
			rl = A.transpose()*xx - b;
			tempProd = mu * sigma;
			rs = xx.cwiseProduct(s);
			rs.array() -= tempProd;
			// 求解稀疏矩阵
			solver.compute(jacobiM);
			delta = solver.solve(right);
			// 更新决策变量
			updataVar();
			// 更新对偶间隙
			checkDualMeasure();
			// 判断是否无解
			if (checkUnbounded())
			{
				isUnbounded = true;
				break;
			}
			// 更新计数器
			curError = delta.norm();
			curCount++;
		}
	};

	// 求解
	init();
	newton();
	// 返回最优解
	if (isUnbounded)
		return INFINITY;
	x = var.head(n);
	auto f = c.transpose()*x;
	return f.value();
}

/*!
	共轭梯度法的核心函数
	\li 只能用于求解G矩阵正定对称的无约束二次优化问题

	\f[
	\begin{aligned}
	\min_{x,y}\quad&\frac{1}{2}x^{T}Gx+p^{T}x\\
	s.t.\quad&\quad
	\end{aligned}
	\f]

	\param G 目标函数二次项及交叉乘项系数矩阵(Eigen矩阵)(对称正定)
	\param p 目标函数一次项系数向量(Eigen向量)
	\param x 决策变量优化结果向量(Eigen向量)

	\return 目标函数优化结果

	\todo 暂无
*/
double QPConjugateGradientKernel(const Eigen::MatrixXd& G, const Eigen::VectorXd& p, Eigen::VectorXd& x)
{
	// Min: (1 / 2)*x'*G*x + p'*x
	// s.t NULL
	// 共轭梯度求解器(稀疏)
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
	// 方程求解
	x = solver.compute(G.sparseView()).solve(-p);
	// 返回最优解
	auto f = p.transpose()*x + 0.5*x.transpose()*G*x;
	return f.value();
}

/*!
	有效约束集法的核心函数
	\li 用于求解不等式约束二次优化问题

	\f[
	\begin{aligned}
	\min_{x,y}\quad&\frac{1}{2}x^{T}Gx+p^{T}x\\
	s.t.\quad&A^{T}x\geq b
	\end{aligned}
	\f]

	\param G 目标函数二次项及交叉乘项系数矩阵(Eigen矩阵)
	\param A 约束条件系数矩阵(Eigen矩阵)
	\param p 目标函数一次项系数矩阵(Eigen向量)
	\param b 约束条件常数向量(Eigen向量)
	\param x 决策变量优化结果向量(Eigen向量)

	\return 目标函数优化结果

	\todo 暂无
*/
double QPASMKernel(const Eigen::MatrixXd& G, const Eigen::MatrixXd& A, const Eigen::VectorXd& p, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{
	// Min: (1 / 2)*x'*G*x + p'*x
	// s.t A'*x >= b
	assert(G.rows() == G.cols());
	assert(G.rows() == A.rows());
	assert(G.rows() == p.rows());
	assert(A.cols() == b.rows());
	assert(x.rows() == G.rows());
	// 两个主要维度
	const auto n = G.rows();
	const auto m = A.cols();
	// 当前的约束集
	std::set<int> S;
	// 当前子问题的解
	Eigen::VectorXd tempx;
	Eigen::VectorXd delta, lambda;
	// 当前更新步长
	double alpha;
	std::map<int, double> candidate;
	// 当前的p，A和b
	Eigen::MatrixXd tempA;
	Eigen::VectorXd tempp, tempb;
	// 退出标志
	bool quit = false;
	// 迭代计数器
	int iterCount = 0;
	// lambda对象，用途：初始化约束集
	auto init = [&]() {
		S.clear();
		Eigen::VectorXd temp = A.transpose()*x - b;
		for (int i = 0; i < m; i++)
		{
			if (temp[i] >= 0 && std::abs(temp[i]) < zeroAct)
				S.emplace(i);
		}
	};
	// lambda对象，用途：根据当前约束集更新当前的p，A和b
	auto updateSubP = [&]() {
		int count = 0;
		tempA.resize(n, S.size());
		tempb.resize(S.size());
		for (auto& i : S)
		{
			tempA.col(count) = A.col(i);
			tempb[count] = 0.0;
			count++;
		}
		tempp = G * x + p;
	};
	// lambda对象，用途：判断当前x+delta是否为原问题可行解
	auto isFeasible = [&]() {
		Eigen::VectorXd temp = A.transpose()*(x + delta) - b;
		double minVal = temp.minCoeff();
		return minVal >= 0;
	};
	// lambda对象，用途：更新alpha和S
	auto updateAlpha = [&]() {
		candidate.clear();
		if (isFeasible())
		{
			alpha = 1.0;
		}
		else
		{
			for (int i = 0; i < m; i++)
			{
				bool okay = S.find(i) == S.end();
				double temp = A.col(i).transpose()*delta;
				okay = okay && temp < 0.0;
				if (okay)
					candidate.emplace(i, (double)((b[i] - A.col(i).transpose()*x) / temp));
			}
			int result = std::min_element(candidate.begin(), candidate.end(), [](auto& p1, auto& p2) {return p1.second <= p2.second; })->first;
			alpha = candidate.at(result);
			S.emplace(result);
		}
	};

	// 初始化
	init();
	// 开始迭代
	while (!quit)
	{
		updateSubP();
		QPLagrangeKernel(G, tempA, tempp, tempb, tempx);
		delta = tempx.head(n);
		lambda = tempx.tail(tempx.rows() - n);
		if (delta.norm() < zeroAct)
		{
			// 说明deleta为0
			// 更新S
			int minRow, maxCol;
			double minVal = lambda.minCoeff(&minRow, &maxCol);
			if (minVal > 0)
			{
				quit = true;
			}
			else
			{
				auto iter = S.begin();
				std::advance(iter, minRow);
				S.erase(iter);
			}
		}
		else
		{
			// 说明deleta不为0
			// 更新alpha
			updateAlpha();
			// 更新x
			x += alpha * delta;
		}
	}
	// 返回最优解
	auto f = p.transpose()*x + 0.5*x.transpose()*G*x;
	return f.value();
}

/*!
	拉格朗日法的核心函数
	\li 用于求解等式约束二次优化问题

	\f[
	\begin{aligned}
	\min_{x,y}\quad&\frac{1}{2}x^{T}Gx+p^{T}x\\
	s.t.\quad&A^{T}x=b
	\end{aligned}
	\f]

	\param G 目标函数二次项及交叉乘项系数矩阵(Eigen矩阵)
	\param A 约束条件系数矩阵(Eigen矩阵)
	\param p 目标函数一次项系数矩阵(Eigen向量)
	\param b 约束条件常数向量(Eigen向量)
	\param x 决策变量优化结果向量(Eigen向量)

	\return 目标函数优化结果

	\todo 暂无
*/
double QPLagrangeKernel(const Eigen::MatrixXd& G, const Eigen::MatrixXd& A, const Eigen::VectorXd& p, const Eigen::VectorXd& b, Eigen::VectorXd& x)
{
	// Min: (1 / 2)*x'*G*x + p'*x
	// s.t A'*x = b

	assert(G.rows() == G.cols());
	assert(G.rows() == A.rows());
	assert(G.rows() == p.rows());
	assert(A.cols() == b.rows());
	auto n = G.rows();
	auto m = A.cols();
	// 初始化拉格朗日矩阵K
	Eigen::MatrixXd K(n + m, n + m);
	K.setZero();
	K.topLeftCorner(n, n) = G;
	K.topRightCorner(n, m) = -A;
	K.bottomLeftCorner(m, n) = -A.transpose();
	// 初始化右端向量d
	Eigen::VectorXd d(n + m);
	d.head(n) = -p;
	d.tail(m) = -b;
	// 解方程
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	solver.compute(K.sparseView());
	K.resize(0, 0);
	x = solver.solve(d);
	//x = K.colPivHouseholderQr().solve(d);
	// 返回最优解
	auto f = p.transpose()*x + 0.5*x.transpose()*G*x;
	return f.value();
}