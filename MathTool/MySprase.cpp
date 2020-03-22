#include "MySprase.h"

using namespace Utility;

std::vector<double> Sprase::GEKernel(Matrix<double>& A, std::vector<double>& b)
{
	auto size = A.GetSize();
	assert(size.first == size.second);
	std::vector<double> x(size.first);

	// 选列主元函数
	auto selectPivot = [&](TIndex k) {
		// std::cout << k << " " << std::endl;
		auto max = 0.0;
		auto q = 0;
		// 选列主元
		for (TIndex i = k; i < size.first; i++)
		{
			if (std::abs(A[i][k]) > max)
			{
				max = A[i][k];
				q = i;
			}
		}
		assert(std::abs(max) > mini_value);
		// 交换行
		if (q != k)
		{
			std::swap(A[q], A[k]);
			std::swap(b[q], b[k]);
		}
	};
	// 消去函数
	auto elimination = [&]() {
		for (TIndex k = 0; k < size.first-1; k++)
		{
			// 选列主元
			// selectPivot(k);
			for (TIndex i = k + 1; i < size.first; i++)
			{
				A[i][k] = A[i][k] / A[k][k];
				for (TIndex j = k + 1; j < size.first; j++)
				{
					A[i][j] = A[i][j] - A[i][k] * A[k][j];
				}
				b[i] = b[i] - A[i][k] * b[k];
				A[i][k] = 0;
			}
		}
	};
	// 回代函数
	auto back = [&]() {
		x[size.first - 1] = b[size.first - 1] / A[size.first - 1][size.first - 1];
		for (int k = size.first - 2; k > -1; k--)
		{
			auto s = b[k];
			for (TIndex j = k + 1; j < size.first; j++)
			{
				s = s - A[k][j] * x[j];
			}
			x[k] = s / A[k][k];
		}
	};

	// 求解
	elimination();
	back();
	return std::move(x);
}


std::vector<double> Sprase::LUDecomposition(Matrix<double>& A)
{
	auto size = A.GetSize();
	assert(size.first == size.second);
	std::vector<double> pi(size.first);
	for (TIndex i = 0; i < size.first; i++)
	{
		pi[i] = i;
	}

	// 选主元函数
	auto selectPivot = [&](TIndex k) {
		auto max = 0.0;
		auto q = 0;
		// 选列主元
		for (TIndex i = k; i < size.first; i++)
		{
			if (std::abs(A[i][k]) > max)
			{
				max = A[i][k];
				q = i;
			}
		}
		assert(std::abs(max) > mini_value);
		// 交换行
		if (q != k)
		{
			A[q].swap(A[k]);
			std::swap(pi[q], pi[k]);
		}
	};
	// 分解函数
	auto decomposition = [&]() {
		// 处理高阶
		for (TIndex k = 0; k < size.first ; k++)
		{
			// selectPivot(k);
			for (TIndex i = k + 1; i < size.first; i++)
			{
				auto& down = A[i][k];
				if (down)
				{
					down /= A[k][k];
				}
				for (TIndex j = k + 1; j < size.first; j++)
				{
					auto& right = A[k][j];
					if (down && right)
					{
						A[i][j] -= down * right;
					}
				}
			}
		}
	};

	decomposition();
	return std::move(pi);
}


std::vector<double> Sprase::LUSolution(Matrix<double>& A, std::vector<double>& b, std::vector<double>& pi)
{
	auto size = A.GetSize();
	assert(size.first == size.second);
	std::vector<double> x(size.first);
	std::vector<double> y(size.first);

	// 正向替换函数
	auto forwardSolve = [&]() {
		for (TIndex i = 0; i < size.first; i++)
		{
			auto sum = 0.0;
			for (TIndex j = 0; j < i; j++)
			{
				auto& right = A[i][j];
				if (right)
				{
					sum += right * y[j];
				}
			}
			y[i] = b[pi[i]] - sum;
		}
	};
	// 逆向替换函数
	auto backSolve = [&]() {
		for (int i = size.first - 1; i > -1; i--)
		{
			auto sum = 0.0;
			for (TIndex j = i + 1; j < size.first; j++)
			{
				auto& right = A[i][j];
				if (right)
				{
					sum += right * x[j];
				}
			}
			x[i] = (y[i] - sum) / A[i][i];
		}
	};

	// 求解
	forwardSolve();
	backSolve();
	return std::move(x);
}

std::vector<double> Sprase::SORKernel(SpraseCSR<double>&A, std::vector<double>& b, std::vector<double>& xPre, double weight)
{
	assert(weight < 2 && weight>0);
	auto size = A.GetSize();
	assert(size.first == size.second);
	// assert(b.size() >= size.first && xPre.size >= size.first);
	// 第k次x的结果
	std::vector<double> x0(size.first);
	// 第k+1次x的结果
	std::vector<double> x(xPre);
	// 残差
	double residual = 0;
	// 迭代次数
	int count = 0;
	// k次迭代过程
	auto forEachStep = [&]()
	{
		int i = 0, j = 0;
		// k+1次结果参与的累加/k次结果参与的累加
		double xSum = 0.0, x0Sum = 0.0;
		double pivot = 0;
		// 修改x0
		x0 = x;
		for (int row = 0; row < size.first; row++)
		{
			// 选中i行
			auto curItem = A.SelectRow(i, j, row);
			while ((curItem.second != false) && (j<i))
			{
				xSum += curItem.first*(x[j]);
				curItem = A.Right(i, j);
			}
			assert(i == j && curItem.first != 0);
			pivot = curItem.first;
			curItem = A.Right(i, j);
			while (curItem.second != false)
			{
				x0Sum += curItem.first*(x0[j]);
				curItem = A.Right(i, j);
			}
			// 修改x
			x[row] = (1 - weight) * x0[row] + weight * (b[row] - xSum - x0Sum) / pivot;
			// 累加值清零
			x0Sum = 0;
			xSum = 0;
		}
	};
	// 判断收敛
	auto checkConverge = [&]()
	{
		residual = 0;
		for (int i = 0; i < size.first; i++)
		{
			residual += std::abs(x[i] - x0[i]);
		}
		if (residual < SORConvergeLimit)
		{
			return true;
		}
		else
		{
			return false;
		}
	};

	// 执行求解
	do
	{
		forEachStep();
	} while (!checkConverge());
	// 返回结果
	return std::move(x);
}

std::vector<double> Sprase::JACKernel(SpraseCSR<double>&A, std::vector<double>& b)
{
	auto size = A.GetSize();
	assert(size.first == size.second);
	// 第k次x的结果
	std::vector<double> x0(size.first);
	// 第k+1次x的结果
	std::vector<double> x(size.first);
	// 残差
	double residual = 0;
	// 迭代次数
	int count = 0;
	// k次迭代过程
	auto forEachStep = [&]()
	{
		int i = 0, j = 0;
		// k+1次结果参与的累加/k次结果参与的累加
		double x0Sum = 0.0;
		double pivot = 0;
		// 修改x0
		x0 = x;
		for (int row = 0; row < size.first; row++)
		{
			// 选中i行
			auto curItem = A.SelectRow(i, j, row);
			while (curItem.second != false)
			{
				if (j != i)
				{
					x0Sum += curItem.first*(x[j]);					
				}
				else
				{
					assert(curItem.first != 0);
					pivot = curItem.first;
				}
				curItem = A.Right(i, j);
			}			
			// 修改x
			x[row] = (b[row] - x0Sum) / pivot;
			// 累加值清零
			x0Sum = 0;
		}
	};
	// 判断收敛
	auto checkConverge = [&]()
	{
		residual = 0;
		for (int i = 0; i < size.first; i++)
		{
			residual += std::abs(x[i] - x0[i]);
		}
		if (residual < JACConvergeLimit)
		{
			return true;
		}
		else
		{
			return false;
		}
	};

	// 执行求解
	do
	{
		forEachStep();
	} while (!checkConverge());
	// 返回结果
	return std::move(x);
}

std::pair<std::vector<double>, std::vector<double> > Sprase::QRDecomposition(Matrix<double>&A, bool isExtension)
{
	// A的维度
	auto size = A.GetSize();
	// 设置行列
	int m = size.first;
	int n = isExtension ? size.second - 1 : size.second;
	// 初始化对角向量d和数量\alpha
	std::vector<double> d(n);
	std::vector<double> alpha(n);
	// 开始分解
	for (int k = 0; k < n - 1; k++)
	{
		double sigma = 0.0;
		for (int i = k; i < m ; i++)
		{
			sigma += A[i][k] * A[i][k];
		}
		sigma = std::sqrt(sigma);
		// sigma取反
		if (A[k][k] >= 0)
		{
			sigma = -sigma;
		}
		d[k] = sigma;
		alpha[k] = sigma * (sigma - A[k][k]);
		// 将Akk修改为存放Wk的第一个元素
		A[k][k] = A[k][k] - sigma;
		// H分解
		// 若此时分解的是增广矩阵
		auto end = isExtension ? n + 1 : n;
		for (int j = k + 1; j < end; j++)
		{
			// 构造Belta
			double belta = A[k][k] * A[k][j];
			for (int i = k + 1; i < m; i++)
			{
				belta += A[i][k] * A[i][j];
			}
			belta /= alpha[k];
			A[k][j] = A[k][j] - belta * A[k][k];
			for (int i = k + 1; i < m; i++)
			{
				A[i][j] = A[i][j] - belta * A[i][k];
			}
		}
	}
	// 阶数为n时，单独处理
	{
		double sigma = 0.0;
		for (int i = n - 1; i < m; i++)
		{
			sigma += A[i][n - 1] * A[i][n - 1];
		}
		sigma = std::sqrt(sigma);
		if (A[n - 1][n - 1] >= 0)
		{
			sigma = -sigma;
		}
		d[n - 1] = sigma;
		alpha[n - 1] = sigma * (sigma - A[n - 1][n - 1]);
		A[n - 1][n - 1] -= sigma;
	}
	if (isExtension)
	{
		// 若此时分解的是增广矩阵
		double belta = A[n-1][n-1] * A[n-1][n];
		for (int i = n; i < m; i++)
		{
			belta += A[i][n-1] * A[i][n];
		}
		belta /= alpha[n-1];
		A[n - 1][n] -= belta * A[n - 1][n - 1];
		for (int i = n; i < m; i++)
		{
			A[i][n] -= belta * A[i][n - 1];
		}
	}

	// 分别返回对角数组和数量数组，用于构成H=QT
	return { std::move(d), std::move(alpha) };
}


std::vector<double> Sprase::QRKernel(Matrix<double>&A, std::vector<double>& b)
{
	// A的维度
	auto size = A.GetSize();
	int m = size.first, n = size.second;
	std::vector<double> d;
	std::vector<double> x(n);

	// lambda对象，用途：生成增广阵
	auto exten = [&]() {
		for (int i = 0; i < m; i++)
		{
			A[i].push_back(b[i]);
		}
	};
	// lambda对象，用途：执行QR分解
	auto decomposition = [&]() {
		d = QRDecomposition(A, true).first;
	};
	// lambda对象，通途：执行回代
	auto backSolve = [&]() {
		x[n - 1] = A[n - 1][n] / d[n - 1];
		for (int k = n - 2; k > -1; k--)
		{
			x[k] = A[k][n];
			for (int j = k + 1; j < n; j++)
			{
				x[k] -= A[k][j] * x[j];
			}
			x[k] /= d[k];
		}
	};

	// 求解
	exten();
	decomposition();
	backSolve();
	//  返回结果
	return std::move(x);
}

std::vector<double> Sprase::CGKernel(Matrix<double>& A, std::vector<double>& b)
{
	// A的维度
	auto size = A.GetSize();
	assert(size.first == size.second);
	// 同时也是迭代次数理论上限
	auto N = size.first;
	// 结果向量x
	std::vector<double> x(N);
	// 步长alpha
	double alpha = 0.0;
	// 临时向量belta
	double belta;
	// 方向向量d
	std::vector<double> d(N);
	// 残差向量r(第k+1次)
	std::vector<double> r(N);
	// 残差向量内积 r0Dot(第k次)
	double r0Dot = 0.0;
	// 当前迭代次数
	int count = 0;

	// lambda对象，通途：计算向量内积
	auto dotProduct = [&](std::vector<double>& v1, std::vector<double>& v2) {
		assert(v1.size() == v2.size());
		auto vSize = v1.size();
		double dotP = 0.0;
		for (size_t i = 0; i < vSize; i++)
		{
			dotP += v1[i] * v2[i];
		}
		return dotP;		
	};
	// lambda对象，用途：计算矩阵和列向量的乘积
	auto matTimeVec = [&](Matrix<double>& M, std::vector<double>& v) {
		assert(M.GetSize().second == v.size());
		auto vSize = v.size();
		std::vector<double> result(vSize);
		for (size_t i = 0; i < vSize; i++)
		{
			result[i] = dotProduct(M[i], v);
		}
		return std::move(result);
	};
	// lambda对象，用途：生成残差
	auto updateR = [&]() {
		if (count == 0)
		{
			r0Dot = 0.0;
			// 生成初次残差，顺带更新d
			for (size_t i = 0; i < N; i++)
			{
				d[i] = b[i] - dotProduct(A[i], x);
				r0Dot += d[i] * d[i];
			}
		}
		else
		{
			// 生成k次残差
			for (size_t i = 0; i < N; i++)
			{
				r[i]=b[i]- dotProduct(A[i], x);
			}
		}
	};
	// lambda对象，用途：更新x
	auto updateX = [&]() {
		for (size_t i = 0; i < N; i++)
		{
			x[i] = x[i] + alpha * d[i];
		}
	};
	// lambda对象，用途：更新d
	auto updateD = [&]() {
		for (size_t i = 0; i < N; i++)
		{
			d[i] = r[i] + belta * d[i];
		}
	};
	// lambda对象，用途：判断收敛
	auto checkout = [&]()
	{
		double error = 0.0;
		for (size_t i = 0; i < N; i++)
		{
			error += r[i] * r[i];
		}
		error = std::sqrt(error);
		// 收敛条件为：迭代次数大于N或者此时精度达到收敛精度要求
		return (count > N || error < CGConvergeLimit) ? true : false;
	};

	// 开始计算
	updateR();
	for (count = 1; ; count++)
	{
		auto tempAd = matTimeVec(A, d);
		alpha = r0Dot / dotProduct(d, tempAd);
		updateX();
		updateR();
		if (checkout())
		{
			break;
		}
		auto tempRR = dotProduct(r, r);
		belta = tempRR / r0Dot;
		updateD();
		// 更新r0
		r0Dot = tempRR;
	}
	// 返回结果
	return std::move(x);
}

std::vector<double> Sprase::CMKernel(SpraseTD<double>& A, std::vector<double>& d)
{
	// A的维度
	auto size = A.GetSize();
	assert(size.first == size.second);
	auto N = size.first;
	// 存放L和U的临时向量
	std::vector<double> L(N), U(N);
	// 存放x和y的临时向量
	std::vector<double> x(N), y(N);

	// lambda对象，用途：追赶过程
	auto chasing = [&]() {
		U[0] = A.AtDiagBRef(0);
		y[0] = d[0];
		for (size_t i = 1; i < N; i++)
		{
			L[i] = A.AtDiagARef(i-1) / U[i - 1];
			U[i] = A.AtDiagBRef(i) - L[i] * A.AtDiagCRef(i - 1);
			y[i] = d[i] - L[i] * y[i - 1];
		}
	};
	// lambda对象，用途：回代过程
	auto back = [&]() {
		x[N - 1] = y[N - 1] / U[N - 1];
		for (int i = N - 2; i > -1; i--)
		{
			x[i] = (y[i] - A.AtDiagCRef(i)*x[i + 1]) / U[i];
		}
	};

	// 开始计算
	chasing();
	back();
	// 返回结果
	return std::move(x);	
}