#include "MyOptimalize.h"

using namespace Utility;

double Optimalize::DFPKernel(std::vector<double>& x, std::function<double(std::vector<double>&)> f, std::function<std::vector<double>(std::vector<double>&)> gf)
{
	// x作为迭代的初值，f为目标函数，gf为目标函数梯度函数
	// 运行结束后x为最优解，返回值为优化结果

	// 这里定义了一些必要的参数
	
	// 维度
	int N = x.size();
	// 最大迭代次数
	int maxIterCount = 500;
	// 最大精度
	double maxError = 1e-12;
	// 参数rho
	double rho = 0.5;
	// 参数sigma
	double sigma = 0.2;
	// 当前迭代次数
	int count = 0;
	// 临时变量-当前梯度
	std::vector<double> Gk(N);
	// 临时变量-当前海森矩阵
	Sprase::Matrix<double> Hk(N, N);
	// 临时变量-当前搜索方向
	std::vector<double> dk(N);
	// 临时变量-梯度差
	std::vector<double> yk(N);
	// 临时变量-x迭代差
	std::vector<double> sk(N);

	// lamdba对象，用途：矩阵乘法(C=A*B)
	auto matTimeOp = [&](Sprase::Matrix<double>& A, Sprase::Matrix<double>& B, Sprase::Matrix<double>& C) {
		assert(A.GetSize().second == B.GetSize().first);
		Sprase::Matrix<double> result(C.GetSize().first, C.GetSize().second);
		for (int i = 0; i < A.GetSize().first; i++)
		{
			for (int j = 0; j < A.GetSize().second; j++)
			{
				for (int k = 0; k < A.GetSize().second; k++)
				{
					result[i][j] += A[i][k] * B[k][j];
				}
			}
		}
		C.Init(std::move(result));
	};
	// lambda对象，用途：向量乘法(C=a*b)
	auto vectorTimeOp = [&](std::vector<double>& a, std::vector<double>& b, Sprase::Matrix<double>& C) {
		assert(a.size()==b.size());
		for (int i = 0; i < a.size(); i++)
		{
			for (int j = 0; j < b.size(); j++)
			{
				C[i][j] = a[i] * b[j];
			}
		}
	};
	// lambda对象，用途：向量乘法(c=a.*b)
	auto vectorDotOp = [&](std::vector<double>& a, std::vector<double>& b, double& c) {
		assert(a.size() == b.size());
		for (int i = 0; i < a.size(); i++)
		{
			c+= a[i] * b[i];
		}
	};
	// lambda对象，用途：向量矩阵乘法(c=A*b)
	auto matVectorTimeOp = [&](Sprase::Matrix<double>& A, std::vector<double>& b, std::vector<double>& c, bool tran=false) {
		assert(A.GetSize().second == b.size());
		if (tran)
		{
			for (int i = 0; i < A.GetSize().first; i++)
			{
				for (int j = 0; j < b.size(); j++)
				{
					c[i] -= A[i][j] * b[j];
				}
			}
		}
		else
		{
			for (int i = 0; i < A.GetSize().first; i++)
			{
				for (int j = 0; j < b.size(); j++)
				{
					c[i] += A[i][j] * b[j];
				}
			}
		}
	};

	// 初始化
	for (int i = 0; i < N; i++)
	{
		Hk[i][i] = 1;
	}
	int mk = 0;
	double f0 = 0, f1 = 0;
	double tempDot=0;
	auto& tempx=yk;
	auto& tempv1 = dk;
	Sprase::Matrix<double> tempH(N, N);
	// 开始迭代
	while (count < maxIterCount)
	{
		// 计算梯度
		Gk = gf(x);
		// 更新搜索方向
		dk = std::vector<double>(N);
		matVectorTimeOp(Hk, Gk, dk, true);
		// 更新mk，Armijo准则线性搜索
		mk = 0;
		tempDot = 0;
		vectorDotOp(Gk, dk, tempDot);
		while (mk < 200)
		{
			for (int i = 0; i < N; i++)
			{
				tempx[i] = x[i] + std::pow(rho, mk)*dk[i];
			}
			f1 = f(tempx);
			f0 = f(x);	
			if (f1 < f0 + sigma * std::pow(rho, mk)*tempDot)
			{
				break;
			}
			else
			{
				mk += 1;
			}
		}
		// 执行DFP校正
		// 更新x和sk
		tempDot = std::pow(rho, mk);
		for (int i = 0; i < N; i++)
		{
			sk[i] = tempDot *dk[i];
		}
		tempDot = 0;
		for (int i = 0; i < N; i++)
		{
			tempDot += sk[i] * sk[i];
		}
		if (tempDot < maxError)
		{
			break;
		}
		for (int i = 0; i < N; i++)
		{
			x[i] += sk[i];
		}
		// 更新Gk和yk
		yk = std::move(Gk);
		Gk = gf(x);
		for (int i = 0; i < N; i++)
		{
			yk[i] = Gk[i] - yk[i];
		}
		// 更新Hk
		tempDot = 0;
		vectorDotOp(sk, yk, tempDot);
		if (tempDot > 0)
		{
			tempv1 = std::vector<double>(N);
			matVectorTimeOp(Hk, yk, tempv1, false);
			tempDot = 0;
			vectorDotOp(yk, tempv1, tempDot);
			vectorTimeOp(tempv1, yk, tempH);
			matTimeOp(tempH, Hk, tempH);
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					Hk[i][j] -= tempH[i][j] / tempDot;
				}
			}
			tempDot = 0;
			vectorDotOp(sk, yk, tempDot);
			vectorTimeOp(sk, sk, tempH);
			for (int i = 0; i < N; i++)
			{
				for (int j = 0; j < N; j++)
				{
					Hk[i][j] += tempH[i][j] / tempDot;
				}
			}
		}
		count += 1;
	}
	return f(x);
}