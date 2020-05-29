/*! \file MyNumerical.cpp */

#include "MyNumerical.h"
//矩阵运算包
#include <Eigen/Eigen>

using namespace Utility;

// 最优一致逼近中用于搜索最大值所抽取的点数(临时)
constexpr int MaxPointNum = 1000;

/*!
	龙贝格积分法
	\li 计算原理如下

	\f[T_1 = \frac{b-a}{2}[f(a)-f(b)]\f]
	\f[T_{2^{k+1}} = \frac{T_{2^{k}}}{2} + \frac{b-a}{2^{k+1}}\sum_{i=1}^{2^{k}}f(a+(2i-1)\frac{b-a}{2^{k+1}})\f]
	\f[S_{2^{k}} = T_{2^{k+1}} + \frac{T_{2^{k+1}}-T_{2^{k}}}{3}\f]
	\f[C_{2^{k}} = S_{2^{k+1}} + \frac{S_{2^{k+1}}-S_{2^{k}}}{15}\f]
	\f[R_{2^{k}} = C_{2^{k+1}} + \frac{C_{2^{k+1}}-C_{2^{k}}}{63}\f]

	\param begin 积分区间首端
	\param end 积分区间末端
	\param fun 被积函数
	\param converge 收敛判据

	\return 积分结果

	\todo 暂无
*/
double Numerical::RombergIntegration(const double begin, const double end, std::function<double(double)> fun, const double converge)
{
	assert(begin < end);
	assert(converge > 0 && converge < 1);
	// 第k次计算结果
	double T0 = 0, S0 = 0, C0 = 0, R0 = 0;
	// 第k+1次计算结果
	double Tk = 0, Sk = 0, Ck = 0, Rk = 0;
	// 计数器
	int k = 0;

	// lambda对象，用途：梯形积分递归公式
	auto updateT = [&]() {
		// 更新T0
		T0 = Tk;
		// 更新Tk
		Tk = 0.5*T0;
		double sum = 0.0;
		double point = 0.0;
		double power = std::pow(2, k);
		double doublePower = (end - begin) / power / 2;
		for (int i = 0; i < (int)power; i++)
		{
			point = begin + (2 * i + 1) * doublePower;
			sum += fun(point);
		}
		Tk += sum * doublePower;
	};
	// lambda对象，用途：辛普森积分递归公式
	auto updateS = [&]() {
		// 更新S0
		S0 = Sk;
		// 更新Sk
		Sk = Tk;
		Sk += (Tk - T0) / 3;
	};
	// lambda对象，用途：科兹积分递归公式
	auto updateC = [&]() {
		// 更新C0
		C0 = Ck;
		// 更新Ck
		Ck = Sk;
		Ck += (Sk - S0) / 15;
	};
	// lambda对象，用途：龙贝格积分递归公式
	auto updateR = [&]() {
		// 更新R0
		R0 = Rk;
		// 更新Sk
		Rk = Ck;
		Rk += (Ck - C0) / 63;
	};
	// lambda对象，用途：初始化
	auto prepare = [&]() {
		// 对前四行进行一个计算，得到第一个R
		// 下列注释中下标均为k的值
		// 积分表第一行，计算T0
		Tk = (end - begin) * (fun(begin) + fun(end)) / 2;
		// 积分表第二行，计算T1、S0
		updateT();
		updateS();
		k++;
		// 积分表第三行，计算T2、S1、C0
		updateT();
		updateS();
		updateC();
		k++;
		// 积分表第四行，计算T3、S2、C1、R0
		updateT();
		updateS();
		updateC();
		updateR();
		k++;
	};
	// lambda对象，用途：判断收敛
	auto check = [&]() {
		return std::abs(Rk - R0) < converge;
	};

	// 初始化积分表前四行
	prepare();
	// 迭代计算龙贝格积分值，直至收敛
	while (!check())
	{
		updateT();
		updateS();
		updateC();
		updateR();
		k++;
	}
	// 返回积分结果
	return Rk;
}

/*!
	默认为多项式形式的一元最小二乘拟合
	\li 升幂排列，从0次幂开始

	\param v 点值序列
	\param b 函数值序列
	\param n 最大的次幂

	\return 多项式系数表(升幂排序)

	\todo 暂无
*/
std::vector<double> Numerical::LeastSquareFitting(std::vector<double>& v, std::vector<double>& b, int n)
{
	// 构造多项式函数序列
	std::vector<std::function<double(double)> > f;
	f.reserve(n + 1);
	for (int i = 0; i < n + 1; i++)
		f.emplace_back([i](double x)->double {return std::pow(x, i); });
	// 构造权值序列，默认为1
	std::vector<double> w(v.size(), 1.0);
	// 调用普通版本
	return LeastSquareFitting(v, b, f, w);
}

/*!
	允许自定义正交函数和权值的一元最小二乘拟合
	\li 正交函数定义为两两内积为0的函数族

	\param v 点值序列
	\param b 函数值序列
	\param f 正交函数序列
	\param w 权重值序列

	\return 多项式系数表(与f的顺序相对应)

	\todo 暂无
*/
std::vector<double> Numerical::LeastSquareFitting(std::vector<double>& v, std::vector<double>& b,
	std::vector<std::function<double(double)> >& f, std::vector<double>& w)
{
	// 检查对齐
	assert(v.size() == b.size());
	assert(v.size() == w.size());
	// lambda对象，用途：计算两个正交函数的内积
	auto funDotProduct1 = [&](int i, int j)->double {
		double sum = 0.0;
		for (int k = 0; k < v.size(); k++)
			sum += f[i](v[k])*f[j](v[k])*w[k];
		return sum;
	};
	// lambda对象，用途：计算一个正交函数和原函数的内积
	auto funDotProduct2 = [&](int i)->double {
		double sum = 0.0;
		for (int k = 0; k < v.size(); k++)
			sum += f[i](v[k])*b[k] * w[k];
		return sum;
	};
	// 初始化Eigen对象
	Eigen::VectorXd eb;
	Eigen::MatrixXd eA;
	eb.resize(f.size());
	eA.resize(f.size(), f.size());
	for (int i = 0; i < f.size(); i++)
		for (int j = 0; j < f.size(); j++)
			eA.coeffRef(i, j) = funDotProduct1(i, j);
	for (int i = 0; i < f.size(); i++)
		eb.coeffRef(i) = funDotProduct2(i);
	// 使用QR分解求解方程组Ac=b
	Eigen::VectorXd ec = eA.colPivHouseholderQr().solve(eb);
	// 取回ec并返回结果
	std::vector<double> result(f.size());
	for (int i = 0; i < f.size(); i++)
		result[i] = ec.coeff(i);
	return std::move(result);
}

/*!
	多元线性最小二乘拟合

	\param A 点值矩阵(行优先，线性化存储)，第一维为点值序列编号，第二维为序列值编号
	\param b 函数值序列
	\param m 点值序列的总数
	\param m 各点值序列的长度

	\return 多项式系数表(与点值序列顺序相对应(与输入数据相对应，如果希望常数项则输入数据应包含单位向量))

	\todo 返回值信息应当修改，添加stat信息，而原本的返回值替换到输入参数中
*/
std::vector<double> Numerical::MutliLinearLeastSquareFitting(std::vector<double>& A, std::vector<double>& b, int n, int m)
{
	// 检查对齐
	assert(A.size() == n * m);
	assert(b.size() == m);
	// 初始化Eigen对象
	Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, 1> >
		eb(b.data(), b.size(), 1, Eigen::Stride<Eigen::Dynamic, 1>(b.size(), 1));
	// 注意这里的eA相当于A^T
	Eigen::Map<Eigen::Matrix<double, -1, -1, 1>, 0, Eigen::Stride<-1, -1> >
		eA(A.data(), n, m, Eigen::Stride<-1, -1>(m, 1));
	// 直接计算ec=(A^T*A)^-1*A^T*b
	Eigen::VectorXd ec = (eA*eA.transpose()).inverse()*(eA*eb);
	// 取回ec并返回
	std::vector<double> result(n);
	for (int i = 0; i < n; i++)
		result[i] = ec.coeff(i);
	return std::move(result);
}


/*!
	默认为多项式形式的一元最优一致逼近(Remes算法)

	\param f 需要逼近的函数
	\param begin 区间前端点
	\param end 区间末端点
	\param n 逼近多项式阶数
	\param maxCount 最大迭代次数
	\param error 最大误差

	\return 多项式系数表(升幂排序，最后一个值为误差En)

	\todo 寻找当前的逼近多项式与实际函数误差最大的位置的算法需要优化
*/
std::vector<double> Numerical::BestUniformApproximation(std::function<double(double)> f,
	double begin, double end, int n, int maxCount, double error)
{
	// 交错点组
	std::vector<double> x;
	// 交错点组对应的函数值
	std::vector<double> y;
	// 逼近多项式的系数
	std::vector<double> d;
	// 当前最大误差
	double curError = 0.0;
	// lambda对象，用途：返回多项式的值
	auto getPolyval = [&](double point)->double {
		double sum = 0.0;
		for (int i = 0; i < n + 1; i++)
			sum += d[i] * std::pow(point, i);
		return sum;
	};
	// lambda对象，用途：初始化
	auto init = [&]() {
		// 生成初始交错点组
		double h = (double)(end - begin) / (n + 1);
		x.resize(n + 2);
		y.resize(n + 2);
		for (int i = 0; i < n + 2; i++)
		{
			x[i] = begin + i * h;
			y[i] = f(x[i]);
		}
		// 生成系数
		d.resize(n + 2);
	};
	// lambda对象，用途：更新交错点组
	auto updateX = [&]() -> bool {
		// 首先寻找当前的逼近多项式与实际函数误差最大的位置，这部分需要优化
		double maxErrorx = begin;
		double maxErrory = 0.0;
		{
			std::default_random_engine randBase(clock());
			std::uniform_real_distribution<double> rand(begin, end);
			std::vector<double> tempx(MaxPointNum);
			std::vector<double> tempy(MaxPointNum);
			for (int i = 0; i < MaxPointNum; i++)
			{
				tempx[i] = rand(randBase);
				tempy[i] = std::abs(f(tempx[i]) - getPolyval(tempx[i]));
			}
			auto index = std::max_element(tempy.begin(), tempy.end()) - tempy.begin();
			maxErrorx = tempx[index];
			maxErrory = tempy[index];
			if (maxErrory < error)
				return true;
		}
		// 利用已经计算得到的最大误差处x更新交错点
		{
			int k = 0;
			for (k = 0; k < n + 2; k++)
			{
				if (x[k] > maxErrorx)
					break;
			}
			x[k] = maxErrorx;
		}
		// 更新交错点组处函数值
		for (int i = 0; i < n + 2; i++)
			y[i] = f(x[i]);
		return false;
	};
	// lambda对象，用途：根据切比雪夫定理求解线性方程组得到多项式
	auto chebyshev = [&]() {
		Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, 1> >
			eb(y.data(), y.size(), 1, Eigen::Stride<Eigen::Dynamic, 1>(y.size(), 1));
		Eigen::MatrixXd eA;	
		eA.resize(n + 2, n + 2);
		for (int i = 0; i < n + 2; i++)
		{
			for (int j = 0; j < n + 1; j++)
				eA.coeffRef(i, j) = std::pow(x[i], j);
			eA.coeffRef(i, n + 2) = std::pow(-1, i + 1);
		}
		// 使用QR分解求解方程组Ac=b
		Eigen::VectorXd ec = eA.colPivHouseholderQr().solve(eb);
		// 更新
		for (int i = 0; i < n + 2; i++)
			d[i] = ec.coeff(i);
	};

	// 执行计算
	init();
	int count = 0;
	while (count++ < maxCount)
	{
		chebyshev();
		if (updateX())
			break;
	}
	// 计算完毕返回结果
	return std::move(d);
}
