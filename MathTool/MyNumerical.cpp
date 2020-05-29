/*! \file MyNumerical.cpp */

#include "MyNumerical.h"
//���������
#include <Eigen/Eigen>

using namespace Utility;

// ����һ�±ƽ��������������ֵ����ȡ�ĵ���(��ʱ)
constexpr int MaxPointNum = 1000;

/*!
	��������ַ�
	\li ����ԭ������

	\f[T_1 = \frac{b-a}{2}[f(a)-f(b)]\f]
	\f[T_{2^{k+1}} = \frac{T_{2^{k}}}{2} + \frac{b-a}{2^{k+1}}\sum_{i=1}^{2^{k}}f(a+(2i-1)\frac{b-a}{2^{k+1}})\f]
	\f[S_{2^{k}} = T_{2^{k+1}} + \frac{T_{2^{k+1}}-T_{2^{k}}}{3}\f]
	\f[C_{2^{k}} = S_{2^{k+1}} + \frac{S_{2^{k+1}}-S_{2^{k}}}{15}\f]
	\f[R_{2^{k}} = C_{2^{k+1}} + \frac{C_{2^{k+1}}-C_{2^{k}}}{63}\f]

	\param begin ���������׶�
	\param end ��������ĩ��
	\param fun ��������
	\param converge �����о�

	\return ���ֽ��

	\todo ����
*/
double Numerical::RombergIntegration(const double begin, const double end, std::function<double(double)> fun, const double converge)
{
	assert(begin < end);
	assert(converge > 0 && converge < 1);
	// ��k�μ�����
	double T0 = 0, S0 = 0, C0 = 0, R0 = 0;
	// ��k+1�μ�����
	double Tk = 0, Sk = 0, Ck = 0, Rk = 0;
	// ������
	int k = 0;

	// lambda������;�����λ��ֵݹ鹫ʽ
	auto updateT = [&]() {
		// ����T0
		T0 = Tk;
		// ����Tk
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
	// lambda������;������ɭ���ֵݹ鹫ʽ
	auto updateS = [&]() {
		// ����S0
		S0 = Sk;
		// ����Sk
		Sk = Tk;
		Sk += (Tk - T0) / 3;
	};
	// lambda������;�����Ȼ��ֵݹ鹫ʽ
	auto updateC = [&]() {
		// ����C0
		C0 = Ck;
		// ����Ck
		Ck = Sk;
		Ck += (Sk - S0) / 15;
	};
	// lambda������;����������ֵݹ鹫ʽ
	auto updateR = [&]() {
		// ����R0
		R0 = Rk;
		// ����Sk
		Rk = Ck;
		Rk += (Ck - C0) / 63;
	};
	// lambda������;����ʼ��
	auto prepare = [&]() {
		// ��ǰ���н���һ�����㣬�õ���һ��R
		// ����ע�����±��Ϊk��ֵ
		// ���ֱ��һ�У�����T0
		Tk = (end - begin) * (fun(begin) + fun(end)) / 2;
		// ���ֱ�ڶ��У�����T1��S0
		updateT();
		updateS();
		k++;
		// ���ֱ�����У�����T2��S1��C0
		updateT();
		updateS();
		updateC();
		k++;
		// ���ֱ�����У�����T3��S2��C1��R0
		updateT();
		updateS();
		updateC();
		updateR();
		k++;
	};
	// lambda������;���ж�����
	auto check = [&]() {
		return std::abs(Rk - R0) < converge;
	};

	// ��ʼ�����ֱ�ǰ����
	prepare();
	// �����������������ֵ��ֱ������
	while (!check())
	{
		updateT();
		updateS();
		updateC();
		updateR();
		k++;
	}
	// ���ػ��ֽ��
	return Rk;
}

/*!
	Ĭ��Ϊ����ʽ��ʽ��һԪ��С�������
	\li �������У���0���ݿ�ʼ

	\param v ��ֵ����
	\param b ����ֵ����
	\param n ���Ĵ���

	\return ����ʽϵ����(��������)

	\todo ����
*/
std::vector<double> Numerical::LeastSquareFitting(std::vector<double>& v, std::vector<double>& b, int n)
{
	// �������ʽ��������
	std::vector<std::function<double(double)> > f;
	f.reserve(n + 1);
	for (int i = 0; i < n + 1; i++)
		f.emplace_back([i](double x)->double {return std::pow(x, i); });
	// ����Ȩֵ���У�Ĭ��Ϊ1
	std::vector<double> w(v.size(), 1.0);
	// ������ͨ�汾
	return LeastSquareFitting(v, b, f, w);
}

/*!
	�����Զ�������������Ȩֵ��һԪ��С�������
	\li ������������Ϊ�����ڻ�Ϊ0�ĺ�����

	\param v ��ֵ����
	\param b ����ֵ����
	\param f ������������
	\param w Ȩ��ֵ����

	\return ����ʽϵ����(��f��˳�����Ӧ)

	\todo ����
*/
std::vector<double> Numerical::LeastSquareFitting(std::vector<double>& v, std::vector<double>& b,
	std::vector<std::function<double(double)> >& f, std::vector<double>& w)
{
	// ������
	assert(v.size() == b.size());
	assert(v.size() == w.size());
	// lambda������;���������������������ڻ�
	auto funDotProduct1 = [&](int i, int j)->double {
		double sum = 0.0;
		for (int k = 0; k < v.size(); k++)
			sum += f[i](v[k])*f[j](v[k])*w[k];
		return sum;
	};
	// lambda������;������һ������������ԭ�������ڻ�
	auto funDotProduct2 = [&](int i)->double {
		double sum = 0.0;
		for (int k = 0; k < v.size(); k++)
			sum += f[i](v[k])*b[k] * w[k];
		return sum;
	};
	// ��ʼ��Eigen����
	Eigen::VectorXd eb;
	Eigen::MatrixXd eA;
	eb.resize(f.size());
	eA.resize(f.size(), f.size());
	for (int i = 0; i < f.size(); i++)
		for (int j = 0; j < f.size(); j++)
			eA.coeffRef(i, j) = funDotProduct1(i, j);
	for (int i = 0; i < f.size(); i++)
		eb.coeffRef(i) = funDotProduct2(i);
	// ʹ��QR�ֽ���ⷽ����Ac=b
	Eigen::VectorXd ec = eA.colPivHouseholderQr().solve(eb);
	// ȡ��ec�����ؽ��
	std::vector<double> result(f.size());
	for (int i = 0; i < f.size(); i++)
		result[i] = ec.coeff(i);
	return std::move(result);
}

/*!
	��Ԫ������С�������

	\param A ��ֵ����(�����ȣ����Ի��洢)����һάΪ��ֵ���б�ţ��ڶ�άΪ����ֵ���
	\param b ����ֵ����
	\param m ��ֵ���е�����
	\param m ����ֵ���еĳ���

	\return ����ʽϵ����(���ֵ����˳�����Ӧ(�������������Ӧ�����ϣ������������������Ӧ������λ����))

	\todo ����ֵ��ϢӦ���޸ģ����stat��Ϣ����ԭ���ķ���ֵ�滻�����������
*/
std::vector<double> Numerical::MutliLinearLeastSquareFitting(std::vector<double>& A, std::vector<double>& b, int n, int m)
{
	// ������
	assert(A.size() == n * m);
	assert(b.size() == m);
	// ��ʼ��Eigen����
	Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, 1> >
		eb(b.data(), b.size(), 1, Eigen::Stride<Eigen::Dynamic, 1>(b.size(), 1));
	// ע�������eA�൱��A^T
	Eigen::Map<Eigen::Matrix<double, -1, -1, 1>, 0, Eigen::Stride<-1, -1> >
		eA(A.data(), n, m, Eigen::Stride<-1, -1>(m, 1));
	// ֱ�Ӽ���ec=(A^T*A)^-1*A^T*b
	Eigen::VectorXd ec = (eA*eA.transpose()).inverse()*(eA*eb);
	// ȡ��ec������
	std::vector<double> result(n);
	for (int i = 0; i < n; i++)
		result[i] = ec.coeff(i);
	return std::move(result);
}


/*!
	Ĭ��Ϊ����ʽ��ʽ��һԪ����һ�±ƽ�(Remes�㷨)

	\param f ��Ҫ�ƽ��ĺ���
	\param begin ����ǰ�˵�
	\param end ����ĩ�˵�
	\param n �ƽ�����ʽ����
	\param maxCount ����������
	\param error ������

	\return ����ʽϵ����(�����������һ��ֵΪ���En)

	\todo Ѱ�ҵ�ǰ�ıƽ�����ʽ��ʵ�ʺ����������λ�õ��㷨��Ҫ�Ż�
*/
std::vector<double> Numerical::BestUniformApproximation(std::function<double(double)> f,
	double begin, double end, int n, int maxCount, double error)
{
	// �������
	std::vector<double> x;
	// ��������Ӧ�ĺ���ֵ
	std::vector<double> y;
	// �ƽ�����ʽ��ϵ��
	std::vector<double> d;
	// ��ǰ������
	double curError = 0.0;
	// lambda������;�����ض���ʽ��ֵ
	auto getPolyval = [&](double point)->double {
		double sum = 0.0;
		for (int i = 0; i < n + 1; i++)
			sum += d[i] * std::pow(point, i);
		return sum;
	};
	// lambda������;����ʼ��
	auto init = [&]() {
		// ���ɳ�ʼ�������
		double h = (double)(end - begin) / (n + 1);
		x.resize(n + 2);
		y.resize(n + 2);
		for (int i = 0; i < n + 2; i++)
		{
			x[i] = begin + i * h;
			y[i] = f(x[i]);
		}
		// ����ϵ��
		d.resize(n + 2);
	};
	// lambda������;�����½������
	auto updateX = [&]() -> bool {
		// ����Ѱ�ҵ�ǰ�ıƽ�����ʽ��ʵ�ʺ����������λ�ã��ⲿ����Ҫ�Ż�
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
		// �����Ѿ�����õ��������x���½����
		{
			int k = 0;
			for (k = 0; k < n + 2; k++)
			{
				if (x[k] > maxErrorx)
					break;
			}
			x[k] = maxErrorx;
		}
		// ���½�����鴦����ֵ
		for (int i = 0; i < n + 2; i++)
			y[i] = f(x[i]);
		return false;
	};
	// lambda������;�������б�ѩ����������Է�����õ�����ʽ
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
		// ʹ��QR�ֽ���ⷽ����Ac=b
		Eigen::VectorXd ec = eA.colPivHouseholderQr().solve(eb);
		// ����
		for (int i = 0; i < n + 2; i++)
			d[i] = ec.coeff(i);
	};

	// ִ�м���
	init();
	int count = 0;
	while (count++ < maxCount)
	{
		chebyshev();
		if (updateX())
			break;
	}
	// ������Ϸ��ؽ��
	return std::move(d);
}
