#include "MyNumerical.h"

using namespace Utility;

/*
	�������������ʽ

	T_1 = \frac{b-a}{2}[f(a)-f(b)]

	for k =0 to ...
		T_{2^{k+1}} = \frac{T_{2^{k}}}{2} + \frac{b-a}{2^{k+1}}\sum_{i=1}^{2^{k}}f(a+(2i-1)\frac{b-a}{2^{k+1}})
		S_{2^{k}} = T_{2^{k+1}} + \frac{T_{2^{k+1}}-T_{2^{k}}}{3}
		C_{2^{k}} = S_{2^{k+1}} + \frac{S_{2^{k+1}}-S_{2^{k}}}{15}
		R_{2^{k}} = C_{2^{k+1}} + \frac{C_{2^{k+1}}-C_{2^{k}}}{63}
	end
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


Numerical::SegmentFunction<double, double> Numerical::CubicSplineInterpolation(std::vector<std::pair<double, double> > points, std::pair<double, double> conditions)
{
	// ��ֵ����ʽ��������
	assert(points.size() > 0);
	auto N = points.size()-1;
	// һ�ײ��h
	// h[i]=x[i]-x[i-1]
	std::vector<double> h;
	// ϵ��M
	std::vector<double> M;
	// �ֶκ���
	SegmentFunction<double, double> result(N);

	// lambda����ͨ;������һ�ײ��h
	auto calculateH = [&]() {
		h= std::vector<double>(N + 1);
		h[0] = 0;
		for (size_t i = 1; i < N + 1; i++)
		{
			h[i] = points[i].first - points[i - 1].first;
		}
	};
	// lambda������;������Mϵ��
	auto calculateM = [&]() {
		std::vector<double> vLambda(N);
		std::vector<double> vMu(N);
		std::vector<double> vConst(N + 1, 2);
		std::vector<double> d(N + 1);
		// �γ�vMu
		for (size_t i = 0; i < N-1; i++)
		{
			vMu[i] = h[i + 1] / (h[i + 1] + h[i + 2]);
		}
		vMu[N - 1] = 1;
		// �γ�vLambda
		for (size_t i = 1; i < N; i++)
		{
			vLambda[i] = h[i + 1] / (h[i] + h[i + 1]);
		}
		vLambda[0] = 1;
		// �γ�d
		for (size_t i = 1; i < N; i++)
		{
			d[i] = 6 / (h[i] + h[i + 1]);
			d[i] *= ((points[i + 1].second - points[i].second) / h[i + 1] - (points[i].second - points[i - 1].second) / h[i]);
		}
		d[0] = 6 * ((points[1].second - points[0].second) / h[1] - conditions.first) / h[1];
		d[N] = 6 * (conditions.second - (points[N].second - points[N - 1].second) / h[N]) / h[N];
		// ���ԽǾ���
		Sprase::SpraseTD<double> matrixA(std::move(vMu), std::move(vConst), std::move(vLambda));
		// ���M
		M = Sprase::CMSolver(matrixA, d);	
	};
	// lambda������;�����ɷֶκ���
	auto createSegFunc = [&]() {
		// ����ʽ������/һ/��/������ϵ��A0/A1/A2/A3
		double A0 = 0, A1 = 0, A2 = 0, A3 = 0;
		for (size_t i = 1; i < N + 1; i++)
		{
			// ����A0/A1/A2/A3
			A0 = (std::pow(points[i].first, 3)*M[i - 1] - std::pow(points[i - 1].first, 3)*M[i]) / h[i] / 6;
			A0 += (points[i - 1].second*points[i].first - points[i].second*points[i - 1].first) / h[i];
			A0 += h[i] * (M[i] * points[i - 1].first - M[i - 1] * points[i].first) / 6;
			A1 = (std::pow(points[i - 1].first, 2)*M[i] - std::pow(points[i].first, 2)*M[i - 1]) / h[i] / 2;
			A1 += (points[i].second - points[i - 1].second) / h[i];
			A1 += h[i] * (M[i - 1] - M[i]) / 6;
			A2 = (points[i].first*M[i - 1] - points[i - 1].first*M[i]) / h[i] / 2;
			A3 = (M[i] - M[i - 1]) / h[i] / 6;
			// ��ӷֶκ���
			auto func = [A0, A1, A2, A3](double x)->double {
				// return A0 + A1 * x + A2 * x*x + A3 * x*x*x;
				// ���ؾ����㷨��д
				return ((A3*x + A2)*x + A1)*x + A0;
			};
			// ��ӷֶκ�������
			std::string des = std::to_string(A0) + " + " + std::to_string(A1) + " * x + " + std::to_string(A2) + " * x^2 +" + std::to_string(A3) + " * x^3";
			// �������
			Interval<double> curInt({ points[i - 1].first, true }, { points[i].first, false });
			// ���ɸöηֶκ���
			result.AddSegment(i - 1, curInt, func, des);
		}
	};

	calculateH();
	calculateM();
	createSegFunc();
	return std::move(result);
}