#include "MyOptimalize.h"

using namespace Utility;

double Optimalize::DFPKernel(std::vector<double>& x, std::function<double(std::vector<double>&)> f, std::function<std::vector<double>(std::vector<double>&)> gf)
{
	// x��Ϊ�����ĳ�ֵ��fΪĿ�꺯����gfΪĿ�꺯���ݶȺ���
	// ���н�����xΪ���Ž⣬����ֵΪ�Ż����

	// ���ﶨ����һЩ��Ҫ�Ĳ���
	
	// ά��
	int N = x.size();
	// ����������
	int maxIterCount = 500;
	// ��󾫶�
	double maxError = 1e-12;
	// ����rho
	double rho = 0.5;
	// ����sigma
	double sigma = 0.2;
	// ��ǰ��������
	int count = 0;
	// ��ʱ����-��ǰ�ݶ�
	std::vector<double> Gk(N);
	// ��ʱ����-��ǰ��ɭ����
	Sprase::Matrix<double> Hk(N, N);
	// ��ʱ����-��ǰ��������
	std::vector<double> dk(N);
	// ��ʱ����-�ݶȲ�
	std::vector<double> yk(N);
	// ��ʱ����-x������
	std::vector<double> sk(N);

	// lamdba������;������˷�(C=A*B)
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
	// lambda������;�������˷�(C=a*b)
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
	// lambda������;�������˷�(c=a.*b)
	auto vectorDotOp = [&](std::vector<double>& a, std::vector<double>& b, double& c) {
		assert(a.size() == b.size());
		for (int i = 0; i < a.size(); i++)
		{
			c+= a[i] * b[i];
		}
	};
	// lambda������;����������˷�(c=A*b)
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

	// ��ʼ��
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
	// ��ʼ����
	while (count < maxIterCount)
	{
		// �����ݶ�
		Gk = gf(x);
		// ������������
		dk = std::vector<double>(N);
		matVectorTimeOp(Hk, Gk, dk, true);
		// ����mk��Armijo׼����������
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
		// ִ��DFPУ��
		// ����x��sk
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
		// ����Gk��yk
		yk = std::move(Gk);
		Gk = gf(x);
		for (int i = 0; i < N; i++)
		{
			yk[i] = Gk[i] - yk[i];
		}
		// ����Hk
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