#include "MySprase.h"

using namespace Utility;

std::vector<double> Sprase::GEKernel(Matrix<double>& A, std::vector<double>& b)
{
	auto size = A.GetSize();
	assert(size.first == size.second);
	std::vector<double> x(size.first);

	// ѡ����Ԫ����
	auto selectPivot = [&](TIndex k) {
		// std::cout << k << " " << std::endl;
		auto max = 0.0;
		auto q = 0;
		// ѡ����Ԫ
		for (TIndex i = k; i < size.first; i++)
		{
			if (std::abs(A[i][k]) > max)
			{
				max = A[i][k];
				q = i;
			}
		}
		assert(std::abs(max) > mini_value);
		// ������
		if (q != k)
		{
			std::swap(A[q], A[k]);
			std::swap(b[q], b[k]);
		}
	};
	// ��ȥ����
	auto elimination = [&]() {
		for (TIndex k = 0; k < size.first-1; k++)
		{
			// ѡ����Ԫ
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
	// �ش�����
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

	// ���
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

	// ѡ��Ԫ����
	auto selectPivot = [&](TIndex k) {
		auto max = 0.0;
		auto q = 0;
		// ѡ����Ԫ
		for (TIndex i = k; i < size.first; i++)
		{
			if (std::abs(A[i][k]) > max)
			{
				max = A[i][k];
				q = i;
			}
		}
		assert(std::abs(max) > mini_value);
		// ������
		if (q != k)
		{
			A[q].swap(A[k]);
			std::swap(pi[q], pi[k]);
		}
	};
	// �ֽ⺯��
	auto decomposition = [&]() {
		// ����߽�
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

	// �����滻����
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
	// �����滻����
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

	// ���
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
	// ��k��x�Ľ��
	std::vector<double> x0(size.first);
	// ��k+1��x�Ľ��
	std::vector<double> x(xPre);
	// �в�
	double residual = 0;
	// ��������
	int count = 0;
	// k�ε�������
	auto forEachStep = [&]()
	{
		int i = 0, j = 0;
		// k+1�ν��������ۼ�/k�ν��������ۼ�
		double xSum = 0.0, x0Sum = 0.0;
		double pivot = 0;
		// �޸�x0
		x0 = x;
		for (int row = 0; row < size.first; row++)
		{
			// ѡ��i��
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
			// �޸�x
			x[row] = (1 - weight) * x0[row] + weight * (b[row] - xSum - x0Sum) / pivot;
			// �ۼ�ֵ����
			x0Sum = 0;
			xSum = 0;
		}
	};
	// �ж�����
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

	// ִ�����
	do
	{
		forEachStep();
	} while (!checkConverge());
	// ���ؽ��
	return std::move(x);
}

std::vector<double> Sprase::JACKernel(SpraseCSR<double>&A, std::vector<double>& b)
{
	auto size = A.GetSize();
	assert(size.first == size.second);
	// ��k��x�Ľ��
	std::vector<double> x0(size.first);
	// ��k+1��x�Ľ��
	std::vector<double> x(size.first);
	// �в�
	double residual = 0;
	// ��������
	int count = 0;
	// k�ε�������
	auto forEachStep = [&]()
	{
		int i = 0, j = 0;
		// k+1�ν��������ۼ�/k�ν��������ۼ�
		double x0Sum = 0.0;
		double pivot = 0;
		// �޸�x0
		x0 = x;
		for (int row = 0; row < size.first; row++)
		{
			// ѡ��i��
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
			// �޸�x
			x[row] = (b[row] - x0Sum) / pivot;
			// �ۼ�ֵ����
			x0Sum = 0;
		}
	};
	// �ж�����
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

	// ִ�����
	do
	{
		forEachStep();
	} while (!checkConverge());
	// ���ؽ��
	return std::move(x);
}

std::pair<std::vector<double>, std::vector<double> > Sprase::QRDecomposition(Matrix<double>&A, bool isExtension)
{
	// A��ά��
	auto size = A.GetSize();
	// ��������
	int m = size.first;
	int n = isExtension ? size.second - 1 : size.second;
	// ��ʼ���Խ�����d������\alpha
	std::vector<double> d(n);
	std::vector<double> alpha(n);
	// ��ʼ�ֽ�
	for (int k = 0; k < n - 1; k++)
	{
		double sigma = 0.0;
		for (int i = k; i < m ; i++)
		{
			sigma += A[i][k] * A[i][k];
		}
		sigma = std::sqrt(sigma);
		// sigmaȡ��
		if (A[k][k] >= 0)
		{
			sigma = -sigma;
		}
		d[k] = sigma;
		alpha[k] = sigma * (sigma - A[k][k]);
		// ��Akk�޸�Ϊ���Wk�ĵ�һ��Ԫ��
		A[k][k] = A[k][k] - sigma;
		// H�ֽ�
		// ����ʱ�ֽ�����������
		auto end = isExtension ? n + 1 : n;
		for (int j = k + 1; j < end; j++)
		{
			// ����Belta
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
	// ����Ϊnʱ����������
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
		// ����ʱ�ֽ�����������
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

	// �ֱ𷵻ضԽ�������������飬���ڹ���H=QT
	return { std::move(d), std::move(alpha) };
}


std::vector<double> Sprase::QRKernel(Matrix<double>&A, std::vector<double>& b)
{
	// A��ά��
	auto size = A.GetSize();
	int m = size.first, n = size.second;
	std::vector<double> d;
	std::vector<double> x(n);

	// lambda������;������������
	auto exten = [&]() {
		for (int i = 0; i < m; i++)
		{
			A[i].push_back(b[i]);
		}
	};
	// lambda������;��ִ��QR�ֽ�
	auto decomposition = [&]() {
		d = QRDecomposition(A, true).first;
	};
	// lambda����ͨ;��ִ�лش�
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

	// ���
	exten();
	decomposition();
	backSolve();
	//  ���ؽ��
	return std::move(x);
}

std::vector<double> Sprase::CGKernel(Matrix<double>& A, std::vector<double>& b)
{
	// A��ά��
	auto size = A.GetSize();
	assert(size.first == size.second);
	// ͬʱҲ�ǵ���������������
	auto N = size.first;
	// �������x
	std::vector<double> x(N);
	// ����alpha
	double alpha = 0.0;
	// ��ʱ����belta
	double belta;
	// ��������d
	std::vector<double> d(N);
	// �в�����r(��k+1��)
	std::vector<double> r(N);
	// �в������ڻ� r0Dot(��k��)
	double r0Dot = 0.0;
	// ��ǰ��������
	int count = 0;

	// lambda����ͨ;�����������ڻ�
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
	// lambda������;�����������������ĳ˻�
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
	// lambda������;�����ɲв�
	auto updateR = [&]() {
		if (count == 0)
		{
			r0Dot = 0.0;
			// ���ɳ��βв˳������d
			for (size_t i = 0; i < N; i++)
			{
				d[i] = b[i] - dotProduct(A[i], x);
				r0Dot += d[i] * d[i];
			}
		}
		else
		{
			// ����k�βв�
			for (size_t i = 0; i < N; i++)
			{
				r[i]=b[i]- dotProduct(A[i], x);
			}
		}
	};
	// lambda������;������x
	auto updateX = [&]() {
		for (size_t i = 0; i < N; i++)
		{
			x[i] = x[i] + alpha * d[i];
		}
	};
	// lambda������;������d
	auto updateD = [&]() {
		for (size_t i = 0; i < N; i++)
		{
			d[i] = r[i] + belta * d[i];
		}
	};
	// lambda������;���ж�����
	auto checkout = [&]()
	{
		double error = 0.0;
		for (size_t i = 0; i < N; i++)
		{
			error += r[i] * r[i];
		}
		error = std::sqrt(error);
		// ��������Ϊ��������������N���ߴ�ʱ���ȴﵽ��������Ҫ��
		return (count > N || error < CGConvergeLimit) ? true : false;
	};

	// ��ʼ����
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
		// ����r0
		r0Dot = tempRR;
	}
	// ���ؽ��
	return std::move(x);
}

std::vector<double> Sprase::CMKernel(SpraseTD<double>& A, std::vector<double>& d)
{
	// A��ά��
	auto size = A.GetSize();
	assert(size.first == size.second);
	auto N = size.first;
	// ���L��U����ʱ����
	std::vector<double> L(N), U(N);
	// ���x��y����ʱ����
	std::vector<double> x(N), y(N);

	// lambda������;��׷�Ϲ���
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
	// lambda������;���ش�����
	auto back = [&]() {
		x[N - 1] = y[N - 1] / U[N - 1];
		for (int i = N - 2; i > -1; i--)
		{
			x[i] = (y[i] - A.AtDiagCRef(i)*x[i + 1]) / U[i];
		}
	};

	// ��ʼ����
	chasing();
	back();
	// ���ؽ��
	return std::move(x);	
}