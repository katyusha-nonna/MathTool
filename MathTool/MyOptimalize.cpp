/*! \file MyOptimalize.cpp */

#include "MyOptimalize.h"
//���������
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
	�ϰ��ڵ㷨�ĺ��ĺ���
	\li ������ⲻ��ʽԼ�������Ż�����

	\f[
	\begin{aligned}
	\min_{x,y}\quad&c^{T}x\\
	s.t.\quad&A^{T}x\geq b
	\end{aligned}
	\f]

	\param c Ŀ�꺯��ϵ������(Eigen����)
	\param A Լ������ϵ������(Eigen����)
	\param b Լ��������������(Eigen����)
	\param x ���߱����Ż��������(Eigen����)

	\return Ŀ�꺯���Ż����

	\todo ����
*/
double LPBarrierInteriorPointKernel(
	const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, double t, double mu, double error)
{
	// Min: c'*x
	// s.t A'*x >= b
	assert(A.cols() == b.rows());
	assert(A.rows() == x.rows());
	// ������Ҫ��ά��(������Ŀ/Լ����Ŀ)
	const auto n = A.rows();
	const auto m = A.cols();
	// ��ʱ����-�ݶ�����
	Eigen::VectorXd gradV;
	// ��ʱ����-��ɭ����
	Eigen::MatrixXd hessianM;
	// ��ʱ����-���߱�������
	Eigen::VectorXd delta;
	// ��ǰ������
	double curGlobalError = 1.0;
	// lambda������;����ʼ��
	auto init = [&]() {
		// ���ﲻ��x���м�飬ֱ����Ϊ���뺯����x���ȵ���n�Ҿ��к��ʵĳ�ֵ
	};
	// lambda������;��ţ�ٵ����������Է�����
	auto newton = [&]() {
		int maxCount = 8, curCount = 0;
		double maxError = 1e-4, curError = 1.0;
		while (curError > maxError && curCount < maxCount)
		{
			// �����ݶ�����
			gradV = c * t - A * (A.transpose()*x - b).cwiseInverse();
			// ���º�ɭ����
			hessianM = A * (A.transpose()*x - b).cwiseAbs2().cwiseInverse().asDiagonal()*A.transpose();
			// ���µ�ǰ��
			delta = hessianM.colPivHouseholderQr().solve(gradV);
			x -= delta;
			// ���¼�������������
			curError = delta.norm() / x.norm();
			curCount++;
		}
	};

	init();
	// ��ʼ����
	while (curGlobalError > error)
	{
		newton();
		curGlobalError = m / t;
		t *= mu;
	}
	// �������Ž�
	auto f = c.transpose()*x;
	return f.value();
}

/*!
	ԭʼ��ż�ڵ㷨�ĺ��ĺ���
	\li ��������ʽԼ�������Ż�����

	\f[
	\begin{aligned}
	\min_{x,y}\quad&c^{T}x\\
	s.t.\quad&A^{T}x=b
	\end{aligned}
	\f]

	\param c Ŀ�꺯��ϵ������(Eigen����)
	\param A Լ������ϵ������(Eigen����)
	\param b Լ��������������(Eigen����)
	\param x ���߱����Ż��������(Eigen����)
	\param alpha ţ�ٷ�����(0<alpha<1����alpha<=0���Զ�ѡ�񲽳�)
	\param sigma ţ�ٷ��½�����(0<=sigma<1����sigma<0���Զ�ѡ���½�����)
	\param error ţ�ٷ����������(ԭ����-��ż����Ķ�ż��϶)

	\return Ŀ�꺯���Ż����(���޽�ʱ����INFITY)

	\todo ����
*/
double LPPrimalDualInteriorPointKernel(
	const Eigen::VectorXd& c, const Eigen::MatrixXd& A, const Eigen::VectorXd& b, Eigen::VectorXd& x, double alpha, double sigma, double error)
{
	// Min: c'*x
	// s.t A'*x = b
	assert(A.cols() == b.rows());
	assert(A.rows() == x.rows());
	// ������Ҫ��ά��(������Ŀ/Լ����Ŀ)
	const auto n = A.rows();
	const auto m = A.cols();
	// һЩ��־λ
	const bool autoAlpha = alpha <= 0;
	const bool autoSigma = sigma < 0;
	bool isUnbounded = false;
	// ��ʱ����-���߱��� var=[x; lambda; s]
	Eigen::VectorXd var;
	// ��ʱ����-�ſ˱Ⱦ���(ϡ��)
	Eigen::SparseMatrix<double, Eigen::ColMajor> jacobiM(2 * n + m, 2 * n + m);
	// ��ʱ����-KKT�������Ҷ�����
	Eigen::VectorXd right;
	// ��ʱ�����������������(��ֵ) delta=[xDelta; lambdaDelta; sDelta]
	Eigen::VectorXd delta;
	// ��ǰ�Ķ�ż��϶
	double mu;
	// lambda������;����ʼ��
	auto init = [&]() {
		// ��ʼ�����߱���var
		var.setOnes(2 * n + m);
		var.head(n) = x;
		// ��ʼ���ſɱȾ���
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
		// ��ʼ���Ҷ�����
		right.setZero(2 * n + m);
		// ��ʼ����������
		delta.setZero(2 * n + m);
		// ��ʼ����ż��϶
		mu = xx.dot(s) / n;
	};
	// lambda������;�����¾��߱���
	auto updataVar = [&]() {
		// �ж�alpha�������²���
		if (autoAlpha)
		{
			// �Զ������������¾��߱���
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
	// lambda������;�������½�����
	auto updateSigma = [&]() {
		if (autoSigma)
		{
			// ��ʱ�̶�Ϊ1.0
			sigma = 1.0;
		}
	};
	// lambda������;�����¶�ż��϶
	auto checkDualMeasure = [&]() {
		// mu=x^{T}s/n
		auto xx = Eigen::VectorXd::Map(var.data(), n);
		auto s = Eigen::VectorXd::Map(var.data() + n + m, n);
		mu = xx.dot(s) / n;
	};
	// lambda������;���ж��Ƿ��޽�
	auto checkUnbounded = [&]() {
		auto xNorm = Eigen::VectorXd::Map(var.data(), n).norm();
		auto lNorm = Eigen::VectorXd::Map(var.data() + n, m).norm();
		// �ж�x��lambda�Ķ������Ƿ���NAN����INFITY(�ж��޽�)
		return isinf(xNorm) || isinf(lNorm) || isnan(xNorm) || isnan(lNorm);
	};
	// lambda������;��ţ�ٷ�����KKT�����������ԭʼ-��ż����
	auto newton = [&]() {
		// �ſ˱Ⱦ���Ϊһ��ϡ�����ÿ�ε���ʱ�������½Ǻ����½�����
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
			// �����½�����
			updateSigma();
			// �����ſ˱Ⱦ���
			for (int k = 0; k < n; k++)
				jacobiM.coeffRef(n + m + k, k) = s.coeff(k);
			for (int k = 0; k < n; k++)
				jacobiM.coeffRef(n + m + k, n + m + k) = xx.coeff(k);
			// �����Ҷ�����
			rx = A * lambda + s - c;
			rl = A.transpose()*xx - b;
			tempProd = mu * sigma;
			rs = xx.cwiseProduct(s);
			rs.array() -= tempProd;
			// ���ϡ�����
			solver.compute(jacobiM);
			delta = solver.solve(right);
			// ���¾��߱���
			updataVar();
			// ���¶�ż��϶
			checkDualMeasure();
			// �ж��Ƿ��޽�
			if (checkUnbounded())
			{
				isUnbounded = true;
				break;
			}
			// ���¼�����
			curError = delta.norm();
			curCount++;
		}
	};

	// ���
	init();
	newton();
	// �������Ž�
	if (isUnbounded)
		return INFINITY;
	x = var.head(n);
	auto f = c.transpose()*x;
	return f.value();
}

/*!
	�����ݶȷ��ĺ��ĺ���
	\li ֻ���������G���������ԳƵ���Լ�������Ż�����

	\f[
	\begin{aligned}
	\min_{x,y}\quad&\frac{1}{2}x^{T}Gx+p^{T}x\\
	s.t.\quad&\quad
	\end{aligned}
	\f]

	\param G Ŀ�꺯��������������ϵ������(Eigen����)(�Գ�����)
	\param p Ŀ�꺯��һ����ϵ������(Eigen����)
	\param x ���߱����Ż��������(Eigen����)

	\return Ŀ�꺯���Ż����

	\todo ����
*/
double QPConjugateGradientKernel(const Eigen::MatrixXd& G, const Eigen::VectorXd& p, Eigen::VectorXd& x)
{
	// Min: (1 / 2)*x'*G*x + p'*x
	// s.t NULL
	// �����ݶ������(ϡ��)
	Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Upper> solver;
	// �������
	x = solver.compute(G.sparseView()).solve(-p);
	// �������Ž�
	auto f = p.transpose()*x + 0.5*x.transpose()*G*x;
	return f.value();
}

/*!
	��ЧԼ�������ĺ��ĺ���
	\li ������ⲻ��ʽԼ�������Ż�����

	\f[
	\begin{aligned}
	\min_{x,y}\quad&\frac{1}{2}x^{T}Gx+p^{T}x\\
	s.t.\quad&A^{T}x\geq b
	\end{aligned}
	\f]

	\param G Ŀ�꺯��������������ϵ������(Eigen����)
	\param A Լ������ϵ������(Eigen����)
	\param p Ŀ�꺯��һ����ϵ������(Eigen����)
	\param b Լ��������������(Eigen����)
	\param x ���߱����Ż��������(Eigen����)

	\return Ŀ�꺯���Ż����

	\todo ����
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
	// ������Ҫά��
	const auto n = G.rows();
	const auto m = A.cols();
	// ��ǰ��Լ����
	std::set<int> S;
	// ��ǰ������Ľ�
	Eigen::VectorXd tempx;
	Eigen::VectorXd delta, lambda;
	// ��ǰ���²���
	double alpha;
	std::map<int, double> candidate;
	// ��ǰ��p��A��b
	Eigen::MatrixXd tempA;
	Eigen::VectorXd tempp, tempb;
	// �˳���־
	bool quit = false;
	// ����������
	int iterCount = 0;
	// lambda������;����ʼ��Լ����
	auto init = [&]() {
		S.clear();
		Eigen::VectorXd temp = A.transpose()*x - b;
		for (int i = 0; i < m; i++)
		{
			if (temp[i] >= 0 && std::abs(temp[i]) < zeroAct)
				S.emplace(i);
		}
	};
	// lambda������;�����ݵ�ǰԼ�������µ�ǰ��p��A��b
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
	// lambda������;���жϵ�ǰx+delta�Ƿ�Ϊԭ������н�
	auto isFeasible = [&]() {
		Eigen::VectorXd temp = A.transpose()*(x + delta) - b;
		double minVal = temp.minCoeff();
		return minVal >= 0;
	};
	// lambda������;������alpha��S
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

	// ��ʼ��
	init();
	// ��ʼ����
	while (!quit)
	{
		updateSubP();
		QPLagrangeKernel(G, tempA, tempp, tempb, tempx);
		delta = tempx.head(n);
		lambda = tempx.tail(tempx.rows() - n);
		if (delta.norm() < zeroAct)
		{
			// ˵��deletaΪ0
			// ����S
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
			// ˵��deleta��Ϊ0
			// ����alpha
			updateAlpha();
			// ����x
			x += alpha * delta;
		}
	}
	// �������Ž�
	auto f = p.transpose()*x + 0.5*x.transpose()*G*x;
	return f.value();
}

/*!
	�������շ��ĺ��ĺ���
	\li ��������ʽԼ�������Ż�����

	\f[
	\begin{aligned}
	\min_{x,y}\quad&\frac{1}{2}x^{T}Gx+p^{T}x\\
	s.t.\quad&A^{T}x=b
	\end{aligned}
	\f]

	\param G Ŀ�꺯��������������ϵ������(Eigen����)
	\param A Լ������ϵ������(Eigen����)
	\param p Ŀ�꺯��һ����ϵ������(Eigen����)
	\param b Լ��������������(Eigen����)
	\param x ���߱����Ż��������(Eigen����)

	\return Ŀ�꺯���Ż����

	\todo ����
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
	// ��ʼ���������վ���K
	Eigen::MatrixXd K(n + m, n + m);
	K.setZero();
	K.topLeftCorner(n, n) = G;
	K.topRightCorner(n, m) = -A;
	K.bottomLeftCorner(m, n) = -A.transpose();
	// ��ʼ���Ҷ�����d
	Eigen::VectorXd d(n + m);
	d.head(n) = -p;
	d.tail(m) = -b;
	// �ⷽ��
	Eigen::SparseQR<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int> > solver;
	solver.compute(K.sparseView());
	K.resize(0, 0);
	x = solver.solve(d);
	//x = K.colPivHouseholderQr().solve(d);
	// �������Ž�
	auto f = p.transpose()*x + 0.5*x.transpose()*G*x;
	return f.value();
}