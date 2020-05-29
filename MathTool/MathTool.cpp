#include "MyOptimalize.h"
#include "MyNumerical.h"
#include <Eigen/Eigen>

int main()
{
	// 测试共轭梯度求解器


	/*
	// 测试拉格朗日法
	Eigen::MatrixXd G(2, 2);
	G << 2, -1, -1, 2;
	Eigen::MatrixXd A(2, 1);
	A << 0, 1;
	Eigen::VectorXd p(2);
	p << 0, -1.5;
	Eigen::VectorXd b(1);
	b << 0;
	Eigen::VectorXd x;
	std::cout << "op: " << LagrangeKernel(G, A, p, b, x) << std::endl;
	std::cout << "x: " << x;
	*/

	/*
	// 测试ACM
	Eigen::MatrixXd G(2, 2);
	G << 2, -1, -1, 2;
	Eigen::MatrixXd A(2, 3);
	A << -1, 1, 0, -1, 0, 1;
	Eigen::VectorXd p(2);
	p << -3, 0;
	Eigen::VectorXd b(3);
	b << -2, 0, 0;
	Eigen::VectorXd x(2);
	x << 0, 0;
	std::cout << "op: " << ACMKernel(G, A, p, b, x) << std::endl;
	std::cout << "x: " << x;
	*/


	/*
	// 测试map
	std::vector<double> y(10);
	for (int i = 0; i < y.size(); i++)
		y[i] = i;
	Eigen::Map<Eigen::VectorXd, Eigen::Unaligned, Eigen::Stride<Eigen::Dynamic, 1> >
		eb(y.data(), y.size(), 1, Eigen::Stride<Eigen::Dynamic, 1>(y.size(), 1));
	std::cout << eb;
	for (int i = 0; i < y.size(); i++)
		y[i] = i * i;
	std::cout << eb.transpose();
	*/

	/*
	// 测试Map
	int n = 5, m = 5;
	std::vector<double> A(n*m);
	for (int i = 0; i < n*m; i++)
		A[i] = i;
	Eigen::Map<Eigen::Matrix<double, -1, -1, Eigen::RowMajor>, 0, Eigen::Stride<-1, -1> >
		eA(A.data(), n, m, Eigen::Stride<-1, -1>(m, 1));
	std::vector<double> tempOne(m, 1.0);
	std::cout << eA << std::endl;
	for (int i = 0; i < n*m; i++)
		A[i] = i*i;
	std::cout << eA << std::endl;
	*/

	/*
	// 测试最小二乘
	int n = 4, m = 7;
	std::vector<double> A{ 1,1,2,2,0,4,2,
	1,2,1,3,3,4,0, 
	1,2,2,3,0,4,3, 
	1,1,1,1,1,1,1};

	std::vector<double> b{ 4,8,7,12,7,18,6 };
	auto result = Utility::Numerical::MutliLinearLeastSquareFitting(A, b, n, m);
	for (auto& s : result)
		std::cout << s << "    ";
	*/

	/*
	// 测试障碍内点法
	Eigen::MatrixXd A(2, 4);
	A << -1, -2, 1, 0, -2, -1, 0, 1;
	//std::cout << "A: " << A << std::endl;
	Eigen::VectorXd c(2);
	c << 1, 1;
	Eigen::VectorXd b(4);
	b << -1, -1, 0, 0;
	Eigen::VectorXd x(2);
	x << 0.3, 0.3;
	std::cout << "op: " << BarrierInteriorPointKernel(c, A, b, x, 1, 1.2, 1e-8) << std::endl;
	std::cout << "x: " << x;
	*/

	/*
	// 测试原始对偶内点法
	// Min x1+x2
	// s.t. x1+2*x2<=1
	//      2*x1+x2<=1
	//      2*x1+4*x2>=1
	Eigen::MatrixXd A(5, 3);
	A << 1, 2, 2, 2, 1, 4, 1, 0, 0, 0, 1, 0, 0, 0, 1;
	//std::cout << "A: " << A << std::endl;
	Eigen::VectorXd c(5);
	c << 1, 1, 0, 0, 0;
	Eigen::VectorXd b(3);
	b << 1, 1, 1;
	Eigen::VectorXd x(5);
	x << 0.2, 0.2, 0, 0, 0;
	std::cout << "op: " << LPPrimalDualInteriorPointKernel(c, A, b, x, -1, -1, 1e-8) << std::endl;
	std::cout << "x: " << x;
	*/

	/*
	// 测试三种向量范数 
	Eigen::VectorXd c(5);
	c.setZero();
	c.coeffRef(0) = NAN;
	c.coeffRef(1) = INFINITY;
	c.coeffRef(2) = 3.1415;
	c.coeffRef(3) = -2.71828;
	c.coeffRef(4) = 0;
	std::cout << "L1:  " << (isinf(c.lpNorm<1>()) || isnan(c.lpNorm<1>())) << std::endl;
	std::cout << "L2:  " << (isinf(c.norm()) || isnan(c.norm())) << std::endl;
	std::cout << "Infty:  " << (isinf(c.lpNorm<Eigen::Infinity>()) || isnan(c.lpNorm<Eigen::Infinity>())) << std::endl;
	*/


	system("pause");
	return 0;
}