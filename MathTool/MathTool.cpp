#include "MyOptimalize.h"

int main()
{
	// 测试DFP算法

	// 定义目标函数
	auto f = [](std::vector<double>& x)->double {
		return x[0] * (x[0] - x[1] + 2) + x[1] * (x[1] - 4);
	};
	// 定义梯度函数
	auto gf = [](std::vector<double>& x)->std::vector<double> {
		std::vector<double> result(2);
		result[0] = 2 * x[0] - x[1] + 2;
		result[1] = 2 * x[1] - x[0] - 4;
		return std::move(result);
	};

	//// 定义目标函数
	//auto f = [](std::vector<double>& x)->double {
	//	return x[0] * x[0] + 4 * x[1] * x[1];
	//};
	//// 定义梯度函数
	//auto gf = [](std::vector<double>& x)->std::vector<double> {
	//	std::vector<double> result(2);
	//	result[0] = 2 * x[0];
	//	result[1] = 8 * x[1];
	//	return std::move(result);
	//};

	std::vector<double> x{2, 2};
	// 目标函数
	std::cout << "Object: " << Utility::Optimalize::DFPKernel(x, f, gf) << std::endl;
	// 最优解
	std::cout << "x: ";
	for (auto& e : x)
	{
		std::cout << "  " << e;
	}
	std::cout << std::endl;
	system("pause");
	return 0;
}