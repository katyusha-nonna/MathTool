#include "MyOptimalize.h"

int main()
{
	// ����DFP�㷨

	// ����Ŀ�꺯��
	auto f = [](std::vector<double>& x)->double {
		return x[0] * (x[0] - x[1] + 2) + x[1] * (x[1] - 4);
	};
	// �����ݶȺ���
	auto gf = [](std::vector<double>& x)->std::vector<double> {
		std::vector<double> result(2);
		result[0] = 2 * x[0] - x[1] + 2;
		result[1] = 2 * x[1] - x[0] - 4;
		return std::move(result);
	};

	//// ����Ŀ�꺯��
	//auto f = [](std::vector<double>& x)->double {
	//	return x[0] * x[0] + 4 * x[1] * x[1];
	//};
	//// �����ݶȺ���
	//auto gf = [](std::vector<double>& x)->std::vector<double> {
	//	std::vector<double> result(2);
	//	result[0] = 2 * x[0];
	//	result[1] = 8 * x[1];
	//	return std::move(result);
	//};

	std::vector<double> x{2, 2};
	// Ŀ�꺯��
	std::cout << "Object: " << Utility::Optimalize::DFPKernel(x, f, gf) << std::endl;
	// ���Ž�
	std::cout << "x: ";
	for (auto& e : x)
	{
		std::cout << "  " << e;
	}
	std::cout << std::endl;
	system("pause");
	return 0;
}