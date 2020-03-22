#pragma once

#include <functional>
#include "MySprase.h"

#ifndef MYOPYIMALIZE_H
#define MYOPYIMALIZE_H

namespace Utility
{
	namespace Optimalize
	{
		// 拟牛顿法之一-DFP算法的核心函数
		double DFPKernel(std::vector<double>& x, std::function<double(std::vector<double>&)> f, std::function<std::vector<double>(std::vector<double>&)> gf);
	}
}

#endif


