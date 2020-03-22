#pragma once

#include <functional>
#include "MySprase.h"

#ifndef MYOPYIMALIZE_H
#define MYOPYIMALIZE_H

namespace Utility
{
	namespace Optimalize
	{
		// ��ţ�ٷ�֮һ-DFP�㷨�ĺ��ĺ���
		double DFPKernel(std::vector<double>& x, std::function<double(std::vector<double>&)> f, std::function<std::vector<double>(std::vector<double>&)> gf);
	}
}

#endif


