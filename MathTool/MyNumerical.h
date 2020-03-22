#pragma once

#include <string>
#include <math.h>
#include <assert.h>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include "MySprase.h"

// author: Katyusha
// date: 2019/11/22

// 包含如下基础数据结构
// 1 区间Interval(TX 只允许是整型/浮点类型)
// 2 分段函数SegmentFunction(定义域连续)


// 包含以下数值算法
// 1 插值类
// 1-1 三次样条插值

// 2 最优逼近类

// 3 数值微积分类
// 3-1 龙贝格积分法


#ifndef MYNUMERICAL_H
#define MYNUMERICAL_H

namespace Utility
{
	namespace Numerical
	{
		// 基础数据结构：区间Interval(TX 只允许是整型/浮点类型)
		template<typename TX, std::enable_if_t<std::is_floating_point_v<TX> || std::is_integral_v<TX>, int> = 0 >
		class Interval
		{
		private:
			// 左端点(true表示闭区间，false表示开区间)
			std::pair<TX, bool> _begin;
			// 右端点(true表示闭区间，false表示开区间)
			std::pair<TX, bool> _end;
		public:
			// 返回区间左端点
			std::pair<TX, bool> Begin() { return _begin; }
			// 设置区间左端点
			void Begin(TX x1, bool isClosed) { _begin = { x1, isClosed }; }
			// 返回区间右端点
			std::pair<TX, bool> End() { return _end; }
			// 设置区间右端点
			void End(TX x2, bool isClosed) { _end = { x2, isClosed }; }
			// 判断是否在区间内
			bool Contain(TX x)
			{
				bool result = true;
				// 判断是否在左端点的右侧
				_begin.second ?
					result = result && x >= _begin.first : result = result && x > _begin.first;
				// 判断是否在右端点左侧
				_end.second ?
					result = result && x <= _end.first : result = result && x < _end.first;
				return result;
			}
			// 打印区间
			std::string Print()
			{
				std::string result;
				// 判断左端点开闭
				_begin.second ?
					result += "[" : result += "(";
				// 添加区间
				result += std::to_string(_begin.first) + ", " + std::to_string(_end.first);
				// 判断右端点开闭
				_end.second ?
					result += "]" : result += ")";
				return std::move(result);
			}
		public:
			Interval()
			{

			}
			Interval(TX x1, TX x2)
				: _begin(x1, false), _end(x2, false)
			{

			}
			Interval(std::pair<TX, bool> x1, std::pair<TX, bool> x2)
				: _begin(std::move(x1)), _end(std::move(x2))
			{

			}
			//Interval(std::initializer_list<std::pair<TX, bool> > list)
			//{
			//	
			//}
		};


		// 基础数据结构：分段函数SegmentFunction
		// 注意，不保证分段的顺序
		template<typename TX, typename TY>
		class SegmentFunction
		{
		private:
			// 存放分段定义域 [a, b]
			std::vector<Interval<TX> > _Segments;
			// 存放分段函数 Fun: [a, b] -> R
			std::vector<std::function<TY(TX)> > _Functions;
			// 存放分段函数的描述
			std::vector<std::string> _Description;
			// 总分段数
			unsigned int _NumSegment;
		public:
			// 返回时所在分段的Offset
			size_t GetSegment(TX x)
			{
				size_t result = 0;
				while (!_Segments[result].Contain(x))
				{
					result++;
				}
				return result;
			}
			// 添加新分段
			void AddSegment(Interval<TX> inter, std::function<TY(TX)> fun, std::string des="")
			{
				_Segments.emplace_back(std::move(inter));
				_Functions.emplace_back(std::move(fun));
				_Description.emplace_back(std::move(des));
				_NumSegment++;
			}
			void AddSegment(size_t offSet, Interval<TX> inter, std::function<TY(TX)> fun, std::string des = "")
			{
				if (offSet >= _NumSegment)
				{
					AddSegment(std::move(inter), std::move(fun), std::move(des));
				}
				else
				{
					_Segments[offSet] = std::move(inter);
					_Functions[offSet] = std::move(fun);
					_Description[offSet] = std::move(des);
				}
			}
			// 清除分段函数
			void Clear()
			{
				_Segments.clear();
				_Functions.clear();
				_Description.clear();
				_NumSegment = 0;
			}
			// 重载operator()
			// 警告，不检查定义域
			TY operator() (TX x)
			{
				return _Functions[GetSegment(x)](x);
			}
			// 打印分段函数
			std::string Print()
			{
				std::string result("Function:\n");
				for (size_t seg = 0; seg < _NumSegment; seg++)
				{
					result += "\t" + _Description[seg] + " : " + _Segments[seg].Print() + "\n";
				}
				return std::move(result);
			}
		public:
			// 构造函数
			SegmentFunction()
			{

			}
			SegmentFunction(int numSegment)
				: _NumSegment(numSegment > 0 ? numSegment : 0)
			{
				_Segments = std::vector<Interval<TX> >(_NumSegment);
				_Functions = std::vector<std::function<TY(TX)> >(_NumSegment);
				_Description = std::vector<std::string>(_NumSegment);
			}
			// 禁止列表初始化
			SegmentFunction(std::initializer_list<std::pair<Interval<TX>, std::function<TY(TX)> > >) = delete;
		};

		// 三次样条插值法
		// 参考《计算方法》(李乃成、梅立泉著，科学出版社)
		SegmentFunction<double, double> CubicSplineInterpolation(std::vector<std::pair<double, double> > points, std::pair<double, double> conditions);

		// 龙贝格积分法收敛精度
		constexpr double RombergConvergeLimit = 1e-12;

		// 龙贝格积分法
		// 参考《计算方法》(李乃成、梅立泉著，科学出版社)
		double RombergIntegration(const double begin, const double end, std::function<double(double)> fun, const double converge = RombergConvergeLimit);
	}
}



#endif
