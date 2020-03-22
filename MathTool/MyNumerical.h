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

// �������»������ݽṹ
// 1 ����Interval(TX ֻ����������/��������)
// 2 �ֶκ���SegmentFunction(����������)


// ����������ֵ�㷨
// 1 ��ֵ��
// 1-1 ����������ֵ

// 2 ���űƽ���

// 3 ��ֵ΢������
// 3-1 ��������ַ�


#ifndef MYNUMERICAL_H
#define MYNUMERICAL_H

namespace Utility
{
	namespace Numerical
	{
		// �������ݽṹ������Interval(TX ֻ����������/��������)
		template<typename TX, std::enable_if_t<std::is_floating_point_v<TX> || std::is_integral_v<TX>, int> = 0 >
		class Interval
		{
		private:
			// ��˵�(true��ʾ�����䣬false��ʾ������)
			std::pair<TX, bool> _begin;
			// �Ҷ˵�(true��ʾ�����䣬false��ʾ������)
			std::pair<TX, bool> _end;
		public:
			// ����������˵�
			std::pair<TX, bool> Begin() { return _begin; }
			// ����������˵�
			void Begin(TX x1, bool isClosed) { _begin = { x1, isClosed }; }
			// ���������Ҷ˵�
			std::pair<TX, bool> End() { return _end; }
			// ���������Ҷ˵�
			void End(TX x2, bool isClosed) { _end = { x2, isClosed }; }
			// �ж��Ƿ���������
			bool Contain(TX x)
			{
				bool result = true;
				// �ж��Ƿ�����˵���Ҳ�
				_begin.second ?
					result = result && x >= _begin.first : result = result && x > _begin.first;
				// �ж��Ƿ����Ҷ˵����
				_end.second ?
					result = result && x <= _end.first : result = result && x < _end.first;
				return result;
			}
			// ��ӡ����
			std::string Print()
			{
				std::string result;
				// �ж���˵㿪��
				_begin.second ?
					result += "[" : result += "(";
				// �������
				result += std::to_string(_begin.first) + ", " + std::to_string(_end.first);
				// �ж��Ҷ˵㿪��
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


		// �������ݽṹ���ֶκ���SegmentFunction
		// ע�⣬����֤�ֶε�˳��
		template<typename TX, typename TY>
		class SegmentFunction
		{
		private:
			// ��ŷֶζ����� [a, b]
			std::vector<Interval<TX> > _Segments;
			// ��ŷֶκ��� Fun: [a, b] -> R
			std::vector<std::function<TY(TX)> > _Functions;
			// ��ŷֶκ���������
			std::vector<std::string> _Description;
			// �ֶܷ���
			unsigned int _NumSegment;
		public:
			// ����ʱ���ڷֶε�Offset
			size_t GetSegment(TX x)
			{
				size_t result = 0;
				while (!_Segments[result].Contain(x))
				{
					result++;
				}
				return result;
			}
			// ����·ֶ�
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
			// ����ֶκ���
			void Clear()
			{
				_Segments.clear();
				_Functions.clear();
				_Description.clear();
				_NumSegment = 0;
			}
			// ����operator()
			// ���棬����鶨����
			TY operator() (TX x)
			{
				return _Functions[GetSegment(x)](x);
			}
			// ��ӡ�ֶκ���
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
			// ���캯��
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
			// ��ֹ�б��ʼ��
			SegmentFunction(std::initializer_list<std::pair<Interval<TX>, std::function<TY(TX)> > >) = delete;
		};

		// ����������ֵ��
		// �ο������㷽����(���˳ɡ�÷��Ȫ������ѧ������)
		SegmentFunction<double, double> CubicSplineInterpolation(std::vector<std::pair<double, double> > points, std::pair<double, double> conditions);

		// ��������ַ���������
		constexpr double RombergConvergeLimit = 1e-12;

		// ��������ַ�
		// �ο������㷽����(���˳ɡ�÷��Ȫ������ѧ������)
		double RombergIntegration(const double begin, const double end, std::function<double(double)> fun, const double converge = RombergConvergeLimit);
	}
}



#endif
