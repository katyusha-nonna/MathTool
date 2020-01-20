#pragma once

#include <assert.h>
#include <vector>
#include <type_traits>
#include <utility>
#include <iostream>
#include <list>
#include <algorithm>
#include <complex>

// author: Katyusha
// date: 2019/11/09

// 包含如下基础数据结构
// 1 三元组Tuple/复数Complex
// 2 稠密矩阵Matrix
// 3 稀疏矩阵(普通，行压缩)SpraseCSR
// 4 稀疏矩阵(对称，行压缩)SpraseSM
// 5 稀疏矩阵(普通，十字链表)SpraseOL
// 6 稀疏矩阵(三对角，对角压缩)SpraseTD

// 包含如下辅助数据结构
// 1 零值类型萃取zero_traits

// 包含如下几处算法
// 1 线性方程组求解(高斯消元，列主元)
// 2 线性方程组求解(LU分解，杜立特尔，列主元)
// 3 线性方程组求解(逐次超松弛迭代法/高斯-赛德尔迭代法，稀疏)
// 4 线性方程组求解(QR分解，不选主元)
// 5 线性方程组求解(共轭梯度法，对称正定，不选主元)
// 6 线性方程组求解(追赶法，三对角矩阵，不选主元)

# ifndef MYSPRASE_H
# define MYSPRASE_H

namespace Unility
{
	namespace Sprase
	{
		// 基础数据结构：三元组
		template<typename TA, typename TB, typename TC>
		class Tuple
		{
		private:
			TA _first;
			TB _second;
			TC _third;
		public:
			const TA& First() const { return std::cref(_first); }
			const TB& Second() const { return std::cref(_second); }
			const TC& Third() const { return std::cref(_third); }
		public:
			Tuple(TA first, TB second, TC third) :
				_first(first), _second(second), _third(third)
			{

			}
		};

		// 基础数据结构：复数
		template<typename TN>
		class Complex
		{
		public:
			TN real;
			TN imag;

			template<typename TS, typename TN>
			friend TS& operator << (TS& stream, const Complex<TN>& complex);

			template<typename TN2>
			friend bool operator == (const Complex<TN2>& complex1, const Complex<TN2>& complex2);

			template<typename TN2>
			friend bool operator != (const Complex<TN2>& complex1, const Complex<TN2>& complex2);
		};

		template<typename TS, typename TN>
		TS& operator << (TS& stream, const Complex<TN>& complex)
		{
			stream << "(" << complex.real << ", " << complex.imag << ")" << std::endl;
			return stream;
		}

		template<typename TN2>
		bool operator == (const Complex<TN2>& complex1, const Complex<TN2>& complex2)
		{
			return complex1.real == complex2.real&&complex1.imag == complex2.imag;
		}

		template<typename TN2>
		bool operator != (const Complex<TN2>& complex1, const Complex<TN2>& complex2)
		{
			return !(complex1 == complex2);
		}

		typedef unsigned int TIndex;
		constexpr double mini_value = 1e-10;

		// 基础数据结构：稠密矩阵
		template<typename TN = double>
		class Matrix
		{
		private:
			std::vector<std::vector<TN> > _mat;
			// 稠密矩阵实例间互为友元
			friend class Matrix<TN>;
		public:
			// 获取稠密矩阵大小
			std::pair<TIndex, TIndex> GetSize() const { return std::pair<TIndex, TIndex>(_mat.size(), _mat.begin()->size()); }
			// 访问元素(Ref)
			TN& AtRef(TIndex row, TIndex col) const { return std::ref(_mat[row][col]); }
			// 访问元素(Val)
			TN AtVal(TIndex row, TIndex col) const { return _mat[row][col]; }
			// 访问行(Ref)
			std::vector<TN>& operator [] (TIndex row) { return std::ref(_mat[row]); }
			// (重新)初始化矩阵
			void Init(TIndex row, TIndex col)
			{
				_mat = std::vector<std::vector<TN> >(row);
				for (int i = 0; i < row; ++i)
				{
					_mat[i] = std::vector<TN>(col);
				}
			}
			void Init(std::vector<std::vector<TN> >&& mat) { _mat = mat; }
			void Init(const std::vector<std::vector<TN> >& mat) { _mat = mat; }
			void Init(Matrix&& mat) { _mat = std::move(mat._mat); }
			void Init(const Matrix& mat) { _mat = mat; }
			// 清除矩阵
			void Clear() { _mat.clear(); }
		public:
			// 向指定流对象输出
			template<typename TS, typename TN2>
			friend TS& operator << (TS& stream, const Matrix<TN2>& mat);
		public:
			template<typename TM>
			Matrix(TM&& mat)
			{
				Init(std::forward<TM>(mat));
			}

			Matrix(TIndex row, TIndex col)
			{
				Init(row, col);
			}

			Matrix()
			{
				
			}

			template<typename TM>
			Matrix& operator = (TM&& mat)
			{
				Init(std::forward<TM>(mat));
				return *this;
			}
		};

		template<typename TS, typename TN2>
		TS& operator << (TS& stream, const Matrix<TN2>& mat)
		{
			for (auto& r : mat._mat)
			{
				for (auto& e : r)
				{
					stream << e << " ";
				}
				stream << std::endl;
			}
			return stream;
		}

		// 定义非零类型的type_traits
		template<typename _Tp>
		struct zero_constant
		{
			typedef _Tp zero_type;
			typedef zero_constant<_Tp> type;
		};
		// 泛化
		template<typename TN>
		struct zero_traits :public zero_constant<TN>
		{
			static constexpr typename TN::zero_type zero_value = TN::zero_value;
		};
		// 指针偏特化版
		template<typename TN>
		struct zero_traits<TN*> :public zero_constant<TN*>
		{
			static constexpr zero_type zero_value = nullptr;
		};
		// const指针偏特化版
		template<typename TN>
		struct zero_traits<const TN*> :public zero_constant<const TN*>
		{
			static constexpr zero_type zero_value = nullptr;
		};
		// double特化版
		template<>
		struct zero_traits<double> :public zero_constant<double>
		{
			static constexpr zero_type zero_value = 0.0;
		};
		// float特化版
		template<>
		struct zero_traits<float> :public zero_constant<float>
		{
			static constexpr zero_type zero_value = 0.0;
		};
		// Int特化版
		template<>
		struct zero_traits<int> :public zero_constant<int>
		{
			static constexpr zero_type zero_value = 0;
		};
		// bool特化版
		template<>
		struct zero_traits<bool> :public zero_constant<bool>
		{
			static constexpr zero_type zero_value = false;
		};
		// unsigned int特化版
		template<>
		struct zero_traits<unsigned int> :public zero_constant<unsigned int>
		{
			static constexpr zero_type zero_value = 0;
		};
		// complex特化版
		template<>
		struct zero_traits< Complex<double> > :public zero_constant< Complex<double> >
		{
			static constexpr zero_type zero_value = Complex<double>();
		};

		// C++17下编译成功
		/*
		template<typename TN>
		using zero_traits_t = typename zero_traits<TN>::zero_type;

		template<typename TN>
		inline constexpr typename zero_traits<TN>::zero_type zero_traits_v = zero_traits<TN>::zero_value;
		*/

		class SpraseBase
		{

		};

		// 基础数据结构：稀疏矩阵(行压缩)
		// 使用vector作为容器，在随机访问时表现较优，在插入新元素时表现较差
		// 可以从三元组链表/二维矩阵/稠密矩阵构造
		template<typename TN>
		class SpraseCSR : SpraseBase
		{
		private:
			// 非零元素
			static constexpr typename zero_traits<TN>::zero_type zeroVal = zero_traits<TN>::zero_value;

			// 存储非零元素
			std::vector<TN> _values;
			// 存储列号
			std::vector<TIndex> _columnIndices;
			// 存储行偏移
			std::vector<TIndex> _rowOffsets;

			// 元素个数
			TIndex _nItems;
			// 矩阵真实大小
			std::pair<TIndex, TIndex> _mSize;
			// 有效行 [_vaildBegin, _vaildEnd)
			// 有效行起始行
			TIndex _vaildBegin;
			// 有效行结束行(的下一行)
			TIndex _vaildEnd;
			// 定位索引(用于Right函数)
			int _currentIndex;
			// 定位索引(用于Right函数)
			int _currentRow;

			// 实例间互为友元
			friend class SpraseCSR<TN>;
		private:
			void clearMat()
			{
				_values.clear();
				_columnIndices.clear();
				_rowOffsets.clear();
			}
			void newMat(TIndex nItems, TIndex rows)
			{
				_values = std::vector<TN>(nItems);
				_columnIndices = std::vector<TIndex>(nItems);
				_rowOffsets = std::vector<TIndex>(rows + 1);
			}
			auto getItem(TIndex row, TIndex col) const
			{
				// 不在有效行内直接返回
				if (row >= _vaildEnd || row < _vaildBegin)
					return _columnIndices.end() - _columnIndices.begin();
				// 在有效行范围内
				const auto& iter0 = _columnIndices.begin();
				const auto& iter1 = iter0 + _rowOffsets[row];
				const auto& iter2 = iter0 + _rowOffsets[row + 1];
				// 使用效率更高的二分查找
				// const auto& iter = std::find(iter1, iter2, col);
				// return (iter == iter2) ? _columnIndices.end() - iter0 : iter - iter0;
				const auto& iter = std::lower_bound(iter1, iter2, col);
				return (iter == iter2 || (*iter) != col) ? _columnIndices.end() - iter0 : iter - iter0;
			}
		public:
			// 向指定流对象输出
			template<typename TS, typename TN2, int z = 0>
			friend TS& operator << (TS& stream, const SpraseCSR<TN2>& mat);
		public:
			// 获取稀疏矩阵大小
			std::pair<TIndex, TIndex> GetSize() const { return _mSize; }
			// 访问元素(Ref)，若元素不存在，将抛出异常；插入新的元素将导致之前的引用失效
			TN& AtRef(TIndex row, TIndex col) const
			{
				return std::ref(_values[getItem(row, col)]);
			}
			// 访问元素(Val)，若元素不存在，将返回zeroVal
			TN AtVal(TIndex row, TIndex col) const
			{
				auto index = getItem(row, col);
				return (index - _nItems == 0) ? zeroVal : _values[index];
			}
			// 选中行
			std::pair<TN, bool> SelectRow(int& I, int& J, TIndex row)
			{
				if (_vaildBegin <= row && row < _vaildEnd&&_rowOffsets[row+1]- _rowOffsets[row] > 0)
				{
					_currentIndex = _rowOffsets[row];
					_currentRow = row;
					I = row;
					J = _columnIndices[_currentIndex];
					return std::pair<TN, bool>(_values[_currentIndex], true);
				}
				else
				{
					return std::pair<TN, bool>(zeroVal, false);
				}
			}
			// (在选中行的基础上)返回右边的元素
			std::pair<TN, bool> Right(int& I, int& J)
			{
				if (_currentIndex < _rowOffsets[_currentRow + 1])
				{
					_currentIndex++;
				}
				if (_currentIndex < _rowOffsets[_currentRow + 1])
				{
					I = _currentRow;
					J = _columnIndices[_currentIndex];
					return std::pair<TN, bool>(_values[_currentIndex], true);
				}
				else
				{
					return std::pair<TN, bool>(zeroVal, false);
				}
			}
			// (重新)初始化稀疏矩阵
			void Init(std::list<Tuple<TIndex, TIndex, TN> >& mat, TIndex nRow, TIndex nCol)
			{
				// 更新矩阵大小
				_mSize.first = nRow;
				_mSize.second = nCol;
				// 对list进行排序
				mat.sort(
					[](const Tuple<TIndex, TIndex, TN>& tuple1, const Tuple<TIndex, TIndex, TN>& tuple2) ->bool {
					if (tuple1.First() != tuple2.First())
					{
						return tuple1.First() < tuple2.First();
					}
					else
					{
						return tuple1.Second() < tuple2.Second();
					}
				});
				// 统计矩阵元素，并重新分配存储空间
				_nItems = mat.size();
				newMat(_nItems, nRow);
				// 初始化矩阵
				TIndex itemCount = 0;
				TIndex row = 0, col = 0, maxCol = 0;
				_vaildBegin = mat.begin()->First();
				for (TIndex i = 0; i < _vaildBegin; i++)
				{
					_rowOffsets[i] = 0;
					row = i;
				}
				_rowOffsets[_vaildBegin] = 0;
				for (auto& item : mat)
				{
					col = item.Second();
					_columnIndices[itemCount] = col;
					if (col > maxCol)
					{
						maxCol = col;
					}
					if (item.First() > row)
					{
						for (TIndex k = row + 1; k < item.First(); k++)
						{
							_rowOffsets[k] = itemCount;
						}
						row = item.First();
						_rowOffsets[row] = itemCount;
					}
					_values[itemCount++] = std::move(item.Third());
				}
				assert(maxCol <= _mSize.second);
				_vaildEnd = row + 1;
				_rowOffsets[_vaildEnd] = itemCount;
				for (TIndex i = _vaildEnd + 1; i < nRow + 1; i++)
				{
					_rowOffsets[i] = _rowOffsets[_vaildEnd];
				}
				_currentIndex = -1;
				_currentRow = -1;
			}
			void Init(Matrix<TN>& mat, TIndex nRow, TIndex nCol)
			{
				auto size = mat.GetSize();
				std::list<Tuple<TIndex, TIndex, TN> > Mat;
				for (TIndex row = 0; row < size.first; row++)
				{
					for (TIndex col = 0; col < size.second; col++)
					{
						if (mat[row][col] != zeroVal)
						{
							Mat.emplace_back(row, col, mat[row][col]);
						}
					}
				}
				Init(Mat, nRow, nCol);
			}
			void Init(const std::vector<std::vector<TN> >& mat, TIndex nRow, TIndex nCol)
			{

			}
			// 清除矩阵
			void Clear()
			{
				clearMat();
				_nItems = 0;
				_mSize.first = 0;
				_mSize.second = 0;
				_vaildBegin = 0;
				_vaildEnd = 0;
				_currentIndex = -1;
				_currentRow = -1;
			}
		public:
			template<typename TM>
			SpraseCSR(TM&& mat, TIndex nRow, TIndex nCol)
				:_mSize(nRow, nCol), _nItems(0), _vaildBegin(0), _vaildEnd(0), _currentIndex(-1), _currentRow(-1)
			{
				Init(std::forward<TM>(mat), nRow, nCol);
			}

			SpraseCSR()
				:_mSize(0, 0), _nItems(0), _vaildBegin(0), _vaildEnd(0), _currentIndex(-1), _currentRow(-1)
			{

			}
		};

		template<typename TS, typename TN2, int z = 0>
		TS& operator << (TS& stream, const SpraseCSR<TN2>& mat)
		{
			auto size = mat.GetSize();
			for (TIndex i = 0; i < size.first; i++)
			{
				for (TIndex j = 0; j < size.second; j++)
				{
					stream << mat.AtVal(i, j) << " ";
				}
				stream << std::endl;
			}
			return stream;
		}

		// 基础数据结构：稀疏矩阵(对称，行压缩)
		// 使用vector作为容器，支持随机访问，插入新元素时性能较差
		template<typename TN>
		class SpraseSM : SpraseBase
		{
		private:
			// 非零元素
			static constexpr typename zero_traits<TN>::zero_type zeroVal = zero_traits<TN>::zero_value;

			// 对角向量
			std::vector<TN> _diag;
			// 上三角向量
			std::vector<TN> _upTria;
			// 存储上三角列号
			std::vector<TIndex> _columnIndices;
			// 存储上三角行偏移
			std::vector<TIndex> _rowOffsets;

			// 元素个数
			TIndex _nItems;
			// 上三角元素个数
			TIndex _nUpItems;
			// 矩阵真实大小
			std::pair<TIndex, TIndex> _mSize;
			// 有效行 [_vaildBegin, _vaildEnd)
			// 有效行起始行
			TIndex _vaildBegin;
			// 有效行结束行(的下一行)
			TIndex _vaildEnd;

			// 实例间互为友元
			friend class SpraseSM<TN>;
		private:
			void clearMat()
			{
				_diag.clear();
				_upTria.clear();
				_columnIndices.clear();
				_rowOffsets.clear();
			}
			void newMat(TIndex nItems, TIndex nUpItems, TIndex rows)
			{
				// 会在对角线上造成少量的空间浪费
				_diag = std::vector<TN>(rows);
				_upTria = std::vector<TN>(nUpItems);
				_columnIndices = std::vector<TIndex>(nUpItems);
				_rowOffsets = std::vector<TIndex>(rows + 1);
			}
			auto getUpItem(TIndex row, TIndex col) const
			{
				// 得到正确的行列索引
				if (row > col)
				{
					std::swap(row, col);
				}
				assert(row != col);
				// 不在有效行内直接返回
				if (row >= _vaildEnd || row < _vaildBegin)
					return _columnIndices.end() - _columnIndices.begin();
				// 在有效行范围内
				const auto& iter0 = _columnIndices.begin();
				const auto& iter1 = iter0 + _rowOffsets[row];
				const auto& iter2 = iter0 + _rowOffsets[row + 1];
				// 使用效率更高的二分查找
				const auto& iter = std::lower_bound(iter1, iter2, col);
				return (iter == iter2 || (*iter) != col) ? _columnIndices.end() - iter0 : iter - iter0;
			}
		public:
			// 向指定流对象输出矩阵
			template<typename TS, typename TN2>
			friend TS& operator << (TS& stream, const SpraseSM<TN2>& mat);
		public:
			// 获取稀疏矩阵大小
			std::pair<TIndex, TIndex> GetSize() const { return _mSize; }
			// 访问元素(Ref)，若元素不存在，将抛出异常；插入新的元素将导致之前的引用失效
			TN& AtRef(TIndex row, TIndex col) const
			{
				if (row == col)
				{
					// 对角元素
					return std::ref(_diag[row]);
				}
				else
				{
					// 上三角元素
					return std::ref(_upTria[getUpItem(row, col)]);
				}
			}
			// 访问元素(Val)，若元素不存在，将返回zeroVal
			TN AtVal(TIndex row, TIndex col) const
			{
				if (row == col)
				{
					// 对角元素
					return _diag[row];
				}
				else
				{
					// 上三角元素
					auto index = getUpItem(row, col);
					return (index - _nUpItems == 0) ? zeroVal : _upTria[index];
				}
			}
			// 初始化稀疏矩阵-由三元组列表
			void Init(std::list<Tuple<TIndex, TIndex, TN> >& diagMat, std::list<Tuple<TIndex, TIndex, TN> >& upMat, TIndex nRow, TIndex nCol)
			{
				assert(nRow == nCol);
				// 更新矩阵大小
				_mSize.first = nRow;
				_mSize.second = nCol;
				// 对list进行排序
				upMat.sort(
					[](const Tuple<TIndex, TIndex, TN>& tuple1, const Tuple<TIndex, TIndex, TN>& tuple2) ->bool {
					if (tuple1.First() != tuple2.First())
					{
						return tuple1.First() < tuple2.First();
					}
					else
					{
						return tuple1.Second() < tuple2.Second();
					}
				});
				diagMat.sort(
					[](const Tuple<TIndex, TIndex, TN>& tuple1, const Tuple<TIndex, TIndex, TN>& tuple2) ->bool {
					return tuple1.First() < tuple2.First();
				});
				// 统计矩阵元素，并重新分配存储空间
				_nUpItems = upMat.size();
				_nItems = upMat.size() + diagMat.size();
				newMat(_nItems, _nUpItems, nRow);
				// 初始化
				// 对角元素
				TIndex begin = 0;
				for (auto& d : diagMat)
				{
					for (TIndex i = begin; i < d.First(); i++)
					{
						_diag[i] = zeroVal;
						begin++;
					}
					begin++;
					_diag[d.First()] = d.Third();
				}
				// 非对角元素
				TIndex itemCount = 0;
				TIndex row = 0, col = 0, maxCol = 0;
				_vaildBegin = upMat.begin()->First();
				for (TIndex i = 0; i < _vaildBegin; i++)
				{
					_rowOffsets[i] = 0;
					row = i;
				}
				_rowOffsets[_vaildBegin] = 0;
				for (auto& up : upMat)
				{
					col = up.Second();
					_columnIndices[itemCount] = col;
					if (col > maxCol)
					{
						maxCol = col;
					}
					if (up.First() > row)
					{
						for (TIndex k = row + 1; k < up.First(); k++)
						{
							_rowOffsets[k] = itemCount;
						}
						row = up.First();
						_rowOffsets[row] = itemCount;
					}
					_upTria[itemCount++] = std::move(up.Third());
				}
				assert(maxCol <= _mSize.second);
				_vaildEnd = row + 1;
				_rowOffsets[_vaildEnd] = itemCount;
				for (TIndex i = _vaildEnd + 1; i < nRow + 1; i++)
				{
					_rowOffsets[i] = _rowOffsets[_vaildEnd];
				}
			}
			// 初始化稀疏矩阵-由稠密矩阵(对称方阵)
			void Init(Matrix<TN>& mat, TIndex nRow, TIndex nCol)
			{
				auto size = mat.GetSize();
				assert(size.first == size.second);
				nRow = std::max(size.first, nRow);
				nCol = std::max(size.second, nCol);
				assert(nRow == nCol);
				// 建立list
				std::list<Tuple<TIndex, TIndex, TN> > diagMat;
				std::list<Tuple<TIndex, TIndex, TN> > upMat;
				for (TIndex row = 0; row < size.first; row++)
				{
					for (TIndex col = row; col < size.second; col++)
					{
						if (row == col) {
							diagMat.emplace_back(row, col, mat[row][col]);
						}
						else
						{
							if (mat[row][col] == zeroVal)
							{
								continue;
							}
							else
							{
								upMat.emplace_back(row, col, mat[row][col]);
							}
						}
					}
				}
				// 初始化
				Init(diagMat, upMat, nRow, nCol);
			}

			// 清除矩阵
			void Clear()
			{
				clearMat();
				_nItems = 0;
				_nUpItems = 0;
				_mSize.first = 0;
				_mSize.second = 0;
				_vaildBegin = 0;
				_vaildEnd = 0;
			}
		public:
			template<typename TM>
			SpraseSM(TM&& diagMat, TM&& upMat, TIndex nRow, TIndex nCol)
				:_mSize(nRow, nCol), _nItems(0), _nUpItems(0), _vaildBegin(0), _vaildEnd(0)
			{
				Init(std::forward<TM>(diagMat), std::forward<TM>(upMat), nRow, nCol);
			}

			template<typename TM>
			SpraseSM(TM&& mat, TIndex nRow, TIndex nCol)
				: _mSize(nRow, nCol), _nItems(0), _nUpItems(0), _vaildBegin(0), _vaildEnd(0)
			{
				Init(std::forward<TM>(mat), nRow, nCol);
			}

			SpraseSM()
				:_mSize(0, 0), _nItems(0), _nUpItems(0), _vaildBegin(0), _vaildEnd(0)
			{

			}
		};

		template<typename TS, typename TN2>
		TS& operator << (TS& stream, const SpraseSM<TN2>& mat)
		{
			auto size = mat.GetSize();
			for (TIndex i = 0; i < size.first; i++)
			{
				for (TIndex j = 0; j < size.second; j++)
				{
					stream << mat.AtVal(i, j) << " ";
				}
				stream << std::endl;
			}
			return stream;
		}


		// 十字链表支路节点(矩阵对应元素)
		template<typename TL>
		struct ListNode
		{
			// 有向图支路结束顶点、起始顶点
			const TIndex row, col;
			// 指向弧头相同的下一条支路、指向弧尾相同的下一条支路
			std::shared_ptr<ListNode<TL> > down, right;
			// 支路数据
			TL info;

			ListNode(const TIndex& r, const TIndex& c, std::shared_ptr<ListNode<TL> > d, std::shared_ptr<ListNode<TL> > rt, const TL& i)
				: row(r), col(c), down(d), right(rt), info(i)
			{

			}

			ListNode()
				: row(0), col(0), down(nullptr), right(nullptr), info(zero_traits<TN>::zero_value)
			{

			}
		};

		// 十字链表顶点节点(矩阵行列索引)
		template<typename TL>
		struct VexNode
		{
			// 指向以该顶点为弧头、弧尾的第一条支路
			std::shared_ptr<ListNode<TL> > asCol, asRow;

			VexNode(std::shared_ptr<ListNode<TL> > r, std::shared_ptr<ListNode<TL> > c)
				: asRow(r), asCol(c)
			{

			}

			VexNode()
				: asRow(nullptr), asCol(nullptr)
			{

			}
		};

		// 对应关系说明：
		// traiVex -> row
		// headVex -> col
		// hLink -> down
		// tLink -> right
		// firstIn -> cols[i]
		// firstOut -> rows[i]

		// 基础数据结构：稀疏矩阵(普通，十字链表)
		// 参考https://blog.csdn.net/bible_reader/article/details/71214096
		// 支持插入和删除元素，适用于存储稀疏方阵或者有向图
		template<typename TN>
		class SpraseOL : SpraseBase
		{
		private:
			// 非零元素
			static constexpr typename zero_traits<TN>::zero_type zeroVal = zero_traits<TN>::zero_value;

			// 顶点列表
			std::vector<VexNode<TN> > _vexList;

			// 矩阵真实大小
			std::pair<TIndex, TIndex> _mSize;
		private:
			// 添加顶点
			void insertVex()
			{
				_vexList.emplace_back(nullptr, nullptr);
			}
			// 删除节点
			void deleteVex(TIndex vex)
			{

			}
		public:
			// 获取稀疏矩阵大小
			std::pair<TIndex, TIndex> GetSize() const { return _mSize; }
			// 访问元素(指针)
			std::shared_ptr<ListNode<TN> > AtPtr(TIndex row, TIndex col) const
			{
				auto cur = this->_vexList[col].asCol;
				if (cur)
				{
					while (cur)
					{
						if (cur->row == row)
						{
							break;
						}
						cur = cur->down;
					}
				}
				return cur;
			}
			// 访问元素(引用)
			TN& AtRef(TIndex row, TIndex col) const
			{
				auto cur = AtPtr(row, col);
				return std::ref(cur->info);
			}
			// 访问元素(值)
			TN AtVal(TIndex row, TIndex col) const
			{
				auto cur = AtPtr(row, col);
				return cur ? cur->info : zeroVal;
			}
			// 插入元素(不保证符合矩阵元素顺序，覆盖原值)
			bool Insert(TIndex row, TIndex col, const TN& value)
			{
				auto pRow = _vexList[row].asRow;
				auto pCol = _vexList[col].asCol;
				auto pItem = std::make_shared<ListNode<TN> >(row, col, nullptr, nullptr, value);
				if (pRow)
				{
					pItem->right = pRow;
				}
				if (pCol)
				{
					pItem->down = pCol;
				}
				_vexList[row].asRow = pItem;
				_vexList[col].asCol = pItem;
				return true;
			}
			// 更新元素(保证矩阵元素顺序，覆盖)
			bool Update(TIndex row, TIndex col, const TN& value)
			{
				// 暂时没完成
				return true;
			}
			// 删除元素
			bool Delete(TIndex row, TIndex col)
			{
				// 断开行链表
				auto deleteOfRow = [&]() -> bool {
					auto cur = this->_vexList[row].asRow;
					auto pre = cur;
					TIndex count = 0;
					if (cur)
					{
						while (cur)
						{
							count++;
							if (cur->col == col)
							{
								break;
							}
							pre = cur;
							cur = cur->right;
						}
					}
					else
					{
						return false;
					}
					if (!cur)
					{
						return false;
					}
					else if (count <= 1)
					{
						this->_vexList[row].asRow
					}
					else
					{
						pre->right = cur->right;
					}
					return true;
				};
				// 断开列链表
				auto deleteOfCol = [&]() -> bool {
					auto cur = this->_vexList[col].asCol;
					auto pre = cur;
					TIndex count = 0;
					if (cur)
					{
						while (cur)
						{
							count++;
							if (cur->row == row)
							{
								break;
							}
							pre = cur;
							cur = cur->down;
						}
					}
					else
					{
						return false;
					}
					if (!cur)
					{
						return false;
					}
					else if (count <= 1)
					{
						this->_vexList[col].asCol = pre->down;
					}
					else
					{
						pre->down = cur->down;
					}
					return true;
				};
				// 两个链表都断开，智能指针自动析构
				auto result1 = deleteOfRow();
				auto result2 = deleteOfCol();
				assert(result1 == result2);
				return result1 && result2;
			}
			// 清除矩阵
			void Clear()
			{
				_vexList.clear();
				_mSize.first = 0;
				_mSize.second = 0;
			}
			// 对系数矩阵进行排序，使其符合矩阵的存储顺序
			void Sort()
			{
				// 保存非零元素(指针)
				std::list<std::shared_ptr<ListNode<TN> > > itemSave;
				for (auto& curNode : _vexList)
				{
					auto curItem = curNode.asRow;
					while (curItem)
					{
						itemSave.emplace_back(curItem);
						curItem = curItem->right;
					}
					curNode.asRow = nullptr;
					curNode.asCol = nullptr;
				}
				// 重新排序
				itemSave.sort(
					[](std::shared_ptr<ListNode<TN> > item1, std::shared_ptr<ListNode<TN> > item2) 
				{
					if (item1->row != item2->row)
					{
						return item1->row > item2->row;
					}
					else
					{
						return item1->col > item2->col;
					}
				});
				// 重新插入
				for (auto pItem : itemSave)
				{
					auto pRow = _vexList[pItem->row].asRow;
					auto pCol = _vexList[pItem->col].asCol;
					pItem->right = pRow;
					pItem->down = pCol;
					_vexList[pItem->row].asRow = pItem;
					_vexList[pItem->col].asCol = pItem;
				}
			}
			// 初始化十字链表稀疏矩阵
			void Init(TIndex nRow, TIndex nCol)
			{
				auto nVex = std::max(nRow, nCol);
				_mSize.first = nRow;
				_mSize.second = nCol;
				_vexList = std::vector<VexNode<TN> >(nVex);
			}
			// 初始化十字链表稀疏矩阵(由三元组列表，同时检查重复元素)
			void Init(std::list<Tuple<TIndex, TIndex, TN> >& mat, TIndex nRow, TIndex nCol)
			{
				// 对list进行排序(大->小)
				mat.sort(
					[](const Tuple<TIndex, TIndex, TN>& tuple1, const Tuple<TIndex, TIndex, TN>& tuple2) ->bool {
					if (tuple1.First() != tuple2.First())
					{
						return tuple1.First() > tuple2.First();
					}
					else
					{
						return tuple1.Second() > tuple2.Second();
					}
				});
				// 删除重复元素
				mat.unique(
					[](const Tuple<TIndex, TIndex, TN>& tuple1, const Tuple<TIndex, TIndex, TN>& tuple2) ->bool {
					return (tuple1.First() == tuple2.First()) && (tuple1.Second() == tuple2.Second());
				});
				// 初始化
				Init(nRow, nCol);
				// 插入
				for (auto& item : mat)
				{
					Insert(item.First(), item.Second(), item.Third());
				}
			}
			// 初始化十字链表稀疏矩阵(由稠密矩阵)
			void Init(const Matrix<TN>& mat, TIndex nRow, TIndex nCol)
			{
				// 初始化
				Init(nRow, nCol);
				// 插入
				for (TIndex i = 0; i < mat.GetSize().first; i++)
				{
					for (TIndex j = 0; j < mat.GetSize().second; j++)
					{
						auto value = mat.AtVal(i, j);
						if (value != zeroVal)
						{
							Insert(i, j, value);
						}
					}
				}
			}
			// 初始化十字链表稀疏矩阵(由二维数组)
			void Init(const std::vector<std::vector<TN> >& mat, TIndex nRow, TIndex nCol)
			{

			}
		public:
			// 向指定流对象输出矩阵
			template<typename TS, typename TN3>
			friend TS& operator << (TS& stream, const SpraseOL<TN3>& mat);
		public:
			template<typename TM>
			SpraseOL(TM&& mat, TIndex nRow, TIndex nCol)
			{
				Init(std::forward<TM>(mat), nRow, nCol);
			}
			
			SpraseOL(TIndex nRow, TIndex nCol)
			{
				Init(nRow, nCol);
			}

			SpraseOL()
			{

			}
		};

		template<typename TS, typename TN3>
		TS& operator << (TS& stream, const SpraseOL<TN3>& mat)
		{
			auto size = mat.GetSize();
			for (TIndex i = 0; i < size.first; i++)
			{
				for (TIndex j = 0; j < size.second; j++)
				{
					stream << mat.AtVal(i, j) << " ";
				}
				stream << std::endl;
			}
			return stream;
		}

		// 基础数据结构：稀疏矩阵(三对角，对角压缩)
		// 主要用于表示弯矩方程组或微分方程组
		// 不支持插入
		template<typename TN>
		class SpraseTD : SpraseBase
		{
		private:
			// 非零元素
			static constexpr typename zero_traits<TN>::zero_type zeroVal = zero_traits<TN>::zero_value;

			// 三对角主对角向量
			std::vector<double> _diagB;
			// 三对角下次主对角线
			std::vector<double> _diagA;
			// 三对角上次主对角线
			std::vector<double> _diagC;

			// 矩阵真实大小
			std::pair<TIndex, TIndex> _mSize;

			// 实例间互为友元
			friend class SpraseTD<TN>;
		public:
			// 向指定流对象输出
			template<typename TS, typename TN4, int z = 0>
			friend TS& operator << (TS& stream, const SpraseTD<TN4>& mat);
		public:
			// 获取稀疏矩阵大小
			std::pair<TIndex, TIndex> GetSize() const { return _mSize; }
			// 访问元素(Ref)，若元素不存在，将抛出异常
			TN& AtRef(TIndex row, TIndex col) const
			{
				if (row == col)
				{
					// 主对角
					return std::ref(_diagB[col]);
				}
				else if(row-col == 1)
				{
					// 下次主对角
					return std::ref(_diagA[col]);
				}
				else if (col - row == 1)
				{
					// 上次主对角
					return std::ref(_diagC[row]);
				}
				else
				{
					// 越界！
					return std::ref(_diagB[_mSize.first]);
				}
			}
			// 访问元素(Val)，若元素不存在，将返回zeroVal
			TN AtVal(TIndex row, TIndex col) const
			{
				if (row == col)
				{
					// 主对角
					return _diagB[col];
				}
				else if (row - col == 1)
				{
					// 下次主对角
					return _diagA[col];
				}
				else if (col - row == 1)
				{
					// 上次主对角
					return _diagC[row];
				}
				else
				{
					// 其他位置，应为0元素
					return zeroVal;
				}
			}
			// 访问主对角(Ref, Offset)，若元素不存在，将抛出异常
			TN& AtDiagBRef(TIndex offset)
			{
				return std::ref(_diagB[offset]);
			}
			// 访问下次对角(Ref, Offset)，若元素不存在，将抛出异常
			TN& AtDiagARef(TIndex offset)
			{
				return std::ref(_diagA[offset]);
			}
			// 访问上次对角(Ref, Offset)，若元素不存在，将抛出异常
			TN& AtDiagCRef(TIndex offset)
			{
				return std::ref(_diagC[offset]);
			}
		public:
			// 清空矩阵
			void Clear()
			{
				_diagB.clear();
				_diagA.clear();
				_diagC.clear();
				_mSize.first = 0;
				_mSize.second = 0;
			}
			// (重新)初始化矩阵
			void Init(TIndex mSize)
			{
				_diagB = std::vector<double>(mSize);
				_diagA = std::vector<double>(mSize - 1);
				_diagC = std::vector<double>(mSize - 1);
				_mSize.first = mSize;
				_mSize.second = mSize;
			}
			void Init(Matrix<TN>& mat)
			{
				auto mSize = mat.GetSize();
				assert(mSize.first == mSize.second);
				Init(mSize.first);
				for (int k = 0; k < mSize.first - 1; k++)
				{
					_diagB[k] = mat[k][k];
					_diagA[k] = mat[k + 1][k];
					_diagC[k] = mat[k][k + 1];
				}
				_diagB[mSize.first-1]= mat[mSize.first - 1][mSize.first - 1];
				_mSize = mSize;
			}
			void Init(std::vector<TN>& diagA, std::vector<TN>& diagB, std::vector<TN>& diagC)
			{
				auto mSize = diagB.size();
				assert(diagA.size() == mSize - 1);
				assert(diagC.size() == mSize - 1);
				_diagB = diagB;
				_diagA = diagA;
				_diagC = diagC;
				_mSize.first = mSize;
				_mSize.second = mSize;
			}
			void Init(std::vector<TN>&& diagA, std::vector<TN>&& diagB, std::vector<TN>&& diagC)
			{
				auto mSize = diagB.size();
				assert(diagA.size() == mSize - 1);
				assert(diagC.size() == mSize - 1);
				_diagB = diagB;
				_diagA = diagA;
				_diagC = diagC;
				_mSize.first = mSize;
				_mSize.second = mSize;
			}
		public:
			// 构造函数
			template<typename TM>
			SpraseTD(TM&& mat, TIndex mSize)
				:_mSize(mSize, mSize)
			{
				Init(std::forward<TM>(mat));
			}

			SpraseTD(TIndex mSize)
				:_mSize(mSize, mSize)
			{
				Init(mSize);
			}

			template<typename TV>
			SpraseTD(TV&& diagA, TV&& diagB, TV&& diagC)
				:_mSize(0, 0)
			{
				Init(std::forward<TV>(diagA), std::forward<TV>(diagB), std::forward<TV>(diagC));
			}

			SpraseTD()
				:_mSize(0, 0)
			{

			}
		};



		// 高斯消去法
		// 参考《计算方法》(李乃成、梅立泉著，科学出版社)
		template<typename TM>
		std::vector<double> GESolver(TM&& A, std::vector<double>& b)
		{
			return GEKernel(std::forward<TM>(A), b);
		}
		// 高斯消去法核心函数，返回结果向量
		std::vector<double> GEKernel(Matrix<double>& A, std::vector<double>& b);

		// LU/LUP分解法
		// 参考《算法导论 第二版》(Thomas H.Cormen等著，潘金贵等译，机械工业出版社)
		template<typename TM>
		std::vector<double> LUSolver(TM&& A, std::vector<double>& b)
		{
			auto&& pi = LUDecomposition(std::forward<TM>(A));
			return LUSolution(std::forward<TM>(A), b, pi);
		}

		// LU分解，返回转换矩阵
		std::vector<double> LUDecomposition(Matrix<double>& A);
		std::vector<double> LUDecomposition(SpraseOL<double>& A);
		// LU回代，返回结果向量
		std::vector<double> LUSolution(Matrix<double>& A, std::vector<double>& b, std::vector<double>& pi);
		std::vector<double> LUSolution(SpraseOL<double>& A, std::vector<double>& b, std::vector<double>& pi);

		// 逐次超松弛迭代法/高斯赛德尔迭代法
		// 参考《计算方法》(李乃成、梅立泉著，科学出版社)
		template<typename TM>
		std::vector<double> SORSolver(TM&& A, std::vector<double>& b, std::vector<double>& xPre, double weight = 1)
		{
			return SORKernel(std::forward<TM>(A), b, xPre, weight);
		}
		// SOR法收敛精度
		constexpr double SORConvergeLimit = 1e-6;
		// 逐次超松弛迭代法核心函数
		std::vector<double> SORKernel(SpraseCSR<double>&A, std::vector<double>& b, std::vector<double>& xPre, double weight = 1.0);

		// 雅克比迭代法
		// 参考《计算方法》(李乃成、梅立泉著，科学出版社)
		template<typename TM>
		std::vector<double> JACSolver(TM&& A, std::vector<double>& b)
		{
			return JACKernel(std::forward<TM>(A), b);
		}
		// 雅克比迭代法收敛精度
		constexpr double JACConvergeLimit = 1e-6;
		// 雅克比迭代法核心函数
		std::vector<double> JACKernel(SpraseCSR<double>&A, std::vector<double>& b);

		// QR分解法
		// 参考《计算方法》(李乃成、梅立泉著，科学出版社)
		template<typename TM>
		std::vector<double> QRSolver(TM&& A, std::vector<double>& b)
		{
			return QRKernel(std::forward<TM>(A), b);
		}
		// QR分解函数
		std::pair<std::vector<double>, std::vector<double> > QRDecomposition(Matrix<double>&A, bool isExtension = false);
		// QR求解核心函数
		std::vector<double> QRKernel(Matrix<double>&A, std::vector<double>& b);

		// 共轭梯度法
		// 参考《计算方法》(李乃成、梅立泉著，科学出版社)
		template<typename TM>
		std::vector<double> CGSolver(TM&& A, std::vector<double>& b)
		{
			return CGKernel(std::forward<TM>(A), b);
		}
		// 共轭梯度法收敛精度
		constexpr double CGConvergeLimit = 1e-6;
		// 共轭梯度法核心函数
		std::vector<double> CGKernel(Matrix<double>& A, std::vector<double>& b);

		// 追赶法
		// 参考《计算方法》(李乃成、梅立泉著，科学出版社)
		template<typename TM>
		std::vector<double> CMSolver(TM&& A, std::vector<double>& d)
		{
			return CMKernel(std::forward<TM>(A), d);
		}
		// 追赶法核心函数
		std::vector<double> CMKernel(SpraseTD<double>& A, std::vector<double>& d);
	}
}

# endif