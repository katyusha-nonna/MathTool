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

// �������»������ݽṹ
// 1 ��Ԫ��Tuple/����Complex
// 2 ���ܾ���Matrix
// 3 ϡ�����(��ͨ����ѹ��)SpraseCSR
// 4 ϡ�����(�Գƣ���ѹ��)SpraseSM
// 5 ϡ�����(��ͨ��ʮ������)SpraseOL
// 6 ϡ�����(���Խǣ��Խ�ѹ��)SpraseTD

// �������¸������ݽṹ
// 1 ��ֵ������ȡzero_traits

// �������¼����㷨
// 1 ���Է��������(��˹��Ԫ������Ԫ)
// 2 ���Է��������(LU�ֽ⣬�����ض�������Ԫ)
// 3 ���Է��������(��γ��ɳڵ�����/��˹-���¶���������ϡ��)
// 4 ���Է��������(QR�ֽ⣬��ѡ��Ԫ)
// 5 ���Է��������(�����ݶȷ����Գ���������ѡ��Ԫ)
// 6 ���Է��������(׷�Ϸ������ԽǾ��󣬲�ѡ��Ԫ)

# ifndef MYSPRASE_H
# define MYSPRASE_H

namespace Utility
{
	namespace Sprase
	{
		// �������ݽṹ����Ԫ��
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

		// �������ݽṹ������
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

		// �������ݽṹ�����ܾ���
		template<typename TN = double>
		class Matrix
		{
		private:
			std::vector<std::vector<TN> > _mat;
			// ���ܾ���ʵ���以Ϊ��Ԫ
			friend class Matrix<TN>;
		public:
			// ��ȡ���ܾ����С
			std::pair<TIndex, TIndex> GetSize() const { return std::pair<TIndex, TIndex>(_mat.size(), _mat.begin()->size()); }
			// ����Ԫ��(Ref)
			TN& AtRef(TIndex row, TIndex col) const { return std::ref(_mat[row][col]); }
			// ����Ԫ��(Val)
			TN AtVal(TIndex row, TIndex col) const { return _mat[row][col]; }
			// ������(Ref)
			std::vector<TN>& operator [] (TIndex row) { return std::ref(_mat[row]); }
			// (����)��ʼ������
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
			// �������
			void Clear() { _mat.clear(); }
		public:
			// ��ָ�����������
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

		// ����������͵�type_traits
		template<typename _Tp>
		struct zero_constant
		{
			typedef _Tp zero_type;
			typedef zero_constant<_Tp> type;
		};
		// ����
		template<typename TN>
		struct zero_traits :public zero_constant<TN>
		{
			static constexpr typename TN::zero_type zero_value = TN::zero_value;
		};
		// ָ��ƫ�ػ���
		template<typename TN>
		struct zero_traits<TN*> :public zero_constant<TN*>
		{
			static constexpr zero_type zero_value = nullptr;
		};
		// constָ��ƫ�ػ���
		template<typename TN>
		struct zero_traits<const TN*> :public zero_constant<const TN*>
		{
			static constexpr zero_type zero_value = nullptr;
		};
		// double�ػ���
		template<>
		struct zero_traits<double> :public zero_constant<double>
		{
			static constexpr zero_type zero_value = 0.0;
		};
		// float�ػ���
		template<>
		struct zero_traits<float> :public zero_constant<float>
		{
			static constexpr zero_type zero_value = 0.0;
		};
		// Int�ػ���
		template<>
		struct zero_traits<int> :public zero_constant<int>
		{
			static constexpr zero_type zero_value = 0;
		};
		// bool�ػ���
		template<>
		struct zero_traits<bool> :public zero_constant<bool>
		{
			static constexpr zero_type zero_value = false;
		};
		// unsigned int�ػ���
		template<>
		struct zero_traits<unsigned int> :public zero_constant<unsigned int>
		{
			static constexpr zero_type zero_value = 0;
		};
		// complex�ػ���
		template<>
		struct zero_traits< Complex<double> > :public zero_constant< Complex<double> >
		{
			static constexpr zero_type zero_value = Complex<double>();
		};

		// C++17�±���ɹ�
		/*
		template<typename TN>
		using zero_traits_t = typename zero_traits<TN>::zero_type;

		template<typename TN>
		inline constexpr typename zero_traits<TN>::zero_type zero_traits_v = zero_traits<TN>::zero_value;
		*/

		class SpraseBase
		{

		};

		// �������ݽṹ��ϡ�����(��ѹ��)
		// ʹ��vector��Ϊ���������������ʱ���ֽ��ţ��ڲ�����Ԫ��ʱ���ֽϲ�
		// ���Դ���Ԫ������/��ά����/���ܾ�����
		template<typename TN>
		class SpraseCSR : SpraseBase
		{
		private:
			// ����Ԫ��
			static constexpr typename zero_traits<TN>::zero_type zeroVal = zero_traits<TN>::zero_value;

			// �洢����Ԫ��
			std::vector<TN> _values;
			// �洢�к�
			std::vector<TIndex> _columnIndices;
			// �洢��ƫ��
			std::vector<TIndex> _rowOffsets;

			// Ԫ�ظ���
			TIndex _nItems;
			// ������ʵ��С
			std::pair<TIndex, TIndex> _mSize;
			// ��Ч�� [_vaildBegin, _vaildEnd)
			// ��Ч����ʼ��
			TIndex _vaildBegin;
			// ��Ч�н�����(����һ��)
			TIndex _vaildEnd;
			// ��λ����(����Right����)
			int _currentIndex;
			// ��λ����(����Right����)
			int _currentRow;

			// ʵ���以Ϊ��Ԫ
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
				// ������Ч����ֱ�ӷ���
				if (row >= _vaildEnd || row < _vaildBegin)
					return _columnIndices.end() - _columnIndices.begin();
				// ����Ч�з�Χ��
				const auto& iter0 = _columnIndices.begin();
				const auto& iter1 = iter0 + _rowOffsets[row];
				const auto& iter2 = iter0 + _rowOffsets[row + 1];
				// ʹ��Ч�ʸ��ߵĶ��ֲ���
				// const auto& iter = std::find(iter1, iter2, col);
				// return (iter == iter2) ? _columnIndices.end() - iter0 : iter - iter0;
				const auto& iter = std::lower_bound(iter1, iter2, col);
				return (iter == iter2 || (*iter) != col) ? _columnIndices.end() - iter0 : iter - iter0;
			}
		public:
			// ��ָ�����������
			template<typename TS, typename TN2, int z = 0>
			friend TS& operator << (TS& stream, const SpraseCSR<TN2>& mat);
		public:
			// ��ȡϡ������С
			std::pair<TIndex, TIndex> GetSize() const { return _mSize; }
			// ����Ԫ��(Ref)����Ԫ�ز����ڣ����׳��쳣�������µ�Ԫ�ؽ�����֮ǰ������ʧЧ
			TN& AtRef(TIndex row, TIndex col) const
			{
				return std::ref(_values[getItem(row, col)]);
			}
			// ����Ԫ��(Val)����Ԫ�ز����ڣ�������zeroVal
			TN AtVal(TIndex row, TIndex col) const
			{
				auto index = getItem(row, col);
				return (index - _nItems == 0) ? zeroVal : _values[index];
			}
			// ѡ����
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
			// (��ѡ���еĻ�����)�����ұߵ�Ԫ��
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
			// (����)��ʼ��ϡ�����
			void Init(std::list<Tuple<TIndex, TIndex, TN> >& mat, TIndex nRow, TIndex nCol)
			{
				// ���¾����С
				_mSize.first = nRow;
				_mSize.second = nCol;
				// ��list��������
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
				// ͳ�ƾ���Ԫ�أ������·���洢�ռ�
				_nItems = mat.size();
				newMat(_nItems, nRow);
				// ��ʼ������
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
			// �������
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

		// �������ݽṹ��ϡ�����(�Գƣ���ѹ��)
		// ʹ��vector��Ϊ������֧��������ʣ�������Ԫ��ʱ���ܽϲ�
		template<typename TN>
		class SpraseSM : SpraseBase
		{
		private:
			// ����Ԫ��
			static constexpr typename zero_traits<TN>::zero_type zeroVal = zero_traits<TN>::zero_value;

			// �Խ�����
			std::vector<TN> _diag;
			// ����������
			std::vector<TN> _upTria;
			// �洢�������к�
			std::vector<TIndex> _columnIndices;
			// �洢��������ƫ��
			std::vector<TIndex> _rowOffsets;

			// Ԫ�ظ���
			TIndex _nItems;
			// ������Ԫ�ظ���
			TIndex _nUpItems;
			// ������ʵ��С
			std::pair<TIndex, TIndex> _mSize;
			// ��Ч�� [_vaildBegin, _vaildEnd)
			// ��Ч����ʼ��
			TIndex _vaildBegin;
			// ��Ч�н�����(����һ��)
			TIndex _vaildEnd;

			// ʵ���以Ϊ��Ԫ
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
				// ���ڶԽ�������������Ŀռ��˷�
				_diag = std::vector<TN>(rows);
				_upTria = std::vector<TN>(nUpItems);
				_columnIndices = std::vector<TIndex>(nUpItems);
				_rowOffsets = std::vector<TIndex>(rows + 1);
			}
			auto getUpItem(TIndex row, TIndex col) const
			{
				// �õ���ȷ����������
				if (row > col)
				{
					std::swap(row, col);
				}
				assert(row != col);
				// ������Ч����ֱ�ӷ���
				if (row >= _vaildEnd || row < _vaildBegin)
					return _columnIndices.end() - _columnIndices.begin();
				// ����Ч�з�Χ��
				const auto& iter0 = _columnIndices.begin();
				const auto& iter1 = iter0 + _rowOffsets[row];
				const auto& iter2 = iter0 + _rowOffsets[row + 1];
				// ʹ��Ч�ʸ��ߵĶ��ֲ���
				const auto& iter = std::lower_bound(iter1, iter2, col);
				return (iter == iter2 || (*iter) != col) ? _columnIndices.end() - iter0 : iter - iter0;
			}
		public:
			// ��ָ���������������
			template<typename TS, typename TN2>
			friend TS& operator << (TS& stream, const SpraseSM<TN2>& mat);
		public:
			// ��ȡϡ������С
			std::pair<TIndex, TIndex> GetSize() const { return _mSize; }
			// ����Ԫ��(Ref)����Ԫ�ز����ڣ����׳��쳣�������µ�Ԫ�ؽ�����֮ǰ������ʧЧ
			TN& AtRef(TIndex row, TIndex col) const
			{
				if (row == col)
				{
					// �Խ�Ԫ��
					return std::ref(_diag[row]);
				}
				else
				{
					// ������Ԫ��
					return std::ref(_upTria[getUpItem(row, col)]);
				}
			}
			// ����Ԫ��(Val)����Ԫ�ز����ڣ�������zeroVal
			TN AtVal(TIndex row, TIndex col) const
			{
				if (row == col)
				{
					// �Խ�Ԫ��
					return _diag[row];
				}
				else
				{
					// ������Ԫ��
					auto index = getUpItem(row, col);
					return (index - _nUpItems == 0) ? zeroVal : _upTria[index];
				}
			}
			// ��ʼ��ϡ�����-����Ԫ���б�
			void Init(std::list<Tuple<TIndex, TIndex, TN> >& diagMat, std::list<Tuple<TIndex, TIndex, TN> >& upMat, TIndex nRow, TIndex nCol)
			{
				assert(nRow == nCol);
				// ���¾����С
				_mSize.first = nRow;
				_mSize.second = nCol;
				// ��list��������
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
				// ͳ�ƾ���Ԫ�أ������·���洢�ռ�
				_nUpItems = upMat.size();
				_nItems = upMat.size() + diagMat.size();
				newMat(_nItems, _nUpItems, nRow);
				// ��ʼ��
				// �Խ�Ԫ��
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
				// �ǶԽ�Ԫ��
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
			// ��ʼ��ϡ�����-�ɳ��ܾ���(�ԳƷ���)
			void Init(Matrix<TN>& mat, TIndex nRow, TIndex nCol)
			{
				auto size = mat.GetSize();
				assert(size.first == size.second);
				nRow = std::max(size.first, nRow);
				nCol = std::max(size.second, nCol);
				assert(nRow == nCol);
				// ����list
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
				// ��ʼ��
				Init(diagMat, upMat, nRow, nCol);
			}

			// �������
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


		// ʮ������֧·�ڵ�(�����ӦԪ��)
		template<typename TL>
		struct ListNode
		{
			// ����ͼ֧·�������㡢��ʼ����
			const TIndex row, col;
			// ָ��ͷ��ͬ����һ��֧·��ָ��β��ͬ����һ��֧·
			std::shared_ptr<ListNode<TL> > down, right;
			// ֧·����
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

		// ʮ��������ڵ�(������������)
		template<typename TL>
		struct VexNode
		{
			// ָ���Ըö���Ϊ��ͷ����β�ĵ�һ��֧·
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

		// ��Ӧ��ϵ˵����
		// traiVex -> row
		// headVex -> col
		// hLink -> down
		// tLink -> right
		// firstIn -> cols[i]
		// firstOut -> rows[i]

		// �������ݽṹ��ϡ�����(��ͨ��ʮ������)
		// �ο�https://blog.csdn.net/bible_reader/article/details/71214096
		// ֧�ֲ����ɾ��Ԫ�أ������ڴ洢ϡ�跽���������ͼ
		template<typename TN>
		class SpraseOL : SpraseBase
		{
		private:
			// ����Ԫ��
			static constexpr typename zero_traits<TN>::zero_type zeroVal = zero_traits<TN>::zero_value;

			// �����б�
			std::vector<VexNode<TN> > _vexList;

			// ������ʵ��С
			std::pair<TIndex, TIndex> _mSize;
		private:
			// ��Ӷ���
			void insertVex()
			{
				_vexList.emplace_back(nullptr, nullptr);
			}
			// ɾ���ڵ�
			void deleteVex(TIndex vex)
			{

			}
		public:
			// ��ȡϡ������С
			std::pair<TIndex, TIndex> GetSize() const { return _mSize; }
			// ����Ԫ��(ָ��)
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
			// ����Ԫ��(����)
			TN& AtRef(TIndex row, TIndex col) const
			{
				auto cur = AtPtr(row, col);
				return std::ref(cur->info);
			}
			// ����Ԫ��(ֵ)
			TN AtVal(TIndex row, TIndex col) const
			{
				auto cur = AtPtr(row, col);
				return cur ? cur->info : zeroVal;
			}
			// ����Ԫ��(����֤���Ͼ���Ԫ��˳�򣬸���ԭֵ)
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
			// ����Ԫ��(��֤����Ԫ��˳�򣬸���)
			bool Update(TIndex row, TIndex col, const TN& value)
			{
				// ��ʱû���
				return true;
			}
			// ɾ��Ԫ��
			bool Delete(TIndex row, TIndex col)
			{
				// �Ͽ�������
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
				// �Ͽ�������
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
				// ���������Ͽ�������ָ���Զ�����
				auto result1 = deleteOfRow();
				auto result2 = deleteOfCol();
				assert(result1 == result2);
				return result1 && result2;
			}
			// �������
			void Clear()
			{
				_vexList.clear();
				_mSize.first = 0;
				_mSize.second = 0;
			}
			// ��ϵ�������������ʹ����Ͼ���Ĵ洢˳��
			void Sort()
			{
				// �������Ԫ��(ָ��)
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
				// ��������
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
				// ���²���
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
			// ��ʼ��ʮ������ϡ�����
			void Init(TIndex nRow, TIndex nCol)
			{
				auto nVex = std::max(nRow, nCol);
				_mSize.first = nRow;
				_mSize.second = nCol;
				_vexList = std::vector<VexNode<TN> >(nVex);
			}
			// ��ʼ��ʮ������ϡ�����(����Ԫ���б�ͬʱ����ظ�Ԫ��)
			void Init(std::list<Tuple<TIndex, TIndex, TN> >& mat, TIndex nRow, TIndex nCol)
			{
				// ��list��������(��->С)
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
				// ɾ���ظ�Ԫ��
				mat.unique(
					[](const Tuple<TIndex, TIndex, TN>& tuple1, const Tuple<TIndex, TIndex, TN>& tuple2) ->bool {
					return (tuple1.First() == tuple2.First()) && (tuple1.Second() == tuple2.Second());
				});
				// ��ʼ��
				Init(nRow, nCol);
				// ����
				for (auto& item : mat)
				{
					Insert(item.First(), item.Second(), item.Third());
				}
			}
			// ��ʼ��ʮ������ϡ�����(�ɳ��ܾ���)
			void Init(const Matrix<TN>& mat, TIndex nRow, TIndex nCol)
			{
				// ��ʼ��
				Init(nRow, nCol);
				// ����
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
			// ��ʼ��ʮ������ϡ�����(�ɶ�ά����)
			void Init(const std::vector<std::vector<TN> >& mat, TIndex nRow, TIndex nCol)
			{

			}
		public:
			// ��ָ���������������
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

		// �������ݽṹ��ϡ�����(���Խǣ��Խ�ѹ��)
		// ��Ҫ���ڱ�ʾ��ط������΢�ַ�����
		// ��֧�ֲ���
		template<typename TN>
		class SpraseTD : SpraseBase
		{
		private:
			// ����Ԫ��
			static constexpr typename zero_traits<TN>::zero_type zeroVal = zero_traits<TN>::zero_value;

			// ���Խ����Խ�����
			std::vector<double> _diagB;
			// ���Խ��´����Խ���
			std::vector<double> _diagA;
			// ���Խ��ϴ����Խ���
			std::vector<double> _diagC;

			// ������ʵ��С
			std::pair<TIndex, TIndex> _mSize;

			// ʵ���以Ϊ��Ԫ
			friend class SpraseTD<TN>;
		public:
			// ��ָ�����������
			template<typename TS, typename TN4, int z = 0>
			friend TS& operator << (TS& stream, const SpraseTD<TN4>& mat);
		public:
			// ��ȡϡ������С
			std::pair<TIndex, TIndex> GetSize() const { return _mSize; }
			// ����Ԫ��(Ref)����Ԫ�ز����ڣ����׳��쳣
			TN& AtRef(TIndex row, TIndex col) const
			{
				if (row == col)
				{
					// ���Խ�
					return std::ref(_diagB[col]);
				}
				else if(row-col == 1)
				{
					// �´����Խ�
					return std::ref(_diagA[col]);
				}
				else if (col - row == 1)
				{
					// �ϴ����Խ�
					return std::ref(_diagC[row]);
				}
				else
				{
					// Խ�磡
					return std::ref(_diagB[_mSize.first]);
				}
			}
			// ����Ԫ��(Val)����Ԫ�ز����ڣ�������zeroVal
			TN AtVal(TIndex row, TIndex col) const
			{
				if (row == col)
				{
					// ���Խ�
					return _diagB[col];
				}
				else if (row - col == 1)
				{
					// �´����Խ�
					return _diagA[col];
				}
				else if (col - row == 1)
				{
					// �ϴ����Խ�
					return _diagC[row];
				}
				else
				{
					// ����λ�ã�ӦΪ0Ԫ��
					return zeroVal;
				}
			}
			// �������Խ�(Ref, Offset)����Ԫ�ز����ڣ����׳��쳣
			TN& AtDiagBRef(TIndex offset)
			{
				return std::ref(_diagB[offset]);
			}
			// �����´ζԽ�(Ref, Offset)����Ԫ�ز����ڣ����׳��쳣
			TN& AtDiagARef(TIndex offset)
			{
				return std::ref(_diagA[offset]);
			}
			// �����ϴζԽ�(Ref, Offset)����Ԫ�ز����ڣ����׳��쳣
			TN& AtDiagCRef(TIndex offset)
			{
				return std::ref(_diagC[offset]);
			}
		public:
			// ��վ���
			void Clear()
			{
				_diagB.clear();
				_diagA.clear();
				_diagC.clear();
				_mSize.first = 0;
				_mSize.second = 0;
			}
			// (����)��ʼ������
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
			// ���캯��
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



		// ��˹��ȥ��
		// �ο������㷽����(���˳ɡ�÷��Ȫ������ѧ������)
		template<typename TM>
		std::vector<double> GESolver(TM&& A, std::vector<double>& b)
		{
			return GEKernel(std::forward<TM>(A), b);
		}
		// ��˹��ȥ�����ĺ��������ؽ������
		std::vector<double> GEKernel(Matrix<double>& A, std::vector<double>& b);

		// LU/LUP�ֽⷨ
		// �ο����㷨���� �ڶ��桷(Thomas H.Cormen�������˽����룬��е��ҵ������)
		template<typename TM>
		std::vector<double> LUSolver(TM&& A, std::vector<double>& b)
		{
			auto&& pi = LUDecomposition(std::forward<TM>(A));
			return LUSolution(std::forward<TM>(A), b, pi);
		}

		// LU�ֽ⣬����ת������
		std::vector<double> LUDecomposition(Matrix<double>& A);
		std::vector<double> LUDecomposition(SpraseOL<double>& A);
		// LU�ش������ؽ������
		std::vector<double> LUSolution(Matrix<double>& A, std::vector<double>& b, std::vector<double>& pi);
		std::vector<double> LUSolution(SpraseOL<double>& A, std::vector<double>& b, std::vector<double>& pi);

		// ��γ��ɳڵ�����/��˹���¶�������
		// �ο������㷽����(���˳ɡ�÷��Ȫ������ѧ������)
		template<typename TM>
		std::vector<double> SORSolver(TM&& A, std::vector<double>& b, std::vector<double>& xPre, double weight = 1)
		{
			return SORKernel(std::forward<TM>(A), b, xPre, weight);
		}
		// SOR����������
		constexpr double SORConvergeLimit = 1e-6;
		// ��γ��ɳڵ��������ĺ���
		std::vector<double> SORKernel(SpraseCSR<double>&A, std::vector<double>& b, std::vector<double>& xPre, double weight = 1.0);

		// �ſ˱ȵ�����
		// �ο������㷽����(���˳ɡ�÷��Ȫ������ѧ������)
		template<typename TM>
		std::vector<double> JACSolver(TM&& A, std::vector<double>& b)
		{
			return JACKernel(std::forward<TM>(A), b);
		}
		// �ſ˱ȵ�������������
		constexpr double JACConvergeLimit = 1e-6;
		// �ſ˱ȵ��������ĺ���
		std::vector<double> JACKernel(SpraseCSR<double>&A, std::vector<double>& b);

		// QR�ֽⷨ
		// �ο������㷽����(���˳ɡ�÷��Ȫ������ѧ������)
		template<typename TM>
		std::vector<double> QRSolver(TM&& A, std::vector<double>& b)
		{
			return QRKernel(std::forward<TM>(A), b);
		}
		// QR�ֽ⺯��
		std::pair<std::vector<double>, std::vector<double> > QRDecomposition(Matrix<double>&A, bool isExtension = false);
		// QR�����ĺ���
		std::vector<double> QRKernel(Matrix<double>&A, std::vector<double>& b);

		// �����ݶȷ�
		// �ο������㷽����(���˳ɡ�÷��Ȫ������ѧ������)
		template<typename TM>
		std::vector<double> CGSolver(TM&& A, std::vector<double>& b)
		{
			return CGKernel(std::forward<TM>(A), b);
		}
		// �����ݶȷ���������
		constexpr double CGConvergeLimit = 1e-6;
		// �����ݶȷ����ĺ���
		std::vector<double> CGKernel(Matrix<double>& A, std::vector<double>& b);

		// ׷�Ϸ�
		// �ο������㷽����(���˳ɡ�÷��Ȫ������ѧ������)
		template<typename TM>
		std::vector<double> CMSolver(TM&& A, std::vector<double>& d)
		{
			return CMKernel(std::forward<TM>(A), d);
		}
		// ׷�Ϸ����ĺ���
		std::vector<double> CMKernel(SpraseTD<double>& A, std::vector<double>& d);
	}
}

# endif