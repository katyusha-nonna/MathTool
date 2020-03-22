#pragma once
#include <algorithm>
#include <memory>
#include <stack>
#include <map>
#include <list>
#include <functional>
#include <string>
#include <sstream>
#include <iostream>

// author: Katyusha
// date: 2020/01/17

// �������»������ݽṹ
// 1 ���ʽ���ڵ�
// 2 ���ʽ��


#ifndef MYSYMBOLIC_H
#define MYSYMBOLIC_H

namespace Utility
{
	namespace Symbolic
	{
		// �������ݽṹ-���ʽ���ڵ�
		template<typename ET>
		struct TreeNode
		{
			// Ԫ������
			std::string name;
			// Ԫ�ض���ָ��(������Ϊ��ָ��)
			ET* element;
			// ���������ڵ�ָ��(Ҷ�ڵ�Ϊ��ָ��)
			TreeNode<ET>* leftNode;
			// ���������ڵ�ָ��(Ҷ�ڵ�Ϊ��ָ��)
			TreeNode<ET>* rightNode;

			// �ж��Ƿ�Ϊ��Ҷ
			bool IsLeave()
			{
				return (leftNode == nullptr) && (rightNode == nullptr);
			}

			// ���캯��
			TreeNode(std::string n, ET* e, TreeNode<ET>* l, TreeNode<ET>* r)
			{
				name = n;
				element = e;
				leftNode = l;
				rightNode = r;
			}
			TreeNode(std::string n, ET* e)
			{
				name = n;
				element = e;
				leftNode = nullptr;
				rightNode = nullptr;
			}
			TreeNode()
			{
				element = nullptr;
				leftNode = nullptr;
				rightNode = nullptr;
			}
			// ��������
			~TreeNode()
			{
				delete element;
			}
		};

		// �������ݽṹ-���ʽ�������ڴ洢�򵥷���ϵͳ/����ϵͳ���������˽ṹ
		// ���У�ET��ʾ���ʽ�в�������Ԫ�ص�����
		// Ϊ��ͳһ����������ַ������ͣ����㺯��ͨ����ֵ������
		// ����ʹ������(���Ų��������ƽ��)����������������Զ��������ȼ������κ��Զ��������
		// ��׺���ʽ�Ĺ���������ο���https://www.cnblogs.com/jinks/archive/2013/04/28/3048990.html��ʵ�� (���ߣ����������䣺jinksw@vip.qq.com)
		template<typename ET>
		class ExprTree
		{
		private:
			// �������ĸ��ڵ�
			TreeNode<ET>* root;
			// �����������в�����
			std::map<std::string, ET*> data;
			// �����������в�����ʵ��
			std::map<std::string, std::function<ET(ET, ET)> > operators;
			// �������ʽ����ν�ʺ������ȼ��ȽϺ���(�ж����ȼ�A>=B)
			std::function<bool(std::string, std::string)> comparer;
		private:
			// �ж��Ƿ�Ϊ������
			bool IsOperator(std::string& curName)
			{
				return operators.find(curName) != operators.end();
			}
			// �������ƽ��
			bool IsBracketBalance(std::stringstream& input)
			{
				std::string save;
				save = input.str();
				std::string temps;
				std::stack<bool> bracketStack;
				while (input >> temps)
				{
					if (temps == "(")
					{
						bracketStack.push(true);
					}
					else if (temps == ")")
					{
						bracketStack.pop();
					}
					else { ; }
				}
				input.clear();
				input.str("");
				input << save;
				return bracketStack.empty();
			}
			// ��ĳ���ڵ㿪ʼ�ݹ����
			ET calcuteFromNode(TreeNode<ET>* curRoot)
			{
				if (IsOperator(curRoot->name))
				{
					return std::move(operators.at(curRoot->name)(calcuteFromNode(curRoot->leftNode), calcuteFromNode(curRoot->rightNode)));
				}
				else
				{
					return std::move(*(curRoot->element));
				}
			}
		public:
			// ��������ʽ��
			void ClearTree()
			{
				if (root)
				{
					RecursionDelete(root);
					root = nullptr;
				}
			}
			// ��ȫ������ʽ��
			void AllClear(bool clearNode = true, bool clearOp = true)
			{
				if (clearNode)
				{
					ClearTree();
					data.clear();
				}
				if (clearOp)
				{
					operators.clear();
				}
			}
			// �ݹ�ɾ����ĳ���ڵ�Ϊroot������
			void RecursionDelete(TreeNode<ET>* rootForDel)
			{
				if (rootForDel->leftNode)
				{
					RecursionDelete(rootForDel->leftNode);
				}
				if (rootForDel->rightNode)
				{
					RecursionDelete(rootForDel->rightNode);
				}
				delete rootForDel;
			}
			// ����ν�ʺ���
			void AddOperator(std::string name, std::function<ET(ET, ET)> op)
			{
				operators.try_emplace(name, op);
			}
			void AddOperator(const std::map<std::string, std::function<ET(ET, ET)> >& ops)
			{
				operators = ops;
			}
			void AddOperator(std::map<std::string, std::function<ET(ET, ET)> >&& ops)
			{
				operators = ops;
			}
			// ɾ��ν�ʺ���
			void DeleteOperator(std::string name)
			{
				operators.erase(name);
			}
			void DeleteOperator()
			{
				operators.clear();
			}
			// �������ȼ��ȽϺ���
			void AddComp(std::function<bool(std::string, std::string)> comp)
			{
				comparer = comp;
			}
			// ���Ӳ�����
			void AddData(std::string name, ET* value)
			{
				data.try_emplace(name, value);
			}
			void AddData(const std::map<std::string, ET*>& values)
			{
				data = values;
			}
			void AddData(std::map<std::string, ET*>&& values)
			{
				data = values;
			}
			// ɾ��������
			void DeleteData(std::string name)
			{
				data.erase(name);
			}
			void DeleteData()
			{
				data.clear();
			}
			// ����׺���ʽ������ʽ��(�����������)
			void BuildFromInfixExpr(std::stringstream& input)
			{
				// ���ȼ������ƽ��
				if (!IsBracketBalance(input))
				{
					return;
				}
				else { ; }

				// ���ԭ�ȵı��ʽ��
				ClearTree();
				// ���ڹ�����ʽ���Ĳ�����ջ(���������Ų�����)
				std::stack<std::string> opStack;
				// ���ڹ�����ʽ���Ĳ�����ջ
				std::stack<TreeNode<ET>* > dataStack;
				// ��ǰ�Ĳ�����/������
				std::string curString;
				// ѭ������
				while (input >> curString)
				{
					if (operators.find(curString) != operators.end() || curString == "(" || curString == ")")
					{
						// ˵���ǲ�����
						if (curString == "(")
						{
							// ������ֱ����ջ
							opStack.push(curString);
						}
						else if (curString == ")")
						{
							// �����ź���һ�������Ŷ�Ӧ
							// Ϊ�˷�ֹ��׸���ţ�֮ǰ�Ѽ��ǿ�Ʊ�֤����ƽ��
							while (!dataStack.empty() && opStack.top() != "(")
							{
								// �Ӳ�����ջȡ������������
								TreeNode<ET>* secondOpd = dataStack.top();
								dataStack.pop();
								TreeNode<ET>* firstOpd = dataStack.top();
								dataStack.pop();
								// ���������Ͳ��������һ���½�����ջ��
								dataStack.push(new TreeNode<ET>(opStack.top(), nullptr, firstOpd, secondOpd));
								opStack.pop();
							}
							// �������ų�ջ
							opStack.pop();
						}
						else
						{
							//���ջ�����������ȼ����ڶ�����������ȼ��������Ӧ���ȼ���ջ��������
							while (!opStack.empty() && comparer(opStack.top(), curString))
							{
								TreeNode<ET>* secondOpd = dataStack.top();
								dataStack.pop();
								TreeNode<ET>* firstOpd = dataStack.top();
								// �Ӳ�����ջȡ������������
								dataStack.pop();
								// ���������Ͳ��������һ���½�����ջ��
								dataStack.push(new TreeNode<ET>(opStack.top(), nullptr, firstOpd, secondOpd));
								opStack.pop();
							}
							// �������������ջ
							opStack.push(curString);
						}
					}
					else
					{
						// ˵���ǲ�����
						auto element = data.at(curString);
						dataStack.push(new TreeNode<ET>(curString, element));
					}
				}
				curString = "";
				// ����������ڵ�
				while (!opStack.empty() && comparer(opStack.top(), curString))
				{
					TreeNode<ET>* secondOpd = dataStack.top();
					dataStack.pop();
					TreeNode<ET>* firstOpd = dataStack.top();
					dataStack.pop();
					dataStack.push(new TreeNode<ET>(opStack.top(), nullptr, firstOpd, secondOpd));
					opStack.pop();
				}
				// ��Ψһ��ʣ��Ľڵ���Ϊ���ڵ�
				root = dataStack.top();
				dataStack.pop();
			}
			// �����׺���ʽ
			template<typename TS>
			void PrintToInfixExpr(TS& out)
			{
				//���ջ������ʹ��
				std::stack<TreeNode<ET>* > nodeStack;
				// ��ǰ�����ڵ��ʼ��Ϊ���ʽ�����ڵ�
				TreeNode<ET>* pointer = root;
				// ���ڼ�¼����ЩԪ�ر����֮��Ҫ���������
				std::list<TreeNode<ET>* > nodeList;

				// lambda�����ж��Ƿ���Ҫ��������ű�֤����˳�����ȷ��
				auto shouldPrintLeftBracket = [&](bool isLeft)->bool
				{
					if (nodeStack.empty())
						return false;
					if (pointer == nullptr)
						return false;
					auto a = nodeStack.top()->name;
					auto b = pointer->name;
					//��������֣����ô�����
					if (operators.find(b) == operators.end())
						return false;
					if (isLeft)
					{
						// ���pointer������
						if (this->comparer(b, a))
						{
							return false;
						}
						else
						{
							return true;
						}
					}
					else
					{
						// ���pointer���ҽ��
						if (this->comparer(a, b))
						{
							return true;
						}
						else
						{
							return false;
						}
					}
					return false;
				};

				while (!nodeStack.empty() || pointer != nullptr)
				{
					// һֱ�����ӽ���ߣ��ҵ����ӽ��ʱ�������Ľ����ȫ����ջ
					while (pointer != nullptr)
					{
						// ���Ӧ�����������,Ϊ���ӽ��������
						if (shouldPrintLeftBracket(true))
						{
							// �ҵ�Ӧ��������ĸ��������������
							auto temp = pointer->rightNode;
							while (temp->rightNode != nullptr)
							{
								temp = temp->rightNode;
							}
							// �����������ŵĽڵ�λ����ջ
							nodeList.push_back(temp);
							out << " ( ";
						}
						nodeStack.push(pointer);
						pointer = pointer->leftNode;
					}
					// ������ 
					out << " " + nodeStack.top()->name + " ";
					auto it = std::find(nodeList.begin(), nodeList.end(), nodeStack.top());
					// ��ջ�����(����ǰ������)Ϊǰ���¼��Ӧ����������ŵĽ��ʱ����������ţ������ж����������Ҫ���������ѭ��
					while (it != nodeList.end())
					{
						out << " ) ";
						nodeList.erase(it);
						it = std::find(nodeList.begin(), nodeList.end(), nodeStack.top());
					}
					// ���������Ѿ����꣬���ҽڵ�
					pointer = nodeStack.top()->rightNode;
					// ���Ӧ�����������,Ϊ���ӽ��������
					if (shouldPrintLeftBracket(false))
					{
						auto temp = pointer->rightNode;
						while (temp->rightNode != nullptr)
						{
							temp = temp->rightNode;
						}
						nodeList.push_back(temp);
						out << " ( ";
					}
					nodeStack.pop();
				}
				out << std::endl;
			}
			// �����ʽ
			ET Calcute()
			{
				// �ݹ����
				return std::move(calcuteFromNode(root));
			}
		public:
			// ��������
			~ExprTree()
			{
				AllClear(true, true);
			}
		};
	}
}



#endif