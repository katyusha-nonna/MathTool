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

// 包含如下基础数据结构
// 1 表达式树节点
// 2 表达式树


#ifndef MYSYMBOLIC_H
#define MYSYMBOLIC_H

namespace Utility
{
	namespace Symbolic
	{
		// 基础数据结构-表达式树节点
		template<typename ET>
		struct TreeNode
		{
			// 元素名称
			std::string name;
			// 元素对象指针(操作符为空指针)
			ET* element;
			// 左子树根节点指针(叶节点为空指针)
			TreeNode<ET>* leftNode;
			// 右子树根节点指针(叶节点为空指针)
			TreeNode<ET>* rightNode;

			// 判断是否为树叶
			bool IsLeave()
			{
				return (leftNode == nullptr) && (rightNode == nullptr);
			}

			// 构造函数
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
			// 析构函数
			~TreeNode()
			{
				delete element;
			}
		};

		// 基础数据结构-表达式树：用于存储简单发电系统/互联系统的网络拓扑结构
		// 其中：ET表示表达式中参与计算的元素的类型
		// 为了统一，输入采用字符串类型，运算函数通过键值来索引
		// 允许使用括号(括号不检测括号平衡)，括号运算符不用自定义且优先级高于任何自定义操作符
		// 中缀表达式的构建和输出参考了https://www.cnblogs.com/jinks/archive/2013/04/28/3048990.html的实现 (作者：王锦，邮箱：jinksw@vip.qq.com)
		template<typename ET>
		class ExprTree
		{
		private:
			// 整个树的根节点
			TreeNode<ET>* root;
			// 整个树中所有操作数
			std::map<std::string, ET*> data;
			// 整个树中所有操作符实现
			std::map<std::string, std::function<ET(ET, ET)> > operators;
			// 整个表达式树的谓词函数优先级比较函数(判断优先级A>=B)
			std::function<bool(std::string, std::string)> comparer;
		private:
			// 判断是否为操作符
			bool IsOperator(std::string& curName)
			{
				return operators.find(curName) != operators.end();
			}
			// 检查括号平衡
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
			// 从某个节点开始递归计算
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
			// 仅清理表达式树
			void ClearTree()
			{
				if (root)
				{
					RecursionDelete(root);
					root = nullptr;
				}
			}
			// 完全清理表达式树
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
			// 递归删除以某个节点为root的子树
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
			// 增加谓词函数
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
			// 删除谓词函数
			void DeleteOperator(std::string name)
			{
				operators.erase(name);
			}
			void DeleteOperator()
			{
				operators.clear();
			}
			// 增加优先级比较函数
			void AddComp(std::function<bool(std::string, std::string)> comp)
			{
				comparer = comp;
			}
			// 增加操作数
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
			// 删除操作数
			void DeleteData(std::string name)
			{
				data.erase(name);
			}
			void DeleteData()
			{
				data.clear();
			}
			// 由中缀表达式构造表达式树(允许采用括号)
			void BuildFromInfixExpr(std::stringstream& input)
			{
				// 首先检查括号平衡
				if (!IsBracketBalance(input))
				{
					return;
				}
				else { ; }

				// 清除原先的表达式树
				ClearTree();
				// 用于构造表达式树的操作符栈(包含了括号操作符)
				std::stack<std::string> opStack;
				// 用于构造表达式树的操作数栈
				std::stack<TreeNode<ET>* > dataStack;
				// 当前的操作数/操作符
				std::string curString;
				// 循环遍历
				while (input >> curString)
				{
					if (operators.find(curString) != operators.end() || curString == "(" || curString == ")")
					{
						// 说明是操作符
						if (curString == "(")
						{
							// 左括号直接入栈
							opStack.push(curString);
						}
						else if (curString == ")")
						{
							// 右括号和上一个左括号对应
							// 为了防止冗赘括号，之前已检查强制保证括号平衡
							while (!dataStack.empty() && opStack.top() != "(")
							{
								// 从操作数栈取出两个操作数
								TreeNode<ET>* secondOpd = dataStack.top();
								dataStack.pop();
								TreeNode<ET>* firstOpd = dataStack.top();
								dataStack.pop();
								// 将操作数和操作符组成一个新结点存入栈中
								dataStack.push(new TreeNode<ET>(opStack.top(), nullptr, firstOpd, secondOpd));
								opStack.pop();
							}
							// 将左括号出栈
							opStack.pop();
						}
						else
						{
							//如果栈顶操作符优先级高于读入操作符优先级，则表名应该先计算栈顶操作符
							while (!opStack.empty() && comparer(opStack.top(), curString))
							{
								TreeNode<ET>* secondOpd = dataStack.top();
								dataStack.pop();
								TreeNode<ET>* firstOpd = dataStack.top();
								// 从操作数栈取出两个操作数
								dataStack.pop();
								// 将操作数和操作符组成一个新结点存入栈中
								dataStack.push(new TreeNode<ET>(opStack.top(), nullptr, firstOpd, secondOpd));
								opStack.pop();
							}
							// 将读入操作符入栈
							opStack.push(curString);
						}
					}
					else
					{
						// 说明是操作数
						auto element = data.at(curString);
						dataStack.push(new TreeNode<ET>(curString, element));
					}
				}
				curString = "";
				// 单独处理根节点
				while (!opStack.empty() && comparer(opStack.top(), curString))
				{
					TreeNode<ET>* secondOpd = dataStack.top();
					dataStack.pop();
					TreeNode<ET>* firstOpd = dataStack.top();
					dataStack.pop();
					dataStack.push(new TreeNode<ET>(opStack.top(), nullptr, firstOpd, secondOpd));
					opStack.pop();
				}
				// 将唯一的剩余的节点设为根节点
				root = dataStack.top();
				dataStack.pop();
			}
			// 输出中缀表达式
			template<typename TS>
			void PrintToInfixExpr(TS& out)
			{
				//结点栈，遍历使用
				std::stack<TreeNode<ET>* > nodeStack;
				// 当前遍历节点初始化为表达式树根节点
				TreeNode<ET>* pointer = root;
				// 用于记录在哪些元素被输出之后要输出反括号
				std::list<TreeNode<ET>* > nodeList;

				// lambda对象：判断是否需要添加左括号保证运算顺序的正确性
				auto shouldPrintLeftBracket = [&](bool isLeft)->bool
				{
					if (nodeStack.empty())
						return false;
					if (pointer == nullptr)
						return false;
					auto a = nodeStack.top()->name;
					auto b = pointer->name;
					//如果是数字，则不用打括号
					if (operators.find(b) == operators.end())
						return false;
					if (isLeft)
					{
						// 如果pointer是左结点
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
						// 如果pointer是右结点
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
					// 一直向左子结点走，找到左子结点时，经过的结点已全部入栈
					while (pointer != nullptr)
					{
						// 如果应该添加左括号,为左子结点的情况下
						if (shouldPrintLeftBracket(true))
						{
							// 找到应该在输出哪个结点后输出右括号
							auto temp = pointer->rightNode;
							while (temp->rightNode != nullptr)
							{
								temp = temp->rightNode;
							}
							// 待插入右括号的节点位置入栈
							nodeList.push_back(temp);
							out << " ( ";
						}
						nodeStack.push(pointer);
						pointer = pointer->leftNode;
					}
					// 输出结点 
					out << " " + nodeStack.top()->name + " ";
					auto it = std::find(nodeList.begin(), nodeList.end(), nodeStack.top());
					// 若栈顶结点(即当前输出结点)为前面记录的应该输出右括号的结点时，输出右括号，可能有多个右括号需要输出，所以循环
					while (it != nodeList.end())
					{
						out << " ) ";
						nodeList.erase(it);
						it = std::find(nodeList.begin(), nodeList.end(), nodeStack.top());
					}
					// 所有左结点已经走完，走右节点
					pointer = nodeStack.top()->rightNode;
					// 如果应该添加左括号,为右子结点的情况下
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
			// 求解表达式
			ET Calcute()
			{
				// 递归求解
				return std::move(calcuteFromNode(root));
			}
		public:
			// 析构函数
			~ExprTree()
			{
				AllClear(true, true);
			}
		};
	}
}



#endif