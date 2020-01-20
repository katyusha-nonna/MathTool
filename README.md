## 数值分析
### author: Katyusha

沙雕玩具。。。只是个玩具，帮助理解书上的伪代码，不要尝试使用

结构如下：

- **Sprase**：矩阵运算&线性方程组求解
   - 基础数据结构：
      - 三元组Tuple/复数Complex
      - 稠密矩阵Matrix(对二维数组的std::vector的直接封装)
      - 稀疏矩阵(普通，行压缩)SpraseCSR
      - 稀疏矩阵(对称，行压缩)SpraseSM
      - 稀疏矩阵(普通，十字链表)SpraseOL
      - 稀疏矩阵(三对角，对角压缩)SpraseTD
      - 零值类型萃取zero_traits
   - 算法：
      - 线性方程组求解(高斯消元，列主元)GESolver
      - 线性方程组求解(LU分解，杜立特尔，列主元)LUSolver
      - 线性方程组求解(逐次超松弛迭代法/高斯-赛德尔迭代法，稀疏，不选主元)SORSolver
      - 线性方程组求解(QR分解，豪斯霍尔德变换，不选主元)QRSolver
      - 线性方程组求解(共轭梯度法，对称正定，不选主元)CGSolver
      - 线性方程组求解(追赶法，三对角矩阵，不选主元)CMSolver

- **Numerical**：数值积分、数值微分与非线性方程求解
   - 基础数据结构：
      - 区间(只允许是整型/浮点类型)Interval
      - 分段函数SegmentFunction
   - 算法：
      - 三次样条插值法CubicSplineInterpolation
      - 龙贝格积分法RombergIntegration

- **Symbolic**：符号数学工具
   - 基础数据结构：
      - 表达式树ExprTree   





