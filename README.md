## 数值分析
### author: Katyusha

帮助理解算法的一套数值分析工具库，代码未经优化，仅供学习

已经移除了自定义的矩阵运算工具Sprase，转为使用Eigen
Eigen的使用详见官网：http://eigen.tuxfamily.org/index.php?title=Main_Page

结构如下：

- **Numerical**：数值运算
   - 基础数据结构：
      - 区间(只允许是整型/浮点类型)Interval
      - 分段函数SegmentFunction
   - 算法：
      - 三次样条插值法CubicSplineInterpolation
      - 龙贝格积分法RombergIntegration
	  - 最优平方逼近(最小二乘法)LeastSquareFitting/MutliLinearLeastSquareFitting
	  - 最优一致逼近(Remes算法)BestUniformApproximation

- **Symbolic**：符号数学工具
   - 基础数据结构：
      - 表达式树ExprTree   

- **Optimalize**：优化求解工具
   - 基础数据结构：
      - NULL  
   - 算法：
      - LP(等式约束)-原始对偶内点法
	  - LP(不等式约束)-障碍函数内点法
	  - QP(无约束)-共轭梯度法
	  - QP(等式约束)-拉格朗日法
	  - QP(不等式约束)-有效集法





