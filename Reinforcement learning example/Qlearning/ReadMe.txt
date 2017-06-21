Qdemo演示程序

Qlearning Q学习主程序
调用 drnd（随机变量生成函数）

robot
alpha gamma 学习参数
best 最优路径
state 状态

task
initialState  初始状态
terminalState 终止状态

Rew 环境给予的回报
N 学习次数

注意：任务改变时，要设当改变execut子函数和一些脚标变换函数。
      用于打印状态的statements也要改一下。