# 研究与历史资料

这里存放实验配置、初值、分析脚本及历史证据；求解仍统一调用 `ipm.solve(opts)`。
运行模块不依赖本目录。普通使用者从根部 [README](../README.md) 的短算例开始即可。

| 目录 | 用途 |
|---|---|
| [experiments](experiments/README.md) | 重网格、CCF 尾部、Bessel 双尺度等实验配置与初值 |
| [analysis](analysis) | 对结果作拟合、比较、可视化，以及保留的推导和报告源码 |
| [notes](notes) | 保存原版技术文档/变更历史，并以页首注记标明当前契约边界 |
| [output/pdf](output/pdf) | 原版已跟踪的研究报告，不是新计算的输出目录 |

在已经添加本副本根目录的 MATLAB 会话中，需要运行研究代码时显式添加：

```matlab
addpath(fullfile(projectRoot,'research','experiments'));
addpath(fullfile(projectRoot,'research','analysis'));
opts = strong_mesh_case('production');
% 先检查 opts 的计算规模、终止时间及输出路径，再调用 ipm.solve(opts)。
```

`projectRoot` 是本副本根目录的绝对路径。研究算例可能启用长时间计算、保存和视频，
不是安装检查；相对输出路径仍相对于 MATLAB 当前工作目录。
新增试验产物应明确放入本副本 `result/` 下，不能覆盖支撑已有结论的数据。

原版技术文档入口：

- [原使用说明](notes/README_IPM.md)、[原架构说明](notes/ARCHITECTURE.md)、[原变更历史](notes/CHANGELOG.md)。
- [四阶证据及限制](notes/HIGH_ORDER_RESEARCH.md)、[六阶证据及限制](notes/SIXTH_ORDER_RESEARCH.md)。

上述原文和 `analysis/` 中的旧笔记保留历史函数名、目录名和结论日期；
页首或局部 retired 注记是后加的当前契约界标。当前结构以
[STRUCTURE.md](../STRUCTURE.md) 为准，旧函数名按 [迁移清单](../tests/source_map.json) 查找。
历史记录不能代替本副本的重新验证，研究计划也不等于已实现的能力或新的执行指令。
归档报告依赖历史图表；批量计算产物没有复制，因此不保证直接编译历史报告。
