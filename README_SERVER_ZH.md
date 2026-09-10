# IPM 长时间服务器运行说明

这个发布包面向“从原始物理时间 `t=0` 自动运行到长时间”的 Profile 实验。数值求解仍只有一个入口 `ipm.solve`；服务器脚本只是把少量常用设置翻译成完整且经过校验的配置。

## 最快启动

需要 MATLAB，建议使用与本机验证相同的 R2026a。解压后进入发布目录：

```bash
cd ipm_long_time_server_20260910
chmod +x server/launch.sh
nohup server/launch.sh > launcher.out 2>&1 &
```

如果 MATLAB 不在 `PATH`：

```bash
MATLAB_BIN=/opt/MATLAB/R2026a/bin/matlab nohup server/launch.sh > launcher.out 2>&1 &
```

正式长跑建议放在 `tmux` 或批处理系统中。启动前请确认磁盘空间和内存；当前分层策略允许最多 310000 个二维网格节点，稀疏 LU 的实际内存还取决于节点形状和平台。

## 用户接口

新算例只需编辑 `server/profile_settings.m`。默认设置是小箱 H8、原始 `321×161` 均匀网格、目标 canonical 时间 `τ=16`、自动 checkpoint、无中途人工换网格：

```matlab
settings.boxHalfWidth = 8;
settings.boxHeight = 4;
settings.canonicalFinalTime = 16;
settings.maximumSteps = 60000;
settings.maximumAdjacentGridRatio = 2;
settings.outputDirectory = fullfile(projectRoot,'runs','profile_H8_tau16');
settings.restartCheckpoint = '';
```

`maximumAdjacentGridRatio=2` 对应内部的 `remeshMaximumCellRatio=2`。在每次记录状态上，若 X、Y 两轴的最大相邻步长比都不超过 `2*(1+1e-10)`，程序不会因为 `grid_smoothness_failure` 停止。质量守恒、最大值原理、振荡、分辨率、时间终点、最大步数、自动网格候选耗尽等独立停止条件仍然有效；这样不会用一个网格比选项掩盖数值失效。

默认自动网格是 version 4：先在原节点规模内作分层参数搜索，仍不足时按 X/Y 方向分别增加节点，节点因子为 `[1,1]`、`[2,1]`、`[1,2]`、`[2,2]`、`[3,1]`、`[1,3]`、`[3,2]`、`[2,3]`。候选必须通过原来的核心分辨率、前沿、局部平滑、求积、迁移峰跳、质量和值域门才会提交。

也可以直接在 MATLAB 中使用配置接口：

```matlab
addpath('/absolute/path/to/ipm_long_time_server_20260910');
settings = struct( ...
    'canonicalFinalTime',16, ...
    'maximumAdjacentGridRatio',2, ...
    'outputDirectory','/data/ipm/run01');
opts = ipm.config.longTimeProfile(settings);
config = ipm.config.resolve(opts);  % 只检查，不启动计算
result = ipm.solve(opts);           % 从物理 t=0 启动
```

## 断点续算

程序按 `checkpointEvery` 保存不可覆盖的原生 checkpoint，并在正常或安全停止时保存末态 checkpoint。把 `profile_settings.m` 中的 `restartCheckpoint` 设为该文件的绝对路径，同时只延长 `canonicalFinalTime` 或 `maximumSteps`，再运行同一启动命令。

恢复会冻结原来的网格、物理、离散、C 规则和自动网格策略。尤其不能在恢复时把网格比从 1.5 改成 2；若需要 `2`，应在新的 `t=0` 算例启动前设置。

## 输出与检查

输出目录包含：

- `console.log`：完整控制台记录；
- `launch_configuration.mat`：用户设置、平铺选项和最终冻结配置；
- `checkpoint_*_step*.mat` 与 `checkpoint_manifest.jsonl`：可信接受步；
- `result*.mat` 与 `manifest.jsonl`：末态结果；
- `launcher_result_pointer.mat`：启动器返回的结果结构。

查看结果：

```matlab
addpath('/absolute/path/to/ipm_long_time_server_20260910');
r = ipm.output.validate('/data/ipm/run01/result_CASE_ID.mat');
ipm.output.plotResult(r);
```

安装后的接口检查：

```matlab
addpath('/absolute/path/to/ipm_long_time_server_20260910');
addpath('/absolute/path/to/ipm_long_time_server_20260910/tests');
r = ipmtests.baseline.serverInterface();
```

完整验证使用 `ipm.verify('all','quick')`。新旧一致性检查需要另行提供旧版目录，服务器发布包本身不依赖旧版。

## 已验证范围

version 4 的分层/方向增点实现已通过 baseline、四阶、六阶和新旧一致性完整回归；真实旧轨道的冻结场迁移也通过原生配对检查。此前 version 2 从 `t=0` 的 H8 已运行到 `τ=11.6731`、H64 已运行到 `τ=9.9086`，均因有限候选族耗尽停止；version 4 正是针对这种“仍有可用网格但搜索族未覆盖”的故障。

version 4 尚未完成从 `t=0` 到 `τ=16` 的整段服务器实跑，因此本包是可运行的研究版本，不代表已经获得无穷时间极限、空间收敛或奇性证明。运行中不要编辑源码或 `profile_settings.m`；需要新配置时启动新的输出目录。
