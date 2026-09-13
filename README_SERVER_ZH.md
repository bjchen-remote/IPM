# IPM 长时间服务器运行说明

这个发布包面向“从原始物理时间 `t=0` 自动运行到长时间”的 Profile 实验。数值求解仍只有一个入口 `ipm.solve`；服务器脚本只是把少量常用设置翻译成完整且经过校验的配置。

冻结态的方向网格成本前沿和已知停机原因见 [网格建议](MESH_GRID_RECOMMENDATION_ZH.md)。
本次实验包还把连续 X 峰值的竖向核心影子观测写入结果/断点元数据，
用于评估正式网格触发的误差；它不改变正式网格决策或 C 规则。
网格问题和新证据见 [研究说明](research/longtime_lab/AUTONOMOUS_GRID_DESIGN_FROM_FROZEN_DATA_20260913.md)。

## 最快启动

需要 MATLAB，建议使用与本机验证相同的 R2026a。解压后进入发布目录：

```bash
cd ipm_long_time_server_v5_tau14_growth_20260913
chmod +x server/launch.sh
nohup server/launch.sh > launcher.out 2>&1 &
```

如果 MATLAB 不在 `PATH`：

```bash
MATLAB_BIN=/opt/MATLAB/R2026a/bin/matlab nohup server/launch.sh > launcher.out 2>&1 &
```

正式长跑建议放在 `tmux` 或批处理系统中。启动前请确认磁盘空间和内存；`maximumTotalNodes` 是启动时登记的二维节点预算，并不是稀疏 LU 的内存保证。更高预算必须按服务器实际内存试验。

## 用户接口

新算例只需编辑 `server/profile_settings.m`。默认设置是小箱 H8、启动前自动选出的原始 `321×161` 均匀网格、自动 checkpoint、无中途人工换网格。时间和步数设成较远的**运行上限**，避免新服务器作业仅因达到旧 `τ=16` 里程碑而停止：

```matlab
settings.boxHalfWidth = 8;
settings.boxHeight = 4;
settings.initialNodeCount = 'auto';
settings.canonicalFinalTime = 1000;
settings.maximumSteps = 10000000;
settings.maximumAdjacentGridRatio = 2;
settings.autonomousMeshVersion = 5;
settings.maximumTotalNodes = 310000;
settings.outputDirectory = fullfile(projectRoot,'runs','profile_H8_long');
settings.restartCheckpoint = '';
```

`maximumAdjacentGridRatio=2` 对应内部的 `remeshMaximumCellRatio=2`。在每次记录状态上，若 X、Y 两轴的最大相邻步长比都不超过 `2*(1+1e-10)`，程序不会因为 `grid_smoothness_failure` 停止。质量守恒、最大值原理、振荡、分辨率、时间终点、最大步数、自动网格候选耗尽等独立停止条件仍然有效；这样不会用一个网格比选项掩盖数值失效。

`τ=1000` 和一千万步不是已验证可达的时间或奇异性判据，只是避免**人为的短时间上限**；作业可能更早因严格的数值安全门或注册节点预算停止。请根据服务器的磁盘、内存和预计运行时间在启动前调整这些上限，并用 checkpoint 检查进展。当前 H8 的远边界告警意味着长跑主要检验网格和时间推进，最终 Profile 还必须做大箱与时空误差比较。

`initialNodeCount='auto'` 对默认 H8 使用已验证的 `321×161` 起点；对更高的计算箱，在构造任何求解器 LU 之前，用原始解析 `t=0` 数据和原质量门按已登记的二维节点成本搜索初始 X/Y 数量，并把选中的数量冻结进配置和 checkpoint。实测 H64/H128 选 `321×161`；H256/Y128 的固定 `321×161` 在初始 X 候选族耗尽，自动模式选 `401×161`，随后原生初始化及短程推进通过。可改为 `[Nx,Ny]` 明确固定初始节点；旧 `ipm.config.longTimeProfile` 未指定此项时保持原 numeric 默认。自动选择不保证稀疏 LU 内存、未来长时容量或计算域收敛。

服务器设置文件显式选 version 5。它依据初始 `Nx×Ny` 和启动时的 `maximumTotalNodes`，一次登记预算内的方向节点级；每轴先用整数因子 1–6，之后按至多约 4/3 的因子增长，直到预算边界。触发重网格时按节点总成本由小到大尝试，每个节点级内从已注册候选中选择最大相邻格宽比最小的合格轴对，并保留质量裕量作为次级排序。相同参考 X/Y 轴的几何试探在单次规划内复用。实际迁移仍须通过核心、前沿、模板、求积、峰跳、质量和值域门；失败不能改变已接受状态。预算 310000、初始 `321×161` 时登记 14 个节点级；改为 1000000 时登记 38 个。选更大预算只需在**新** `t=0` 算例启动前修改设置，不能在 checkpoint 续算时修改。旧 H8 轨道曾在约 τ=11.67 因固定节点族耗尽；若服务器内存允许，建议在新长跑开始前把预算提高到 1000000 并先做短程资源试跑。

`ipm.config.longTimeProfile` 在未指定 `autonomousMeshVersion` 时仍默认 version 4，以便已有脚本配置逐值不变；它的注册节点族固定到 310000。version 5 是新增的自动扩容实验路径，尚需完成原生从零长跑和不同计算域验证。

也可以直接在 MATLAB 中使用配置接口：

```matlab
addpath('/absolute/path/to/ipm_long_time_server_v5_tau14_growth_20260913');
settings = struct( ...
    'canonicalFinalTime',1000, ...
    'maximumSteps',10000000, ...
    'autonomousMeshVersion',5, ...
    'initialNodeCount','auto', ...
    'maximumTotalNodes',1000000, ...
    'maximumAdjacentGridRatio',2, ...
    'outputDirectory','/data/ipm/run01');
opts = ipm.config.longTimeProfile(settings);
config = ipm.config.resolve(opts);  % 只检查，不启动计算
result = ipm.solve(opts);           % 从物理 t=0 启动
```

## 断点续算

程序按 `checkpointEvery` 保存不可覆盖的原生 checkpoint，并在正常或安全停止时保存末态 checkpoint。把 `profile_settings.m` 中的 `restartCheckpoint` 设为该文件的绝对路径，同时只延长 `canonicalFinalTime` 或 `maximumSteps`，再运行同一启动命令。

恢复会冻结原来的网格、物理、离散、C 规则和自动网格策略。尤其不能在恢复时把网格比从 1.5 改成 2；若需要 `2`，应在新的 `t=0` 算例启动前设置。跨机器复制旧 checkpoint 时，先用目标机器的 `profile_status` 校验：旧数值签名含浮点归约，线程数或平台差异可能使严格校验失败。新 `t=0` 算例没有这个迁移限制。

## 输出与检查

输出目录包含：

- `console.log`：完整控制台记录；
- `launch_configuration.mat`：用户设置、平铺选项和最终冻结配置；
- `checkpoint_*_step*.mat` 与 `checkpoint_manifest.jsonl`：可信接受步；
- `result*.mat` 与 `manifest.jsonl`：末态结果；
- `launcher_result_pointer.mat`：启动器返回的结果结构。

运行中只读查看最新可信 checkpoint（不建 LU、不改变轨道）：

```bash
matlab -batch "addpath('/absolute/path/to/ipm_long_time_server_v5_tau14_growth_20260913/server'); profile_status('/data/ipm/run01');"
```

`profile_status` 返回时间、步数、节点数、累计自动重布次数、核心格数、两轴相邻网格比、物理梯度以及远边界的源与速度指标。它验证 checkpoint，但 checkpoint 的存在本身不能证明求解器进程仍在运行；进程状态需由作业系统或 `ps` 单独确认。
若 checkpoint 来自本实验包，状态还会显示连续峰值影子 Y 核心格数及累计请求分歧。

查看结果：

```matlab
addpath('/absolute/path/to/ipm_long_time_server_v5_tau14_growth_20260913');
r = ipm.output.validate('/data/ipm/run01/result_CASE_ID.mat');
ipm.output.plotResult(r);
```

安装后的接口检查：

```matlab
addpath('/absolute/path/to/ipm_long_time_server_v5_tau14_growth_20260913');
addpath('/absolute/path/to/ipm_long_time_server_v5_tau14_growth_20260913/tests');
r = ipmtests.baseline.serverInterface();
```

完整验证使用 `ipm.verify('all','quick')`。新旧一致性检查需要另行提供旧版目录，服务器发布包本身不依赖旧版。

## 已验证范围

2026-09-13 按用户指令暂停本机两条原始 `t=0` H8 长跑，进程已退出。
暂停后只读验签的最新完整 checkpoint 如下；周期保存意味着中断前最后几个
接受步可能未保存。重启必须使用各自原来的源码、配置和线程条件，不要把
研究分支的设置强加到旧 checkpoint；本发布包不包含这些大型 checkpoint。

| 轨道 | 最新可信时间、步数 | 网格、重布 | 核心 X/Y、相邻格比 X/Y | 本机 checkpoint |
| --- | --- | --- | --- | --- |
| version 4 主线 | `τ=14.2002012595`，step 19135 | `961×321`，42 次 | 28.139/33.731，1.055504/1.064449 | `../ipm_structured/runs/profile_H8_tau16/checkpoint_*_step0000019135.mat` |
| version 5 | `τ=12.4002900296`，step 12890 | `641×321`，37 次 | 32.656/35.149，1.073693/1.059257 | `../ipm_grid_v5/runs/v5_H8_tau16_continuation_20260913/checkpoint_*_step0000012890.mat` |

两者均未达到 `τ=16`，远边界指标仍高于 0.01 告警阈值。
继续前以 `profile_status` 重新验签具体文件；从各自运行目录的
`checkpoint_manifest.jsonl` 获取完整文件名。下述数值是各阶段的历史验收记录。

2026-09-13 更新：原始 `t=0` version-4 H8 主轨道已在 τ11.932175 自动完成第 35 次
重布，一级分层 X 回退候选在真实场迁移中通过；τ12.00025 的验签 checkpoint 为
641×321、核心 X/Y=30.780/36.208、格宽比 1.072844/1.057942。
旧固定候选族的 τ≈11.67 边界已由同一原生轨道跨越。第五版原始 `t=0` 轨道已
原生验证至 τ6.5，之后继续运行。影子观察器在八个冻结 checkpoint 上与整轴 PCHIP
逐值一致，并通过早期 12 步、晚期 47 步、跨第 35 次自动重布的 123 步隔离原生
续跑与签名结果/断点验收；三个窗口没有影子/原生网格请求分歧。
跨重布轨道与主轨道共同的自由步上，180 个数值历史字段和重布事务账本逐值相等。
H8 的远边界告警仍在，尚无大箱与时空收敛证据。

最新的请求时刻配对实验（τ11.932175）显示，X90 备选网格通过真实场迁移和
安全审计；单次隔离续跑到下一次重网格请求的间隔比主轨道 X82 长 57.8%，
而最大相邻比从 1.072844 增至 1.075200。证据和复现脚本随包附在
`research/longtime_lab/`。这是研究分支的一次候选排序实验，正式 version-5
网格策略、服务器接口和 C 规则没有因此改变；仍需后续独立触发与长期精度验证。
第37、38次独立请求上的同成本寿命排序均保留原候选1；第37次的绝对寿命
预测明显偏长，具体限制见 `research/longtime_lab/AMORTIZED_MESH_LIFETIME_HOLDOUT_20260913.md`。
第38次的前瞻预测间隔为 0.41017，随后验签的第39次实际请求间隔为
0.24994；拟合寿命目前不能当作准确的重布时钟。
同一次请求的跨级 `961×321` 真实迁移也通过；相邻比降到 1.060356，
但比原选网格多 49.92% 节点，预测重布间隔只多 7.13%，本次原 RHS 计时
较慢。原始方程的墙面归一化形状导数在两候选间仍相差 16.67%，
因此暂不把跨级最小网格比改成正式选择规则。证据见
`research/longtime_lab/CROSS_LEVEL_MESH_COST_AND_RHS_20260913.md`。
该跨级候选的隔离原方程续跑随后到 `τ=13.039977` 才再次请求网格，
间隔比主轨道同轮长 32.23%，明显高于冻结态预测的 7.13%；单位时间步数
接近，而更大网格单次 RHS 成本明显更高。该收益不能直接换算为计算加速，
也没有独立精度认证。
主轨道随后在第42次请求 `τ=13.935876` 自主将 `641×321` 增至
`961×321`，真实迁移通过；见
`research/longtime_lab/evidence/native_growth_tx42_20260913.json`。
第42次增长的同时间墙面形状跳变仅 0.01896%，原 RHS 形状导数却几乎反向，
见 `research/acceleration_lab/ORIGINAL_T0_TX42_GROWTH_SHAPE_20260913.md`。

原始 `t=0` H8 到 τ13 的实际剖面图及研究审计也随包提供，见
`research/acceleration_lab/ORIGINAL_T0_TAU13_SHAPE_AND_REMESH_20260913.md`。
同一轨道更晚的 τ13.0–13.9 五帧显示物理墙面峰继续增长至 156.948，
墙核宽收缩到 2.8062e−5；归一化墙面形状首尾仍相差 0.46094%，见
`research/acceleration_lab/ORIGINAL_T0_TAU139_PROFILE_20260913.md`。
第38、39次真实重网格的同一时刻归一化剖面跳变很小，但缓存原 RHS
给出的形状导数变化显著；这提示必须先做空间/时间与更大箱的误差检验，
不能把光滑的有限时曲线直接外推为极限 Profile。
迁移两侧原 RHS 与连续 C1 峰值的只读分解见
`research/acceleration_lab/REMESH_RHS_C1_DECOMPOSITION_20260913.md`：
第38、39次的墙面 `Omega_tau` 相对差约 1.40%/2.54%，却对应更大的
归一化形状导数差；仅改变 C1 幅度率不能解决网格导数敏感性。

version 4 的分层/方向增点实现已通过 baseline、四阶、六阶和新旧一致性完整回归；真实旧轨道的冻结场迁移也通过原生配对检查。此前 version 2 从 `t=0` 的 H8 已运行到 `τ=11.6731`、H64 已运行到 `τ=9.9086`，均因有限候选族耗尽停止；version 4 正是针对这种“仍有可用网格但搜索族未覆盖”的故障。version 5 在 τ≈7.4 和 τ≈9.3 的真实冻结场上完成 38 级方向配置的无 LU 网格筛选：较晚时 `961×321` 用约 30.8 万节点把最佳相邻比从 `1.062837` 降到 `1.048762`；`961×481` 用约 46.2 万节点降到 `1.036248`；单独加到 `1281×321` 没有进一步改善。原生 H8 `t=0→τ=0.5` 完成 5 次自然自动重布，`321×161→641×161` 的真实第五版场/算子迁移及安全审计通过；H64/H128 和 H256 自动初始节点的短程原生运行通过，H256 checkpoint 续算至 τ=.008 也通过。修改后的 `ipm.verify('all','quick')` 与 `ipm.verify('equivalence','quick')` 均 exit0。更晚的大节点原生迁移和整段长跑仍待验证。

旧 H8 失败请求的 70 个基础 X 候选全部被拒；相邻步长比、局部求积权重和核心格数是主要限制，Y 方向仍有可用候选。新 version 4 原始 `t=0` H8 实跑已验证到至少 `τ≈9`：同节点自动重布至少 17 次，并在 `τ≈4.30` 自主从 `321×161` 增至 `641×161`，之后达到 `641×321`。远边界源/速度指标超过 `0.01` 告警阈值；H8 长跑可检验网格和时间推进，但最终 Profile 必须另作更大计算域的匹配比较。

已暂停的 version 4 H8 原始 `t=0` 轨道还未到 `τ=16`，version 5 也未完成该整段原生验证。本包只支持继续实验，不代表已经获得无穷时间极限、空间收敛或奇性证明。运行中不要编辑源码或 `profile_settings.m`；需要新配置时启动新的输出目录。
