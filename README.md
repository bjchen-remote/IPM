# IPM：结构化独立副本

本目录保留原求解器的数值方法，按职责拆分为 MATLAB 包。原版
[`../sixth_order_integration`](../sixth_order_integration) 原地保留；本副本可独立运行，
不调用原版文件。来源为提交 `65f999b7e205dd4decef06c937d1aa003cf8b744`；
原版 140 个已跟踪文件的 SHA-256 与新旧函数对应见 [source_map.json](tests/source_map.json)。

求解二维上半平面截断域上的不可压多孔介质方程：

```text
rho_t + u · grad(rho) = 0,
u = grad^perp psi,    -Delta psi = d_x rho.
```

## 从一个短算例开始

在 MATLAB 中将当前目录设为本目录，再运行：

```matlab
projectRoot = pwd;
addpath(projectRoot);  % 只添加包的父目录，不递归添加 +ipm 子目录
run(fullfile(projectRoot,'examples','quick_start.m'));
```

[quick_start.m](examples/quick_start.m)、[fourth_order.m](examples/fourth_order.m) 和
[sixth_order.m](examples/sixth_order.m) 均为短算例，关闭结果保存、绘图和视频。
正式计算统一调用 `result = ipm.solve(opts)`，其中 `opts` 是平铺的标量结构体。
可调用 `ipm.config.resolve(opts)` 检查配置；运行时仍向 `ipm.solve` 传入原始 `opts`。

**不要用 `ipm.solve()` 做安装检查。** 无参数或空输入会加载
[activeCase.m](+ipm/+config/activeCase.m)：`513 × 257` 大域动态重标度算例，
启用解析初始重网格、自适应重网格及常规输出。`ipm.solve(struct())` 则采用配置表默认值，
同样不是短算例。实际选项与默认值以 [schema.m](+ipm/+config/schema.m) 为准。

长时间 Profile 的服务器接口使用
`opts = ipm.config.longTimeProfile(settings)`；可编辑的启动文件、断点续算和
`maximumAdjacentGridRatio=2` 的准确语义见
[README_SERVER_ZH.md](README_SERVER_ZH.md)。该辅助函数只构造并校验配置，数值求解仍统一调用
`ipm.solve`。

本次实验包附有连续峰值 X 处的竖向核心**影子观测**。它在 `result.metadata.continuousVerticalShadow`
中记录与原网格请求的分歧，不修改正式触发、硬停机或 C 规则。冻结数据、当前自适应网格的容量
问题和验证范围见 [自动网格研究说明](research/longtime_lab/AUTONOMOUS_GRID_DESIGN_FROM_FROZEN_DATA_20260913.md)。
R2.2 新算例的服务器默认幅值规范改为 `(2,0)` 外壁面密度窗；旧峰值规范仍可显式选择。
这项 C 规则变更和短测试范围见 [R2.2 发行说明](RELEASE_NOTES_R2_ZH.md)。
R2.3 增加了只读的外区截断 `L^p` 收敛诊断入口
`server/analyze_outer_profile.m`，并修复新 C 下影子纵向核心观测缺失；
正式 C 和网格决策不变。初步数据和未验证范围见
[R2.3 发行说明](RELEASE_NOTES_R2_3_ZH.md)。

## 数值选择与边界

| 路径 | 空间离散 | 输运 | 时间推进 |
|---|---|---|---|
| 配置表基线 | `legacy_second_order` | `muscl_minmod` | `ssprk3` |
| 无参数算例 | `legacy_second_order` | `weno5_nonuniform` | `ssprk3` |
| 光滑四阶 | `high_order` | `weno5_fd` | `ssprk54` |
| 光滑六阶 | `sixth_order` | `weno7_fd` | `rk6` |

六阶使用 `opts = ipm.config.sixthOrder(overrides)`，同时指定
`remeshTransferScheme='sixth_order'`；四个选择必须配套。辅助函数不会绕过配置校验，
也不会替调用者选择合适的区域、初值、分辨率和终止时间。

## 动态重标度契约

新配置是 schema 4，并固定
`scalingContract='exact_gauge_no_feedback_v1'`。动态规范只解瞬时精确约束：
`lengthScaleGain=1`，`omegaGaugeGain=widthGaugeGain=travelingWaveGain=0`。
`maxDynamicRate`已删除；`adaptiveGain` 运行组件与
`adaptiveLengthScaling`、`widthExpansionStrength`、`widthContractionOnset`、
`widthContractionStrength`、`maxWidthRateCorrection` 这些 width-controller 选项已退役。
规范的参考量可用于定义、
符号、正性和条件数安全门，但不得以 restoring 项回馈到变化率。
峰值、层级和连通宽度的网格单元计数只是诊断、重网格触发/提案和硬停机输入，
不参与 `c_l`、`c_omega` 或 `c_r` 的方程。

四、六阶的整体阶数证据限于相应测试覆盖的光滑固定网格问题；六阶仅覆盖均匀网格及
通过质量检查的光滑 sinh 网格。动态规范、触发式重网格、任意非均匀网格均不因此获得
整体高阶保证。RK6 不具有 SSP/TVD/保正保证。`weno5_nonuniform` 在边界层退回 MUSCL，
不能视为整体四阶方法。原始证据和限制保留在 [research/notes](research/README.md)。

## 验证与结果

```matlab
ipm.verify();                    % 基线；不运行生产算例
ipm.verify('fourth');            % 四阶，默认 quick 稳定性扫描
ipm.verify('sixth');             % 六阶
ipm.verify('all');               % 独立结构检查 + 基线、四阶、六阶
ipm.verify('equivalence');       % 另与保留的原版比较，需要原版目录
```

完整说明、`heavy` 模式与指定原版路径见 [tests/README.md](tests/README.md)。
本次副本的验证记录见 [VALIDATION.md](tests/VALIDATION.md)，变更见 [CHANGELOG.md](CHANGELOG.md)；
历史测试记录不等于本次验证结果。

输出保持 `schemaVersion=2`：末态在 `result.state`，物理量在 `result.physical`，
网格在 `result.grid`，记录在 `result.history`，实际输出路径在 `result.metadata`。
用 `ipm.output.validate(resultOrMatFile)` 读取和校验，
用 `ipm.output.plotResult(result)` 绘图；不提供旧的扁平字段兼容层。
开启保存时，默认写入当前工作目录下的 `result/production/`，已有结果使用新名字保留；
从其他目录运行时，建议显式设置绝对 `resultFile`、`videoFile`。

长算例可独立于 `storeSnapshots` 开启可信接受步 checkpoint；`every` 使用 canonical
时间，实际文件是带 case ID 和步数的不可覆盖产物：

```matlab
opts.checkpoint = struct('file','result/production/run_checkpoint.mat', ...
    'every',0.05,'atExit',true);
result = ipm.solve(opts);

resume = struct('maxSteps',300000,'physicalFinalTime',2.1, ...
    'checkpoint',opts.checkpoint);
result = ipm.solve(resume,result.metadata.latestCheckpointFile);
```

恢复只允许延长终止界及修改保存、绘图、视频、元数据和 checkpoint 控制；网格、物理、
离散与重标度选择保持冻结。原生 checkpoint 保存完整内部尺度、参考网格、重标度参考量、
标量历史与输出游标，可作 same-grid exact resume。当前 checkpoint 为 schema v4，
只接受 config schema v4 且要求
`scalingContract='exact_gauge_no_feedback_v1'`。旧 schema-v3 checkpoint 仍属可读的
分析证据，但它保存的是允许 restoring/width feedback 的旧运行期契约；
schema-v1/2 还属于 `maxDynamicRate` 契约。这些旧 checkpoint 都不会被隐式升级或续算。
可信 terminal v2 结果仅在其冻结配置为 schema v4 时可用
`ipm.output.checkpointFromResult(resultOrFile,targetFile)` 受控转换；带 config schema 1--3 的
历史 result 仍可读取分析，但不能桥接到当前续算。

## 去哪里修改

[STRUCTURE.md](STRUCTURE.md) 说明模块边界与数据契约；每个模块的 `INTERNAL.md`
说明局部实现；[research/README.md](research/README.md) 管理实验入口和历史笔记。
结构借鉴 `iclm` 可见的包目录、就近文档及测试/研究分离方式，不增加不必要的类或包装层。
