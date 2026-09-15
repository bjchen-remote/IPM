# IPM 第一象限求解器

本仓库维护二维不可压多孔介质方程的 MATLAB 求解器，生产路径只存储
`X>=0, Y>=0` 第一象限，并在 `(1,0)` 附近通过 level-set 直接重网格捕捉尖峰。
唯一求解入口是 `result = ipm.solve(opts)`。

## 快速运行

在 MATLAB 中把仓库根目录加入路径，然后运行短例：

```matlab
projectRoot = pwd;
addpath(projectRoot);
run(fullfile(projectRoot,'examples','quadrant_level_set.m'));
```

生产配置必须显式设置：

```matlab
opts.quadrantOnly = true;
opts.xlim = [0,8];
opts.adaptiveRemesh = true;
opts.finalTime = Inf;
opts.physicalFinalTime = Inf;
opts.maxSteps = Inf;
opts.maxRemeshes = Inf;
opts.rhoXStop = 1000;
result = ipm.solve(opts);
```

`nx` 是实际存储的正半轴节点数，不再记录负半轴镜像。

## REMESH 规则

重网格只使用 90% 壁面核心、左侧前沿及对应纵向核心的 level-set 位置和尺度，
一次确定一对新坐标轴。每次最多构造一组几何、执行一次原生迁移；候选未通过
网格、守恒、值域或有限性验收时保留原网格。

候选失败后不会逐步重复反演。程序冻结失败时的 90% 核心中心、横向核心宽度和
纵向核心宽度；只有这些 level-set 几何累计变化至少一个目标网格单元时才再次尝试。
重试门不使用步数或时间。`max|rho_x|` 直接复用 RHS 特征跟踪已经计算的最大值，
不会为终止门再次扫描整个二维场。

历史上的梯度、值域、振荡、相邻格比和核心分辨率停机阈值只在成功 REMESH 后
检查。越界项输出 `WARNING` 并写入 `result.metadata.remeshWarnings`，但不回滚、
不改写 `stopReason`、不终止 PDE。

生产服务器入口没有 canonical/physical 时间终点、步数上限或重网格次数上限；唯一
配置的正常终止条件是物理 `max|rho_x| >= 1000`。每次输出都会播报该值。重网格构造
或迁移异常也会保留当前已接受网格并记录警告。非有限算术或时间步下溢到机器零属于
无法继续推进的数值故障，程序会保存最后有限态后退出。

完整设计见[架构记录](ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)，模块关系见
[结构说明](STRUCTURE.md)，包内接口从[包导航](+ipm/README.md)进入。

## 服务器发布

服务器采用复制 `main.m` 并按绝对路径提交的方式。生产入口位于
[server/main_quadrant_release.m](server/main_quadrant_release.m)，只从自身路径取得
作业输出目录；求解器始终使用代码中固定的服务器发布绝对路径。

运行 [server/package_quadrant_release.sh](server/package_quadrant_release.sh) 生成最小
发行包。发行包只包含 `+ipm`、根 `main.m`、必要说明和校验文件，不再携带研究脚本、
测试套件或旧版入口。部署与预飞见
[服务器说明](server/README_QUADRANT_ABSOLUTE.md)。

## 验证

```matlab
ipm.verify('quadrant');
ipm.verify('baseline');
ipm.verify('all','quick');
```

当前验证记录见 [tests/VALIDATION.md](tests/VALIDATION.md)。本地计算产物放在
`runs/` 或 `result/`，两者均不进入 Git。
