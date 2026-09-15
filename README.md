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
opts.maxSteps = Inf;
result = ipm.solve(opts);
```

`nx` 是实际存储的正半轴节点数，不再记录负半轴镜像。

## REMESH 规则

重网格只使用 90% 壁面核心、左侧前沿及对应纵向核心的 level-set 位置和尺度，
一次确定一对新坐标轴。每次最多构造一组几何、执行一次原生迁移；候选未通过
网格、守恒、值域或有限性验收时保留原网格。

历史上的梯度、值域、振荡、相邻格比和核心分辨率停机阈值只在成功 REMESH 后
检查。越界项输出 `WARNING` 并写入 `result.metadata.remeshWarnings`，但不回滚、
不改写 `stopReason`、不终止 PDE。

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
