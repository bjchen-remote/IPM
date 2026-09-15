# 第一象限服务器发布

当前大版本名：`ipm_long_time_server_20260915_v3`。

完整发行包部署到固定位置：

```text
/data/user/hd58131/ipm/ipm_long_time_server_20260915_v3
```

随后把包根目录的 `main.m` 复制到独立作业目录，并按复制件的绝对路径运行。
`main.m` 只用 `mfilename('fullpath')` 确定作业目录；求解源码位置始终使用上述固定
绝对路径，不从复制件位置推断。

## 本地打包

```bash
./server/package_quadrant_release.sh
```

脚本生成同名目录、ZIP、包内 `SHA256SUMS` 和 ZIP SHA-256。发行包只包含：

- `+ipm/`
- `main.m`
- `README.md`
- `STRUCTURE.md`
- 第一象限架构记录
- `CHANGELOG.md`
- `PACKAGE_INFO.txt` 与校验文件

研究脚本、测试套件、本地结果和旧服务器入口不会进入发行包。

## 部署预飞

先在服务器设置预飞开关：

```bash
IPM_PREFLIGHT_ONLY=1 /opt/MATLAB/R2026a/bin/matlab -batch \
  "run('/绝对路径/作业目录/main.m')"
```

预飞只检查固定发布路径和配置，不建立 Poisson/LU，也不创建输出目录。去掉环境变量
才开始 PDE。新作业必须使用空输出目录；旧全域 checkpoint 不能改成第一象限续算。

入口默认实际正象限网格为 `1025 x 513`。canonical/physical 终点、最大步数及最大
REMESH 次数均为 `Inf`，唯一正常终止门是物理 `max|rho_x|>=1000`。每隔 `0.01`
canonical time 的播报包含当前 `max|rho_x|`，checkpoint 间隔为 `0.1`。REMESH 保持
开启；质量越界、提案拒绝和 REMESH 异常只发布警告并保留当前网格。

如需恢复同版本原生 checkpoint，设置 `IPM_RESTART_CHECKPOINT` 为其绝对路径。
恢复时冻结网格、物理和数值配置，只延长原终止界并生成不覆盖旧文件的新产物。

详细数学和网格约束见
[第一象限架构记录](../ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)。
