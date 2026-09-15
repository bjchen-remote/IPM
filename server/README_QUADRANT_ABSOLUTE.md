# 第一象限：复制主文件、固定绝对发布路径

新发行包根目录直接含可用的 `main.m`，由
[main_quadrant_release.m](main_quadrant_release.m) 安装。将**完整包**部署到
`/data/user/hd58131/ipm/ipm_long_time_server_20260915`，再把该包中的
`main.m` 复制到独立服务器作业目录，并按**复制件的绝对路径**运行。
`main.m` 只用 `mfilename('fullpath')` 定位 `jobRoot`；源码
`releaseRoot='/data/user/hd58131/ipm/ipm_long_time_server_20260915'` 是固定字面量，
相对于用户提供的 `/data/user/hd58131/ipm/ipm_long_time_server_20260910`
仅改版本号。随后 `restoredefaultpath; addpath(releaseRoot)`，结果与日志都写入
`jobRoot/output_quadrant_NxXNy`。不从复制件位置推断源码发布路径。
旧发布包和本仓库旧全域 [main.m](../main.m) 均未修改。

本地用 [package_quadrant_release.sh](package_quadrant_release.sh) 生成新的不可覆盖
`ipm_long_time_server_20260915` 目录、ZIP 与 SHA-256 清单。
部署后先以 `IPM_PREFLIGHT_ONLY=1` 运行作业目录的复制件：它只解析配置，
打印固定发布路径、作业目录、实际节点数及输出位置，不建 Poisson/LU、不创建输出。

```bash
IPM_PREFLIGHT_ONLY=1 /opt/MATLAB/R2026a/bin/matlab -batch "run('/绝对路径/作业目录/main.m')"
```

去掉 `IPM_PREFLIGHT_ONLY=1` 才启动全新 PDE 作业。已有非空输出目录会被拒绝，
不得把旧全域 checkpoint 改成象限续算。

9.14 象限作业若因旧版 rounded-axis 可行性断言退出，可从最后一个原生 checkpoint
恢复到本修复版。设置 `IPM_RESTART_CHECKPOINT` 为该 checkpoint 的绝对路径后运行
新的主文件复制件；程序保留原冻结网格与数值配置，只延长既有时间上限，并为恢复日志、
结果和指针选择不覆盖旧文件的名称。不要设置它来迁移旧全域 checkpoint。

`main.m` 直接构造第一象限选项，不调用再删去 `longTimeProfile` 的旧自动网格政策。
`adaptiveRemesh=true` 与 `initialAnalyticRemesh=true` 使用新的 level-set 单提案方法；
不设置 `autonomousMesh`，因此旧 `maximumTotalNodes=310000` 不适用，但这也不构成
大网格内存保证。`1025×513=525825` 是实际存储节点；与旧全域 `1025×513`
相比，正半轴 X 的单元数从 512 增至 1024，而总场节点数并未减半。
若只想保持旧正半轴的 X 分辨率，可用 `513×513=263169`；若要实际
`1025×1025` 象限节点，则是 1,050,625 个场节点，须先做目标服务器的
单步 LU 峰值内存与耗时试跑。入口中的 `4.5/4000` 是提交上限，不是已验收长跑。

数学/网格限制见 [第一象限架构记录](../ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)，
旧全域服务器入口见 [README_SERVER_ZH.md](../README_SERVER_ZH.md)。
