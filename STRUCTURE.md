# 项目结构

本仓库只保留生产求解器、验证、一个示例和一个服务器发布入口。历史研究材料已从
当前分支移除，可由 Git 历史提交 `3b58d93` 恢复。

```text
+ipm/                 MATLAB 运行包
  +config/            配置解析与冻结
  +mesh/              第一象限坐标、差分、求积和 Poisson 算子
  +field/             初值、速度、Green/Poisson 与输运
  +diagnostics/       level set、质量指标和 REMESH 后警告
  +remesh/            单提案网格、迁移和事务验收
  +evolve/            状态、时间推进和 REMESH 调度
  +output/            历史、结果、checkpoint 和绘图
examples/             第一象限短例
server/               固定绝对路径入口与最小打包脚本
tests/                数值和接口回归
```

## 运行路径

```text
opts
  -> config.resolve
  -> mesh.build -> field.initialDensity -> evolve.flow
  -> evolve.advance
       -> selectTimestep -> Runge-Kutta -> finite check
       -> remeshIfNeeded -> one proposal -> transfer -> accept/fallback
       -> accepted REMESH only: publishRemeshWarnings
  -> physical max|rho_x| threshold
  -> output.record / maybeCheckpoint
  -> output.finalize / write / report
```

`state` 是唯一可变运行对象，持有 `rho/ops/flow/scale/config/runMetadata`。
REMESH 是事务：候选状态只有在几何、迁移和有限性全部通过后才替换当前状态。
质量阈值警告在事务提交后生成，只进入控制台和 `runMetadata.remeshWarnings`。

## 文件与数据契约

- `result.schemaVersion = 2`：末态、物理场、网格、历史和元数据分组保存。
- checkpoint schema 4：保存完整接受态、当前非均匀轴、冻结配置和输出游标。
- `quadrantOnly=true`：所有运行场、快照和 checkpoint 都是 `ny x nx` 正象限数组。
- 服务器 v3 将 canonical/physical 时间、步数和 REMESH 次数上限全部设为 `Inf`，
  每步只以物理 `max|rho_x|=1000` 作为正常终止条件。
- REMESH 提案拒绝或异常保留当前接受网格并持久化警告。非有限状态和机器时间步为零
  仍是不可继续的数值故障，不冒充可恢复警告。

局部实现见各子包 `INTERNAL.md`；第一象限数学、网格及内存边界见
[ARCHITECTURE_QUADRANT_LEVELSET_20260914.md](ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)。
