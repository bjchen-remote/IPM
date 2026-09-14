# `ipm` 包导航

`ipm.solve(opts)` 是唯一求解入口；`ipm.verify(...)` 是独立验证入口。
公共数据契约与完整调用路径见 [结构总览](../STRUCTURE.md)，第一象限设计与证据见
[架构记录](../ARCHITECTURE_QUADRANT_LEVELSET_20260914.md)。

| 子包 | 负责什么 | 下游接口 |
|---|---|---|
| [config](+config/INTERNAL.md) | 冻结配置及组合准入 | `config.resolve` |
| [mesh](+mesh/INTERNAL.md) | 坐标、导数、求积、完整运行算子 | `mesh.build` |
| [field](+field/INTERNAL.md) | 密度、Green/Poisson、速度与输运 | `field.velocity`, `field.transport` |
| [diagnostics](+diagnostics/INTERNAL.md) | 峰、level set、分辨率及停止诊断 | `diagnostics.meshFeatureIntervals` |
| [remesh](+remesh/INTERNAL.md) | 候选轴、迁移、验收事务 | `remesh.adapt` |
| [output](+output/INTERNAL.md) | 历史、v2 结果与原生 checkpoint | `output.finalize` |
| [evolve](+evolve/INTERNAL.md) | 初始化、时间推进、重网格调度 | `evolve.initialize`, `evolve.advance` |

第一象限的运行链是
`config → mesh → field → diagnostics → remesh → field → output`，由 `evolve` 协调。
`remesh` 失败返回原场/原轴，不提交计数；旧全域自动候选族与新象限直接方法
有各自的适用域，详见 [remesh 包说明](+remesh/INTERNAL.md)。
验证范围与命令见 [测试说明](../tests/README.md)；历史研究材料不属于运行依赖。
服务器按新发行包附带的 `main.m` 复制件绝对路径提交时，见
[第一象限部署说明](../server/README_QUADRANT_ABSOLUTE.md)。
