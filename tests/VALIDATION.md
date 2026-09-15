# 验证记录

日期：2026-09-15；MATLAB R2026a。

## REMESH 后警告

- `ipm.verify('quadrant')`：通过。
- 强制 `gradientStop=realmin` 后，第 1 步接受第 1 次 REMESH，产生并持久化
  `gradient_threshold`，但求解仍以 `final_time` 结束。
- 同阈值固定网格对照没有 `metadata.remeshWarnings`，证明只在成功 REMESH 后检查。
- checkpoint 可保存并恢复第一象限非均匀网格及告警元数据。
- `rhoXStop=realmin` 回归以 `rho_x_threshold` 结束，`maxRemeshes=Inf` 保持为无限；
  checkpoint 恢复可显式覆盖 `rhoXStop=1000`。

## 无限时域终止探针

- `65 x 33` 第一象限探针将 canonical/physical 终点、步数与 REMESH 次数上限全部
  设为 `Inf`，临时使用 `rhoXStop=0.84`。
- 第 263、693 步的 REMESH 拒绝均只发布警告并继续；最终在第 1172 步以
  `rho_x_threshold` 结束，末值 `max|rho_x|=0.840078666`。

## 完整回归

- v3 `ipm.verify('baseline')`：通过。布局报告 131 个运行文件、
  4 个公开入口和 123 个 `ipm.solve` 依赖；新第一象限服务器接口测试通过。
- `ipm.verify('all','quick')`：通过，覆盖 baseline、四阶、六阶和第一象限。
- 变更 MATLAB 文件 `checkcode`：0 条。
- `git diff --check`：通过。
- v3 最小发行包的压缩结构、包内 SHA-256 和 ZIP SHA-256 在提交后生成并复核。

这些测试验证代码与数值接口，不替代服务器 1024 级内存预飞、长时间网格独立性或
爆破结论的数学证明。
