# 从零 A/B/C/D/E 的唯一只读比较器

入口：

```matlab
report = mesh_compare_registered_zero_cases(factorialRegistration, timeRegistration, newOutputDirectory);
```

只加载登记的真实结果与其原生checkpoint，默认10线程；无mesh.build、flow、transfer、advance或solve。
A的expectedConfig由原始A登记的opts重新正规化；B/C/D/E使用各自事先冻结的完整expectedConfig。
缺result或配对原生CP是pending；已有数据未过原生签名、精确result/CP/config/时钟/历史/lineage或终点门是rejected。
不能从登记文件、求解开始或最大步数退出推断算例已完成。

每项明确读取原历史字段：

- `history.common.physicalGradInf`：max hypot(physical rhoX, physical rhoY)，全梯度模长最大值。
- `history.common.physicalRhoXInf`：max abs(physical rhoX)，横向导数最大值。
- `history.common.physicalQuadraticPeak`：二次峰值。

这三个量分列，旧helper的同名`physicalGradientMaximum`不直接互相比。
对旧linear comparator只显式映射rhoXInf。物理时间log尺度率严格取真正末态：

```matlab
amplitude = exp(state.scale.logC_l - state.scale.logC_omega);
physicalRateCL = history.common.c_l(end) * amplitude;
physicalRateCOmega = history.common.c_omega(end) * amplitude;
```

绝不使用`timeStep.physicalClockSpeed`；后者是源步的时间步遥测，不能代替当前端点尺度。
canonicalCL/canonicalCOmega及其双向门另外保存。物理率门对A作分母和对B作分母分别执行，不能只取较宽的一侧。

## 原门与观察协议

固定共同物理终点1.9140859724939365、core [0,.8]×[0,.25]、holdout [0,2]×[0,1]；
各项实际物理箱都要覆盖窗口，不外推。原线性场/峰/legacy宽/率门调用原
`mesh_compare_fresh_box_data`并反向再次调用，其原判定完整保存。

Hermite场比较和同连续峰90%宽以相同预注册阈值另作辅助判定。
continuous peak-window [.1,.6]先映回每项真实规范轴；唯一非退化峰与横纵交点必须存在。
连续观察失败不会把原legacy值改写成通过，也不作为修改native gauge的依据。
原continuous helper绝对误差和协变资格可复用，但逐相位每层误差下降≥8没有全部通过。

每条真实边 B−A、D−C、C−A、D−B、E−A 分两个窗口，双向linear与Hermite并列，
保留原legacy总判定和独立Hermite/continuous辅助总判定。
full physicalGradInf是另列观察量；rhoXInf/二次峰才进入原峰门。
`allScientificGatesPassed`只有全部五项真实native资格和所有这些原/辅助门满足才true，
即使true也不是无穷域、整体收敛阶或奇性结论。

D−C−B+A只在四项真实终态齐全且覆盖共同窗后计算。
每窗口先将各自真实场观察到同一321×161物理均匀query网格，分别保存linear/Hermite数组，
对场做有符号`D-C-B+A`后再取相对A范数，禁止先把四个范数相减。
实际观察数组在各项`*_observations.mat`，有符号交互场在`signed_interaction_fields.mat`。

321×161是预注册观察网格，不是求解网格，也不保证窄核心的交互场被充分采样。
例如A实际纵向90%宽约.00323，而holdout query的dy=.00625。
因此共同query交互仅作固定观察诊断；双向原生节点场门另行执行，不能拿共同query加密/缺采样冒充PDE精度。
每项同时保存自己的linear−Hermite query观察差，以识别这层观察误差。

## 实际A自比较验证

冻结目录：`result/verification/autonomous_runtime_20260909/zero_case_comparator_v1/source`。
生产+ipm来自实际成功v2冻源，研究依赖完整冻结；静态闭包均在该目录，CodeAnalyzer=0。
实跑session26158 exit0，profile无求解路径。结果`zero_case_comparator_v1/run/report.json`：

- A真正step1653的原生/expectedConfig/coverage全通过；两窗口双向linear与Hermite节点自比较**精确0**，
  所有峰、率、legacy与连续宽自比较也为0。
- A的全梯度4.8904365810、rhoXInf4.0590845991、二次峰4.0593749138分列。
- 末态物理率CL=.94421165023、COmega=−1.64125247802。
- A自身固定query的Hermite−linear梯度最大相对Inf差：core3.78872508e-4、holdout3.61151165e-4。
  这不是B/C/D/E的场差，也不是PDE误差估计。
- A的legacy宽X/Y=.01216834505/.003165743708；连续同峰宽=.01216547173/.003234089843。
  不同纵向采样列的定义自身即可改变宽度，这个差异单独保存，未改任何.005跨算例门。
- B/C/D/E实际文件缺失，四项均pending；所有真实边与交互项pending，`allScientificGatesPassed=false`。

后续真实终态齐全后，在下述v4冻源path下重调同一helper，使用**新的输出目录**。
已有A自比较及pending报告不可覆盖，也不可把A复制成其他标签来模拟四格通过。

## 控制边界复核（2026-09-09）

独立审阅确认了末态率、双向门和有符号交互公式；另修正两处接口防护。
比较器现在严格要求factorial四标签A/B/C/D及time单标签E的完整顺序，
核对两份注册的peak-window，以及time引用的完整factorial注册；若time显式提供query形状，也必须相等。
每项观察MAT成功保存后才发布供比较的view/query，保存或其他异常会清空两者，避免rejected项进入真实边。

当前完整冻结源码为`result/verification/autonomous_runtime_20260909/zero_case_comparator_v4/source`。
`boundary_report.json`记录六个元数据反例实际通过：错序、缺项、错误E标签、错误peak-window、
错误显式query形状、错误parent注册。session57104 exit0，CodeAnalyzer=0，profile无native读取、LU或PDE；
保存后发布与catch清空另作源序静态检查。`final_source_bytes.json`确认冻结字节未变。
这一轮没有重算A的数值观察；真实A的全部数值证据仍在上述v1目录。
v2的CodeAnalyzer测试脚本风格失败及v3的脚本结果收集结构赋值失败保留，
后续以隔离函数收集六项实际错误标识完成v4回归，没有将测试框架异常记为PDE失败。
