# A/B 的原生零时刻诊断

2026-09-09。**两支的规范率差在真实 t=0 已经存在，不是只能在晚时出现的
差异。** 但它同时包含箱大小和不同自动初选网格，不能据此归因于纯边界
误差，也不能替代已预登记的共同物理终点 A/B 比较。

使用 B 的冻结 runtime `box_space_factorial_v3/source`、默认 10 线程，
严格读取 A 的原生 step1653 checkpoint 与 B 已完整写出的 step501
checkpoint。只取两份可信 history 的首行及初始化元数据；两个 checkpoint
的当前时刻仅作为来源信息保存，没有互相比较。没有 restore、LU 或 PDE。
MATLAB session74495 exit0，两个原 CP 的 SHA 在读取前后均未改变。

两份首行均为 step=normalizedTime=canonicalTau=physicalTime=0，
Cx=Cy=Comega=1、Xshift=0、remeshCount=0。它们记录的是**解析初值在自动
选择的初始网格重新采样之后、任何 PDE 步之前**的状态。此时
dτ/dt=Cx/Comega=1，因此下表中的原 canonical 率与 physical-time 率相同；
没有使用 checkpoint 末态的尺度来换算首行。

| 首行量 | A：[-8,8]×[0,4] | B：[-16,16]×[0,8] | B−A，相对 A |
| --- | ---: | ---: | ---: |
| c_l = d log Cx / dt | .434918543086 | .461731121152 | +6.16497% |
| c_omega = d log Comega / dt | .125936730181 | .125264307368 | −.533937% |
| physicalRhoXInf | .686061449435 | .686061979512 | +7.72637e-7 |
| physicalRhoYInf | .686073237518 | .685259943674 | −.118543% |
| physicalGradInf | .835138794552 | .835112308462 | −3.17146e-5 |
| physicalQuadraticPeak | .686075619254 | .686075574520 | −6.52025e-8 |
| 物理 wall 90% 宽 | .481312446411 | .481311543210 | −1.87654e-6 |
| 物理 vertical 90% 宽 | .983065183093 | .982194070847 | −.0886119% |

`physicalRhoXInf` 是全计算域 max |rho_x|；`physicalGradInf` 是全计算域
max hypot(rho_x,rho_y)。二者不是同一个峰值，不能交换名称。
`physicalQuadraticPeak` 是规范所用墙上三点二次峰；它也不同于前两者。
二者 c_r 均为零。

## 初选网格的实际差异

原解析 physics、scaling、transport、elliptic 完整配置均逐值相同；原始
节点数都是 321×161，箱不同。解析初值均为 degenerate_primitive k=8。
两个实际初始化都选择第一个候选，其全部 x/y 坐标与各自预登记候选
逐值相同。初选后双方 x 轴和 y 轴都不相同。

| 来源/初选量 | A | B |
| --- | ---: | ---: |
| 原 uniform dx / dy | .05 / .025 | .1 / .05 |
| 原 uniform core x/y（纯预检） | 9.6810 / 40.0817 | 4.8295 / 20.0384 |
| 实际初选后 core x/y（原生首行） | 34.85945 / 39.32261 | 34.83374 / 36.36114 |
| 实际 x∈[.8,1.8] 的 dx 范围 | .0135133–.0218206 | .0135470–.0238069 |
| 实际 y∈[0,1] 的 dy 范围 | .025（均匀） | .0165046–.0449436 |
| 初始 safetyFactor | .2348061 | .3157057 |

原 uniform 行来自此前纯解析预检，不是未选轴时的原生 flow 记录；该处
没有保存瞬时速度/规范率。本轮没有为了补这些数据重算 Poisson。
实际初选 audit 的 core 与原生首行一致，解析重采样、零时钟、未创建新
physical epoch 的所有标志均通过严格 CP 检查。

由于 transport_anchor 把 c_l 直接联系到锚点速度，这些数据说明初始
非局部离散速度/规范已经具有明显的箱与网格联合敏感性，即使墙峰高度
非常接近。它没有把连续有限箱效应、Green/椭圆离散误差与不同纵向初选
分开，也没有提供纯 box 因果归因。仍须等待原定 B 终点及 C/D/E 对照。

小产物位于
`result/longtime/20260908_campaign_v1/acceleration_lab/zero_initial_AB_78cda38da37348399790cf469fcc80e5/run/`：
`summary.json` 为紧凑表，`report.json/.mat` 保留完整首行、原生初始化
audit/轴与纯预检证据，`file_audit.json` 保留 CP 不变校验。
没有修改运行文件、终点门或任何 C 规则。
