# 测试

测试代码不进入服务器发行包，只用于本地验证当前生产源码。

```matlab
ipm.verify('baseline');
ipm.verify('quadrant');
ipm.verify('fourth');
ipm.verify('sixth');
ipm.verify('all','quick');
```

目录：

```text
tests/+ipmtests/+baseline/   配置、网格、椭圆、输运、REMESH、输出和 checkpoint
tests/+ipmtests/+fourth/     高阶/WENO5/SSPRK54
tests/+ipmtests/+sixth/      六阶算子、WENO7 和 RK6
tests/+ipmtests/quadrant.m   第一象限、直接重网格、警告和恢复
tests/+ipmtests/layout.m     包入口和运行依赖边界
```

服务器接口测试检查当前第一象限固定绝对路径、无限运行界、`rhoXStop=1000` 及最小发行包内容，
不再验证已删除的旧全域 R2 脚本。

最近执行结果见 [VALIDATION.md](VALIDATION.md)。测试通过不等于服务器大网格性能、
长时间收敛或爆破结论已经得到证明。
