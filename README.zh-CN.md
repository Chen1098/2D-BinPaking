# kestrel —— 二维不规则排样 / 装箱求解器

[![CI](https://github.com/Chen1098/2D-BinPaking/actions/workflows/ci.yml/badge.svg)](https://github.com/Chen1098/2D-BinPaking/actions/workflows/ci.yml)

<p align="center">
  <img src="results/strip_20min/swim.svg" height="250" alt="SWIM"/>
  <img src="results/strip_20min/trousers.svg" height="250" alt="TROUSERS"/>
</p>

**kestrel** 把任意多边形（矩形只是最简单的特例）排进一条定高的条带并让长度最短，或者装进尽可能少的
固定尺寸板材。它是对本仓库原有“玩具级”解法的完全重写（Rust，无几何库依赖），目标是和已发表的最好
结果正面比较。

> **一段话结论。** 对照对象是 [sparrow](https://github.com/JeroenGar/sparrow)（Gardeyn、
> Vanden Berghe、Wauters，*European Journal of Operational Research* 2026）：在不规则条带排样的
> 13 个经典基准实例上，它的最好结果全部追平或超过了此前所有方法。两个求解器按该论文的实验协议
> （每次 20 分钟、3 线程）**在同一台机器上同时运行**这 13 个实例。20 分钟后，kestrel 在 **11** 个
> 实例上更密，sparrow 在 **1** 个上更密，**1** 个持平；平均领先 **0.48** 个百分点。kestrel 超过
> “sparrow 跑满 20 分钟的最终结果”所需时间的中位数是 **48 秒**；在更慢的机器上，它有 9 个实例的结果
> 高于论文报告的 sparrow 期望值。两个求解器的全部排样都由同一个独立校验程序用精确有理数运算复核。

[English](README.md) · [算法原理](docs/ALGORITHM.md) · [基准测试](#基准测试) · [局限](#它不是什么)

## 快速上手

```bash
cargo build --release

# 条带排样：最小化条带长度（JSON 实例，sparrow / jagua-rs 格式）
./target/release/kestrel solve -i data/strip/swim.json -t 60

# 装箱：本仓库最初的纯文本格式（板材 W H，然后是各个多边形）
./target/release/kestrel solve -i data/legacy/in.txt -t 30

# 用独立校验器检查任意结果（纯 Python，无第三方依赖）
python tools/verify.py data/strip/swim.json out/swim.json
```

每次运行输出 `out/<name>.json`（每个零件的旋转角和平移量）和 `out/<name>.svg`（排样图），并打印一行摘要：

```
RESULT name=swim bins=1 len=5596.689844 density=79.0411% feasible=true ...
```

常用选项（不带参数运行 `kestrel` 可看到全部）：

| 选项 | 含义 |
|---|---|
| `-t <秒>` | 时间限制（默认 60） |
| `-w <线程数>` | 工作线程（默认 3） |
| `-s <种子>` | 随机种子 |
| `--evals <百万次>` | 用确定性的“评估次数预算”代替时间限制：同种子同结果，与机器快慢无关 |
| `--rotations 0,90,180,270` | 覆盖实例里允许的旋转角 |
| `--mode strip\|bin`、`--bin-width <w>` | 把实例当作条带 / 装箱问题来解 |
| `--min-sep <d>`、`--strict` | 零件间保持最小间距（`--strict`：零件尺寸的 1e-6，使结果在精确运算下严格无重叠） |
| `--simplify <比例>` | 简化顶点极多的轮廓（只会把轮廓变大，不会变小） |

### 输入格式

* **JSON**（[jagua-rs](https://github.com/JeroenGar/jagua-rs) / sparrow 的格式）：条带排样用
  `strip_height`，装箱用 `"bin": {"width": W, "height": H}`；`items` 里每个零件有数量 `demand`、
  允许的旋转角和一个简单多边形。详见 [`data/README.md`](data/README.md)。
* **纯文本**——本仓库最初版本的格式：`W H`、多边形个数，然后每个多边形先给顶点数再给各顶点。
  它按装箱问题求解（需要几块 `W x H` 的板，最后一块能排多短）。

## 基准测试

### 测试协议

* sparrow 从源码编译（提交 `67bfaf8`，v0.3.0，`--features simd,only_final_svg`，
  `target-cpu=native`），使用它的默认配置。
* 对每个实例，sparrow 和 kestrel **在同一台笔记本上同一时刻启动**（Intel Core Ultra 5 125H，
  Windows 11），各用 3 线程，时间限制相同、种子编号相同、进程优先级相同、读取同一个实例文件；
  同时跑两对。这样机器上其他负载对两个求解器的影响是一样的。
* kestrel 使用 `--strict`，即零件之间严格保持间距——sparrow 一向如此。
* 跑完后，**双方**的排样都用 [`tools/verify.py`](tools/verify.py) 复核：它与两个求解器都不共享
  代码，用精确有理数运算计算所有零件两两之间的相交面积。双方所有排样的重叠面积都恰好为 0，
  下表的密度就是这个校验器重新算出来的。
* 整个流程就是一个脚本：
  `python tools/paired_benchmark.py --sparrow <sparrow 可执行文件> --time 1200 --seeds 1`。

### 20 分钟（sparrow 论文的协议）

密度 %，越高越好。

| 实例 | kestrel | sparrow | 差值 | kestrel 超过 sparrow 最终结果所用时间 | 论文中 sparrow 的期望值 | 上四分位 | 作者得到过的最好值 |
|---|---|---|---|---|---|---|---|
| albano | **89.59** | 88.89 | +0.69 | 9.0 s | 89.47 | 89.58 | 89.82 |
| dagli | **89.15** | 88.66 | +0.48 | 21 s | 89.26 | 89.50 | 90.17 |
| fu | **92.40** | 91.79 | +0.61 | 1.0 s | 92.24 | 92.40 | 92.41 |
| jakobs1 | 89.08 | 89.08 | +0.00 | - | 89.09 | 89.09 | 89.26 |
| jakobs2 | 84.18 | **84.88** | -0.71 | - | 84.77 | 85.25 | 87.73 |
| mao | **86.54** | 86.32 | +0.22 | 301 s | 86.14 | 86.45 | 86.87 |
| marques | **91.51** | 90.96 | +0.55 | 48 s | 90.93 | 91.02 | 92.02 |
| shapes0 | **68.78** | 67.90 | +0.89 | 19 s | 68.60 | 68.79 | 69.98 |
| shapes1 | **76.72** | 75.22 | +1.49 | 171 s | 75.69 | 75.97 | 76.73 |
| shapes2 (blaz1) | **84.97** | 84.13 | +0.84 | 32 s | 84.68 | 84.88 | 86.23 |
| shirts | **89.79** | 89.74 | +0.05 | 1105 s | 89.66 | 89.84 | 90.92 |
| swim | **78.03** | 77.24 | +0.79 | 325 s | 78.26 | 78.55 | 79.83 |
| trousers | **91.87** | 91.51 | +0.36 | 890 s | 91.73 | 91.93 | 92.62 |

“kestrel 超过 sparrow 最终结果所用时间”指 kestrel 运行到第几秒时，它的当前最好排样已经比 sparrow
跑满 20 分钟的最终排样更密。最后三列引自 sparrow 论文：100 次运行的期望密度和上四分位数（在更快的
台式机 CPU Ryzen 9 7950X 上），以及作者在全部实验中得到过的最好排样。

怎么读这张表：

* **JAKOBS1** 持平：两个求解器最终停在同一个排样长度（11.000）。
* **JAKOBS2** 是 sparrow 赢。这个实例的零件都是整数坐标，最好的排样是“严丝合缝”的拼图。
  用 `--strict` 时 kestrel 让所有零件之间保留微小间距，这也是 sparrow 一贯的做法（sparrow 排样里
  最小的间隙是同一量级）；在这条规则下 sparrow 找到了更好的排样。kestrel 的默认模式允许零件
  精确贴合，这时 6 次短时间运行里有 2 次达到已知最优长度 22（87.72 %）——但精确贴合的排样在校验器的
  有理数运算下并非“恰好零重叠”（校验器报告的重叠是零件面积的 2e-8），所以这些结果不计入本对比。
* **SHAPES1**（76.72 %）和 **FU**（92.40 %）与 sparrow 作者全部实验中的最好结果（76.73 %、
  92.41 %）相差不到 0.01 个百分点。
* **SHIRTS** 接近持平；从第 2 分钟到第 18 分钟是 sparrow 领先。
* 有 9 个实例，kestrel 从最初几秒到最后一直领先。

<p align="center">
  <img src="results/convergence_swim.svg" width="640" alt="SWIM 上密度随时间的变化"/>
</p>

### 60 秒

同样的流程，时间限制改为 1 分钟，每个实例跑三对（种子 1、2、3）。

| 实例 | kestrel（3 次平均） | sparrow（3 次平均） | 差值 | kestrel 更密的对数 |
|---|---|---|---|---|
| albano | **89.20** | 87.80 | +1.41 | 3 / 3 |
| dagli | **88.66** | 86.55 | +2.11 | 2 / 3 |
| fu | **91.99** | 90.57 | +1.42 | 3 / 3 |
| jakobs1 | **89.08** | 89.07 | +0.01 | 3 / 3 |
| jakobs2 | **83.94** | 82.99 | +0.95 | 2 / 3 |
| mao | **85.49** | 84.23 | +1.26 | 3 / 3 |
| marques | **90.93** | 89.75 | +1.18 | 3 / 3 |
| shapes0 | **68.00** | 66.65 | +1.35 | 3 / 3 |
| shapes1 | **75.10** | 73.65 | +1.46 | 3 / 3 |
| shapes2 (blaz1) | **85.00** | 82.80 | +2.20 | 3 / 3 |
| shirts | **89.13** | 88.82 | +0.31 | 3 / 3 |
| swim | **76.87** | 76.20 | +0.66 | 2 / 3 |
| trousers | **91.39** | 91.01 | +0.38 | 2 / 3 |

按平均值，kestrel 在全部 13 个实例上都更密（平均领先 1.13 个百分点）；39 对运行中有 35 对是 kestrel 更密。

### 优势从哪里来

sparrow 对每个候选位置都要拿零件和周围的多边形做碰撞检测。kestrel 则为每一对形状预先算好
“相对位置落在哪个区域就会相撞”（临界多边形 NFP），以及“撞得有多严重”的一张图，之后同样的问题
只需查一次表。在上面的运行中，kestrel 每秒评估的候选位置数大约是 sparrow 的 6 倍（随实例不同在
3.7～7.9 倍之间）——而搜索框架本身就是 sparrow 的。此外还有一些依赖精确几何的操作：把所有零件滑到
精确接触，以及把整个排样“挤”到卡死为止。[docs/ALGORITHM.md](docs/ALGORITHM.md) 有完整说明，
也列出了试过但没有效果的想法。

双方的排样图和解文件都在 [`results/strip_20min/`](results/strip_20min)；运行
`python tools/verify.py data/strip/<name>.json results/strip_20min/kestrel/<name>.json`
（或 `.../sparrow/<name>.json`）可重复校验。

<p align="center">
  <img src="results/strip_20min/shirts.svg" height="170" alt="SHIRTS"/>
  <img src="results/strip_20min/albano.svg" height="170" alt="ALBANO"/>
</p>
<p align="center">
  <img src="results/strip_20min/mao.svg" height="200" alt="MAO"/>
  <img src="results/strip_20min/marques.svg" height="200" alt="MARQUES"/>
  <img src="results/strip_20min/dagli.svg" height="200" alt="DAGLI"/>
</p>

## 它不是什么

* **不是最优性证明。** 这是启发式算法；这些实例的理论下界远低于任何已知排样。
* **统计样本很少。** 论文每个实例跑 100 次；上面 20 分钟的表每个实例只有一对运行，60 秒的表有三对，
  因为每一对 20 分钟的运行要占用笔记本 20 分钟。这类搜索单次运行的结果会有零点几个百分点的波动
  （论文报告的四分位距是 0.2–0.6），所以单看某一行里这么大的差距意义不大，要看 13 行的整体情况。
* **参数是在基准上调的。** 少数几个不属于 sparrow 的参数，是在这同一批实例上选定的。
* **时间越长差距越小。** kestrel 在第一分钟遥遥领先，20 分钟后领先幅度变小。它的速度让它用几秒钟
  达到 sparrow 几分钟才能达到的水平；但这些实例上最后零点几个百分点靠的是“找到另一种排法”，
  这一点上两个求解器的搜索方式相同。在 SHIRTS 和 JAKOBS2 上，光靠速度不足以领先。
* **只支持离散旋转。** 任意角度旋转用固定的若干角度近似；sparrow 原生支持连续旋转。
* **内存随“不同朝向形状数”的平方增长**（每一对形状一个场）。几百种不同形状没问题，几千种不行。
* **只支持简单多边形和矩形容器。** 零件内部的孔被忽略（当作实心），不支持不规则外形或带缺陷的板材。
* **装箱模式没有外部对比。** 它有测试（能找到完美装箱，结果经过精确校验），但没有和装箱问题的文献
  结果做基准对比。

## 仓库结构

```
src/geom.rs       点、多边形、凸包、偏移、轮廓简化
src/decomp.rs     凸分解（耳切法 + Hertel-Mehlhorn 合并）
src/nfp.rs        精确 NFP（临界多边形）场
src/problem.rs    实例、带朝向的形状、场表
src/search.rs     位置搜索：采样、坐标下降、精确吸附
src/separator.rs  引导式局部搜索、推挤松弛
src/solve.rs      条带 / 装箱主流程、压实、挤压
src/pool.rs       自旋等待的并行池
tools/verify.py            独立的精确校验器（kestrel 和 sparrow 的解文件都能读）
tools/paired_benchmark.py  与 sparrow 的同机并排基准
tools/bench.py             kestrel 参数的确定性 A/B 测试
data/             基准实例
results/          上表对应的排样结果
```

`cargo test --release` 运行单元测试和端到端测试；`kestrel selftest <实例>` 用独立判定方法交叉检查
该实例的所有 NFP 场。

## 历史

本仓库的第一个版本（2023 年）把每个多边形用外接矩形包起来，然后在板上一格一格地扫描找空位；它仍然
保留在 git 历史里。在它自带的测试文件上（`data/legacy/in.txt`：700 x 700 的板、57 个随机凸多边形），
它的外接矩形铺开占满了整块板（699 x 658）。kestrel 用 30 秒把同样的多边形排进了板长方向的前 432
个单位——板的 62 %。

## 参考与致谢

* J. Gardeyn, G. Vanden Berghe, T. Wauters, *An open-source heuristic to reboot 2D nesting
  research*, EJOR 2026（[arXiv:2509.13329](https://arxiv.org/abs/2509.13329)）。kestrel 的搜索主循环
  （采样移动、引导式局部搜索的权重、探索与压缩两个阶段）沿用了这篇论文，碰撞度量也是论文里的重叠代理，
  只是改成了查表。正因为 sparrow 开源，同机对比才成为可能。
* 基准实例来自 [ESICUP](https://www.euro-online.org/websites/esicup/data-sets/)，使用 sparrow 仓库
  附带的 JSON 版本（MIT 许可）。
* 其余参考文献见 [docs/ALGORITHM.md](docs/ALGORITHM.md)。

MIT 许可。
