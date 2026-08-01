# 18 槽 / 20 极 PMSM 快速筛选与 FEMM 验收

[English Version](README.md)

本仓库提供一套使用同一配置文件的 18 槽 / 20 极表贴式永磁同步电机分析流程：

- `design_pmsm.py` 用于带电流、电压约束的快速稳态 dq 筛选。
- `femm_spm_template.py` 用于建立二维 FEMM 几何，并完成电磁验收和数值收敛检查。
- 两条路径都严格读取采用 SI 单位的 [`motor_config.json`](motor_config.json)。

下图是本次重新生成的中等网格 FEMM 快照，工况为 `Id = 0`、
`Iq = 5 A 峰值`、转子机械角 `0 deg`。`600 rpm` 只作为工况元数据；
转速不会改变磁静态 FEMM 解。

![当前 FEMM 磁密快照](docs/assets/femm_field_density.png)

## 当前状态

转子、18 槽绕组、Park 变换、相序和电角度零点已经使用同一套约定。
当前缩减回归通过了全部 6 项电磁门限和全部 9 项收敛门限。这些结果适合做
代码与模型回归检查，但不能替代工程定型验收。

原 `output_femm/` 目录中的图片由旧模型生成，本 README 已不再引用。旧图中的
`33 mH` 电感、`15 A` 时约 `0.01 N m` 的转矩以及 `20 deg` 齿槽扫描都与
当前模型不一致。README 的全部图片链接现已统一指向 `docs/assets/`。

## 快速开始

创建 Python 环境：

```powershell
conda create -n motor python=3.11 -y
conda activate motor
pip install -r requirements.txt
```

运行默认快速筛选：

```powershell
.\motor.bat quick
```

将整个包络限制在当前已经验证的 `5 A 峰值` 范围：

```powershell
.\motor.bat quick --current-limit 5 --torque 2 --max-rpm 1000
```

运行全部单元测试：

```powershell
.\motor.bat test
```

在 Windows 上运行短时 FEMM 电磁冒烟验收：

```powershell
.\motor.bat femm-validate
```

启动器依次查找：显式指定的 `-PythonExe`、仓库 `.venv`、已有的 `motor`
Conda 环境，最后是 `PATH` 中通过测试的 Python。所有路径都基于仓库目录，
因此可以从其他工作目录调用。

## 已验证配置

| 项目 | 当前值 |
|---|---:|
| 槽数 / 极数 / 极对数 | `18 / 20 / 10` |
| 绕组 | 单层、星形连接、每相 3 个串联线圈 |
| 匝数 | 每槽 `25` 匝，每相串联 `75` 匝 |
| 基波绕组因数 | `0.9452136366` |
| 叠长 | `80 mm` |
| 非磁孔半径 / 转子铁芯外半径 | `25 / 30 mm` |
| 转子铁芯径向范围 / 厚度 | `25–30 mm / 5 mm` |
| 磁钢厚度 / 气隙 | `1.5 / 0.5 mm` |
| 定子内半径 / 槽外半径 | `32 / 40 mm` |
| 定子轭部径向范围 / 厚度 | `40–45 mm / 5 mm` |
| 定子外半径 / 外空气边界半径 | `45 / 70 mm` |
| 永磁磁链 `psi_pm` | `0.03201896826 Wb` |
| 增量 `Ld / Lq` | `0.962518 / 0.939532 mH` |
| 20 deg C 相电阻 | `0.22 ohm` |
| 直流母线 / 调制上限 | `48 V / 0.95` |
| 相电压峰值上限 | `26.3272 V` |
| 配置中的电流上限 | `15 A 峰值` |
| FEMM 参数已支持范围 | 到 `5 A 峰值` |

下图直接由 FEMM 使用的槽位表生成。每相包含三个正向和三个反向线圈边，
三相基波幅值平衡，电角度相差 120 度。

图中的浅灰环是定子铁芯：连续齿顶桥位于 `32–32.2 mm`，绕组槽区位于
`32.2–40 mm`，线圈外侧的
`40–45 mm` 浅灰部分才是 5 mm 定子轭部。深灰环是 `25–30 mm` 的
5 mm 转子铁芯。中心白色 `0–25 mm` 区域在当前 FEMM 模型中设为非磁性孔，
不是铁芯；`45–70 mm` 的外部计算空气没有画在这张绕组示意图中。

共享配置为兼容原有文件仍使用字段名 `shaft_radius_m`，但该字段当前定义的是
非磁孔边界，FEMM 在孔内赋予 `Air`。如果 25 mm 实际要表示实体钢轴，必须先增加
轴材料建模并重新执行全部电磁验收和收敛，不能直接沿用本文结果。

![当前 18 槽 / 20 极绕组](docs/assets/winding_18s20p.png)

## 电气约定

全部 dq 电流和相电流均使用峰值。正弦相电流满足
`I_rms = I_peak / sqrt(2)`。统一约定为：

- 相序为 `A-B-C`。
- 正 `Iq` 产生正转矩。
- 机械角 `theta_m = 0` 时，空载永磁磁链位于正 d 轴。
- 角度以度表示时，`theta_e = -10 * theta_m - 100 deg`。
- Park 逆变换采用
  `i_phase = Id*cos(theta_phase) + Iq*sin(theta_phase)`。

负电角度方向与 FEMM 几何中的转子旋转方向一致。若只单独修改其中一个符号，
本文报告的验收结果将不再成立。

## 当前结果

### 快速 dq 筛选

默认结果采用配置中的 `15 A 峰值`驱动上限，请求轴转矩为 `2 N m`。
超过 `5 A` 的结果属于电磁参数外推，CSV 会用 `load_extrapolated` 和
`envelope_extrapolated` 两列明确标记。

| 指标 | 结果 | 含义 |
|---|---:|---|
| 空载、`Id=0` 电压基速 | `785.18 rpm` | 解析电压极限 |
| 网格估计恒转矩基速 | `650 rpm` | 扫速间隔为 25 rpm |
| `2 N m` 请求转矩最高可行点 | `1225 rpm` | 后续点被约束裁剪 |
| 低速最大轴转矩 | `7.204 N m` | 外推到 15 A |
| 包络峰值轴功率 | `503.71 W` | 外推到 15 A |
| 最大 dq / 系统功率残差 | `1.99e-13 W` | 数值闭合误差 |

![当前快速约束扫速结果](docs/assets/quick_sweep.png)

效率曲线明确表示“仅计铜耗的效率上界”。给定配置中的铁耗、磁钢损耗、
机械损耗和逆变器损耗尚未标定。温度曲线是稳态热阻估算，不是热瞬态。

### 电磁验收

本次报告的回归使用中等网格、`Iq = +/-1 A 峰值`、中心差分
`delta-I = 1 A`，并在一个电周期内取 3 个转子位置。转矩对称性和斜率均先扣除
同一转子角度下的齿槽转矩。

| 验收项 | 测量结果 | 门限 | 状态 |
|---|---:|---:|:---:|
| 空载 `psi_d0` 为正 | `32.0190 mWb` | `> 0` | 通过 |
| `rms(psi_q0) / psi_d0` | `0.179794%` | `<= 5%` | 通过 |
| `T(+Iq)` 与 `-T(-Iq)` 误差 | `0.051590%` | `<= 10%` | 通过 |
| 转矩斜率 | `0.482547 N m/A`，理论值 `0.480285 N m/A`；误差 `0.471050%` | 误差 `<= 20%` | 通过 |
| 增量 `Ld / Lq` 差异 | `2.416970%` | `<= 20%` | 通过 |
| 加载转矩 / 齿槽转矩 | `109.76` | `>= 5` | 通过 |

理论转矩斜率为 `1.5 * p * psi_pm`。下图同时给出空载 dq 对齐、扣除齿槽转矩
后的正负电流转矩幅值，以及使用正负电流中心差分得到的增量电感。

由于转子铁芯和定子轭部都只有 5 mm，另在机械角 `0 deg` 直接比较了 medium/fine
两套全局网格。相对细网格，中网格的 `psi_d0`、转矩斜率、`Ld`、`Lq` 变化分别为
`0.0051%`、`0.0194%`、`0.0666%`、`0.0882%`；最大变化低于 `0.1%`。
这项检查只验证全局网格敏感性，不能替代上面的三转角验收。

![当前 FEMM 电磁验收](docs/assets/electromagnetic_acceptance.png)

### 转矩收敛

仓库中的收敛结果是有意缩减的回归：每档气隙网格取 3 个位置，再比较 3 点和
6 点角度扫描。转矩积分还在三个气隙半径上比较周期性 180/360/720 点 Maxwell
积分，并与 FEMM 加权应力张量转矩交叉检查。

| 验收项 | 测量变化或比值 | 门限 | 状态 |
|---|---:|---:|:---:|
| 加载转矩，气隙网格 | `0.00899%` | `<= 2%` | 通过 |
| 加载转矩，角度步长 | `0.09441%` | `<= 2%` | 通过 |
| 齿槽峰峰值，网格 | `4.0559%` | `<= 10%` | 通过 |
| 齿槽峰峰值，角度步长 | `2.7917%` | `<= 10%` | 通过 |
| 加载转矩 / 齿槽转矩 | `76.59` | `>= 5` | 通过 |
| 加载转矩 / 重网格噪声底 | `398.27` | `>= 20` | 通过 |
| Maxwell 积分点数变化 | `0.19536%` | `<= 2%` | 通过 |
| Maxwell 积分半径离散 | `0.24958%` | `<= 5%` | 通过 |
| Maxwell 与 WST 差异 | `0.13276%` | `<= 10%` | 通过 |

因此，在这次回归中，加载转矩明显高于齿槽转矩和重网格数值噪声。

![当前 FEMM 转矩收敛](docs/assets/torque_convergence.png)

### 可复现 FEMM 场快照

README 中的场图使用一个明确指定的中等网格工况：

| 量 | 值 |
|---|---:|
| 工况标签转速 | `600 rpm` |
| 转子机械角 | `0 deg` |
| `Id / Iq` | `0 / 5 A 峰值` |
| FEMM 转矩 | `2.411687 N m` |
| `psi_d / psi_q` | `0.03202130 / 0.00474991 Wb` |
| 场图采样 | `24 径向 x 120 周向` |
| 气隙采样 | `360 周向` |

加载时的非零 `psi_q` 是 q 轴电感磁链；上面的 `psi_q0` 验收是在零电流下完成。

![当前 FEMM 磁力线](docs/assets/femm_flux_lines.png)

![当前 FEMM 气隙磁密](docs/assets/femm_airgap_flux_density.png)

## FEMM 命令

安装 FEMM 4.2；如 COM 服务尚未注册，请在管理员 CMD 中执行：

```bat
"C:\FEMM42\femm.exe" /regserver
```

生成 README 使用的同一场快照：

```powershell
python .\femm_spm_template.py `
  --config .\motor_config.json `
  --analysis field airgap `
  --rpm-list 600 --iq-list 5 `
  --snapshot-angle-deg 0 `
  --field-radial-points 24 --field-angular-points 120 `
  --airgap-points 360 --mesh-level medium --workers 1 `
  --out .\output_readme_femm
```

复现仓库中的缩减电磁验收与收敛表：

```powershell
python .\femm_spm_template.py `
  --config .\motor_config.json `
  --analysis validate convergence `
  --mesh-level medium `
  --validation-iq 1 --validation-delta-current 1 --validation-steps 3 `
  --convergence-iq 5 --convergence-mesh-points 3 `
  --convergence-angle-points "3,6" `
  --out .\output_femm_fix_smoke
```

对 5 mm 转子铁芯和 5 mm 定子轭部做一个代表转角的细网格交叉检查：

```powershell
python .\femm_spm_template.py `
  --config .\motor_config.json `
  --analysis meshcheck `
  --mesh-check-reference-level medium --mesh-check-level fine `
  --validation-iq 1 --validation-delta-current 1 `
  --out .\output_femm_meshcheck
```

`meshcheck` 在同一转角直接比较 medium/fine 两套全局网格，不能替代至少 3 个
转角的 `validate` 验收。

在已验证电流范围内运行一个完整电周期波形：

```powershell
python .\femm_spm_template.py `
  --analysis basic --rpm-list 600 --iq-list 5 `
  --points 72 --mesh-level medium --out .\output_femm_basic
```

分别运行齿槽转矩或增量电感扫描：

```powershell
python .\femm_spm_template.py --analysis cogging --cogging-steps 72 --out .\output_femm_cogging
python .\femm_spm_template.py --analysis inductance --ind-current 1 --ind-steps 37 --out .\output_femm_inductance
```

新研究应使用新的输出目录。CLI 会写出解析后的配置，但不会删除已有输出目录中的
无关文件。

### 分析类型

| 名称 | 用途 | 主要输出 |
|---|---|---|
| `basic` | 指定电流下一个电周期的转矩/损耗采样 | `femm_waveforms_*.csv/.png`、`femm_summary.csv` |
| `field` | 指定角度的一次磁密与磁力线快照 | `field_snapshot_*.csv`、`field_density_*.png`、`field_lines_*.png` |
| `airgap` | 指定角度的一次气隙 `Bn/Bt` 快照 | `airgap_B_*.csv/.png` |
| `cogging` | 零电流下一个真实齿槽周期（`2 deg`） | `cogging_torque.csv/.png` |
| `inductance` | 中心差分增量 `Ld/Lq` 扫描 | `inductance.csv/.png` |
| `validate` | dq 对齐、转矩对称/斜率、`Ld/Lq`、加载/齿槽门限 | 电磁验收 CSV |
| `meshcheck` | 单转角 medium/fine 全局网格比较，不代替完整验收 | `electromagnetic_mesh_*.json` |
| `convergence` | 气隙网格、角度步长和 Maxwell/WST 收敛 | 转矩收敛 CSV |
| `emap` | 原始指定电流的电磁 Map | 实验功能；README 不展示 |
| `tncurve` | 基于原始 Map 最大值的比较 | 实验功能；README 不展示 |
| `all` | 包含长时间和实验分析的全部项目 | 日常使用不推荐 |

FEMM 的 rpm/Iq 点是直接指定电流的磁静态结果，没有经过 48 V 逆变器电压椭圆
筛选。系统可达性应使用快速 dq 模型判断。

### 运行时间与求解次数

FEMM 时间与计算机有关。本机本次单张 README 快照约 `4 分 7 秒`，中网格
三转角验收约 `19 分 25 秒`，单转角 medium/fine 组合网格检查约 `35 分 53 秒`。
这些时间均包括几何建立；场快照时间还包括绘图采样。规划任务时，求解次数比
估计分钟数更可靠：

| 分析 | 近似 FEM 求解次数 |
|---|---:|
| 场图 + 气隙快照 | `1` |
| 单工况 basic，72 点 | `72` |
| 齿槽转矩，72 点 | `72` |
| 电感，37 个位置 | `185`（`5 x 37`） |
| 验收，3 个位置且测试/扰动电流相同 | `15`（`5 x 3`） |
| 单转角 `meshcheck` | `10`（两套网格各 `5` 次） |
| 默认原始 emap | `192`（`4 x 4 x 12`） |

仓库中的缩减收敛回归本次约耗时 2 小时 49 分，并可从 `.convergence_cache` 续跑。
生产级收敛研究可能需要数小时以上。

## 重新生成 README 图片

当验收 CSV/JSON、两份收敛 CSV、medium/fine 比较 JSON、全部场图和快照元数据
均已生成，并且各输出目录同时含 `resolved_femm_config.json` 和
`femm_run_manifest.json` 时，使用下面的命令重建全部 README 图片和来源清单：

```powershell
python .\generate_readme_assets.py `
  --femm-image-dir .\output_readme_femm `
  --femm-case-tag rpm600_iq5_a0
```

生成器使用非交互 Matplotlib 后端和稳定文件名，并在
[`docs/assets/manifest.json`](docs/assets/manifest.json) 中记录来源 SHA-256、
图片尺寸和图片哈希。每次成功的 FEMM 调用还会追加
`femm_run_manifest.json`，把分析参数、物理模型指纹和本次实际改写的每个结果文件
绑定起来；图片生成器会逐个核对这些记录。只有在明确不生成场图时才使用
`--skip-femm`；该选项只生成快速模型、绕组、电磁验收和收敛图片。

## 模型边界

- 快速路径只处理正转速、正转矩的稳态筛选，不模拟逆变器开关、闭环 FOC、
  机械动力学或瞬态工况。
- 当前 FEMM 检查只支持降阶电磁参数到 `5 A 峰值`。配置中的 `15 A 峰值`
  上限仍是外推，必须完成更高电流下的饱和与退磁检查后才能用于定型。
- 给定配置中的附加损耗未标定且为零。因此效率只是铜耗上界，热结果也只是筛选值。
- FEMM 采用二维磁静态求解。端部绕组、端部效应、制造公差、斜槽/斜极、瞬态涡流、
  转子机械应力和退磁裕量不在当前模型中。
- 3 位置电磁验收和缩减收敛数据属于回归检查。设计发布前应提高角度分辨率并与实测校核。

将结果用于下游前，必须检查解析后的配置、两列外推标记、电磁验收 CSV 和收敛 CSV。

## 仓库结构

```text
motor design/
  motor_config.json             权威共享配置
  motor_config.py               严格配置结构、校验和适配器
  design_pmsm.py                带约束的快速 dq 筛选
  femm_spm_template.py          FEMM 几何与分析流程
  generate_readme_assets.py     可复现文档图片生成器
  motor.bat / motor.ps1         推荐的 Windows 启动器
  tests/                        纯 Python 回归测试
  docs/assets/                  当前 README 图片和哈希清单
  output_femm_fix_smoke/        当前缩减验收/收敛 CSV
  output_femm_meshcheck/        单转角 medium/fine 敏感性 JSON 和运行清单
  output_readme_femm/           FEMM 场快照源文件和元数据
```

## FEMM 故障排查

出现 `Invalid class string` 或 `femm.ActiveFEMM` 时，注册 FEMM 后检查：

```bat
reg query HKCR\femm.ActiveFEMM
```

在 PowerShell 中传递列表时，将逗号分隔值放入引号：

```powershell
python .\femm_spm_template.py --rpm-list "300,600,900" --iq-list "1,3,5"
```

分析名称现在会被严格校验；`--analysis` 拼写错误会直接退出并报告有效选项，
不再静默跳过。
