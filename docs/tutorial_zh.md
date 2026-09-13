# 用 MultiAgeClock 批量计算生物年龄

MultiAgeClock 将冻结模型应用于已有表格。每行代表一次测量，各指标占一列，另需实际年龄 `chronological_age`，单位为年。计算在本机完成，不上传表格。

代码、文档、模型参数与权重均按 [MultiAgeClock Noncommercial Research License 1.0](../LICENSE) 提供，仅限非商业科研使用。商业产品、收费计算服务及面向商业利益的研发需另行取得权利人的书面授权。

## 1. 安装和选择模型

```r
install.packages("remotes")
remotes::install_github("Misaka-15134/MultiAgeClock@v0.1.0")
library(MultiAgeClock)

list_models()
model_features("clinical_k17")
```

`clinical_k17` 是论文外部验证使用的 17 指标临床模型。`clinical` 为 47 指标完整版，`olink` 为 2,920 蛋白模型，`nmr` 为 168 指标核磁共振代谢组模型，`integrated` 联合三个完整模态。选用哪套模型取决于已有指标，不能用 K17 代替整合模型要求的完整临床面板。

`clinical_k4` 是四指标探索性模型，只需胱抑素 C（mg/L）、收缩压（mmHg）、腰围（cm）、HbA1c（mmol/mol），另加实际年龄（年）。K4 有独立的五种子权重、标准化参数及年龄标定；它在外部结局读取后作为探索性方案评价，年龄映射只使用 UK Biobank 数据拟合。默认模型仍为 K17。

```r
# 使用 K4 时显式选择模型。
# download_models("clinical_k4")
# k4 <- load_model("clinical_k4")
# k4_ages <- predict_age(example_data("clinical_k4"), k4, id_col = "sample_id")
# score_file("k4_measurements.csv", "k4_ages.csv", model = k4, id_col = "sample_id")
```

K4 的输入列名为 `chronological_age`、`cystatin_c`、`systolic_blood_pressure`、`waist_circumference`、`hba1c`，可另加样本标识列。四项指标均需完整，参考年龄范围为 40–70 岁。

模型权重只需下载一次。它们独立于 R 包保存在用户缓存中；离线计算可指定已下载的模型目录。

```r
download_models("clinical_k17")
fit <- load_model("clinical_k17")

# 也可指定存储位置；同一目录用于下载和加载。
# download_models("clinical_k17", model_dir = "D:/MultiAgeClock-models")
# fit <- load_model("clinical_k17", model_dir = "D:/MultiAgeClock-models")
```

## 2. 准备 K17 表格

先用 `example_data()` 查看列名。这九行是人工生成的软件示例，不是受试者，也不代表真实参考人群。

```r
demo <- example_data("clinical_k17")
head(demo)
write.csv(demo, "k17_input_template.csv", row.names = FALSE)
```

| 列名 | 指标 | 单位 |
|---|---|---|
| `chronological_age` | 实际年龄 | 年 |
| `cystatin_c` | 胱抑素 C | mg/L |
| `systolic_blood_pressure` | 收缩压 | mmHg |
| `waist_circumference` | 腰围 | cm |
| `hba1c` | 糖化血红蛋白 | mmol/mol |
| `diastolic_blood_pressure` | 舒张压 | mmHg |
| `cholesterol` | 总胆固醇 | mmol/L |
| `urate` | 尿酸 | µmol/L |
| `urea` | 尿素 | mmol/L |
| `platelet_count` | 血小板计数 | 10^9/L |
| `mean_corpuscular_volume` | 平均红细胞体积 | fL |
| `standing_height` | 身高 | cm |
| `creatinine` | 肌酐 | µmol/L |
| `glucose` | 葡萄糖 | mmol/L |
| `triglycerides` | 甘油三酯 | mmol/L |
| `haematocrit_percentage` | 红细胞压积 | % |
| `hdl` | 高密度脂蛋白胆固醇 | mmol/L |
| `crp` | C 反应蛋白 | mg/L |

输入测量值保持上述原始单位。不要事先做 z 分数标准化，也不要对 CRP 等临床指标自行取对数。模型内部应用冻结的标准化参数。百分比应按百分数输入，例如红细胞压积 42% 填 `42`。

如实验室采用其他单位，先明确转换。以下例子只适用于列名所指的单位：

```r
# 示例转换；请在自己的数据框中使用实际列名。
# d$creatinine <- d$creatinine_mg_dl * 88.4
# d$glucose <- d$glucose_mg_dl / 18.0182
# d$hba1c <- (d$hba1c_ngsp_percent - 2.15) * 10.929
# d$crp <- d$crp_mg_dl * 10
```

HbA1c 的 `%` 与 `mmol/mol` 不可直接混用；尿素氮（BUN）与尿素也不可仅改列名。测量方法、采血和样本处理差异仍可能影响跨队列适用性。

## 3. 计算数据框

```r
ages <- predict_age(demo, fit, id_col = "sample_id")
ages[, c("sample_id", "chronological_age", "ba", "raw_gap", "baa")]
```

返回结果保留输入行顺序：

- `row_id`：原表行号。
- `sample_id`：通过 `id_col` 指定的样本标识。
- `ba`：生物年龄，单位为年。
- `raw_gap`：生物年龄减实际年龄。
- `baa`：生物年龄减冻结参考年龄曲线的期望值，单位为年。
- `age_outside_reference`：实际年龄是否超出该模型的参考年龄范围。

`raw_gap` 与 `baa` 的定义不同。BAA 校正参数始终来自冻结模型，包不会用本次输入的样本重新拟合年龄曲线。单人计算也采用同一规则。K17 的参考年龄范围为 40–70 岁；超出范围会标记，结果不作截断。

如果表格列名与模型字段不同，可以显式映射，不必修改原表：

```r
renamed <- demo
names(renamed)[names(renamed) == "chronological_age"] <- "age_years"
ages <- predict_age(renamed, fit, id_col = "sample_id",
                     column_map = c(chronological_age = "age_years"))
```

首版要求所选面板完整。报错会列出缺失列或含缺失值的列及行号；不会以另一套面板替代，也不会调用训练样本填补缺失值。

## 4. 直接计算 CSV、TSV 和 Excel

```r
score_file("k17_input_template.csv", "k17_ages.csv",
           model = fit, id_col = "sample_id")

# Excel 读写依赖只需安装一次。
install.packages(c("readxl", "writexl"))
# score_file("measurements.xlsx", "ages.xlsx", model = fit,
#            id_col = "sample_id", sheet = 1)
```

原文件保持不变，输出文件包含预测列和样本标识。已存在的输出文件会报错，请另取文件名。CSV/TSV 的文本编号会保留前导零。Excel 中需要保留前导零的编号应预先存成文本，数字单元格的显示格式不能恢复已经丢失的零。

## 5. 其他模态和整合模型

```r
# 先确认指标和单位，再下载所需权重。
model_features("olink")
model_features("nmr")

# download_models("olink")
# olink_ages <- predict_age(example_data("olink"), "olink", id_col = "sample_id")

# download_models("integrated")  # 同时下载三套完整模态的权重
# integrated_ages <- predict_age(example_data("integrated"), "integrated",
#                                id_col = "sample_id", batch_size = 256)
```

Olink 输入为与训练流程一致的 NPX（标准化蛋白表达量）log2 数值，不要再取一次 log2。NMR 输入字典同时给出冻结 UK Biobank 字段名、Nightingale 名称和单位，可据此使用 `column_map`。

整合模型要求三套模态在同一行对应同一个样本。多张表应先按唯一标识核对并合并，不能直接按行号拼接。所有模型均保留五种子集成；`batch_size` 只控制每批计算的行数。

## 6. 用 R 调用 PyTorch

默认 `backend = "r"`，无需 Python。大表格或 GPU 计算可以选择 `backend = "pytorch"`，两者使用相同的模型资产和年龄标定参数。

在已有 Python 环境中安装 NumPy 和 PyTorch，调用时传入该环境的 Python 路径。PyTorch 的 CPU、CUDA 安装命令按[官方安装页面](https://pytorch.org/get-started/locally/)选择。

```r
# 改成已经安装 numpy、torch 的 Python 可执行文件。
# py <- "D:/python-envs/multiageclock/python.exe"
# ages <- predict_age(demo, fit, id_col = "sample_id", backend = "pytorch", python = py)
# GPU 环境已正确配置时：
# ages <- predict_age(demo, fit, backend = "pytorch", device = "cuda", python = py)
# score_file("measurements.csv", "ages_gpu.csv", model = fit,
#            backend = "pytorch", device = "cuda", python = py, id_col = "sample_id")
```

选择 PyTorch 后，缺少依赖或设备不可用会报错。R 启动独立 Python 进程调用随包提供的推理模块，输入通过本机临时二进制文件传递，计算结束后清理。PyTorch 使用 float32 全精度前向计算。论文历史 GPU 批量导出使用混合精度，逐样本值可能存在小幅数值差异，不能要求与历史文件逐位相同。

这些输出用于科研。计算结果的一致性不等于在任意检测平台或人群中完成了外部验证，也不能直接据此判断个人疾病或制定治疗。
