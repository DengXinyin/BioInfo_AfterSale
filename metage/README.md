# metage:v2.87 可视化函数库

本目录整理自 Docker 镜像
`192.168.30.202:23099/metage_megahit/metage:v2.87`，镜像 ID 为
`sha256:8a920c76b8179a855ff56fe1fdbf71a002c9623427296d880c6c361d7bcad698`。
原流程的 39 个 R 脚本（4,482 行）包含大量 tax/func 重复实现和固定目录副作用；这里按图型合并为可直接 `source()` 的函数。

## 目录

```text
metage/
├── R/                 # 可复用函数，输入数据、返回图或模型
├── docs/              # 与功能域对应的中文说明
├── scripts/           # 虚构数据示例生成脚本
├── tests/             # 不依赖客户数据的 smoke test
├── README.md
└── SOURCE_INVENTORY.md
```

## 快速使用

```r
source("metage/R/utils.R")
source("metage/R/02_composition_heatmap.R")

# 行为物种/功能，列为样本
p <- metage_abundance_bar(abundance_matrix, top_n = 15)
metage_save_plot(p, "Output/composition.pdf", width = 10, height = 7)
```

批量载入全部函数：

```r
invisible(lapply(list.files("metage/R", full.names = TRUE), source))
```

## 整理约定

- 丰度矩阵统一为“行=特征、列=样本”；命名分组向量按样本名对齐。
- 绘图函数不写死项目目录，不自动输出 Excel/HTML，也不包含客户数据。
- 已有坐标、P 值、STAMP 和 LEfSe marker 表可直接传入，避免不必要的重算。
- 统计函数明确返回模型；raw P 与 adjusted P 由参数/列名明确区分。
- 默认字体依次尝试 Times New Roman、Arial，再回退系统 `sans`；正式输出使用 PDF/cairo。
- PNG 仅由虚构数据生成用于视觉检查，不作为本目录的正式交付物。

运行 smoke test：

```bash
/root/anaconda3/envs/r/bin/Rscript metage/tests/smoke_test.R
```

生成临时验收图：

```bash
/root/anaconda3/envs/r/bin/Rscript metage/scripts/generate_examples.R /tmp/metage_examples
```

完整来源映射见 [SOURCE_INVENTORY.md](SOURCE_INVENTORY.md)。
