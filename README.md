```bash
bayesian-survival-project/
│
├── data/                  # 原始数据 (小数据可放这里，大数据用链接)
│   ├── lung.csv
│   └── ...
│
├── code/                  # R 脚本和 Rmd
│   ├── 01_cleaning.R
│   ├── 02_exponential_fit.R
│   ├── 03_plots.R
│   └── ...
│
├── results/               # 运行 R 后生成的所有图、表
│   ├── survival_curve.pdf
│   ├── traceplot_lambda.png
│   ├── tables/
│   │   ├── param_summary.csv
│   │   └── ...
│   └── ...
│
├── paper/                 # 你的 Overleaf 主文件（LaTeX）
│   ├── main.tex
│   ├── sections/
│   │   ├── introduction.tex
│   │   ├── methods.tex
│   │   ├── results.tex
│   │   ├── discussion.tex
│   ├── figures/           # 用于 Overleaf 引用的图（从 results 复制过来）
│   │   ├── survival_curve.pdf
│   │   ├── traceplot_lambda.png
│   │   └── ...
│   ├── refs.bib
│   ├── statsmsc.cls
│   ├── supplementary.tex
│   └── ...
│
├── README.md              # 项目描述 + 如何 reproduce
└── .gitignore             # 忽略不必要的中间文件
```
