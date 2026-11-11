# fig5/data/ — 输入对象挂钩目录（仓库不放真身）

Figure 5 的 notebook 只按相对路径 `data/<name>` 读输入。这些对象都太大或是上游产物，
**不进公开树**；clone 之后跑一次挂钩脚本即可：

```bash
export DSCI_WORK_ROOT=<你的分析盘 online 根>       # 见 ../../README.md
bash fig5/data/link_data.sh                          # 已存在则跳过，--force 重建
```

| 文件 | 大小 | md5(前 12) | 真身（由环境变量解析） | 来源 / 再生成方式 |
|---|---|---|---|---|
| `10x_subset_harmony.RDS` | 752 MB | `033bf057d5d4` | `$DSCI_FY_ROOT/tenx/` | 10x Multiome 成年小鼠脑细胞核 Seurat 对象（Harmony 校正后子集）。原始数据 CNGBdb **CNP0010263**；重建走 `00_Preprocessing/` + `fig2/`（图 2 的整合流程写出同一对象） |
| `HT_harmony.RDS` | 1.9 GB | `571d09b51601` | `$DSCI_FY_ROOT/HT/` | HT-scCAT-seq 皮层/海马对象（Harmony 校正后子集），同上 |
| `scplus_obj_final.pkl` | 7.5 GB | `7276a15fda7a` | `$DSCI_FY_ROOT/oligo/3kdsci/GRN/scplus/` | SCENIC+ 最终 pyscenic 对象（regulons / eRegulons / TF 活性）。由 `figS5` 上游 SCENIC+ 流程产出，见 `figS6/upstream/`、`figS5/README.md` |
| `NIHMS70326-supplement-Table_S1.xlsx` | 56 KB | `975354087775` | `$DSCI_WORK_ROOT/` | Yao et al. 全脑图谱的作者补充表（细胞类型命名与 marker 依据，正文引文 42/43）。期刊公开补充材料，非本研究生成 |

说明
----
- `link_data.sh` 只做 `ln -sfn`，不复制数据；找不到真身会逐个报 `MISS` 并以非零码退出，
  方便 CI/自查。
- notebook 里出现 `data/scplus_obj_final.pkl` 等相对路径，就是指望这个目录；
  未挂钩就执行会报 file-not-found，不会静默出错。
- 这些对象的绝对路径**不允许**再写回仓库文本里；若你新增脚本需要它们，请改用
  `DSCI_*` 环境变量（`00_Figure_Tools/port_paths.py` 会把违规写法抓出来）。
