# figS5/data/ — 输入对象挂钩目录（仓库不放真身）

Figure S5 的 notebook 只按相对路径 `data/<name>` 读输入。真身都在分析盘上（太大或属上游产物），
clone 之后跑一次挂钩脚本即可：

```bash
export DSCI_WORK_ROOT=<你的分析盘 online 根>       # 见 ../../README.md
bash figS5/data/link_data.sh                         # 已存在则跳过，--force 重建
```

| 文件 | 大小 | md5(前 12) | 真身（由环境变量解析） | 来源 / 再生成方式 |
|---|---|---|---|---|
| `scplus_obj_final.pkl` | 7.5 GB | `7276a15fda7a` | `$DSCI_FY_ROOT/oligo/3kdsci/GRN/scplus/` | SCENIC+ 最终对象（与 `fig5/data/` 同一个文件）；上游流程见 `figS6/upstream/` |
| `eregulon_stability_table.csv` | 10 KB | `720c2b576734` | `$DSCI_SCENIC_ROOT/stability/` | 降采样稳定性账本（eRegulon 在 1k/5k/10k 细胞下的重现情况），由 `figS5/assembly/` 的稳定性步骤写出 |
| `nrxn3_mef2c_enhancers.fa` | 20 KB | `7a87b9205861` | `$DSCI_SCENIC_ROOT/` | Mef2c eRegulon 命中、位于 Nrxn3 上游 1.05 Mb 窗口内的四个增强子序列（JASPAR MA0052.4） |
| `nrxn3_tracks.npz` | 65 KB | `18cbc3faf1fa` | `$DSCI_FIG_ROOT/scripts/cache/` | 图 S5C 用的按亚型聚合可及性轨迹缓存 |
| `nrxn3_expr.npz` | 56 KB | `556424f352da` | `$DSCI_FIG_ROOT/scripts/cache/` | 图 S5C 用的按亚型聚合表达量缓存 |
| `NIHMS70326-supplement-Table_S1.xlsx` | 56 KB | `975354087775` | 仓内 `fig5/data/` | 指向 `fig5/data/` 的同名文件（同一份作者补充表，两处都要用） |

说明
----
- 最后一条是仓内自链：只要 `fig5/data/link_data.sh` 跑过，`figS5/data/` 这条就自动可用。
- 想核对真身是否被换过：`md5sum data/<name>` 对比上表前 12 位即可。
