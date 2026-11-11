import os
import time, numpy as np, pandas as pd, glob, anndata as ad
from scipy.stats import pearsonr, spearmanr
t0 = time.time(); log = lambda *a: print("[%5ds]" % (time.time() - t0), *a, flush=True)
import scdrs, scdrs.method, scdrs.pp
log("scdrs from %s" % scdrs.__file__)
a = ad.read_h5ad(os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/3k_dsci_oligo_subset.h5ad'))
log("矩阵 %s type=%s X.max=%.2f" % (a.shape, type(a.X).__name__, a.X.max()))
scdrs.pp.preprocess(a); log("preprocess 完成，累计 %.0fs" % (time.time() - t0))
line = [l.rstrip("\n") for l in open(os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/scdrs_ms_out/gs_files/gs_file/magma_10kb_top1000_zscore.74_traits.rv1.gs')) if l.startswith("PASS_Multiple_sclerosis\t")][0]
w = {}
for tok in line.split("\t")[1].split(","):
    if ":" in tok:
        g, z = tok.rsplit(":", 1)
        try: w[g.strip().upper()] = float(z)
        except ValueError: pass
pos = {g.upper(): i for i, g in enumerate(a.var_names)}
gl = [g for g in w if g in pos]
log("MS 基因集命中 %d/%d  权重中位 |z|=%.2f" % (len(gl), len(w), np.median(list(map(abs, w.values())))))
r = scdrs.method.score_cell(a, gene_list=gl, n_ctrl=200, random_seed=20261030, verbose=False)
log("score_cell(全基因集) 返回列=%s 累计 %.0fs" % (list(r.columns)[:8], time.time() - t0))
d = pd.read_csv(glob.glob(os.path.join(os.environ.get('DSCI_FY_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/scdrs_traits_draft/full74/*/*Multiple_sclerosis*full_score.gz'))[0], sep="\t", index_col=0)
inter = r.index.intersection(d.index)
log("交集核 %d（交付文件 %d，本次 %d）" % (len(inter), len(d), len(r)))
for c in r.columns:
    if c in d.columns:
        log("   复现校验 %-12s Pearson=%.4f Spearman=%.4f  交付中位=%.4f 本次中位=%.4f"
            % (c, pearsonr(r.loc[inter, c], d.loc[inter, c])[0],
               spearmanr(r.loc[inter, c], d.loc[inter, c]).correlation,
               d.loc[inter, c].median(), r.loc[inter, c].median()))
t1 = time.time()
r1 = scdrs.method.score_cell(a, gene_list=[gl[0]], n_ctrl=200, random_seed=20261030, verbose=False)
log("单基因 %s 耗时 %.1fs → 890 基因串行 %.0f min（32 核并行约 %.1f min）"
    % (gl[0], time.time() - t1, (time.time() - t1) * 890 / 60, (time.time() - t1) * 890 / 60 / 24))
