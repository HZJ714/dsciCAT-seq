#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
preflight_checks.py — 首次公开推送前的硬闸（在仓库根目录跑）
=============================================================
每条闸输出 PASS / WARN / FAIL；FAIL≠0 时 push_preflight.sh 以非零码退出。
闸的设计原则：只检查"公开后无法撤回"的东西——身份、绝对路径、误传的大文件、
以及 fresh clone 能不能跑起来。

G1  身份 token（用户名 / 同事名 / 集群项目号 / 机器账号）
G2  作者机绝对路径（/data/…、/home/…），豁免名单 dev/PATH_LEAK_ALLOWLIST.txt
G3  被跟踪的符号链接（指向仓外一律 FAIL）
G4  被跟踪的历史垃圾（_superseded 二进制、regen.log、.ipynb_checkpoints、*.bak）
G5  体积（>50 MB FAIL，10–50 MB WARN）
G6  语法（python / notebook JSON+code cell / R / shell / WDL 括号）
G7  本地模块 import 是否都在仓内（fresh clone 的 ModuleNotFoundError）
G8  环境变量契约：默认值里不得出现绝对路径；用到的 DSCI_* 必须在 README 有行
G9  SI 表号是否还是旧编号（Table S9+、"S6a–S6g + Tables"、ten-panel 等）
G10 提交身份 / LICENSE / 提交说明里的机器路径（只报，不改）
用法: python3 dev/preflight_checks.py [--verbose]
"""
import ast, json, os, re, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.chdir(ROOT)
VERBOSE = '--verbose' in sys.argv
results = []


def gate(name, kind, items, hint=''):
    results.append([name, kind, items, hint])


def tracked():
    out = subprocess.run(['git', 'ls-files', '-z'], capture_output=True, text=True).stdout
    return [f for f in out.split('\0') if f]


FILES = tracked()
TEXT_EXT = ('.py', '.R', '.r', '.sh', '.bash', '.ipynb', '.md', '.txt', '.json', '.yaml',
            '.yml', '.csv', '.tsv', '.log', '.wdl')


def read(f):
    try:
        return open(f, encoding='utf-8', errors='ignore').read()
    except Exception:
        return ''


def readable():
    for f in FILES:
        if f.endswith(TEXT_EXT):
            yield f, read(f)


ALLOW = {l.strip() for l in read('dev/PATH_LEAK_ALLOWLIST.txt').split('\n')
         if l.strip() and not l.startswith('#')}
G1_ALLOW = ALLOW | {'dev/identity_tokens.example.txt'}   # 公开树里不再有明文清单，只豁免扫描器自身

# ── G1 身份 token ─────────────────────────────────────────────────────
# 2026-11-25：仓库转 public 并准备 Zenodo 归档，公开树里不能再出现明文身份串，
# 因此 G1 的匹配来源拆成两半：
#   dev/identity_tokens.example.txt  结构型正则（入库，任何人可复查规则本身）
#   dev/identity_tokens.txt          逐字精确清单（本工作环境私有，.gitignore 排除）
PATTERNS = [l.strip() for l in read('dev/identity_tokens.example.txt').split('\n')
            if l.strip() and not l.startswith('#')]
TOKENS = []
if os.path.exists('dev/identity_tokens.txt'):
    TOKENS = [l.strip() for l in read('dev/identity_tokens.txt').split('\n')
              if l.strip() and not l.startswith('#')]
ID_RX = re.compile('|'.join(PATTERNS + [re.escape(t) for t in TOKENS]))
hits = []
for f, t in readable():
    if f in G1_ALLOW:
        continue
    for m in set(ID_RX.findall(t)):
        hits.append('%s: %s' % (f, m))
gate('G1 身份 token 清零', 'FAIL' if hits else 'PASS', hits,
     '%d 条结构模式 + 本地精确清单 %d 条；公开树里不得有用户名/同事名/项目号/机房前缀'
     % (len(PATTERNS), len(TOKENS)))

# ── G2 机器绝对路径 ───────────────────────────────────────────────────
PATH_RX = re.compile(r'/data/(?:work|input|users|temp)\b|/data/users/')
hits = []
for f, t in readable():
    if f in ALLOW:
        continue
    n = len(PATH_RX.findall(t))
    if n:
        hits.append('%s (%d 处)' % (f, n))
gate('G2 作者机绝对路径清零', 'FAIL' if hits else 'PASS', hits,
     '跑 00_Figure_Tools/port_paths.py 自动改写；豁免写进 dev/PATH_LEAK_ALLOWLIST.txt')

# ── G3 符号链接 ───────────────────────────────────────────────────────
hits = []
for f in FILES:
    if os.path.islink(f):
        tgt = os.readlink(f)
        if not tgt.startswith(('fig', '../', './')):
            hits.append('%s -> %s' % (f, tgt))
gate('G3 无仓外符号链接', 'FAIL' if hits else 'PASS', hits,
     '改用 figX/data/link_data.sh + data/README.md')

# ── G4 历史垃圾 ───────────────────────────────────────────────────────
JUNK = re.compile(r'(_superseded[^/]*/(?!README\.md)|\.ipynb_checkpoints/|regen\.log$|\.bak($|\.)|Rplots\.pdf$)')
hits = [f for f in FILES if JUNK.search(f)]
gate('G4 不留历史垃圾/缓存', 'FAIL' if hits else 'PASS', hits,
     '_superseded 只保留 README.md，日志与 .bak 不入库')

# ── G5 体积 ───────────────────────────────────────────────────────────
big, mid = [], []
for f in FILES:
    try:
        sz = os.path.getsize(f)
    except OSError:
        continue
    if sz > 50 * 1024 * 1024:
        big.append('%s (%.1f MB)' % (f, sz / 1048576.0))
    elif sz > 10 * 1024 * 1024:
        mid.append('%s (%.1f MB)' % (f, sz / 1048576.0))
gate('G5 无 >50 MB 文件', 'FAIL' if big else 'PASS', big, 'GitHub 会拒收/变慢')
gate('G5b 10–50 MB 文件盘点', 'WARN' if mid else 'PASS', mid, '确认每个都值得进公开历史')

# ── G6 语法 ───────────────────────────────────────────────────────────
hits = []
HAS_R = subprocess.run(['which', 'Rscript'], capture_output=True, text=True).stdout.strip()
for f in FILES:
    if not os.path.exists(f):
        continue
    if f.endswith('.py') and os.path.getsize(f) < 3_000_000:
        try:
            ast.parse(read(f))
        except SyntaxError as e:
            hits.append('PY %s: %s' % (f, str(e)[:70]))
    elif f.endswith('.ipynb'):
        try:
            nb = json.loads(read(f))
        except Exception as e:
            hits.append('NBJSON %s: %s' % (f, str(e)[:60]))
            continue
        lang = (nb.get('metadata', {}).get('kernelspec', {}) or {}).get('language', 'python').lower()
        for i, c in enumerate(nb.get('cells', [])):
            if c.get('cell_type') != 'code' or lang in ('r', 'ir_kernel'):
                continue
            try:
                ast.parse(''.join(c.get('source', [])))
            except SyntaxError as e:
                hits.append('NBCODE %s cell%d: %s' % (f, i, str(e)[:60]))
    elif f.endswith(('.sh', '.bash')):
        if subprocess.run(['bash', '-n', f], capture_output=True).returncode:
            hits.append('SH %s' % f)
    elif f.endswith(('.R', '.r')) and HAS_R:
        if subprocess.run(['Rscript', '-e', 'invisible(parse(%r))' % f], capture_output=True).returncode:
            hits.append('R %s' % f)
gate('G6 语法/JSON 全通过', 'FAIL' if hits else 'PASS', hits,
     'R 文件在缺少 Rscript 的机器上跳过，推送前请在有 R 的环境补跑')

# ── G7 本地 import 自洽 ───────────────────────────────────────────────
mods = {os.path.splitext(os.path.basename(f))[0] for f in FILES if f.endswith('.py')}
hits = []
IMP = re.compile(r'^\s*(?:from|import)\s+([a-zA-Z_][\w.]*)', re.M)
for f in FILES:
    if f.endswith('.py'):
        texts = [read(f)]
    elif f.endswith('.ipynb'):
        try:
            nb = json.loads(read(f))
        except Exception:
            continue
        texts = [''.join(c.get('source', [])) for c in nb.get('cells', []) if c.get('cell_type') == 'code']
    else:
        continue
    for t in texts:
        for m in IMP.finditer(t):
            base = m.group(1).split('.')[0]
            if base in mods and not os.path.exists(os.path.join(os.path.dirname(os.path.abspath(f)), base + '.py')):
                # 同目录没有 → 必须在 sys.path 铺垫里能找到（近似判断：文件里出现该目录名或 _ROOT 铺垫）
                if not re.search(r'00_Figure_Tools|fig5/assembly|_ROOT|sys\.path', t):
                    hits.append('%s imports %s（仓内有同名模块但不在同目录，未见 sys.path 铺垫）' % (f, base))
gate('G7 fresh clone 不缺模块', 'FAIL' if hits else 'PASS', hits,
     '缺模块 = 读者 clone 后第一步就报错')

# ── G8 环境变量契约 ───────────────────────────────────────────────────
hits = []
DEFAULT_ABS = re.compile(r"(?:environ\.get|Sys\.getenv)\(\s*['\"][A-Z0-9_]+['\"]\s*,\s*['\"]/(?!tmp\b|dev/null)")   # 默认值以 / 开头 = 绝对路径
for f, t in readable():
    if f in ALLOW:
        continue
    for m in DEFAULT_ABS.finditer(t):
        hits.append('%s: %s' % (f, m.group(0)[:70]))
vars_used = set()
for f, t in readable():
    vars_used.update(re.findall(r'\bDSCI_[A-Z_]+|FIG1_DATA|\bPREPROC_[A-Z]+\b', t))
documented = read('README.md') + read('00_Figure_Tools/env_check.py') + read('00_Preprocessing/README.md')
undoc = sorted(v for v in vars_used if v not in documented)
gate('G8a getenv/environ.get 默认值不含绝对路径', 'FAIL' if hits else 'PASS', hits)
gate('G8b 每个环境变量都在文档里有一行', 'WARN' if undoc else 'PASS', undoc,
     '新增变量请同步 README《Data layout》或 env_check.py 的表')

# ── G9 SI 表号口径 ────────────────────────────────────────────────────
STALE = [(r'Table S(?:9|1[0-9])\b', '出现了不存在的 Table S9+'),
         (r'S6a[–-]S6g\s*\+\s*Tables', '旧的 "S6a–S6g + Tables S6/S7" 说法'),
         (r'ten[- ]panel[^\n]{0,30}S6|S6[^\n]{0,30}ten[- ]panel',
          '十格 S6 的旧描述（S6 现行 v16/v18 = 7 格；Figure S7 v30 本就是 ten panel PDFs，不在此列）'),
         (r'eRSS.{0,40}Table S6|Table S6.{0,40}eRSS', 'eRSS 现在叫 Table S5')]
hits = []
SELF = 'dev/preflight_checks.py'        # 本闸的正则字面量本身就是被检串
for f, t in readable():
    if f == SELF or f.startswith(('legends/figure_legends.md', 'legends/README.md', 'dev/TABLES.md',
                                  'dev/TIDY_RECEIPT', 'figS6/STATUS.md')) or '/_superseded' in f:
        continue
    for rx, why in STALE:
        for m in re.finditer(rx, t, re.I):
            line = t[:m.start()].count('\n') + 1
            hits.append('%s:%d  %s  ← %s' % (f, line, m.group(0)[:50], why))
gate('G9 无旧编号残留', 'WARN' if hits else 'PASS', hits, '口径见 dev/TABLES.md')

# ── G10 提交身份 / LICENSE / 提交说明 ─────────────────────────────────
hits = []
try:
    email = subprocess.run(['git', 'config', 'user.email'], capture_output=True, text=True).stdout.strip()
    name = subprocess.run(['git', 'config', 'user.name'], capture_output=True, text=True).stdout.strip()
    hits.append('commit 身份: %s <%s>' % (name, email))
except Exception:
    pass
for l in subprocess.run(['git', 'log', '--format=%an <%ae> | %s'], capture_output=True, text=True).stdout.split('\n'):
    if l and re.search(ID_RX, l):        # 提交说明/作者字段里不得有身份串（按形态匹配）
        hits.append('提交说明含机器路径: ' + l[:100])
lic = read('LICENSE')
m = re.search(r'Copyright[^\n]*', lic)
if m:
    hits.append('LICENSE: ' + m.group(0)[:80])
gate('G10 身份与署名（需作者拍板）', 'WARN', hits,
     '已压平为单一干净提交（2026-11-11）；剩余 2 项 = noreply 身份 + LICENSE 法律口径，作者拍板')

# ── 汇总 ──────────────────────────────────────────────────────────────
nfail = sum(1 for r in results if r[1] == 'FAIL')
nwarn = sum(1 for r in results if r[1] == 'WARN')
for name, kind, items, hint in results:
    print('%-4s %-46s %d 项' % (kind, name, len(items)))
    if (items and VERBOSE) or kind == 'FAIL':
        lim = 12 if kind == 'FAIL' else 20
        for x in items[:lim]:
            print('      ·', x)
        if len(items) > lim:
            print('      · …另有 %d 项' % (len(items) - lim))
    if hint and kind != 'PASS':
        print('      ↳', hint)
print('\n合计 PASS %d / WARN %d / FAIL %d' % (
    sum(1 for r in results if r[1] == 'PASS'), nwarn, nfail))
sys.exit(1 if nfail else 0)
