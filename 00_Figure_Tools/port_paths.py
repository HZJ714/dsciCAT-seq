#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
port_paths.py — 把作者机路径从仓库里彻底摘掉（v3, 2026-11-11）
=============================================================
规矩（v3 起收紧）：**仓库里不允许出现任何机器绝对路径**，连“默认值”也不行。
所有数据根目录一律走环境变量，未设置时退化为 `<DSCI_WORK_ROOT>/…`，由用户按
README《Environment variables》一节 export；仓内文件之间一律用相对路径。

    DSCI_WORK_ROOT    分析机的 online 工作盘根        必须 export（无默认）
    DSCI_REPO_ROOT      本仓库检出                      默认 '.'（在仓库根目录运行）
    DSCI_FY_ROOT        FY 数据根                       <WORK>/FY
    DSCI_FIG_ROOT       figures 工作盘                  <WORK>/figures
    DSCI_PYLIBS         持久 PyMuPDF                    <WORK>/pylibs
    DSCI_FONTS          Arial TTF 目录                  <WORK>/fonts
    DSCI_FONTS_NAMED    按字重命名的字体                 <WORK>/fonts_named
    DSCI_SEURAT_ROOT    Seurat 产物                     <WORK>/pi/plots-seurat
    DSCI_SCENIC_ROOT    SCENIC+ 产物                    <WORK>/pi/plots-scenic+
    DSCI_SCENIC_SCRIPTS SCENIC+ 脚本                    <WORK>/pi/scripts-scenic+
    DSCI_PROBE_ROOT     routeA 探针产物                  <WORK>/routeA_probe
    DSCI_REF_ROOT       referrence                      <WORK>/referrence
    DSCI_MOTIF_ROOT     aertslab motif collection       <WORK>/aertslab_motif_colleciton
    DSCI_SCRIPTS_ROOT   分析脚本                        <WORK>/scripts
    FIG1_DATA           原始数据（ResultData/<user>）    必须 export（无默认）
    DSCI_REFDATA_ROOT   参考注释（ReferenceData）        必须 export（无默认）

改写形式（幂等；三种语言各一套）：
    Python  '/data/work/FY/oligo/x.h5ad'
            -> os.path.join(os.environ.get('DSCI_FY_ROOT',
               os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'FY')), 'oligo/x.h5ad')
    R       -> file.path(Sys.getenv("DSCI_FY_ROOT",
               file.path(Sys.getenv("DSCI_WORK_ROOT", ""), "FY")), "oligo/x.h5ad")
    shell   -> ${DSCI_FY_ROOT:-${DSCI_WORK_ROOT:-}/FY}/oligo/x.h5ad
    已经 port 过但默认值仍是机器路径的 getenv(V,'/data/…') 会被就地换成上面的形式；
    被生成出来的**别的语言**的代码字符串（Python cell 里写 R、写 shell）由兜底规则处理；
    文档 / notebook 输出 / markdown cell → 尖括号占位符 <repo> <work> <raw> <refdata>。

⚠ 三条老坑仍然成立：
  1) 本脚本自身不能被自己 port（MAP 里必须保留机器路径字面量），列在 SKIP；
     这也是推送前检查里唯一被豁免的文件（见 dev/PATH_LEAK_ALLOWLIST.txt）。
  2) 幂等性只靠“文本里已无 /data 绝对路径”这一事实判断；不要改用“含 DSCI_* 变量名就跳过”，
     v1 就是这样漏掉了所有“变量名已加、默认值仍是机器路径”的文件。
  3) 会改写别人源码的脚本（figs5_repanels.py 用 src.replace("OUTDIR = '…'", …) 生成 cell）
     里的字面量是**搜索串**，一起改会破坏匹配 ⇒ 列进 SKIP，由人工处理。
用法: python3 port_paths.py <文件或目录>... [--dry]
"""
import os, re, sys, json

WORK = 'DSCI_WORK_ROOT'
REPO = 'DSCI_REPO_ROOT'

# (机器路径前缀, 环境变量, 相对 DSCI_WORK_ROOT 的子路径 or None=不派生)
MAP = [
    ('/data/work/GitHub',                      REPO,                    None),
    ('/data/work/pi/plots-seurat',             'DSCI_SEURAT_ROOT',      'pi/plots-seurat'),
    ('/data/work/pi/plots-scenic+',            'DSCI_SCENIC_ROOT',      'pi/plots-scenic+'),
    ('/data/work/pi/scripts-scenic+',          'DSCI_SCENIC_SCRIPTS',   'pi/scripts-scenic+'),
    ('/data/work/pi',                          'DSCI_PI_ROOT',          'pi'),
    ('/data/work/pylibs',                      'DSCI_PYLIBS',           'pylibs'),
    ('/data/work/fonts_named',                 'DSCI_FONTS_NAMED',      'fonts_named'),
    ('/data/work/fonts',                       'DSCI_FONTS',            'fonts'),
    ('/data/work/routeA_probe',                'DSCI_PROBE_ROOT',       'routeA_probe'),
    ('/data/work/referrence',                  'DSCI_REF_ROOT',         'referrence'),
    ('/data/work/aertslab_motif_colleciton',   'DSCI_MOTIF_ROOT',       'aertslab_motif_colleciton'),
    ('/data/work/scripts',                     'DSCI_SCRIPTS_ROOT',     'scripts'),
    ('/data/work/figures',                     'DSCI_FIG_ROOT',         'figures'),
    ('/data/work/FY',                          'DSCI_FY_ROOT',          'FY'),
    ('/data/input/Files/ResultData/HZJ',       'FIG1_DATA',             None),
    ('/data/input/Files/ResultData',           'DSCI_RESULTDATA_ROOT',  None),
    ('/data/input/Files/ReferenceData',        'DSCI_REFDATA_ROOT',     None),
    ('/data/temp',                             'TMPDIR',                None),
    ('/data/work',                             WORK,                  None),
]
MAP = sorted(MAP, key=lambda kv: -len(kv[0]))              # 长前缀优先
SKIP = {'port_paths.py', 'figs5_repanels.py'}             # 见头注坑 1 / 坑 3
PLACEHOLDER = [('/data/work/GitHub', '<repo>'),
               ('/data/input/Files/ResultData/HZJ', '<raw>'),
               ('/data/input/Files/ResultData', '<raw>'),
               ('/data/input/Files/ReferenceData', '<refdata>'),
               ('/data/temp', '<tmp>'),
               ('/data/work', '<work>')]
PATH_RE = r'/data/(?:work|input|temp|users)[A-Za-z0-9_./+*@-]*'
PATH_RE_RX = re.compile(PATH_RE)
# 家目录一律按 `/home/<任意用户>` 识别，工具源码里不留任何真实账号名
HOME_ANY = re.compile(r'/home/[A-Za-z0-9_.+-]+(/[A-Za-z0-9_.+/-]*)?')
HOME_LIT = re.compile(r"""(['"])/home/[A-Za-z0-9_.+-]+([^'"]*)\1""")
USERS_RX = re.compile(r'/data/users/[^\s\'"`)\]},;]*')      # 文档里的“另一环境”路径（可能带 …）


def prefix_of(path):
    """返回 (匹配到的前缀, 变量名, online 子路径, 剩余尾串)。"""
    for pref, var, sub in MAP:
        if path == pref or path.startswith(pref + '/'):
            return pref, var, sub, path[len(pref):].lstrip('/')
    m = USERS_RX.match(path)
    if m and '/online' in m.group(0):
        head = m.group(0)[:m.group(0).index('/online') + len('/online')]
        return head, WORK, None, path[len(head):].lstrip('/')
    if m:
        return m.group(0), WORK, None, ''
    return None, None, None, None


def to_placeholder(path):
    m = USERS_RX.match(path)
    if m:
        txt = m.group(0)
        if '/online' in txt:
            return '<work>' + path[len(txt[:txt.index('/online') + 7]):]
        return '<work>'
    for pref, ph in PLACEHOLDER:
        if path.startswith(pref):
            return ph + path[len(pref):]
    return '<path>' + path


def py_base(var, sub):
    if var == REPO:
        return "os.environ.get('%s', '.')" % REPO
    if sub is None:
        return "os.environ.get('%s', '')" % var
    return "os.environ.get('%s', os.path.join(os.environ.get('%s', ''), '%s'))" % (var, WORK, sub)


def r_base(var, sub):
    if var == REPO:
        return 'Sys.getenv("%s", ".")' % REPO
    if sub is None:
        return 'Sys.getenv("%s", "")' % var
    return 'Sys.getenv("%s", file.path(Sys.getenv("%s", ""), "%s"))' % (var, WORK, sub)


def sh_base(var, sub):
    if var == REPO:
        return '${%s:-.}' % REPO
    if sub is None:
        return '${%s:-}' % var
    return '${%s:-${%s:-}/%s}' % (var, WORK, sub)


def expr_for(path, kind):
    _p, var, sub, rest = prefix_of(path)
    if var is None:
        return None
    base = {'py': py_base, 'r': r_base, 'sh': sh_base}[kind](var, sub)
    if not rest:
        return base
    if kind == 'py':
        return "os.path.join(%s, '%s')" % (base, rest)
    if kind == 'r':
        return 'file.path(%s, "%s")' % (base, rest)
    return '%s/%s' % (base, rest)


# ── 兜底：代码字符串里嵌的“另一种语言”的写法 ────────────────────────────
GETENV_ANY = re.compile(r"""((?:os\.)?environ\.get|Sys\.getenv)\(\s*(['"])([A-Z0-9_]+)\2\s*,\s*(['"])((?:%s))\4\s*\)""" % PATH_RE)
SH_DEFAULT = re.compile(r"""\$\{([A-Z0-9_]+):-((?:%s)[^}]*)\}""" % PATH_RE)


def _sweep_getenv(text):
    n = 0

    def f(m):
        nonlocal n
        fn, var, path = m.group(1), m.group(3), m.group(5)
        ispy = fn.endswith('environ.get')
        _p, dvar, dsub, drest = prefix_of(path)
        if dvar is None or dvar == var:
            e = expr_for(path, 'py' if ispy else 'r')
            if e is None:
                return m.group(0)
            n += 1
            return e
        n += 1
        base = (py_base if ispy else r_base)(dvar, dsub)
        call = ("os.environ.get('%s', %s)" % (var, base)) if ispy else 'Sys.getenv("%s", %s)' % (var, base)
        if not drest:
            return call
        return ("os.path.join(%s, '%s')" % (call, drest)) if ispy else ('file.path(%s, "%s")' % (call, drest))
    return GETENV_ANY.sub(f, text), n


def _sweep_sh(text):
    n = 0

    def f(m):
        nonlocal n
        var, val = m.group(1), m.group(2)
        mm = PATH_RE_RX.match(val)
        if not mm:
            return m.group(0)
        _p, dvar, dsub, drest = prefix_of(mm.group(0))
        if dvar is None:
            return m.group(0)
        n += 1
        tail = ('/' + drest if drest else '') + val[len(mm.group(0)):]
        base = sh_base(dvar, dsub)
        if var == dvar and not tail:
            return base
        return '${%s:-%s%s}' % (var, base, tail)
    return SH_DEFAULT.sub(f, text), n


def port_code(text, kind):
    """按语言改写：整串字面量与 getenv 默认值 → 环境变量表达式。"""
    n = 0
    q = r"""(['"])"""
    if kind in ('py', 'r'):
        call = r'(?:os\.)?environ\.get' if kind == 'py' else r'Sys\.getenv'
        rx = re.compile(r"""%s\(%s([A-Z0-9_]+)\1,\s*%s((?:%s))\3\)""" % (call, q, q, PATH_RE))

        def f1(m):
            nonlocal n
            var, path = m.group(2), m.group(4)
            e = expr_for(path, kind)
            if e is None:
                return m.group(0)
            n += 1
            return e if var == prefix_of(path)[1] else (
                "os.environ.get('%s', %s)" % (var, e) if kind == 'py' else 'Sys.getenv("%s", %s)' % (var, e))
        text = rx.sub(f1, text)
        guard = r'(?<![A-Za-z0-9_])' if kind == 'py' else r'(?<![A-Za-z0-9_.])'
        rx2 = re.compile(r"""%s%s((?:%s))\1""" % (guard, q, PATH_RE))

        def f2(m):
            nonlocal n
            e = expr_for(m.group(2), kind)
            if e is None:
                return m.group(0)
            n += 1
            return e
        text = rx2.sub(f2, text)
    else:
        rx = re.compile(r"""\$\{([A-Z0-9_]+):-((?:%s))\}""" % PATH_RE)

        def f1(m):
            nonlocal n
            e = expr_for(m.group(2), 'sh')
            if e is None:
                return m.group(0)
            n += 1
            return e if m.group(1) == prefix_of(m.group(2))[1] else '${%s:-%s}' % (m.group(1), e)
        text = rx.sub(f1, text)
        rx3 = re.compile(r"""(?<![A-Za-z0-9_:$'"])((?:%s))""" % PATH_RE)

        def f3(m):
            nonlocal n
            e = expr_for(m.group(1), 'sh')
            if e is None:
                return m.group(0)
            n += 1
            return e
        text = rx3.sub(f3, text)
    text, k = _sweep_getenv(text)
    n += k
    text, k = _sweep_sh(text)
    n += k
    return text, n


def port_home(text, kind):
    """机器家目录 → 可移植写法（代码）或 ~（文档/输出）。"""
    n = 0
    if kind in ('py', 'r'):
        def f(m):
            nonlocal n
            n += 1
            if kind == 'py':
                return "os.path.expanduser('%s~%s')" % (m.group(1), m.group(2))
            return 'path.expand("%s~%s")' % (m.group(1), m.group(2))
        text = HOME_LIT.sub(f, text)

    def g(m):
        tail = m.group(1) or ''
        if kind == 'sh':
            return '$HOME%s' % tail
        if kind == 'py':
            return "os.path.expanduser('~%s')" % tail
        if kind == 'r':
            return 'path.expand("~%s")' % tail
        return '~%s' % tail
    out, k = HOME_ANY.subn(g, text)
    return out, n + k


def sweep_bare(text):
    """代码里最终剩下的非整串字面量路径（docstring/用法行）→ 占位符；引号内的留给人工。"""
    hits = []

    def f(m):
        if m.start() > 0 and text[m.start() - 1] in ('"', "'"):
            return m.group(0)
        hits.append(m.group(0))
        return to_placeholder(m.group(0))
    return PATH_RE_RX.sub(f, text), len(hits)


def port_code_full(src, kind):
    out, n = port_code(src, kind)
    out, k = port_home(out, kind)
    n += k
    if kind in ('py', 'r'):
        out, k = sweep_bare(out)
        n += k
    return out, n


def leftovers(text):
    return PATH_RE_RX.findall(text) + re.findall(r'/home/[A-Za-z0-9_.+-]+', text)


def port_py_like(src):
    out, n = port_code_full(src, 'py')
    if n and 'os.environ.get(' in out and not re.search(r'\bimport os\b', out):
        out = 'import os\n' + out
    return out, n


def port_nb(path, inplace=True):
    """notebook：code cell 按 kernel 语言改写；markdown/outputs/metadata 用占位符。"""
    nb = json.load(open(path, encoding='utf-8'))
    lang = (nb.get('metadata', {}).get('kernelspec', {}) or {}).get('language', 'python').lower()
    kind = 'r' if lang in ('r', 'ir_kernel', 'rds') else 'py'
    total = 0
    left = []

    def join(lines):
        return ''.join(lines)

    def split(text):
        parts = text.split('\n')
        return [p + '\n' for p in parts[:-1]] + ([parts[-1]] if parts[-1] else [])

    def ph(text):
        out = doc_users(text)
        out = PATH_RE_RX.sub(lambda m: to_placeholder(m.group(0)), out)
        return HOME_ANY.sub(lambda m: '~%s' % (m.group(1) or ''), out)

    for cell in nb.get('cells', []):
        src = join(cell.get('source', []))
        if cell.get('cell_type') == 'markdown':
            new, n = ph(src), len(leftovers(src))
        else:
            new, n = port_code_full(src, kind)
            left += leftovers(new)
        if n:
            cell['source'] = split(new)
            total += n
        for od in cell.get('outputs', []) or []:
            for key in ('text', 'traceback'):
                if key in od:
                    buf = join(od[key])
                    if leftovers(buf):
                        od[key] = split(ph(buf))
                        total += len(leftovers(buf))
            if 'data' in od:
                for mime, payload in od['data'].items():
                    buf = join(payload if isinstance(payload, list) else [payload])
                    if leftovers(buf):
                        od['data'][mime] = split(ph(buf)) if isinstance(payload, list) else ph(buf)
                        total += len(leftovers(buf))
    meta = nb.get('metadata', {})
    if leftovers(json.dumps(meta)):
        nb['metadata'] = json.loads(ph(json.dumps(meta, ensure_ascii=False)))
        total += len(leftovers(json.dumps(meta)))
    if inplace and total:
        with open(path, 'w', encoding='utf-8') as fh:
            json.dump(nb, fh, ensure_ascii=False, indent=1)
            fh.write('\n')
    for x in sorted(set(left))[:6]:
        print('   ! 需人工确认（字符串内嵌路径）:', x[:90])
    return total


def doc_users(text):
    return USERS_RX.sub(lambda m: ('<work>' + m.group(0)[m.group(0).index('/online') + 7:]) if '/online' in m.group(0) else '<work>', text)


def port_doc(src):
    out, k = USERS_RX.subn(lambda m: ('<work>' + m.group(0)[m.group(0).index('/online') + 7:]) if '/online' in m.group(0) else '<work>', src)
    out, n = PATH_RE_RX.subn(lambda m: to_placeholder(m.group(0)), out)
    out, j = HOME_ANY.subn(lambda m: '~%s' % (m.group(1) or ''), out)
    return out, n + k + j


def port(path, inplace=True):
    base = os.path.basename(path)
    if base in SKIP:
        print('%-46s 跳过（列在 SKIP，见头注）' % base)
        return 0
    src = open(path, encoding='utf-8', errors='ignore').read()
    if not leftovers(src):
        return 0
    if path.endswith('.ipynb'):
        if inplace:
            n = port_nb(path, True)
            if n:
                print('%-46s 改写 %d 处（notebook）' % (base, n))
            return n
        nb = json.loads(src)
        lang = (nb.get('metadata', {}).get('kernelspec', {}) or {}).get('language', 'python').lower()
        kind = 'r' if lang in ('r', 'ir_kernel', 'rds') else 'py'
        n = 0
        for cell in nb.get('cells', []):
            txt = ''.join(cell.get('source', []))
            n += len(leftovers(txt)) if cell.get('cell_type') == 'markdown' else port_code_full(txt, kind)[1]
            for od in cell.get('outputs', []) or []:
                for key in ('text', 'traceback'):
                    n += len(leftovers(''.join(od.get(key, []) or [])))
                for payload in (od.get('data', {}) or {}).values():
                    n += len(leftovers(''.join(payload) if isinstance(payload, list) else str(payload)))
        print('%-46s 待改写约 %d 处（notebook, dry）' % (base, n))
        return n
    if path.endswith(('.py',)):
        out, n = port_py_like(src)
    elif path.endswith(('.R', '.r')):
        out, n = port_code_full(src, 'r')
    elif path.endswith(('.sh', '.bash')):
        out, n = port_code_full(src, 'sh')
    elif path.endswith(('.md', '.log', '.txt', '.json', '.yaml', '.yml', '.tsv', '.csv')):
        out, n = port_doc(src)
    else:
        return 0
    left = leftovers(out)
    print('%-46s 改写 %d 处%s' % (base, n, ('   ! 残留 %d 处需人工' % len(left)) if left else ''))
    for x in sorted(set(left))[:6]:
        print('       ·', x[:110])
    if n and inplace:
        open(path, 'w', encoding='utf-8').write(out)
    return n


if __name__ == '__main__':
    dry = '--dry' in sys.argv
    tot = 0
    for a in [x for x in sys.argv[1:] if not x.startswith('--')]:
        if os.path.isdir(a):
            for root, dirs, fs in os.walk(a):
                dirs[:] = [d for d in dirs if d not in ('__pycache__', '.git', '.ipynb_checkpoints')]
                for f in sorted(fs):
                    if f.endswith(('.py', '.R', '.r', '.sh', '.ipynb', '.md', '.log', '.txt', '.json', '.yaml', '.yml', '.csv', '.tsv')):
                        tot += port(os.path.join(root, f), inplace=not dry)
        else:
            tot += port(a, inplace=not dry)
    print('合计 %d 处' % tot)
