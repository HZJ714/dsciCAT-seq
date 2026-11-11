#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""memwatch.py — 长时间作业的资源看护：每 20 s 记一行 cgroup/内存/CPU 水位。

为什么需要：本项目多次被 OOM killer 静默杀掉（旧环境 32 GB 限额），而且 `free -g` 报的是宿主的
503 GB，不是配额；只有 cgroup 的 memory.current/peak 才是真相。CPU 也一样：主机 load 常 >100，
配额只有 32 核，所以并行度必须按 cpu.max 算。

用法（后台）：nohup python3 memwatch.py --log <file> --minutes 240 >/dev/null 2>&1 &
阈值：memory.current 超过 LIMIT 的 80% → 在日志里打 ! 行（先只告警，不自作主张杀作业）。
"""
import os, sys, time, subprocess

LOG = os.path.join(os.environ.get('DSCI_FIG_ROOT', os.path.join(os.environ.get('DSCI_WORK_ROOT', ''), 'figures')), 'routeA_out/GENELEVEL_v1/memwatch.log')
MINS = 240
args = sys.argv[1:]
if "--log" in args:
    LOG = args[args.index("--log") + 1]
if "--minutes" in args:
    MINS = float(args[args.index("--minutes") + 1])


def rd(p):
    try:
        return open(p).read().strip()
    except Exception:
        return "n/a"


def human(x):
    try:
        x = int(x)
    except Exception:
        return str(x)
    return "%.1fGB" % (x / 1024**3) if x > 1024**2 else "%.1fMB" % (x / 1024**2)


def top_rss(n=3):
    try:
        o = subprocess.run(["ps", "-eo", "rss,comm", "--sort=-rss"], capture_output=True, text=True, timeout=10).stdout
        return " | " + " ".join("%s=%s" % (l.split()[1], human(l.split()[0])) for l in o.splitlines()[1:n + 1])
    except Exception:
        return ""


os.makedirs(os.path.dirname(LOG), exist_ok=True)
t_end = time.time() + MINS * 60
with open(LOG, "a") as f:
    f.write("# memwatch 开始 %s | memory.max=%s cpu.max=%s\n" % (time.strftime("%F %T"),
          human(rd("/sys/fs/cgroup/memory.max")), rd("/sys/fs/cgroup/cpu.max")))
    while time.time() < t_end:
        cur = rd("/sys/fs/cgroup/memory.current"); peak = rd("/sys/fs/cgroup/memory.peak")
        mx = rd("/sys/fs/cgroup/memory.max")
        try:
            frac = int(cur) / int(mx)
        except Exception:
            frac = 0.0
        if cur == "n/a":   # cgroup v1 路径名不同，退到 v1 名
            cur = rd("/sys/fs/cgroup/memory/memory.current"); peak = rd("/sys/fs/cgroup/memory/memory.peak")
            mx = rd("/sys/fs/cgroup/memory/memory.limit_in_bytes")
            try:
                frac = int(cur) / int(mx)
            except Exception:
                frac = 0.0
        load = open("/proc/loadavg").read().split()[:3]
        nj = subprocess.run(["pgrep", "-fc", "compute-score"], capture_output=True, text=True).stdout.strip() or "0"
        line = "%s cur=%s peak=%s (%.0f%% of %s) load=%s/%s/%s jobs=%s%s" % (
            time.strftime("%T"), human(cur), human(peak), 100 * frac, human(mx),
            load[0], load[1], load[2], nj, top_rss())
        if frac > 0.80:
            line += "  !! 内存超 80%% 配额，建议降并行"
        try:
            lavg = float(load[0])
        except Exception:
            lavg = 0.0
        if lavg > 130:
            line += "  !! load 过高"
        f.write(line + "\n"); f.flush()
        time.sleep(20)
    f.write("# memwatch 结束 %s\n" % time.strftime("%F %T"))
