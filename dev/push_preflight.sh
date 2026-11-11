#!/usr/bin/env bash
# push_preflight.sh — 首次公开推送前的统一闸门
#   bash dev/push_preflight.sh            # 只看结果
#   bash dev/push_preflight.sh --verbose  # 展开每条明细
# FAIL≠0 时不要 push。闸的内容见 dev/preflight_checks.py 头注。
set -uo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")/.." || exit 2
echo "仓库: $(pwd)"
echo "分支: $(git rev-parse --abbrev-ref HEAD)  HEAD: $(git log -1 --format='%h %s')"
echo "工作树: $(git status --porcelain | wc -l) 处改动   tracked 文件: $(git ls-files | wc -l)"
echo ".git 体积: $(du -sh .git | cut -f1)   工作树: $(du -sh --exclude=.git . | cut -f1)"
echo "remote: $(git remote -v | head -1 | awk '{print $2}')   upstream: $(git rev-parse --abbrev-ref --symbolic-full-name @{u} 2>/dev/null || echo '无（首推需 -u origin main）')"
echo
python3 dev/preflight_checks.py "$@"
rc=$?
echo
if [ $rc -ne 0 ]; then
  echo "结论：FAIL≠0 → 先修再推（不要 --force 绕过）。"
else
  echo "结论：硬闸通过。推送前还差：①推送凭据（GH_TOKEN/SSH，沙箱无）②作者对 G10 两项（noreply 身份 + LICENSE 口径）拍板。"
fi
exit $rc
