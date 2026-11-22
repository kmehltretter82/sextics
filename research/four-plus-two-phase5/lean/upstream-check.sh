#!/usr/bin/env bash
set -euo pipefail
phase_dir="$(cd "$(dirname "$0")" && pwd)"
export PATH="$HOME/.elan/bin:$PATH"
export MATHLIB_NO_CACHE_ON_UPDATE=1
cd "$phase_dir"
if [ ! -e sextics/.git ]; then
    git init -q sextics
    git -C sextics remote add origin https://github.com/hodge-computer/sextics.git
fi
git -C sextics fetch --quiet --depth 1 origin 2388542bd16eb4df849515418e1c2ab291e6dc6a
git -C sextics checkout -q --detach FETCH_HEAD
cd sextics/formal/lean
lean --version
python3 - <<'PY'
from concurrent.futures import ThreadPoolExecutor
import json, pathlib, subprocess
def get_pkg(p):
    dest=pathlib.Path(".lake/packages")/p["name"]
    dest.mkdir(parents=True, exist_ok=True)
    def git(*args):
        return subprocess.run(["git","-C",str(dest),*args],check=True)
    if not (dest/".git").exists():
        git("init","-q")
        git("remote","add","origin",p["url"])
    git("fetch","-q","--depth","1","origin",p["rev"])
    git("checkout","-q","--detach","FETCH_HEAD")
with ThreadPoolExecutor(max_workers=4) as pool:
    list(pool.map(get_pkg,json.loads(pathlib.Path("lake-manifest.json").read_text())["packages"]))
PY
mapfile -t mods < <(awk '/^import Mathlib\./ {print $2}' "$phase_dir"/{FourPlusTwoCertificate,RatioTransport}.lean HodgeComputer/Verify.lean HodgeTestbed/Algebra.lean | sort -u)
lake exe cache get "${mods[@]}"
lake --no-cache build HodgeComputer.Verify
export LEAN_PATH="$(lake env printenv LEAN_PATH):$phase_dir"
export LEAN_SYSROOT="$(lean --print-prefix)"
export PATH="$LEAN_SYSROOT/bin:$PATH"
lean_cmd="$LEAN_SYSROOT/bin/lean"
cd "$phase_dir"
for target in HodgeTarget FourPlusTwoCertificate RatioTransport AxiomCheck; do
    "$lean_cmd" -o "$target.olean" "$target.lean"
done
for target in FourPlusTwoCertificate RatioTransport; do
    leanchecker --verbose "$target"
done
