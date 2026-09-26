"""app.design 的纯标准库测试：python3 tests/test_design.py

重点验证：
- /design-panel 复用 /digest 的联合酶切规则（谱与 run_digest 多酶结果一致）；
- 不能只合并单酶鉴别结果：单酶谱两两不可区分、唯有联合酶切才可区分；
- 最少酶数 + 下标字典序；容差比对；环状多样本；no_solution；各类 400。
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.digest import DigestError, run_digest  # noqa: E402
from app.design import run_design_panel  # noqa: E402

PASS = 0
FAIL = 0


def check(name, cond, extra=""):
    global PASS, FAIL
    if cond:
        PASS += 1
        print(f"  PASS {name}")
    else:
        FAIL += 1
        print(f"  FAIL {name} {extra}")


def expect_error(name, payload, snippet):
    try:
        run_design_panel(payload)
    except DigestError as e:
        check(name, snippet in e.reason, f"reason={e.reason!r}")
        return
    check(name, False, "未抛出 DigestError")


def digest_lengths(seq, topology, enzymes):
    """用现有 /digest 逻辑计算联合酶切片段长度谱（升序、保留重复）。"""
    res = run_digest(
        {"sequence": seq, "topology": topology, "enzymes": enzymes}
    )
    return sorted(f["length"] for f in res["fragments"])


def panel(samples, topology, enzymes, tolerance=0):
    return run_design_panel(
        {
            "samples": [
                {"name": name, "sequence": seq} for name, seq in samples
            ],
            "topology": topology,
            "enzymes": enzymes,
            "length_tolerance": tolerance,
        }
    )


ECO = {"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1}
# 非回文酶：正/反向识别序列不同，每个位点只产生一条链的切点
#   GGTT (rc=AACC), offset=1：正向切点 s+1，反向切点 s+3
#   GAGT (rc=ACTC), offset=2：正向/反向切点均为 s+2
E1 = {"name": "E1", "recognition": "GGTT", "cut_offset": 1}
E2 = {"name": "E2", "recognition": "GAGT", "cut_offset": 2}


def build(n, motifs):
    """长度 n 的 T 填充序列上，在指定起点放置 (motif, start)。"""
    bases = ["T"] * n
    for motif, start in motifs:
        for k, ch in enumerate(motif):
            bases[start + k] = ch
    return "".join(bases)


# ---------------------------------------------------------------------------
# 反作弊用例（线性 n=24）：两样本单酶谱完全相同，只有联合酶切谱不同
#   s1: GGTT@0 -> E1 切点 1；GAGT@6 -> E2 切点 8；ACTC@12 -> E2 切点 14
#   s2: GAGT@6 -> 8；ACTC@12 -> 14；AACC@20 -> E1 反向切点 23
# 单酶 E1：s1 [1,23] == s2 [1,23]（切点 1 与 23 的长度多重集相同）
# 单酶 E2：均为切点 {8,14} -> [6,8,10]
# 联合：s1 {1,8,14}  -> [1,6,7,10]
#       s2 {8,14,23} -> [1,6,8,9]
# ---------------------------------------------------------------------------
S1 = build(24, [("GGTT", 0), ("GAGT", 6), ("ACTC", 12)])
S2 = build(24, [("GAGT", 6), ("ACTC", 12), ("AACC", 20)])
check("anti-cheat lengths", len(S1) == 24 and len(S2) == 24,
      f"{len(S1)},{len(S2)}")
check("anti-cheat seq s1", S1 == "GGTTTTGAGTTTACTCTTTTTTTT", S1)
check("anti-cheat seq s2", S2 == "TTTTTTGAGTTTACTCTTTTAACC", S2)

e1_s1 = digest_lengths(S1, "linear", [E1])
e1_s2 = digest_lengths(S2, "linear", [E1])
e2_s1 = digest_lengths(S1, "linear", [E2])
e2_s2 = digest_lengths(S2, "linear", [E2])
check("anti-cheat E1 singles equal", e1_s1 == e1_s2 == [1, 23],
      f"{e1_s1} vs {e1_s2}")
check("anti-cheat E2 singles equal", e2_s1 == e2_s2 == [6, 8, 10],
      f"{e2_s1} vs {e2_s2}")
joint_s1 = digest_lengths(S1, "linear", [E1, E2])
joint_s2 = digest_lengths(S2, "linear", [E1, E2])
check("anti-cheat joint differ", joint_s1 == [1, 6, 7, 10]
      and joint_s2 == [1, 6, 8, 9]
      and joint_s1 != joint_s2, f"{joint_s1} vs {joint_s2}")

res = panel([("s1", S1), ("s2", S2)], "linear", [E1, E2])
check("anti-cheat result found", res["status"] == "found"
      and res["enzyme_indices"] == [0, 1], str(res))
check("anti-cheat profiles reuse digest rules",
      res["profiles"][0]["fragment_lengths"] == joint_s1
      and res["profiles"][1]["fragment_lengths"] == joint_s2
      and [p["name"] for p in res["profiles"]] == ["s1", "s2"], str(res))

# 单酶组合（更短）确实都不行：直接验证没有 size-1 方案
res0 = panel([("s1", S1), ("s2", S2)], "linear", [E1])
res1 = panel([("s1", S1), ("s2", S2)], "linear", [E2])
check("anti-cheat E1 alone no_solution", res0["status"] == "no_solution", str(res0))
check("anti-cheat E2 alone no_solution", res1["status"] == "no_solution", str(res1))

# 三样本扩展：任何单酶都无法区分全部样本对，唯有双酶组合可以
# s3: GGTT@0 -> E1 切点 1；GAGT@5 -> E2 切点 7；ACTC@16 -> E2 切点 18
S3 = build(24, [("GGTT", 0), ("GAGT", 5), ("ACTC", 16)])
check("s3 length", len(S3) == 24, str(len(S3)))
e1_s3 = digest_lengths(S3, "linear", [E1])
e2_s3 = digest_lengths(S3, "linear", [E2])
check("s3 E1 single same", e1_s3 == [1, 23], str(e1_s3))
# E2 切点 {7,18} -> [7,11,6]
check("s3 E2 single differs", e2_s3 == [6, 7, 11], str(e2_s3))
joint_s3 = digest_lengths(S3, "linear", [E1, E2])
# 切点 {1,7,18} -> [1,6,11,6]
check("s3 joint exact", joint_s3 == [1, 6, 6, 11]
      and joint_s3 not in (joint_s1, joint_s2), str(joint_s3))
res = panel([("a", S1), ("b", S2), ("c", S3)], "linear", [E1, E2])
check("three-sample needs pair", res["status"] == "found"
      and res["enzyme_indices"] == [0, 1], str(res))
check("three-sample profiles",
      [p["fragment_lengths"] for p in res["profiles"]]
      == [joint_s1, joint_s2, joint_s3], str(res))

# ---------------------------------------------------------------------------
# 基础：单酶即可区分（位点数量不同）；返回下标 0
# ---------------------------------------------------------------------------
A = "AAGAATTCTTTGAATTCGG"     # 两个 EcoRI 位点，长度 19
B = "AAGAATTCTTTCCCCCCCCCCC"  # 一个 EcoRI 位点，长度 22
res = panel([("a", A), ("b", B)], "linear", [ECO])
expect_a = digest_lengths(A, "linear", [ECO])
expect_b = digest_lengths(B, "linear", [ECO])
check("basic found", res["status"] == "found"
      and res["enzyme_indices"] == [0], str(res))
check("basic profiles",
      res["profiles"][0]["fragment_lengths"] == expect_a
      and res["profiles"][1]["fragment_lengths"] == expect_b
      and expect_a == [3, 3, 4, 4, 5] and expect_b == [3, 4, 15], str(res))

# 无切割酶也能区分（样本本身长度不同）
res = panel([("a", "AAAA"), ("b", "AAA")], "linear",
            [{"name": "X", "recognition": "GGGG", "cut_offset": 0}])
check("no-cut enzyme distinguishes by length", res["status"] == "found"
      and res["enzyme_indices"] == [0]
      and [p["fragment_lengths"] for p in res["profiles"]] == [[4], [3]], str(res))

# 片段数不同：容差再大也无法混淆（无切点单片段 vs 两个片段）
res = panel([("a", "AAGAATTCTTTGAATTCGG"), ("b", "AAA")], "linear",
            [ECO], 100)
check("different fragment count always distinct",
      res["status"] == "found" and res["enzyme_indices"] == [0], str(res))

# ---------------------------------------------------------------------------
# 容差：谱 [4] 与 [3]，容差 0 可区分，容差 1 不可区分 -> no_solution
# ---------------------------------------------------------------------------
tol_enzymes = [{"name": "X", "recognition": "GGGG", "cut_offset": 0}]
check("tolerance 0 found",
      panel([("a", "AAAA"), ("b", "AAA")], "linear", tol_enzymes, 0)["status"]
      == "found")
check("tolerance 1 no_solution",
      panel([("a", "AAAA"), ("b", "AAA")], "linear", tol_enzymes, 1)["status"]
      == "no_solution")
# [4,4] vs [4,5]：容差 1 下对应差均为 0/1 -> 不可区分；容差 0 可区分
res = panel([("a", "AAAAAA"), ("b", "AAAAA")], "linear", tol_enzymes, 1)
check("tolerance merges close profiles", res["status"] == "no_solution", str(res))
res = panel([("a", "AAAAAA"), ("b", "AAAAA")], "linear", tol_enzymes, 0)
check("tolerance 0 keeps apart", res["status"] == "found", str(res))

# ---------------------------------------------------------------------------
# 环状 + 三样本：GATC 为回文 (rc=自身)，offset=1 -> 每位点正切 s+1、反切 s+3
# ---------------------------------------------------------------------------
GATC = {"name": "E", "recognition": "GATC", "cut_offset": 1}
# GATCGATC (n=8)：位点 0/4 -> 切点 {1,3,5,7} -> [2,2,2,2]
# GATCGATCG (n=9)：位点 0/4 -> 切点 {1,3,5,7}，含跨首尾片段 3 -> [2,2,2,3]
cs1, cs2, cs3 = "GATCGATC", "GATCGATCG", "AAAAAAAA"
res = panel([("c1", cs1), ("c2", cs2), ("c3", cs3)], "circular", [GATC])
check("circular three samples found", res["status"] == "found"
      and res["enzyme_indices"] == [0], str(res))
check("circular profiles",
      [p["fragment_lengths"] for p in res["profiles"]]
      == [[2, 2, 2, 2], [2, 2, 2, 3], [8]], str(res["profiles"]))
res = panel([("c1", cs1), ("c3", cs3)], "circular", [GATC])
check("circular pair", res["status"] == "found"
      and res["profiles"][0]["fragment_lengths"] == [2, 2, 2, 2]
      and res["profiles"][1]["fragment_lengths"] == [8], str(res))

# ---------------------------------------------------------------------------
# 字典序：两个单酶都可行时选下标 0；双酶方案中 [0,1] 优先于 [0,2]
# 两个非回文酶，每个只切一处，避免回文酶一条链两切点干扰判断
# ---------------------------------------------------------------------------
# GGTT@0(切1) 命中 a；CCAA@0(切2) 命中 b；两者都能区分 -> 取字典序下标 0
G_ENZ = {"name": "g", "recognition": "GGTT", "cut_offset": 1}  # rc=AACC 非回文
C_ENZ = {"name": "c", "recognition": "CCAA", "cut_offset": 2}  # rc=TTGG 非回文
res = panel([("a", "GGTTAAAA"), ("b", "CCAAAAAA")], "linear", [G_ENZ, C_ENZ])
check("lex single index 0", res["status"] == "found"
      and res["enzyme_indices"] == [0]
      and res["profiles"][0]["fragment_lengths"] == [1, 7]
      and res["profiles"][1]["fragment_lengths"] == [8], str(res))

# [E1, E2, E2重复]：[0,1] 与 [0,2] 均可行，选字典序更小的 [0,1]
res = panel([("s1", S1), ("s2", S2)], "linear", [E1, E2, dict(E2)])
check("lex pair tie", res["status"] == "found"
      and res["enzyme_indices"] == [0, 1], str(res))
# [不切酶, E1, E2]：唯一可行双酶为 [1,2]
dummy = {"name": "none", "recognition": "GAGAGA", "cut_offset": 1}
res = panel([("s1", S1), ("s2", S2)], "linear", [dummy, E1, E2])
check("pair offset indices", res["status"] == "found"
      and res["enzyme_indices"] == [1, 2]
      and res["profiles"][0]["fragment_lengths"] == joint_s1, str(res))

# 候选酶允许重复定义，下标按输入位置计（两个相同的 GGTT 酶都可行 -> 下标 0）
res = panel([("a", "GGTTAAAA"), ("b", "AAAAAAAA")], "linear",
            [dict(G_ENZ, name="g1"), dict(G_ENZ, name="g2")])
check("duplicate enzymes accepted", res["status"] == "found"
      and res["enzyme_indices"] == [0]
      and res["profiles"][0]["fragment_lengths"] == [1, 7]
      and res["profiles"][1]["fragment_lengths"] == [8], str(res))

# ---------------------------------------------------------------------------
# no_solution：名称唯一但序列完全相同的两样本，任何组合都无法区分
# ---------------------------------------------------------------------------
res = panel([("x", "ACGTACGT"), ("y", "ACGTACGT")], "linear",
            [ECO, E1, E2, dummy])
check("identical sequences no_solution", res["status"] == "no_solution"
      and "enzyme_indices" not in res and "profiles" not in res, str(res))

# ---------------------------------------------------------------------------
# 400 校验
# ---------------------------------------------------------------------------
BASE = {
    "samples": [{"name": "a", "sequence": "ACGT"},
                {"name": "b", "sequence": "ACGT"}],
    "topology": "linear",
    "enzymes": [{"name": "E", "recognition": "AC", "cut_offset": 1}],
    "length_tolerance": 0,
}


def panel_raw(**overrides):
    payload = dict(BASE)
    payload.update(overrides)
    return run_design_panel(payload)


def expect_raw(name, snippet, **overrides):
    try:
        panel_raw(**overrides)
    except DigestError as e:
        check(name, snippet in e.reason, f"reason={e.reason!r}")
        return
    check(name, False, "未抛出 DigestError")


expect_raw("missing samples", "samples", samples=None)
expect_raw("samples not list", "数组", samples={"a": 1})
expect_raw("too few samples", "2..8",
           samples=[{"name": "a", "sequence": "ACGT"}])
expect_raw("too many samples", "2..8",
           samples=[{"name": f"s{i}", "sequence": "ACGT"} for i in range(9)])
expect_raw("sample not object", "对象",
           samples=["ACGT", "ACGT"])
expect_raw("duplicate sample names", "重复",
           samples=[{"name": "a", "sequence": "ACGT"},
                    {"name": " a ", "sequence": "TTTT"}])
expect_raw("empty sample name", "名称不能为空",
           samples=[{"name": "  ", "sequence": "ACGT"},
                    {"name": "b", "sequence": "TTTT"}])
expect_raw("sample name not string", "name 必须是字符串",
           samples=[{"name": 1, "sequence": "ACGT"},
                    {"name": "b", "sequence": "TTTT"}])
expect_raw("sample illegal base", "非法碱基",
           samples=[{"name": "a", "sequence": "ACGN"},
                    {"name": "b", "sequence": "TTTT"}])
expect_raw("sample empty seq", "不能为空",
           samples=[{"name": "a", "sequence": ""},
                    {"name": "b", "sequence": "TTTT"}])
expect_raw("sample too long", "2000",
           samples=[{"name": "a", "sequence": "A" * 2001},
                    {"name": "b", "sequence": "TTTT"}])
expect_raw("bad topology", "拓扑", topology="plasmid")
expect_raw("enzymes empty", "1..10", enzymes=[])
expect_raw("enzymes too many", "1..10",
           enzymes=[{"name": f"e{i}", "recognition": "AC", "cut_offset": 1}
                    for i in range(11)])
expect_raw("enzyme bad offset", "越界",
           enzymes=[{"name": "E", "recognition": "AC", "cut_offset": 5}])
expect_raw("enzyme bad recognition", "非法碱基",
           enzymes=[{"name": "E", "recognition": "AX", "cut_offset": 1}])
expect_raw("tolerance negative", "非负整数", length_tolerance=-1)
expect_raw("tolerance float", "非负整数", length_tolerance=1.5)
expect_raw("tolerance bool", "非负整数", length_tolerance=True)
expect_raw("tolerance string", "非负整数", length_tolerance="0")
# 缺少 length_tolerance 字段（须从 payload 中真正删除该键）
missing_tol = {k: v for k, v in BASE.items() if k != "length_tolerance"}
try:
    run_design_panel(missing_tol)
    check("tolerance missing", False, "未抛出 DigestError")
except DigestError as e:
    check("tolerance missing", "length_tolerance" in e.reason, f"reason={e.reason!r}")
expect_error("payload not dict", ["not", "dict"], "JSON 对象")

# 边界：长度恰为 2000 合法（酶切谱不同即可行）
edge = {
    "samples": [{"name": "a", "sequence": "A" * 2000},
                {"name": "b", "sequence": "A" * 1999 + "G"}],
    "topology": "linear",
    "enzymes": [{"name": "E", "recognition": "AA", "cut_offset": 0}],
    "length_tolerance": 0,
}
# b 中 AAAAA...A(1999)G 含相邻 AA 切点数与 a 不同
res = run_design_panel(edge)
check("length 2000 accepted", res["status"] == "found", str(res))

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
