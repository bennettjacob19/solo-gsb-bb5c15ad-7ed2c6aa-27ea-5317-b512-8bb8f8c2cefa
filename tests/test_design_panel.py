"""app.digest 的 /design-panel 纯标准库测试：python3 tests/test_design_panel.py"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.digest import DigestError, run_digest, run_design_panel  # noqa: E402

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


# ---------------------------------------------------------------------------
# 核心样本对（长度均为 40）：
#   A=AAG/off0、B=GTA/off1 单独酶切两样本的长度谱完全相同，
#   但 A+B 联合酶切后的长度谱不同——任何单酶都无法区分，只有联合可以。
# ---------------------------------------------------------------------------
S1 = "AAAAAAGAACGTAAAAAAAAGAAAAAAAAAACGTAAAAAA"
S2 = "AAAAAAAAAAAAAAAAGAAACGTAAAACGTAAAAAAAAGA"
check("case lengths", len(S1) == 40 and len(S2) == 40)

A = {"name": "A", "recognition": "AAG", "cut_offset": 0}
B = {"name": "B", "recognition": "GTA", "cut_offset": 1}
C = {"name": "C", "recognition": "GGG", "cut_offset": 0}  # 两个样本上均不切割
SAMPLES = [{"name": "s1", "sequence": S1}, {"name": "s2", "sequence": S2}]


def joint_lengths(seq):
    res = run_digest(
        {"sequence": seq, "topology": "linear", "enzymes": [A, B]}
    )
    return sorted(f["length"] for f in res["fragments"])


# --- 先确认前提：两个单酶各自的谱相同，联合谱不同 ---
for enzyme, expected in ((A, [4, 14, 22]), (B, [7, 11, 22])):
    p1 = run_digest({"sequence": S1, "topology": "linear", "enzymes": [enzyme]})
    p2 = run_digest({"sequence": S2, "topology": "linear", "enzymes": [enzyme]})
    l1 = sorted(f["length"] for f in p1["fragments"])
    l2 = sorted(f["length"] for f in p2["fragments"])
    check(f"single {enzyme['name']} profiles equal and fixed",
          l1 == l2 == expected, str((l1, l2)))

check("combined profiles differ", joint_lengths(S1) != joint_lengths(S2),
      str((joint_lengths(S1), joint_lengths(S2))))
check("combined S1 lengths", joint_lengths(S1) == [4, 7, 7, 7, 15],
      str(joint_lengths(S1)))
check("combined S2 lengths", joint_lengths(S2) == [4, 7, 7, 8, 14],
      str(joint_lengths(S2)))

# --- 单酶无方案，联合酶组合 [A, B] 有方案（不能只合并单酶鉴别结果） ---
rA = run_design_panel({"samples": SAMPLES, "topology": "linear",
                       "enzymes": [A], "tolerance": 0})
rB = run_design_panel({"samples": SAMPLES, "topology": "linear",
                       "enzymes": [B], "tolerance": 0})
check("single A no_solution", rA == {"status": "no_solution"}, str(rA))
check("single B no_solution", rB == {"status": "no_solution"}, str(rB))

rAB = run_design_panel({"samples": SAMPLES, "topology": "linear",
                        "enzymes": [A, B], "tolerance": 0})
check("AB status", rAB["status"] == "solution", str(rAB))
check("AB indices", rAB["enzyme_indices"] == [0, 1], str(rAB))
check("AB enzyme_count", rAB["enzyme_count"] == 2)
check("AB profile count", len(rAB["profiles"]) == 2)
check("AB profile s1", rAB["profiles"][0] == {
    "name": "s1", "length": 40, "fragment_count": 5,
    "fragment_lengths": [4, 7, 7, 7, 15]}, str(rAB["profiles"][0]))
check("AB profile s2", rAB["profiles"][1] == {
    "name": "s2", "length": 40, "fragment_count": 5,
    "fragment_lengths": [4, 7, 7, 8, 14]}, str(rAB["profiles"][1]))

# 返回谱与 /digest 对“所选酶组合联合酶切”的片段长度逐样本一致
for i, seq in enumerate((S1, S2)):
    check(f"profile{i} matches /digest joint fragments",
          rAB["profiles"][i]["fragment_lengths"] == joint_lengths(seq))

# 全部片段长度升序排列且保留重复长度
lens = rAB["profiles"][0]["fragment_lengths"]
check("sorted with duplicates kept", lens == sorted(lens) and lens.count(7) == 3,
      str(lens))

# --- 酶数最少：存在单酶区分器时，不返回两酶组合；下标字典序最小 ---
AC = {"name": "AC", "recognition": "AC", "cut_offset": 0}
AG = {"name": "AG", "recognition": "AG", "cut_offset": 0}
r1 = run_design_panel({"samples": SAMPLES, "topology": "linear",
                       "enzymes": [C, AG, AC, A, B], "tolerance": 0})
check("min enzyme count picks singleton index 1",
      r1["status"] == "solution" and r1["enzyme_indices"] == [1]
      and r1["enzyme_count"] == 1, str(r1))
r2 = run_design_panel({"samples": SAMPLES, "topology": "linear",
                       "enzymes": [AC, AG, C], "tolerance": 0})
check("lexicographic AC before AG", r2["enzyme_indices"] == [0], str(r2))
r3 = run_design_panel({"samples": SAMPLES, "topology": "linear",
                       "enzymes": [AG, AC, C], "tolerance": 0})
check("lexicographic AG before AC", r3["enzyme_indices"] == [0], str(r3))

# 重复酶定义合法；B 的重复体与 B 联合等价于单酶 B（无法区分），
# 此时 [A, B, Bdup, C] 中唯一可行二元组是 {0,1}
Bdup = {"name": "Bdup", "recognition": "GTA", "cut_offset": 1}
r4 = run_design_panel({"samples": SAMPLES, "topology": "linear",
                       "enzymes": [A, B, Bdup, C], "tolerance": 0})
check("duplicate enzyme defs, only {0,1} works",
      r4["status"] == "solution" and r4["enzyme_indices"] == [0, 1], str(r4))

# 不切酶 C：两样本谱均为 [40]，不能区分
rC = run_design_panel({"samples": SAMPLES, "topology": "linear",
                       "enzymes": [C], "tolerance": 0})
check("non-cutter no_solution", rC == {"status": "no_solution"}, str(rC))

# --- 三个样本：任意样本对都必须可区分 ---
S3 = "A" * 40
r3s = run_design_panel({
    "samples": [{"name": "s1", "sequence": S1},
                {"name": "s2", "sequence": S2},
                {"name": "s3", "sequence": S3}],
    "topology": "linear", "enzymes": [A, B, C], "tolerance": 0,
})
check("3-sample solution", r3s["status"] == "solution"
      and r3s["enzyme_indices"] == [0, 1], str(r3s))
check("3-sample s3 uncut profile", r3s["profiles"][2] == {
    "name": "s3", "length": 40, "fragment_count": 1,
    "fragment_lengths": [40]}, str(r3s["profiles"][2]))
# s2 与 s3 必须被所选组合区分（s2 谱 5 段 vs s3 谱 1 段）
check("3-sample s2 vs s3 distinguished",
      r3s["profiles"][1]["fragment_count"]
      != r3s["profiles"][2]["fragment_count"])

# --- 容差：片段数相同、逐位长度差均 <= 容差才算不可区分 ---
# AAG/off0：t1 切在 3 -> [3,7]；t2 切在 4 -> [4,6]；逐位差为 1
t1, t2 = "AAAAAGTTTT", "AAAAAAGTTT"
check("tolerance lens", len(t1) == len(t2) == 10)
r_t0 = run_design_panel({"samples": [{"name": "p", "sequence": t1},
                                     {"name": "q", "sequence": t2}],
                         "topology": "linear", "enzymes": [A], "tolerance": 0})
r_t1 = run_design_panel({"samples": [{"name": "p", "sequence": t1},
                                     {"name": "q", "sequence": t2}],
                         "topology": "linear", "enzymes": [A], "tolerance": 1})
check("tolerance 0 distinguishes", r_t0["status"] == "solution", str(r_t0))
check("tolerance 1 indistinguishable", r_t1["status"] == "no_solution", str(r_t1))

# 片段数不同时，无论容差多大都不可视为相同谱（计数不同 => 可区分）
r_t9 = run_design_panel({"samples": [{"name": "p", "sequence": t1},
                                     {"name": "q", "sequence": "AAAAAAAAAA"}],
                         "topology": "linear", "enzymes": [A], "tolerance": 999})
check("different fragment count always distinguishable",
      r_t9["status"] == "solution", str(r_t9))

# --- no_solution：相同序列（即便名称不同）永远不可区分 ---
r_same = run_design_panel({
    "samples": [{"name": "x", "sequence": "ACGT"},
                {"name": "y", "sequence": "acgt"}],
    "topology": "linear", "enzymes": [A, B, C], "tolerance": 5,
})
check("identical sequences no_solution", r_same == {"status": "no_solution"},
      str(r_same))

# --- 环状拓扑：与 /digest 环状联合酶切一致，可区分时正常返回 ---
# AAG/off0 在环上：无位点 -> [8]；两个位点切点 {0,3} -> 升序长度 [3,5]
c1, c2 = "AAAAAAAA", "AAGAAGAA"
rc = run_design_panel({"samples": [{"name": "c1", "sequence": c1},
                                   {"name": "c2", "sequence": c2}],
                       "topology": "circular", "enzymes": [A], "tolerance": 0})
check("circular solution", rc["status"] == "solution"
      and rc["enzyme_indices"] == [0], str(rc))
for i, s in enumerate((c1, c2)):
    d = run_digest({"sequence": s, "topology": "circular", "enzymes": [A]})
    expect = sorted(f["length"] for f in d["fragments"])
    check(f"circular profile {i} matches digest",
          rc["profiles"][i]["fragment_lengths"] == expect,
          str(rc["profiles"][i]["fragment_lengths"]))

# --- 非法输入：全部返回 DigestError（HTTP 层映射为 400） ---
expect_error("payload not dict", ["not", "dict"], "JSON 对象")
expect_error("missing samples", {"topology": "linear", "enzymes": [A],
                                 "tolerance": 0}, "samples")
expect_error("samples not list",
             {"samples": "nope", "topology": "linear", "enzymes": [A],
              "tolerance": 0}, "samples")
expect_error("too few samples",
             {"samples": SAMPLES[:1], "topology": "linear", "enzymes": [A],
              "tolerance": 0}, "数量")
expect_error("duplicate sample name",
             {"samples": [{"name": "s1", "sequence": S1},
                          {"name": "s1", "sequence": S2}],
              "topology": "linear", "enzymes": [A], "tolerance": 0}, "重复")
expect_error("blank sample name",
             {"samples": [{"name": "  ", "sequence": S1},
                          {"name": "s2", "sequence": S2}],
              "topology": "linear", "enzymes": [A], "tolerance": 0}, "name")
expect_error("sample name not string",
             {"samples": [{"name": 1, "sequence": S1},
                          {"name": "s2", "sequence": S2}],
              "topology": "linear", "enzymes": [A], "tolerance": 0}, "name")
expect_error("sample not object",
             {"samples": ["seq", SAMPLES[1]], "topology": "linear",
              "enzymes": [A], "tolerance": 0}, "对象")
expect_error("illegal sample base",
             {"samples": [{"name": "s1", "sequence": "AXGT"},
                          {"name": "s2", "sequence": S2}],
              "topology": "linear", "enzymes": [A], "tolerance": 0}, "非法碱基")
expect_error("sequence too long",
             {"samples": [{"name": "s1", "sequence": "A" * 2001},
                          {"name": "s2", "sequence": S2}],
              "topology": "linear", "enzymes": [A], "tolerance": 0}, "上限")
expect_error("bad topology",
             {"samples": SAMPLES, "topology": "plasmid", "enzymes": [A],
              "tolerance": 0}, "拓扑")
expect_error("no enzymes",
             {"samples": SAMPLES, "topology": "linear", "enzymes": [],
              "tolerance": 0}, "候选酶数量")
expect_error("too many enzymes",
             {"samples": SAMPLES, "topology": "linear",
              "enzymes": [A] * 11, "tolerance": 0}, "候选酶数量")
expect_error("bad enzyme offset",
             {"samples": SAMPLES, "topology": "linear",
              "enzymes": [{"name": "E", "recognition": "AA",
                           "cut_offset": 9}], "tolerance": 0}, "越界")
expect_error("tolerance negative",
             {"samples": SAMPLES, "topology": "linear", "enzymes": [A],
              "tolerance": -1}, "tolerance")
expect_error("tolerance float",
             {"samples": SAMPLES, "topology": "linear", "enzymes": [A],
              "tolerance": 1.0}, "tolerance")
expect_error("tolerance bool",
             {"samples": SAMPLES, "topology": "linear", "enzymes": [A],
              "tolerance": False}, "tolerance")
expect_error("missing tolerance",
             {"samples": SAMPLES, "topology": "linear", "enzymes": [A]},
             "tolerance")

# --- 边界：2000 碱基接受；8 个样本接受 ---
long_seq = "ACGT" * 500  # 正好 2000
check("2000 length accepted", len(long_seq) == 2000)
r_long = run_design_panel({
    "samples": [{"name": "lo", "sequence": long_seq},
                {"name": "hi", "sequence": "TTTT" + long_seq[4:]}],
    "topology": "linear",
    "enzymes": [{"name": "E", "recognition": "AAAA", "cut_offset": 0}],
    "tolerance": 0,
})
check("2000-length panel runs", r_long["status"] in ("solution", "no_solution"),
      str(r_long)[:200])

bases8 = ["ACGTACGT", "TTTTAAAA", "GCGCGCAT", "AATTCCGG",
          "GGGGCCCC", "TATATATA", "CCTTGGAA", "AAACCCGT"]
r8 = run_design_panel({
    "samples": [{"name": f"s{i}", "sequence": bases8[i]} for i in range(8)],
    "topology": "linear",
    "enzymes": [{"name": "E", "recognition": "AC", "cut_offset": 0},
                {"name": "F", "recognition": "GT", "cut_offset": 0},
                {"name": "G", "recognition": "GC", "cut_offset": 0}],
    "tolerance": 0,
})
check("8-sample panel well-formed",
      r8["status"] == "no_solution"
      or (r8["status"] == "solution" and len(r8["profiles"]) == 8),
      str(r8)[:200])

print(f"\n{PASS} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
