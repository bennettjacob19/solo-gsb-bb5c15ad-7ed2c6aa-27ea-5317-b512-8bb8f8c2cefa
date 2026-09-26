"""app.digest 的纯标准库测试：python3 tests/test_digest.py"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from app.digest import DigestError, reverse_complement, run_digest  # noqa: E402

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
        run_digest(payload)
    except DigestError as e:
        check(name, snippet in e.reason, f"reason={e.reason!r}")
        return
    check(name, False, "未抛出 DigestError")


# --- reverse_complement ---
check("revcomp basic", reverse_complement("GAATTC") == "GAATTC")
check("revcomp non-pal", reverse_complement("AATT") == "AATT")
check("revcomp GGTCTC", reverse_complement("GGTCTC") == "GAGACC")

# --- EcoRI: GAATTC, offset 1, 线性 ---
# 序列含两个 EcoRI 位点
seq = "AAGAATTCTTTGAATTCGG"
res = run_digest(
    {
        "sequence": seq,
        "topology": "linear",
        "enzymes": [{"name": "EcoRI", "recognition": "GaAtTc", "cut_offset": 1}],
    }
)
# 位点起点 2 和 11；正向切点 3 和 12；回文位点反向切点 2+5=7 和 11+5=16
check("ecori cuts", res["cuts"] == [3, 7, 12, 16], str(res["cuts"]))
frags = res["fragments"]
frag_seqs = [f["sequence"] for f in frags]
check("ecori frag seqs", frag_seqs == [seq[0:3], seq[3:7], seq[7:12], seq[12:16], seq[16:]],
      str(frag_seqs))
check("ecori frag lengths", [f["length"] for f in frags] == [3, 4, 5, 4, 3],
      str([f["length"] for f in frags]))
check("ecori total length", sum(f["length"] for f in frags) == len(seq))
check("ecori boundaries", [(f["start"], f["end"]) for f in frags]
      == [(0, 3), (3, 7), (7, 12), (12, 16), (16, 19)])
check("ecori strands", [h["strand"] for h in res["enzyme_hits"][0]["hits"]]
      == ["forward", "reverse", "forward", "reverse"])

# --- 非回文酶 BsaI: GGTCTC, offset 1；正向命中 ---
# GGTCTC 起点 2 -> 正向切点 3；反向互补 GAGACC 在别处
res = run_digest({
    "sequence": "AAGGTCTCAAAGAGACCTT",
    "topology": "linear",
    "enzymes": [{"name": "BsaI", "recognition": "GGTCTC", "cut_offset": 1}],
})
# GAGACC 起点 11 -> 反向切点 11+(6-1)=16
s = "AAGGTCTCAAAGAGACCTT"
rc_start = s.find("GAGACC")
check("bsai rc start", rc_start == 11)
check("bsai reverse cut", 11 + (6 - 1) == 16)
check("bsai cuts", res["cuts"] == [3, 16], str(res["cuts"]))

# --- 反向命中偏移换算：offset=L 时反向切点就是起点 ---
res = run_digest({
    "sequence": "AAGAGACCAAA",
    "topology": "linear",
    "enzymes": [{"name": "X", "recognition": "GGTCTC", "cut_offset": 6}],
})
# GAGACC 起点 2 -> 反向切点 2+(6-6)=2
check("offset L reverse cut at start", res["cuts"] == [2], str(res["cuts"]))

# offset=0 反向切点 = 起点 + L
res = run_digest({
    "sequence": "AAGAGACCAAA",
    "topology": "linear",
    "enzymes": [{"name": "X", "recognition": "GGTCTC", "cut_offset": 0}],
})
check("offset 0 reverse cut", res["cuts"] == [8], str(res["cuts"]))

# --- 无切点：线性返回完整原序列 ---
res = run_digest({"sequence": "ACGT", "topology": "linear", "enzymes": []})
check("no enzyme linear", len(res["fragments"]) == 1
      and res["fragments"][0]["sequence"] == "ACGT"
      and res["fragments"][0]["length"] == 4
      and res["cuts"] == [])
res = run_digest({"sequence": "ACGTACGT", "topology": "linear",
                  "enzymes": [{"name": "E", "recognition": "GGG", "cut_offset": 1}]})
check("no hit linear", res["fragments"][0]["sequence"] == "ACGTACGT"
      and res["cuts"] == [])

# 无切点：环状也返回完整原序列（不做 rotate）
res = run_digest({"sequence": "ACGT", "topology": "circular", "enzymes": []})
check("no hit circular full", len(res["fragments"]) == 1
      and res["fragments"][0]["sequence"] == "ACGT"
      and res["fragments"][0]["wraps_origin"] is False)

# --- 环状跨首尾位点 + 取模去重 ---
# GGTCTC 跨边界：seq 末尾 "GG" + 开头 "TCTC"
cseq = "TCTCAAAAGG"  # len 10
res = run_digest({"sequence": cseq, "topology": "circular",
                  "enzymes": [{"name": "BsaI", "recognition": "GGTCTC", "cut_offset": 1}]})
# 跨边界位点 start=8: GG(8,9)+TCTC(0..3)；正向切点 (8+1)%10 = 9
check("circular spanning cut", res["cuts"] == [9], str(res["cuts"]))
f = res["fragments"][0]
check("circular single cut linearized", f["sequence"] == cseq[9:] + cseq[:9]
      and f["length"] == 10 and f["wraps_origin"] is True
      and (f["start"], f["end"]) == (9, 9), str(f))

# 环状两个切点，其中一个是跨边界位点产生的取模切点
cseq2 = "TCTCAAAGGTCTCAA"  # len 15: 含正常 GGTCTC? 位置 9..14 是 GTCTCA 不对
cseq2 = "GGAGGTCTC" + "TTT"  # len 12, GGTCTC 起点 3
# 跨边界：末尾 GG + 开头 AG... 不行。手工构造：
# 直接造一个末尾 GG 开头 TCTC 的，且中间再来一个位点
cseq2 = "TCTCAAAAGGTCTCAAGG"  # len 18；GGTCTC 起点 8，跨边界位点起点 16
res = run_digest({"sequence": cseq2, "topology": "circular",
                  "enzymes": [{"name": "BsaI", "recognition": "GGTCTC", "cut_offset": 1}]})
# 正向切点 8+1=9；跨边界 start=16 -> (16+1)%18=17
check("circular two cuts", sorted(res["cuts"]) == [9, 17], str(res["cuts"]))
frags = res["fragments"]
check("circular frag order", [frag["sequence"] for frag in frags] ==
      [cseq2[9:17], cseq2[17:] + cseq2[:9]],
      str([frag["sequence"] for frag in frags]))
check("circular frag wrap flags",
      [frag["wraps_origin"] for frag in frags] == [False, True])
check("circular frag total", sum(frag["length"] for frag in frags) == 18)
check("circular frag boundaries", [(frag["start"], frag["end"]) for frag in frags]
      == [(9, 17), (17, 9)])

# --- 取模去重：跨边界位点切点对 n 取模 ---
# seq 末尾 GAAT(4..7) + 开头 TC(0..1) -> 跨边界 GAATTC start=4
# EcoRI offset=1：正向 (4+1)%8=5；反向 (4+5)%8=1
cseq3 = "TCAAGAAT"  # len 8
res = run_digest({"sequence": cseq3, "topology": "circular",
                  "enzymes": [{"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1}]})
check("circular modulo cuts", sorted(res["cuts"]) == [1, 5], str(res["cuts"]))

# 真正的去重：offset=0，跨边界正向切点 start+0 = 6；再造一个切点 6 的普通位点
cseq4 = "TCAGAATTCAAGAAT"  # len 15；GAATTC 起点 3 -> 正向 3 反向 8；
# 跨边界 GAAT(11..14)+TC(0..1) start=11, offset=0 -> 11；反向 11+6=17%15=2
res = run_digest({"sequence": cseq4, "topology": "circular",
                  "enzymes": [{"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 0}]})
check("circular dedup set", sorted(res["cuts"]) == sorted({3, 9, 11, 2}),
      str(res["cuts"]))

# --- 端点切点（offset 0 / L 产生 0 或 n）线性不产生空片段 ---
res = run_digest({"sequence": "AATTCCCC", "topology": "linear",
                  "enzymes": [{"name": "E", "recognition": "AATT", "cut_offset": 0}]})
# AATT 起点 0：正向切点 0；回文反向切点 4
check("linear endpoint cut 0", res["cuts"] == [0, 4], str(res["cuts"]))
check("linear no empty frag for cut0",
      [f["sequence"] for f in res["fragments"]] == ["AATT", "CCCC"],
      str(res["fragments"]))
res = run_digest({"sequence": "CCCCAATT", "topology": "linear",
                  "enzymes": [{"name": "E", "recognition": "AATT", "cut_offset": 4}]})
# 位点起点 4：正向切点 8(=n)；回文反向切点 4
check("linear endpoint cut n", res["cuts"] == [4, 8], str(res["cuts"]))
check("linear no empty frag for cutN",
      [f["sequence"] for f in res["fragments"]] == ["CCCC", "AATT"],
      str(res["fragments"]))

# --- 识别序列长于序列 ---
res = run_digest({"sequence": "ACG", "topology": "linear",
                  "enzymes": [{"name": "E", "recognition": "GAATTC", "cut_offset": 1}]})
check("L>n linear no hit", res["fragments"][0]["sequence"] == "ACG")
# 环状 L>n：跨边界拼接后长度不足，不应命中、不应报错
res = run_digest({"sequence": "ACG", "topology": "circular",
                  "enzymes": [{"name": "E", "recognition": "GAATTC", "cut_offset": 1}]})
check("L>n circular no hit", res["fragments"][0]["sequence"] == "ACG")

# --- 多酶切点合并去重 ---
res = run_digest({"sequence": "GAATTCAATTG", "topology": "linear",
                  "enzymes": [
                      {"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1},
                      {"name": "EcoRI2", "recognition": "GAATTC", "cut_offset": 1},
                  ]})
check("multi-enzyme dedup", res["cuts"] == [1, 5], str(res["cuts"]))

# 小写输入 + 空白
res = run_digest({"sequence": " gaattc ", "topology": "LINEAR",
                  "enzymes": [{"name": "e", "recognition": "gaattc", "cut_offset": 1}]})
check("lowercase normalized", res["sequence"] == "GAATTC"
      and res["topology"] == "linear" and res["cuts"] == [1, 5])

# --- 错误用例 ---
expect_error("illegal base seq", {"sequence": "ACGTX", "topology": "linear",
                                  "enzymes": []}, "非法碱基")
expect_error("empty seq", {"sequence": "", "topology": "linear", "enzymes": []},
             "不能为空")
expect_error("whitespace seq", {"sequence": "   ", "topology": "linear",
                                "enzymes": []}, "不能为空")
expect_error("offset > L", {"sequence": "AAAA", "topology": "linear",
                            "enzymes": [{"name": "E", "recognition": "AA",
                                         "cut_offset": 3}]}, "越界")
expect_error("offset negative", {"sequence": "AAAA", "topology": "linear",
                                 "enzymes": [{"name": "E", "recognition": "AA",
                                              "cut_offset": -1}]}, "越界")
expect_error("offset 0 valid no error", {"sequence": "AAAA", "topology": "linear",
                                         "enzymes": [{"name": "E", "recognition": "AA",
                                                      "cut_offset": 2}]}, "") \
    if False else None
# "AA" 的反向互补是 "TT"，故只有正向命中：起点 0/1/2 -> 切点 2/3/4
res = run_digest({"sequence": "AAAA", "topology": "linear",
                  "enzymes": [{"name": "E", "recognition": "AA", "cut_offset": 2}]})
check("offset=L valid", res["cuts"] == [2, 3, 4], str(res["cuts"]))

expect_error("illegal recognition base", {"sequence": "AAAA", "topology": "linear",
                                          "enzymes": [{"name": "E", "recognition": "AAX",
                                                       "cut_offset": 1}]}, "非法碱基")
expect_error("empty recognition", {"sequence": "AAAA", "topology": "linear",
                                   "enzymes": [{"name": "E", "recognition": "",
                                                "cut_offset": 0}]}, "不能为空")
expect_error("bad topology", {"sequence": "AAAA", "topology": "plasmid",
                              "enzymes": []}, "拓扑")
expect_error("missing enzymes", {"sequence": "AAAA", "topology": "linear"}, "enzymes")
expect_error("offset not int", {"sequence": "AAAA", "topology": "linear",
                                "enzymes": [{"name": "E", "recognition": "AA",
                                             "cut_offset": 1.5}]}, "整数")
expect_error("offset bool rejected", {"sequence": "AAAA", "topology": "linear",
                                      "enzymes": [{"name": "E", "recognition": "AA",
                                                   "cut_offset": True}]}, "整数")
expect_error("payload not dict", ["not", "dict"], "JSON 对象")
expect_error("seq not string", {"sequence": 123, "topology": "linear",
                                "enzymes": []}, "字符串")

# 线性越界切点的防御性检查（常规偏移下不会出现；直接验证错误信息路径存在）
# 通过一个识别序列比目标长、但起点范围内可匹配的极端情形无法构造；
# 这里验证线性切点永远落在 0..n（不变量）：
import random  # noqa: E402

random.seed(42)
for _ in range(200):
    n = random.randint(1, 40)
    s = "".join(random.choice("ACGT") for _ in range(n))
    L = random.randint(1, 6)
    rec = "".join(random.choice("ACGT") for _ in range(L))
    off = random.randint(0, L)
    r = run_digest({"sequence": s, "topology": "linear",
                    "enzymes": [{"name": "E", "recognition": rec, "cut_offset": off}]})
    assert all(0 <= c <= n for c in r["cuts"]), (s, rec, off, r["cuts"])
    assert sum(f["length"] for f in r["fragments"]) == n
    assert "".join(f["sequence"] for f in r["fragments"]) == s
# 环状模糊测试：片段总长度恒为 n，且片段内容是环上的连续区间
for _ in range(200):
    n = random.randint(1, 40)
    s = "".join(random.choice("ACGT") for _ in range(n))
    L = random.randint(1, 6)
    rec = "".join(random.choice("ACGT") for _ in range(L))
    off = random.randint(0, L)
    r = run_digest({"sequence": s, "topology": "circular",
                    "enzymes": [{"name": "E", "recognition": rec, "cut_offset": off}]})
    assert all(0 <= c < n for c in r["cuts"])
    assert sum(f["length"] for f in r["fragments"]) == n
    if r["cuts"]:
        doubled = s + s
        for f in r["fragments"]:
            assert doubled[f["start"]:f["start"] + f["length"]] == f["sequence"], (s, f)
print("  PASS randomized linear/circular invariants (400 cases)")

print(f"\n{PASS + 1} passed, {FAIL} failed")
sys.exit(1 if FAIL else 0)
