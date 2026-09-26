"""DNA 限制酶切核心逻辑（纯标准库实现，不依赖 Web 框架，便于独立测试）。"""

from __future__ import annotations

from typing import Any

_BASES = "ACGT"
_TOPOLOGIES = ("linear", "circular")


class DigestError(ValueError):
    """输入非法或切点越界时抛出，reason 为可直接展示给调用方的中文原因。"""

    def __init__(self, reason: str):
        super().__init__(reason)
        self.reason = reason


def reverse_complement(seq: str) -> str:
    """返回 DNA 序列的反向互补序列（仅含 A/C/G/T）。"""
    table = {"A": "T", "T": "A", "C": "G", "G": "C"}
    return "".join(table[b] for b in reversed(seq))


def _is_real_int(value: Any) -> bool:
    # bool 是 int 的子类，需要显式排除
    return isinstance(value, int) and not isinstance(value, bool)


def _validate_sequence(raw: Any) -> str:
    if not isinstance(raw, str):
        raise DigestError("sequence 必须是字符串类型")
    seq = raw.strip().upper()
    if len(seq) == 0:
        raise DigestError("序列不能为空（或仅含空白字符）")
    illegal = sorted({b for b in seq if b not in _BASES})
    if illegal:
        raise DigestError(
            f"序列含非法碱基 {illegal}；仅允许 A/C/G/T（大小写均可）"
        )
    return seq


def _validate_topology(raw: Any) -> str:
    if not isinstance(raw, str):
        raise DigestError("topology 必须是字符串，取值 linear 或 circular")
    topology = raw.strip().lower()
    if topology not in _TOPOLOGIES:
        raise DigestError(
            f"拓扑类型非法：{raw!r}；仅支持 linear（线性）或 circular（环状）"
        )
    return topology


def _validate_enzymes(raw: Any) -> list[dict]:
    if raw is None:
        raise DigestError("缺少 enzymes 字段；无酶时请传空数组 []")
    if not isinstance(raw, list):
        raise DigestError("enzymes 必须是数组")
    enzymes: list[dict] = []
    for i, item in enumerate(raw):
        label = f"第 {i + 1} 个酶"
        if not isinstance(item, dict):
            raise DigestError(f"{label} 必须是对象")

        name = item.get("name", "")
        if name is None:
            name = ""
        if not isinstance(name, str):
            raise DigestError(f"{label} 的 name 必须是字符串")
        label = f"酶 {name!r}" if name else label

        recognition = item.get("recognition")
        if not isinstance(recognition, str):
            raise DigestError(f"{label} 的 recognition 必须是字符串")
        recognition = recognition.strip().upper()
        if len(recognition) == 0:
            raise DigestError(f"{label} 的识别序列不能为空")
        illegal = sorted({b for b in recognition if b not in _BASES})
        if illegal:
            raise DigestError(
                f"{label} 的识别序列含非法碱基 {illegal}；仅允许 A/C/G/T"
            )

        offset = item.get("cut_offset")
        if not _is_real_int(offset):
            raise DigestError(
                f"{label} 的 cut_offset 必须是整数（从识别序列左端起计，含两端）"
            )
        length = len(recognition)
        if offset < 0 or offset > length:
            raise DigestError(
                f"{label} 的切点偏移 {offset} 越界；识别长度为 {length}，"
                f"合法范围为 0..{length}（含两端）"
            )

        enzymes.append(
            {"name": name, "recognition": recognition, "cut_offset": offset}
        )
    return enzymes


def _find_cuts(
    seq: str, enzymes: list[dict], topology: str
) -> tuple[set[int], list[dict]]:
    """返回 (去重后的切点集合, 命中明细)。

    切点使用 0 基索引：k 表示切在第 k 个碱基与第 k+1 个碱基之间
    （k=0 表示最左端，k=len(seq) 表示最右端）。
    """
    n = len(seq)
    circular = topology == "circular"
    cuts: set[int] = set()
    hits: list[dict] = []

    for enzyme in enzymes:
        recognition = enzyme["recognition"]
        offset = enzyme["cut_offset"]
        length = len(recognition)
        rc = reverse_complement(recognition)
        enzyme_hits: list[dict] = []

        if circular:
            # 首尾相接：最多额外拼接 length-1 个碱基即可覆盖所有跨边界位点
            search = seq + seq[: length - 1]
            starts = range(n)  # 环上只有 n 个不同的位点起点
        else:
            search = seq
            starts = range(max(0, n - length + 1))

        for start in starts:
            window = search[start : start + length]
            if len(window) < length:
                continue
            if window == recognition:
                cut = start + offset
                if not circular and not (0 <= cut <= n):
                    raise DigestError(
                        f"酶 {enzyme['name']!r} 在线性序列位置 {start} 处的正向命中"
                        f"产生越界切点 {cut}（合法范围 0..{n}）"
                    )
                if circular:
                    cut %= n
                cuts.add(cut)
                enzyme_hits.append(
                    {"site_start": start, "strand": "forward", "cut": cut}
                )
            if window == rc:
                # 反向命中：切点按 识别长度 - 偏移 换算
                cut = start + (length - offset)
                if not circular and not (0 <= cut <= n):
                    raise DigestError(
                        f"酶 {enzyme['name']!r} 在线性序列位置 {start} 处的反向命中"
                        f"产生越界切点 {cut}（合法范围 0..{n}）"
                    )
                if circular:
                    cut %= n
                cuts.add(cut)
                enzyme_hits.append(
                    {"site_start": start, "strand": "reverse", "cut": cut}
                )

        hits.append(
            {
                "name": enzyme["name"],
                "recognition": recognition,
                "cut_offset": offset,
                "hits": enzyme_hits,
            }
        )

    return cuts, hits


def _build_fragments(seq: str, cuts: set[int], topology: str) -> list[dict]:
    n = len(seq)

    # 无切点：两种拓扑均返回完整原序列
    if not cuts:
        return [
            {
                "index": 0,
                "sequence": seq,
                "length": n,
                "start": 0,
                "end": n,
                "wraps_origin": False,
            }
        ]

    if topology == "linear":
        # 0 和 n 处的切点只切在端点，不产生额外片段
        inner = sorted(c for c in cuts if 0 < c < n)
        bounds = [0, *inner, n]
        fragments = []
        for i, (a, b) in enumerate(zip(bounds, bounds[1:])):
            fragments.append(
                {
                    "index": i,
                    "sequence": seq[a:b],
                    "length": b - a,
                    "start": a,
                    "end": b,
                    "wraps_origin": False,
                }
            )
        return fragments

    # 环状：切点已对 n 取模去重；从最小切点起按原链顺序排列，
    # 最后一个片段跨过索引 0（首尾相接处）
    points = sorted(cuts)
    pairs = list(zip(points, [*points[1:], points[0]]))
    fragments = []
    for i, (a, b) in enumerate(pairs):
        if a == b:  # 仅一个切点：环在此处线性化
            fragment = seq[a:] + seq[:a]
            wraps = True
        elif a < b:
            fragment = seq[a:b]
            wraps = False
        else:  # 跨首尾的片段
            fragment = seq[a:] + seq[:b]
            wraps = True
        length = (b - a) % n or n
        fragments.append(
            {
                "index": i,
                "sequence": fragment,
                "length": length,
                "start": a,
                "end": b,
                "wraps_origin": wraps,
            }
        )
    return fragments


def run_digest(payload: Any) -> dict:
    """校验请求体并执行酶切，返回可直接 JSON 序列化的响应字典。"""
    if not isinstance(payload, dict):
        raise DigestError("请求体必须是 JSON 对象")

    seq = _validate_sequence(payload.get("sequence"))
    topology = _validate_topology(payload.get("topology"))
    enzymes = _validate_enzymes(payload.get("enzymes"))

    cuts, hits = _find_cuts(seq, enzymes, topology)
    fragments = _build_fragments(seq, cuts, topology)

    return {
        "sequence": seq,
        "length": len(seq),
        "topology": topology,
        "cuts": sorted(cuts),
        "fragments": fragments,
        "enzyme_hits": hits,
    }
