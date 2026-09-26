"""DNA 限制酶切核心逻辑（纯标准库实现，不依赖 Web 框架，便于独立测试）。"""

from __future__ import annotations

import itertools
from typing import Any

_BASES = "ACGT"
_TOPOLOGIES = ("linear", "circular")

_PANEL_MIN_SAMPLES = 2
_PANEL_MAX_SAMPLES = 8
_PANEL_MAX_SEQUENCE_LENGTH = 2000
_PANEL_MIN_ENZYMES = 1
_PANEL_MAX_ENZYMES = 10


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


# ---------------------------------------------------------------------------
# /design-panel：为多个 DNA 样本寻找能相互鉴别的最少酶组合
# ---------------------------------------------------------------------------


def _validate_samples(raw: Any) -> list[dict]:
    if raw is None:
        raise DigestError("缺少 samples 字段")
    if not isinstance(raw, list):
        raise DigestError("samples 必须是数组")
    if not (_PANEL_MIN_SAMPLES <= len(raw) <= _PANEL_MAX_SAMPLES):
        raise DigestError(
            f"samples 数量必须在 {_PANEL_MIN_SAMPLES}..{_PANEL_MAX_SAMPLES} 之间，"
            f"当前为 {len(raw)}"
        )

    samples: list[dict] = []
    seen_names: set[str] = set()
    for i, item in enumerate(raw):
        label = f"第 {i + 1} 个样本"
        if not isinstance(item, dict):
            raise DigestError(f"{label} 必须是对象")

        name = item.get("name")
        if not isinstance(name, str):
            raise DigestError(f"{label} 的 name 必须是字符串")
        name = name.strip()
        if not name:
            raise DigestError(f"{label} 的 name 不能为空（或仅含空白字符）")
        if name in seen_names:
            raise DigestError(f"样本名称重复：{name!r}；样本名称必须唯一")
        seen_names.add(name)

        seq = _validate_sequence(item.get("sequence"))
        if len(seq) > _PANEL_MAX_SEQUENCE_LENGTH:
            raise DigestError(
                f"样本 {name!r} 的序列长度为 {len(seq)}，"
                f"超过上限 {_PANEL_MAX_SEQUENCE_LENGTH} 碱基"
            )
        samples.append({"name": name, "sequence": seq})
    return samples


def _validate_tolerance(raw: Any) -> int:
    if not _is_real_int(raw):
        raise DigestError("tolerance 必须是非负整数")
    if raw < 0:
        raise DigestError("tolerance 必须是非负整数（不能为负数）")
    return raw


def _length_profile(seq: str, cuts: set[int], topology: str) -> list[int]:
    """联合酶切后的片段长度谱：全部片段长度升序排列，重复长度保留。"""
    n = len(seq)
    if topology == "linear":
        inner = sorted(c for c in cuts if 0 < c < n)
        bounds = [0, *inner, n]
        return sorted(b - a for a, b in zip(bounds, bounds[1:]))

    # 环状切点均已对 n 取模
    points = sorted(cuts)
    if not points:
        return [n]
    profile = sorted((b - a) % n or n for a, b in zip(points, [*points[1:], points[0]]))
    return profile


def _indistinguishable(profile_a: list[int], profile_b: list[int], tolerance: int) -> bool:
    """两谱片段数相同且对应（升序后）长度差均不超过容差时视为不可区分。"""
    if len(profile_a) != len(profile_b):
        return False
    return all(abs(a - b) <= tolerance for a, b in zip(profile_a, profile_b))


def run_design_panel(payload: Any) -> dict:
    """校验请求体并搜索最少鉴别酶组合，返回可直接 JSON 序列化的响应字典。"""
    if not isinstance(payload, dict):
        raise DigestError("请求体必须是 JSON 对象")

    samples = _validate_samples(payload.get("samples"))
    topology = _validate_topology(payload.get("topology"))
    enzymes = _validate_enzymes(payload.get("enzymes"))
    if not (_PANEL_MIN_ENZYMES <= len(enzymes) <= _PANEL_MAX_ENZYMES):
        raise DigestError(
            f"候选酶数量必须在 {_PANEL_MIN_ENZYMES}..{_PANEL_MAX_ENZYMES} 之间，"
            f"当前为 {len(enzymes)}"
        )
    tolerance = _validate_tolerance(payload.get("tolerance"))

    sequences = [s["sequence"] for s in samples]

    # 对“每个样本 × 每个候选酶”做一次单酶酶切并缓存切点。
    # 任意酶组合的联合酶切等价于各酶切点集合的并（切点位置只依赖序列与酶，
    # 与是否与其他酶同时切无关），因此枚举组合时只需合并切点集合；
    # 这是对联合酶切结果的直接计算，而非合并单酶的鉴别结论。
    cut_table: list[list[set[int]]] = []
    for seq in sequences:
        row: list[set[int]] = []
        for enzyme in enzymes:
            cuts, _hits = _find_cuts(seq, [enzyme], topology)
            row.append(set(cuts))
        cut_table.append(row)

    # 按酶数最少、候选酶输入下标列表字典序最小枚举所有非空组合。
    # 每个组合独立判断：必须同时区分所有样本对（不能只合并各单酶的鉴别结论）。
    m = len(enzymes)
    sample_pairs = list(itertools.combinations(range(len(samples)), 2))
    indices: list[int] | None = None
    for size in range(1, m + 1):
        for combo in itertools.combinations(range(m), size):
            profiles = [
                _length_profile(
                    sequences[j],
                    set().union(*(cut_table[j][e] for e in combo)),
                    topology,
                )
                for j in range(len(samples))
            ]
            if all(
                not _indistinguishable(
                    profiles[a], profiles[b], tolerance
                )
                for a, b in sample_pairs
            ):
                indices = list(combo)
                break
        if indices is not None:
            break

    if indices is None:
        return {"status": "no_solution"}

    selected_profiles = []
    for j, sample in enumerate(samples):
        cuts: set[int] = set().union(*(cut_table[j][e] for e in indices))
        profile = _length_profile(sequences[j], cuts, topology)
        selected_profiles.append(
            {
                "name": sample["name"],
                "length": len(sequences[j]),
                "fragment_count": len(profile),
                "fragment_lengths": profile,
            }
        )

    return {
        "status": "solution",
        "enzyme_indices": indices,
        "enzyme_count": len(indices),
        "profiles": selected_profiles,
    }
