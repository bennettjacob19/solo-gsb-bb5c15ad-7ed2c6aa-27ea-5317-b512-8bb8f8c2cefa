"""多样本酶切鉴别组合（design panel）核心逻辑。

纯标准库实现，复用 digest.py 的酶切规则，不依赖 Web 框架，便于独立测试。

对候选酶的每个非空子集做**真正的联合酶切**（多酶切点取并集后再切片段，
与 /digest 的多酶行为完全一致），取各样本片段长度的升序谱（保留重复长度），
按长度容差两两比对，找出能区分全部样本对的最小组合。
"""

from __future__ import annotations

from itertools import combinations
from typing import Any

from .digest import (
    DigestError,
    _build_fragments,
    _find_cuts,
    _is_real_int,
    _validate_enzymes,
    _validate_sequence,
    _validate_topology,
)

MIN_SAMPLES = 2
MAX_SAMPLES = 8
MAX_SEQUENCE_LENGTH = 2000
MIN_ENZYMES = 1
MAX_ENZYMES = 10


def _validate_samples(raw: Any) -> list[tuple[str, str]]:
    """校验样本数组，返回 [(归一化名称, 归一化序列), ...]（保持输入顺序）。"""
    if raw is None:
        raise DigestError("缺少 samples 字段")
    if not isinstance(raw, list):
        raise DigestError("samples 必须是数组")
    if not (MIN_SAMPLES <= len(raw) <= MAX_SAMPLES):
        raise DigestError(
            f"样本数量为 {len(raw)}；仅支持 {MIN_SAMPLES}..{MAX_SAMPLES} 个唯一样本"
        )

    samples: list[tuple[str, str]] = []
    seen_names: set[str] = set()
    for i, item in enumerate(raw):
        label = f"第 {i + 1} 个样本"
        if not isinstance(item, dict):
            raise DigestError(f"{label} 必须是对象")

        raw_name = item.get("name")
        if not isinstance(raw_name, str):
            raise DigestError(f"{label} 的 name 必须是字符串")
        name = raw_name.strip()
        if not name:
            raise DigestError(f"{label} 的名称不能为空（或仅含空白字符）")
        if name in seen_names:
            raise DigestError(f"样本名称重复：{name!r}；样本名称必须唯一")

        seq = _validate_sequence(item.get("sequence"))
        if len(seq) > MAX_SEQUENCE_LENGTH:
            raise DigestError(
                f"样本 {name!r} 的序列长度为 {len(seq)}，"
                f"超过上限 {MAX_SEQUENCE_LENGTH} 碱基"
            )

        seen_names.add(name)
        samples.append((name, seq))
    return samples


def _validate_tolerance(raw: Any) -> int:
    if not _is_real_int(raw):
        raise DigestError("length_tolerance 必须是非负整数")
    if raw < 0:
        raise DigestError(
            f"length_tolerance 必须是非负整数，不能为负数（收到 {raw}）"
        )
    return raw


def _fragment_length_profile(seq: str, cuts: set[int], topology: str) -> list[int]:
    """联合酶切后的片段长度谱：全部片段长度升序排列，重复长度保留。"""
    fragments = _build_fragments(seq, cuts, topology)
    return sorted(fragment["length"] for fragment in fragments)


def _indistinguishable(
    profile_a: list[int], profile_b: list[int], tolerance: int
) -> bool:
    """片段数相同且升序后对应位置长度差均不超过容差，即视为不可区分。"""
    if len(profile_a) != len(profile_b):
        return False
    return all(abs(a - b) <= tolerance for a, b in zip(profile_a, profile_b))


def run_design_panel(payload: Any) -> dict:
    """校验请求体并搜索最小区分酶组合，返回可直接 JSON 序列化的响应字典。"""
    if not isinstance(payload, dict):
        raise DigestError("请求体必须是 JSON 对象")

    topology = _validate_topology(payload.get("topology"))
    samples = _validate_samples(payload.get("samples"))
    enzymes = _validate_enzymes(payload.get("enzymes"))
    if not (MIN_ENZYMES <= len(enzymes) <= MAX_ENZYMES):
        raise DigestError(
            f"候选酶数量为 {len(enzymes)}；仅支持 {MIN_ENZYMES}..{MAX_ENZYMES} 个"
        )
    tolerance = _validate_tolerance(payload.get("length_tolerance"))

    # 预计算每个 (样本, 候选酶) 的切点集合。/digest 中多酶联合酶切就是把各酶
    # 切点并入同一个集合后统一构建片段，因此子集联合酶切 = 切点并集，
    # 这里在保证规则完全一致的同时避免对每个组合重复扫描序列。
    prepared: list[tuple[str, str, list[set[int]]]] = []
    for name, seq in samples:
        per_enzyme_cuts: list[set[int]] = []
        for enzyme in enzymes:
            cuts, _hits = _find_cuts(seq, [enzyme], topology)
            per_enzyme_cuts.append(cuts)
        prepared.append((name, seq, per_enzyme_cuts))

    pair_indexes = [
        (i, j)
        for i in range(len(prepared))
        for j in range(i + 1, len(prepared))
    ]

    n_enzymes = len(enzymes)
    # 按酶数升序；同酶数下 combinations 按下标列表字典序产出，首个可行组合即最优
    for size in range(1, n_enzymes + 1):
        for combo in combinations(range(n_enzymes), size):
            profiles: list[list[int]] = []
            for _name, seq, per_enzyme_cuts in prepared:
                cuts: set[int] = set()
                for enzyme_index in combo:
                    cuts |= per_enzyme_cuts[enzyme_index]
                profiles.append(_fragment_length_profile(seq, cuts, topology))

            if all(
                not _indistinguishable(profiles[i], profiles[j], tolerance)
                for i, j in pair_indexes
            ):
                return {
                    "status": "found",
                    "enzyme_indices": list(combo),
                    "profiles": [
                        {"name": name, "fragment_lengths": lengths}
                        for (name, _seq, _cuts), lengths in zip(prepared, profiles)
                    ],
                }

    return {"status": "no_solution"}
