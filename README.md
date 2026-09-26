# DNA 限制酶切模拟 API

基于 Python + FastAPI 实现的 DNA 限制性内切酶酶切模拟服务。输入一条仅含 A/C/G/T 的 DNA
序列、其拓扑结构（线性 / 环状）以及若干酶的识别序列与切点偏移，服务会：

- 同时搜索**正向识别序列**和其**反向互补序列**（正反两条链上的命中）；
- 按切点把序列切成片段，返回每个片段的**序列、长度和边界**；
- 环状序列识别**跨越首尾连接处**的位点，切点对序列长度**取模并去重**；
- 线性序列对越界切点**报错拒绝**；
- 无任何切点时返回**完整原序列**作为唯一片段；
- 对非法碱基、空序列、偏移越界等输入返回 **400** 并写明中文原因。

`POST /design-panel` 进一步支持一次提交 2～8 个 DNA 样本与 1～10 个候选酶，
**搜索能用联合酶切谱把任意两个样本相互区分开的最少酶组合**（详见下文）。

## 单条命令启动（一键启动）

需要 Docker（含 Compose 插件 v2）。在项目根目录执行：

```bash
docker compose up --build
```

该命令会自动完成：构建镜像 → 安装 `requirements.txt` 中的全部依赖（FastAPI、uvicorn）
→ 构建期导入校验（应用与 `/digest`、`/design-panel` 路由可正常装配）→ 启动服务。
首次构建需要联网下载依赖，之后的重建会命中 Docker 层缓存。加 `-d` 可后台运行
（`docker compose up --build -d`）。

服务启动后通过宿主机 **8000 端口**访问：

| 用途 | 地址 |
| --- | --- |
| 接口文档（Swagger UI，可直接在线调试） | http://localhost:8000/docs |
| OpenAPI Schema | http://localhost:8000/openapi.json |
| 健康检查 | http://localhost:8000/health |
| 单样本酶切接口 | `POST http://localhost:8000/digest` |
| 多样本鉴别酶组合设计接口 | `POST http://localhost:8000/design-panel` |

停止服务：前台运行按 `Ctrl+C`；后台运行执行：

```bash
docker compose down
```

## 必要配置

默认无需任何配置即可运行，容器内固定监听 `0.0.0.0:8000`。唯一常用项是宿主机映射端口
（默认 `8000`，即上面的访问地址），有两种改法：

```bash
# 方式一：环境变量（docker compose 自动读取同目录 .env）
echo "HOST_PORT=9000" > .env
docker compose up --build -d

# 方式二：命令行临时指定
HOST_PORT=9000 docker compose up --build -d
```

容器内监听地址/端口由 Compose 环境变量 `HOST`、`PORT` 控制（默认 `0.0.0.0:8000`，
一般无需修改），通过端口映射暴露到宿主机。

## API 说明

### `POST /digest`

**请求体（JSON）：**

| 字段 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `sequence` | string | 是 | DNA 序列，仅含 `A/C/G/T`（大小写、首尾空白均可，内部会归一化）；不得为空 |
| `topology` | string | 是 | `linear`（线性）或 `circular`（环状），大小写不敏感 |
| `enzymes` | array | 是 | 酶定义数组；传 `[]` 表示不酶切（直接返回完整序列） |
| `enzymes[].name` | string | 否 | 酶名称，仅用于结果标注，可留空 |
| `enzymes[].recognition` | string | 是 | 识别序列，仅含 `A/C/G/T`，不得为空 |
| `enzymes[].cut_offset` | integer | 是 | 切点偏移，从**正向识别序列左端**起计、含两端，合法范围 `0 .. 识别序列长度` |

**切点偏移语义（0 基索引）：** 若识别序列长度为 L、在序列中的起点为 `s`，则

- 正向命中切点：`s + cut_offset`
- 反向（反向互补链）命中切点：`s + (L − cut_offset)`

例如 EcoRI 识别 `GAATTC`（L=6），切点在 G 与 A 之间，`cut_offset = 1`：
正向命中切点为 `s+1`，反向互补命中切点为 `s+5`。
`cut_offset = 0` 表示切在识别位点左端之外，`cut_offset = L` 表示右端之外。

**线性拓扑：** 切点合法范围为 `0..N`（N 为序列长度）；若任何命中产生越界切点，
整个请求返回 400。位于 `0` 或 `N` 的端点切点只在端点处断开，不产生空片段。

**环状拓扑：** 把首尾相连后搜索，因此可识别跨首尾的位点；所有切点对 N 取模并跨酶
去重。片段按原链方向、从最小切点开始依次排列；跨过索引 0 的片段排在最后，并以
`"wraps_origin": true` 标记。仅有一个（去重后）切点时，环在该处线性化为一条全长片段。

**片段边界：** 采用 0 基半开区间 `[start, end)`。线性序列的片段不环绕；环状序列中
环绕片段满足“内容 = sequence[start:] + sequence[:end]”，单切点片段 `start == end`。

### 请求示例

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "AAGAATTCTTTGAATTCGG",
    "topology": "linear",
    "enzymes": [{"name": "EcoRI", "recognition": "GAATTC", "cut_offset": 1}]
  }'
```

环状序列识别跨首尾位点的示例：

```bash
curl -s -X POST http://localhost:8000/digest \
  -H 'Content-Type: application/json' \
  -d '{
    "sequence": "TCTCAAAAGG",
    "topology": "circular",
    "enzymes": [{"name": "BsaI", "recognition": "GGTCTC", "cut_offset": 1}]
  }'
```

### 响应示例

```json
{
  "sequence": "AAGAATTCTTTGAATTCGG",
  "length": 19,
  "topology": "linear",
  "cuts": [3, 7, 12, 16],
  "fragments": [
    {"index": 0, "sequence": "AAG",   "length": 3, "start": 0,  "end": 3,  "wraps_origin": false},
    {"index": 1, "sequence": "AATT",  "length": 4, "start": 3,  "end": 7,  "wraps_origin": false},
    {"index": 2, "sequence": "CTTTG", "length": 5, "start": 7,  "end": 12, "wraps_origin": false},
    {"index": 3, "sequence": "AATT",  "length": 4, "start": 12, "end": 16, "wraps_origin": false},
    {"index": 4, "sequence": "CGG",   "length": 3, "start": 16, "end": 19, "wraps_origin": false}
  ],
  "enzyme_hits": [
    {
      "name": "EcoRI",
      "recognition": "GAATTC",
      "cut_offset": 1,
      "hits": [
        {"site_start": 2,  "strand": "forward", "cut": 3},
        {"site_start": 2,  "strand": "reverse", "cut": 7},
        {"site_start": 11, "strand": "forward", "cut": 12},
        {"site_start": 11, "strand": "reverse", "cut": 16}
      ]
    }
  ]
}
```

无切点（无酶或没有命中）时，两种拓扑都返回一条完整原序列片段：

```json
{"cuts": [], "fragments": [{"index": 0, "sequence": "原序列", "length": 4,
  "start": 0, "end": 4, "wraps_origin": false}]}
```

### `POST /design-panel`

为多个 DNA 样本寻找**能相互鉴别的最少酶组合**。对每个候选酶组合，服务会把组合内的酶
**联合酶切每一个样本**（复用与 `/digest` 完全相同的酶切规则：正反向识别、环状跨首尾、
切点取模去重等），把每个样本得到的全部片段长度**升序排列、保留重复长度**，作为该样本
在该组合下的片段长度谱。

**两样本不可区分的判定：** 两谱的**片段数相同**，且升序后**对应位置长度之差的绝对值都
不超过** `tolerance`（长度容差）。任一条件不满足（片段数不同，或存在某一对长度差超过
容差）即可区分。所选组合必须让**任意两个样本都能相互区分**。

> 注意：判定针对的是**所选酶联合酶切后的完整长度谱**，而不是“各单酶能不能区分”的简单
> 合并。确实存在“每个单酶单独都无法区分、但联合酶切后可以区分”的样本对，本接口会正确
> 找出这种组合。

**选择规则（确定性）：** 在所有满足条件的**非空**候选酶组合中，

1. 先取**酶数最少**的；
2. 酶数相同则取候选酶**输入下标列表字典序最小**的（下标从 0 开始，按输入顺序）。

不存在任何可行组合时返回 `{"status": "no_solution"}`。

**请求体（JSON）：**

| 字段 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `samples` | array | 是 | 样本数组，数量 `2..8` |
| `samples[].name` | string | 是 | 样本唯一名称，去首尾空白后不得为空且不得与其他样本重名 |
| `samples[].sequence` | string | 是 | DNA 序列，规则同 `/digest`，长度不得超过 `2000` 碱基 |
| `topology` | string | 是 | 所有样本**共同**的拓扑：`linear` 或 `circular` |
| `enzymes` | array | 是 | 候选酶数组，数量 `1..10`，每项格式与 `/digest` 的 `enzymes[]` 完全一致（允许重复定义） |
| `tolerance` | integer | 是 | 非负整数长度容差；逐位长度差 `<= tolerance` 视为相等 |

**响应（找到方案）：**

| 字段 | 说明 |
| --- | --- |
| `status` | 固定 `"solution"` |
| `enzyme_indices` | 所选酶在候选酶输入数组中的下标列表（升序、去重），即酶数最少且字典序最小者 |
| `enzyme_count` | 所选酶数量 |
| `profiles` | 各样本在所选组合**联合酶切**下的片段谱，顺序与输入 `samples` 一致 |
| `profiles[].name` / `length` | 样本名 / 序列长度 |
| `profiles[].fragment_count` | 片段数量 |
| `profiles[].fragment_lengths` | 全部片段长度，**升序排列并保留重复长度** |

找不到方案时仅返回：

```json
{"status": "no_solution"}
```

### 请求示例

```bash
curl -s -X POST http://localhost:8000/design-panel \
  -H 'Content-Type: application/json' \
  -d '{
    "samples": [
      {"name": "sample-1", "sequence": "AAAAAAGAACGTAAAAAAAAGAAAAAAAAAACGTAAAAAA"},
      {"name": "sample-2", "sequence": "AAAAAAAAAAAAAAAAGAAACGTAAAACGTAAAAAAAAGA"}
    ],
    "topology": "linear",
    "enzymes": [
      {"name": "A", "recognition": "AAG", "cut_offset": 0},
      {"name": "B", "recognition": "GTA", "cut_offset": 1}
    ],
    "tolerance": 0
  }'
```

上例中两个酶**单独酶切**两样本的长度谱分别相同（均无法区分），但两酶**联合酶切**后
长度谱分别为 `[4,7,7,7,15]` 与 `[4,7,7,8,14]`，因此唯一可行（也是最少）的组合是下标
`[0, 1]`。响应：

```json
{
  "status": "solution",
  "enzyme_indices": [0, 1],
  "enzyme_count": 2,
  "profiles": [
    {"name": "sample-1", "length": 40, "fragment_count": 5,
     "fragment_lengths": [4, 7, 7, 7, 15]},
    {"name": "sample-2", "length": 40, "fragment_count": 5,
     "fragment_lengths": [4, 7, 7, 8, 14]}
  ]
}
```

容差示例：两样本某酶下长度谱为 `[3, 7]` 和 `[4, 6]`：`tolerance=0` 时可区分（有解），
`tolerance>=1` 时逐位差均不超过容差，视为不可区分（若无其他酶可用则返回
`no_solution`）。片段数不同时（如 `[8]` 与 `[3, 5]`），无论容差多大都可区分。

## 错误响应

所有输入类错误返回 HTTP `400`，响应体形如 `{"error": "中文原因说明"}`，例如：

- `序列含非法碱基 ['N']；仅允许 A/C/G/T（大小写均可）`
- `序列不能为空（或仅含空白字符）`
- `酶 'EcoRI' 的切点偏移 7 越界；识别长度为 6，合法范围为 0..6（含两端）`
- `拓扑类型非法：'plasmid'；仅支持 linear（线性）或 circular（环状）`
- `酶 'X' 的识别序列含非法碱基 ['Y']；仅允许 A/C/G/T`
- `酶 'X' 在线性序列位置 ... 产生越界切点 ...（合法范围 0..N）`
- `samples 数量必须在 2..8 之间，当前为 1`
- `样本名称重复：'s1'；样本名称必须唯一`
- `样本 's1' 的序列长度为 2001，超过上限 2000 碱基`
- `候选酶数量必须在 1..10 之间，当前为 0`
- `tolerance 必须是非负整数`

`/design-panel` 的无方案结果（`{"status": "no_solution"}`）属于正常业务结果，
返回 HTTP `200`，不是错误。

## 本地开发（不使用 Docker）

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

## 测试

核心酶切逻辑不依赖 Web 框架，可用系统 Python 直接运行（无需安装第三方依赖）：

```bash
python3 tests/test_digest.py
python3 tests/test_design_panel.py
```

`test_digest.py` 覆盖：正向/反向命中与切点换算、回文酶的成对切点、环状跨首尾位点、切点
取模去重、跨首尾片段顺序与边界、端点切点不产生空片段、识别序列长于序列、多酶去重、
各类非法输入报错，以及 400 组线性/环状随机不变量测试（片段长度之和恒等于 N、
片段按原链顺序可拼回原序列等）。

`test_design_panel.py` 覆盖：单酶谱相同但**联合酶切谱不同**的样本对（证明结果来自联合
酶切而非合并单酶鉴别结论）、最少酶数与下标字典序择优、重复酶定义、非切割酶、三样本
两两可区分、容差边界（片段数不同恒可区分、逐位差与容差比较）、相同序列 `no_solution`、
环状联合酶切、返回谱与 `/digest` 逐样本一致，以及样本数量/重名/长度上限/酶数量/
容差类型等各类非法输入 400。

## 项目结构

```
.
├── app/
│   ├── __init__.py
│   ├── digest.py          # 酶切与鉴别酶组合核心逻辑（纯函数，可独立测试）
│   └── main.py            # FastAPI 入口与错误处理（/digest、/design-panel）
├── tests/
│   ├── test_digest.py         # /digest 纯标准库测试脚本
│   └── test_design_panel.py   # /design-panel 纯标准库测试脚本
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
└── README.md
```
