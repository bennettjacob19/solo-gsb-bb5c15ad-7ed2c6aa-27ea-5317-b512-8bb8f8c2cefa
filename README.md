# DNA 限制酶切模拟 API

基于 Python + FastAPI 实现的 DNA 限制性内切酶酶切模拟服务。提供两个接口：

- `POST /digest`：对**单条**序列执行酶切模拟；
- `POST /design-panel`：对 **2..8 个 DNA 样本**与 **1..10 个候选酶**做联合酶切，
  搜索能让**任意两个样本相互鉴别**的**最少酶组合**（酶切面板设计）。

输入序列仅含 A/C/G/T，可指定线性 / 环状拓扑以及酶的识别序列与切点偏移，服务会：

- 同时搜索**正向识别序列**和其**反向互补序列**（正反两条链上的命中）；
- 按切点把序列切成片段，返回每个片段的**序列、长度和边界**；
- 环状序列识别**跨越首尾连接处**的位点，切点对序列长度**取模并去重**；
- 线性序列对越界切点**报错拒绝**；
- 无任何切点时返回**完整原序列**作为唯一片段；
- 对非法碱基、空序列、偏移越界等输入返回 **400** 并写明中文原因。

## 单条命令启动

需要 Docker（含 Compose 插件 v2）。在项目根目录执行：

```bash
docker compose up --build
```

加 `-d` 可在后台运行：`docker compose up --build -d`。首次执行会构建镜像（安装
`requirements.txt` 中的 FastAPI/uvicorn 依赖），并在构建阶段运行一次初始化自检，
确认应用可导入、`/digest` 与 `/design-panel` 路由均已注册。

服务启动后监听 **8000 端口**：

| 用途 | 地址 |
| --- | --- |
| 接口文档（Swagger UI，可直接在线调试） | http://localhost:8000/docs |
| OpenAPI Schema | http://localhost:8000/openapi.json |
| 健康检查 | http://localhost:8000/health |
| 单序列酶切接口 | `POST http://localhost:8000/digest` |
| 酶组合鉴别面板接口 | `POST http://localhost:8000/design-panel` |

停止服务：

```bash
docker compose down
```

## 必要配置

默认无需任何配置即可运行，开箱即在宿主机 `8000` 端口提供服务。唯一常用项是宿主机
端口（默认 `8000`，被占用时需要修改），有两种改法：

```bash
# 方式一：环境变量（docker compose 自动读取同目录 .env）
echo "HOST_PORT=9000" > .env
docker compose up --build -d

# 方式二：命令行临时指定
HOST_PORT=9000 docker compose up --build -d
```

容器内固定监听 `0.0.0.0:8000`（由 Dockerfile / Compose 中的 `HOST`、`PORT` 控制，
一般无需修改），通过端口映射暴露到宿主机。Compose 还开启了 `init: true`，使容器以
轻量 init 进程托管 uvicorn，停止时能干净转发信号、回收子进程。


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

为多个 DNA 样本寻找能**相互鉴别**的**最少候选酶组合**。每个组合对每个样本执行的是
**真正的多酶联合酶切**（把组合内各酶的切点取并集后再切片段，与 `/digest` 传多个酶
完全同规则），而不是把各酶单独鉴别的结果简单合并。

**请求体（JSON）：**

| 字段 | 类型 | 必填 | 说明 |
| --- | --- | --- | --- |
| `samples` | array | 是 | 样本数组，长度 **2..8** |
| `samples[].name` | string | 是 | 样本名称，去首尾空白后不得为空，且**必须唯一** |
| `samples[].sequence` | string | 是 | DNA 序列，仅含 `A/C/G/T`，非空，长度 **≤ 2000** 碱基（归一化后计长） |
| `topology` | string | 是 | 所有样本共用的拓扑：`linear` / `circular`，大小写不敏感 |
| `enzymes` | array | 是 | 候选酶数组，长度 **1..10**，字段格式与 `/digest` 的 `enzymes` 完全一致 |
| `length_tolerance` | integer | 是 | 长度容差，**非负整数**（不接受布尔值/浮点/字符串） |

**鉴别判据：** 对某个酶组合，联合酶切每个样本后把**全部片段长度升序排列**（重复长度
保留，即长度多重集）。两个样本当且仅当满足以下条件时视为**不可区分**：

1. 片段**数量相同**；
2. 升序后**对应位置**的长度差均 **≤ `length_tolerance`**。

只有所选组合让**任意一对**样本都可区分时才算可行。可行组合的择优规则为：

1. **酶数最少**（非空组合，酶数从 1 起递增搜索）；
2. 酶数相同时，选候选酶**输入下标列表字典序最小**的一个
   （下标从 0 起，按请求中 `enzymes` 的顺序）。

**成功响应（`200`）：**

```json
{
  "status": "found",
  "enzyme_indices": [0, 1],
  "profiles": [
    {"name": "s1", "fragment_lengths": [1, 6, 7, 10]},
    {"name": "s2", "fragment_lengths": [1, 6, 8, 9]}
  ]
}
```

- `enzyme_indices`：所选候选酶在请求 `enzymes` 中的下标（已排序）；
- `profiles`：按请求 `samples` 顺序给出每个样本的片段长度谱（升序、保留重复）。

**无解响应（`200`）：** 任何非空组合都不能区分全部样本对时：

```json
{"status": "no_solution"}
```

> 注意：单酶分别可区分的样本对取并集并不保证联合可区分性——鉴别基于联合酶切后的
> 完整长度谱。例如两个样本可能在酶 A、酶 B **各自**酶切下谱都相同，但 **A+B 联合**
> 酶切后切点交错产生不同的长度谱；本接口会直接枚举组合做联合酶切来处理这种情况。

**请求示例：**

```bash
curl -s -X POST http://localhost:8000/design-panel \
  -H 'Content-Type: application/json' \
  -d '{
    "samples": [
      {"name": "s1", "sequence": "GGTTTTGAGTTTACTCTTTTTTTT"},
      {"name": "s2", "sequence": "TTTTTTGAGTTTACTCTTTTAACC"}
    ],
    "topology": "linear",
    "enzymes": [
      {"name": "E1", "recognition": "GGTT", "cut_offset": 1},
      {"name": "E2", "recognition": "GAGT", "cut_offset": 2}
    ],
    "length_tolerance": 0
  }'
```

上例中两样本在 E1 单独酶切下谱均为 `[1,23]`、E2 单独酶切下均为 `[6,8,10]`，
只有 `E1+E2` 联合酶切才分别得到 `[1,6,7,10]` 与 `[1,6,8,9]`，因此返回
`"enzyme_indices": [0, 1]`。

## 错误响应

所有输入类错误返回 HTTP `400`，响应体形如 `{"error": "中文原因说明"}`，例如：

- `序列含非法碱基 ['N']；仅允许 A/C/G/T（大小写均可）`
- `序列不能为空（或仅含空白字符）`
- `酶 'EcoRI' 的切点偏移 7 越界；识别长度为 6，合法范围为 0..6（含两端）`
- `拓扑类型非法：'plasmid'；仅支持 linear（线性）或 circular（环状）`
- `酶 'X' 的识别序列含非法碱基 ['Y']；仅允许 A/C/G/T`
- `酶 'X' 在线性序列位置 ... 产生越界切点 ...（合法范围 0..N）`
- `样本数量为 1；仅支持 2..8 个唯一样本`
- `样本名称重复：'a'；样本名称必须唯一`
- `样本 'long' 的序列长度为 2001，超过上限 2000 碱基`
- `候选酶数量为 0；仅支持 1..10 个`
- `length_tolerance 必须是非负整数，不能为负数（收到 -1）`

## 本地开发（不使用 Docker）

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
uvicorn app.main:app --reload --port 8000
```

## 测试

核心酶切逻辑不依赖 Web 框架，可用系统 Python 直接运行（无需安装第三方依赖）：

```bash
python3 tests/test_digest.py    # /digest 核心规则
python3 tests/test_design.py    # /design-panel 组合搜索
```

覆盖内容包括：正向/反向命中与切点换算、回文酶的成对切点、环状跨首尾位点、切点取模
去重、跨首尾片段顺序与边界、端点切点不产生空片段、识别序列长于序列、多酶去重、
各类非法输入报错，以及 400 组线性/环状随机不变量测试（片段长度之和恒等于 N、
片段按原链顺序可拼回原序列等）。

`/design-panel` 测试额外覆盖：

- **联合酶切不可由单酶结果合并替代**：两样本在每个单酶下谱均相同，唯有联合酶切可区分；
- 三样本两两可区分、最少酶数与候选酶下标字典序择优（含候选酶重复定义）；
- 容差比对（片段数不同再大容差也不可混淆、逐位置长度差判定）；
- 环状多样本、无切点酶、无解 `no_solution`；
- 样本数量/重名/空名/序列长度上限、候选酶数量、`length_tolerance` 类型与负数等 400。

## 项目结构

```
.
├── app/
│   ├── __init__.py
│   ├── digest.py          # 酶切核心逻辑（纯函数，可独立测试）
│   ├── design.py          # /design-panel 组合搜索（复用 digest 规则）
│   └── main.py            # FastAPI 入口与错误处理
├── tests/
│   ├── test_digest.py     # /digest 纯标准库测试脚本
│   └── test_design.py     # /design-panel 纯标准库测试脚本
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
└── README.md
```
