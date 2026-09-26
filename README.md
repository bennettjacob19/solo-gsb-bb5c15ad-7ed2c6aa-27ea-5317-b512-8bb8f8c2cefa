# DNA 限制酶切模拟 API

基于 Python + FastAPI 实现的 DNA 限制性内切酶酶切模拟服务。输入一条仅含 A/C/G/T 的 DNA
序列、其拓扑结构（线性 / 环状）以及若干酶的识别序列与切点偏移，服务会：

- 同时搜索**正向识别序列**和其**反向互补序列**（正反两条链上的命中）；
- 按切点把序列切成片段，返回每个片段的**序列、长度和边界**；
- 环状序列识别**跨越首尾连接处**的位点，切点对序列长度**取模并去重**；
- 线性序列对越界切点**报错拒绝**；
- 无任何切点时返回**完整原序列**作为唯一片段；
- 对非法碱基、空序列、偏移越界等输入返回 **400** 并写明中文原因。

## 单条命令启动

需要 Docker（含 Compose 插件 v2）。在项目根目录执行：

```bash
docker compose up --build -d
```

服务启动后监听：

| 用途 | 地址 |
| --- | --- |
| 接口文档（Swagger UI，可直接在线调试） | http://localhost:8000/docs |
| OpenAPI Schema | http://localhost:8000/openapi.json |
| 健康检查 | http://localhost:8000/health |
| 酶切接口 | `POST http://localhost:8000/digest` |

停止服务：

```bash
docker compose down
```

## 必要配置

默认无需任何配置即可运行。唯一常用项是宿主机端口（默认 `8000`），有两种改法：

```bash
# 方式一：环境变量（docker compose 自动读取同目录 .env）
echo "HOST_PORT=9000" > .env
docker compose up --build -d

# 方式二：命令行临时指定
HOST_PORT=9000 docker compose up --build -d
```

容器内固定监听 `0.0.0.0:8000`（由 Dockerfile / Compose 中的 `HOST`、`PORT` 控制，
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

## 错误响应

所有输入类错误返回 HTTP `400`，响应体形如 `{"error": "中文原因说明"}`，例如：

- `序列含非法碱基 ['N']；仅允许 A/C/G/T（大小写均可）`
- `序列不能为空（或仅含空白字符）`
- `酶 'EcoRI' 的切点偏移 7 越界；识别长度为 6，合法范围为 0..6（含两端）`
- `拓扑类型非法：'plasmid'；仅支持 linear（线性）或 circular（环状）`
- `酶 'X' 的识别序列含非法碱基 ['Y']；仅允许 A/C/G/T`
- `酶 'X' 在线性序列位置 ... 产生越界切点 ...（合法范围 0..N）`

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
```

覆盖内容包括：正向/反向命中与切点换算、回文酶的成对切点、环状跨首尾位点、切点取模
去重、跨首尾片段顺序与边界、端点切点不产生空片段、识别序列长于序列、多酶去重、
各类非法输入报错，以及 400 组线性/环状随机不变量测试（片段长度之和恒等于 N、
片段按原链顺序可拼回原序列等）。

## 项目结构

```
.
├── app/
│   ├── __init__.py
│   ├── digest.py          # 酶切核心逻辑（纯函数，可独立测试）
│   └── main.py            # FastAPI 入口与错误处理
├── tests/
│   └── test_digest.py     # 纯标准库测试脚本
├── Dockerfile
├── docker-compose.yml
├── requirements.txt
└── README.md
```
