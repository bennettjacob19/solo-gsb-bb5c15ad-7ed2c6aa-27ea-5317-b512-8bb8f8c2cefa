"""FastAPI 入口：DNA 限制酶切模拟 API。"""

from __future__ import annotations

import os
from typing import Any

from fastapi import FastAPI, Request
from fastapi.responses import JSONResponse, RedirectResponse

from .digest import DigestError, run_design_panel, run_digest

app = FastAPI(
    title="DNA 限制酶切模拟 API",
    description=(
        "输入仅含 A/C/G/T 的序列、线性/环状拓扑与若干酶的识别序列和切点偏移，"
        "按切点返回原链顺序的片段。同时搜索正向识别序列与其反向互补序列。"
    ),
    version="1.0.0",
)


@app.exception_handler(DigestError)
async def digest_error_handler(_request: Request, exc: DigestError) -> JSONResponse:
    return JSONResponse(status_code=400, content={"error": exc.reason})


@app.get("/", include_in_schema=False)
async def root() -> RedirectResponse:
    return RedirectResponse(url="/docs")


@app.get("/health")
async def health() -> dict[str, str]:
    return {"status": "ok"}


@app.post(
    "/digest",
    summary="执行限制酶切模拟",
    responses={
        400: {"description": "非法碱基、空序列、偏移越界或线性序列切点越界"},
    },
)
async def digest(request: Request) -> dict[str, Any]:
    try:
        payload = await request.json()
    except Exception:
        raise DigestError("请求体必须是合法的 JSON 对象")
    return run_digest(payload)


@app.post(
    "/design-panel",
    summary="为多个 DNA 样本寻找可相互鉴别的最少酶组合",
    responses={
        400: {"description": "样本/酶/拓扑/容差非法，或线性序列切点越界"},
    },
)
async def design_panel(request: Request) -> dict[str, Any]:
    try:
        payload = await request.json()
    except Exception:
        raise DigestError("请求体必须是合法的 JSON 对象")
    return run_design_panel(payload)


if __name__ == "__main__":
    import uvicorn

    uvicorn.run(
        "app.main:app",
        host=os.environ.get("HOST", "0.0.0.0"),
        port=int(os.environ.get("PORT", "8000")),
        reload=False,
    )
