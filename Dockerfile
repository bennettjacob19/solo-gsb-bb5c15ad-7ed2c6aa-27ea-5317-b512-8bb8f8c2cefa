FROM python:3.11-slim

ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1

WORKDIR /app

# 先单独拷贝依赖清单，利用 Docker 层缓存
COPY requirements.txt .
RUN pip install --no-cache-dir -r requirements.txt

COPY app ./app

# 构建期初始化自检：确认应用模块可导入、/digest 与 /design-panel 均已注册
RUN python -c "from app.main import app; \
paths = sorted(getattr(r, 'path', None) for r in app.routes); \
assert '/digest' in paths and '/design-panel' in paths, paths; \
print('init check ok:', paths)"

EXPOSE 8000

HEALTHCHECK --interval=30s --timeout=5s --start-period=5s --retries=3 \
    CMD python -c "import json,urllib.request,sys; \
sys.exit(0 if json.load(urllib.request.urlopen('http://127.0.0.1:8000/health'))['status']=='ok' else 1)"

CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]
