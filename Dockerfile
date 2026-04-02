FROM python:3.10-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
    git \
    build-essential \
    procps \
  && rm -rf /var/lib/apt/lists/*

ARG PORECHOP_GIT_URL="https://github.com/trangdo-hsc/Porechop.git"
ARG PORECHOP_REF="e0e2f34"

RUN pip install --no-cache-dir "git+${PORECHOP_GIT_URL}@${PORECHOP_REF}"

ENTRYPOINT ["porechop"]
