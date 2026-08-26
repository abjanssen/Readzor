FROM python:3.11-slim
ARG VERSION=latest
LABEL maintainer="Axel B. Janssen" \
      version="${VERSION}" \
      description="Readzor"
ENV PYTHONUNBUFFERED=1 \
    PYTHONDONTWRITEBYTECODE=1
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir --upgrade pip \
    && VERSION_CLEAN=$(echo "${VERSION}" | sed 's/^v//') \
    && if [ "${VERSION}" = "latest" ] || [ "${VERSION}" = "main" ]; then \
         pip install --no-cache-dir readzor; \
       else \
         pip install --no-cache-dir readzor=="${VERSION_CLEAN}"; \
       fi
WORKDIR /data
ENTRYPOINT ["readzor"]
CMD ["--help"]