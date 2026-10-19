# Dockerfile para Illumina Manta v1.6.0 SV Calling Container
# Compatible con GCP Cloud Batch y ejecución Docker local x86_64

FROM python:3.10-slim-bullseye

LABEL maintainer="Bio Stelith / HSelith Team"
LABEL description="Contenedor optimizado para Illumina Manta v1.6.0 y filtrado de variantes MVA"

ENV DEBIAN_FRONTEND=noninteractive
ENV MANTA_VERSION=1.6.0
ENV MANTA_HOME=/opt/manta

# 1. Dependencias del sistema y herramientas genómicas
RUN apt-get update && apt-get install -y --no-install-recommends \
    bzip2 \
    ca-certificates \
    curl \
    gcc \
    g++ \
    libbz2-dev \
    liblzma-dev \
    make \
    samtools \
    tabix \
    zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# 2. Descarga e instalación de binarios oficiales de Manta v1.6.0
RUN mkdir -p ${MANTA_HOME} && \
    curl -sSL "https://github.com/Illumina/manta/releases/download/v${MANTA_VERSION}/manta-${MANTA_VERSION}.centos6_x86_64.tar.bz2" \
    | tar -xj -C ${MANTA_HOME} --strip-components=1

ENV PATH="${MANTA_HOME}/bin:${PATH}"

# 3. Directorio de la aplicación y dependencias Python
WORKDIR /app

COPY requirements.txt /app/
RUN pip install --no-cache-dir -r requirements.txt

# 4. Copia del código fuente del módulo y shared
COPY src/ /app/src/

ENV PYTHONPATH=/app/src:/app

# 5. Punto de entrada
ENTRYPOINT ["python", "-m", "src.modules.sv_calling.call_sv"]

