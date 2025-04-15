FROM python:3.12-bookworm

# Set environment variables to make apt-get non-interactive
ENV DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

# The installer requires curl (and certificates) to download the release archive
RUN apt-get update && apt-get install -y --no-install-recommends curl ca-certificates
RUN apt-get install libexpat1 libexpat1-dev ffmpeg libsm6 libxext6 -y

# Download the latest installer
ADD https://astral.sh/uv/install.sh /uv-installer.sh

# Run the installer then remove it
RUN sh /uv-installer.sh && rm /uv-installer.sh

# Ensure the installed binary is on the `PATH`
ENV PATH="/root/.local/bin/:$PATH"

RUN mkdir -p /boston_twin
WORKDIR /boston_twin
RUN uv init

# COPY src .
# COPY generate_dataset.py .
# COPY pyproject.toml .
CMD ["tail", "-f", "/dev/null"]
# CMD ["uv", "run", "generate_dataset.py"]
# CMD ["jupyter", "lab", "--port=8888", "--no-browser", "--NotebookApp.token=''", "--ip=0.0.0.0", "--allow-root", "--NotebookApp.allow_origin='*'"]
