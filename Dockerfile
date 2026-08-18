# NESTbench25 — full environment: C++ engines, Python (PyTorch) engines, and the test suite.
#
#   docker build -t nestbench25 .
#   docker run --rm nestbench25                                  # run the quick test suite
#   docker run --rm -w /nestbench/erasure/cpp_code nestbench25 ./sim -p input_control_parameters_learned.dat
#   docker run --rm -w /nestbench/abp/python_code nestbench25 python engine_abp.py -n 10000
#
# For a small C++-only image, see Dockerfile.cpp.
#
# BASE_IMAGE can be overridden (e.g. with a mirror of python:3.12-slim, or a base
# augmented with a corporate CA certificate on networks with TLS-intercepting proxies):
#   docker build --build-arg BASE_IMAGE=<your-base> -t nestbench25 .

ARG BASE_IMAGE=python:3.12-slim
FROM ${BASE_IMAGE}

# g++/make for the C++ engines; ffmpeg only for movie generation (-v)
RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ make ffmpeg \
    && rm -rf /var/lib/apt/lists/*

# The defaults fetch CPU-only PyTorch, which keeps the image small (the engines fall
# back to CPU automatically). On networks with a package-repository mirror or a
# TLS-intercepting proxy, override at build time, e.g.:
#   docker build --build-arg PIP_INDEX_URL=https://mirror.example.com/pypi/simple \
#                --build-arg PIP_EXTRA_INDEX_URL= \
#                --build-arg PIP_OPTS="--trusted-host mirror.example.com" -t nestbench25 .
ARG PIP_INDEX_URL=https://download.pytorch.org/whl/cpu
ARG PIP_EXTRA_INDEX_URL=https://pypi.org/simple
ARG PIP_OPTS=""

RUN pip install --no-cache-dir $PIP_OPTS \
    --index-url $PIP_INDEX_URL \
    ${PIP_EXTRA_INDEX_URL:+--extra-index-url $PIP_EXTRA_INDEX_URL} \
    torch numpy matplotlib pytest

WORKDIR /nestbench
COPY . .

# build every C++ engine
RUN set -e; for d in trap_overdamped trap_underdamped erasure ising abp; do \
        make -C $d/cpp_code standalone; \
    done

CMD ["pytest", "tests/", "-v"]
