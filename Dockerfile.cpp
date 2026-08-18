# NESTbench25 — minimal C++-only image (no Python/PyTorch; visualization movies unavailable).
#
#   docker build -f Dockerfile.cpp -t nestbench25-cpp .
#   docker run --rm -w /nestbench/trap_overdamped/cpp_code nestbench25-cpp ./sim -n 1000 -s 4
#
# Engines write results to report_answer.dat inside the container; mount a volume
# or add `cat report_answer.dat` to retrieve them, e.g.:
#   docker run --rm -w /nestbench/ising/cpp_code nestbench25-cpp sh -c "./sim -n 100 -s 10 && cat report_answer.dat"

FROM debian:stable-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends g++ make \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /nestbench
COPY . .

RUN set -e; for d in trap_overdamped trap_underdamped erasure ising abp; do \
        make -C $d/cpp_code standalone; \
    done

CMD ["/bin/bash"]
