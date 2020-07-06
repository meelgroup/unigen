FROM ubuntu:16.04 as builder

LABEL maintainer="Mate Soos"
LABEL version="1.0"
LABEL Description="Approxmc"

# get curl, etc
RUN apt-get update && apt-get install --no-install-recommends -y software-properties-common
# RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
# RUN apt-get update
RUN apt-get install --no-install-recommends -y libboost-program-options-dev gcc g++ make cmake zlib1g-dev wget make libgmp-dev

# get M4RI
WORKDIR /
RUN wget https://bitbucket.org/malb/m4ri/downloads/m4ri-20200125.tar.gz
RUN tar -xvf m4ri-20200125.tar.gz
WORKDIR m4ri-20200125
RUN ./configure
RUN make \
    && make install

# build CMS
WORKDIR /
RUN wget https://github.com/msoos/cryptominisat/archive/5.8.0.tar.gz
RUN tar -xvf 5.8.0.tar.gz
WORKDIR /cryptominisat-5.8.0
RUN mkdir build
WORKDIR /cryptominisat-5.8.0/build
RUN cmake -DSTATICCOMPILE=ON ..
RUN make -j6 \
    && make install

# build approxmc
WORKDIR /
RUN wget https://github.com/meelgroup/approxmc/archive/4.0.1.tar.gz
RUN tar -xvf 4.0.1.tar.gz
WORKDIR /approxmc-4.0.1
RUN mkdir build
WORKDIR /approxmc-4.0.1/build
RUN cmake -DSTATICCOMPILE=ON ..
RUN make -j6 \
    && make install

# build unigen
USER root
COPY . /home/solver/unigen
WORKDIR /home/solver/unigen
RUN mkdir build
WORKDIR /home/solver/unigen/build
RUN cmake -DSTATICCOMPILE=ON ..
RUN make -j6 \
    && make install

# set up for running
FROM alpine:latest
COPY --from=builder /usr/local/bin/unigen /usr/local/bin/
ENTRYPOINT ["/usr/local/bin/unigen"]

# --------------------
# HOW TO USE
# --------------------
# on file through STDIN:
#    zcat mizh-md5-47-3.cnf.gz | docker run --rm -i -a stdin -a stdout msoos/approxmc

# on a file:
#    docker run --rm -v `pwd`/myfile.cnf.gz:/in msoos/approxmc in

# echo through STDIN:
#    echo "1 2 0" | docker run --rm -i -a stdin -a stdout msoos/approxmc

# hand-written CNF:
#    docker run --rm -ti -a stdin -a stdout msoos/approxmc

