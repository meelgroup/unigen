FROM ubuntu:20.04 as builder

LABEL maintainer="Mate Soos"
LABEL version="1.0"
LABEL Description="Approxmc"

# get curl, etc
RUN apt-get update && apt-get install --no-install-recommends -y software-properties-common
# RUN add-apt-repository -y ppa:ubuntu-toolchain-r/test
# RUN apt-get update
RUN apt-get install --no-install-recommends -y libboost-program-options-dev gcc g++ make cmake zlib1g-dev wget make libgmp-dev git libmpfr-dev

# get M4RI
WORKDIR /
RUN wget https://bitbucket.org/malb/m4ri/downloads/m4ri-20200125.tar.gz
RUN tar -xvf m4ri-20200125.tar.gz
WORKDIR m4ri-20200125
RUN ./configure
RUN make \
    && make install

WORKDIR /
RUN git clone https://github.com/meelgroup/cadical
WORKDIR /cadical
RUN git checkout mate-only-libraries-1.8.0
RUN ./configure
RUN make


WORKDIR /
RUN git clone https://github.com/meelgroup/cadiback
WORKDIR /cadiback
RUN git checkout mate
RUN ./configure
RUN make

WORKDIR /
RUN git clone https://github.com/msoos/cryptominisat
WORKDIR /cryptominisat
RUN mkdir build
WORKDIR /cryptominisat/build
RUN cmake -DSTATICCOMPILE=ON ..
RUN make -j6 \
    && make install

#build SBVA
WORKDIR /
RUN git clone https://github.com/meelgroup/SBVA.git
WORKDIR /SBVA
RUN mkdir build
WORKDIR /SBVA/build
RUN cmake -DSTATICCOMPILE=ON ..
RUN make -j6 \
    && make install

#build arjun
WORKDIR /
RUN git clone https://github.com/meelgroup/arjun
WORKDIR /arjun
RUN mkdir build
WORKDIR /arjun/build
RUN cmake -DSTATICCOMPILE=ON ..
RUN make -j6 \
    && make install


# build approxmc
WORKDIR /
RUN git clone https://github.com/meelgroup/approxmc
WORKDIR /approxmc
RUN mkdir build
WORKDIR /approxmc/build
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

