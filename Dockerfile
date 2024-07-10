# Use the ubuntu:16.04 image for linux/386 platform
FROM --platform=linux/386 ubuntu:16.04

# Update package lists and install gfortran
# (perhaps also libxt6 libccp4-dev libgpp4f-dev
RUN apt-get update && apt-get install -y \
    libgfortran3 \
    command-not-found \
    x11-apps \
    wget \
    less \
    gzip \
    git \
    python \
    vim \
    libtool \
    sharutils \
    libxerces-c-dev \
    libmotif-dev \
    gcc \
    gfortran \
    m4 \
    curl \
    libxm4 && \
    apt autoremove --yes &&\
    apt-get clean && rm -rf /var/lib/apt/lists/*

RUN apt update && \
    apt install libcbf-dev/xenial

# RUN wget http://downloads.sf.net/cbflib/CBFlib-0.9.6.tar.gz && \
#     tar xvf CBFlib-0.9.6.tar.gz && \
#     cd CBFlib-0.9.6 && \
#     make tests

# RUN git clone https://github.com/aaronfinke/d-star-trek.git && \
#     cd d-star-trek && \
#     git checkout modernize && \
#     make && \
#     make install


ENTRYPOINT ["/bin/bash", "-c", "exec /bin/bash"]

