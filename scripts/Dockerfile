FROM ubuntu:latest
  
RUN apt update
RUN DEBIAN_FRONTEND=noninteractive apt install -y g++ cmake python3 python3-pip python3-venv pkg-config libboost-dev libopenmpi-dev libsuitesparse-dev gmsh git git-lfs

ENV VIRTUAL_ENV=~/dune-venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN CXXFLAGS='-Wno-narrowing' pip install dune-mmesh matplotlib

RUN git clone https://github.com/samuelburbulla/dune-mmesh.git

WORKDIR dune-mmesh/doc/examples
RUN python coupling.py

RUN echo 'export PATH="$VIRTUAL_ENV/bin:$PATH"' >> ~/.bashrc
