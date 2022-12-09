FROM ubuntu:latest

RUN apt update
ENV DEBIAN_FRONTEND=noninteractive
RUN apt install -yq g++ cmake python3 python3-pip python3-venv pkg-config libboost-dev libopenmpi-dev libsuitesparse-dev gmsh git git-lfs > /dev/null

RUN python3 -m venv venv

RUN venv/bin/pip install dune-common
RUN venv/bin/pip install dune-geometry
RUN venv/bin/pip install dune-grid
RUN venv/bin/pip install dune-istl
RUN venv/bin/pip install dune-localfunctions
RUN venv/bin/pip install dune-alugrid
RUN venv/bin/pip install dune-fem
RUN venv/bin/pip install dune-mmesh
RUN venv/bin/pip install matplotlib
RUN venv/bin/pip install mpi4py

RUN git clone https://github.com/samuelburbulla/dune-mmesh.git
WORKDIR dune-mmesh/doc/examples
RUN . ../../../venv/bin/activate && python coupling.py

WORKDIR $HOME
RUN echo '. venv/bin/activate' >> ~/.bashrc
