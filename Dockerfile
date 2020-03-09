FROM ubuntu:18.04 as builder
USER root

# Version of tools that are going to be installed.
ENV PARABAM_VER '2.3.0'

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
g++ \
make \
gcc \
pkg-config \
python3 python3-dev python3-pip python3-setuptools \
zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev \
git
# zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev are for building pysam
# git is for parabam installation from source

ENV CGP_OPT /opt/wtsi-cgp
RUN mkdir $CGP_OPT
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.6/site-packages

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN pip3 install wheel cython

RUN pip3 install --install-option="--prefix=$CGP_OPT/python-lib" https://github.com/cancerit/parabam/releases/download/${PARABAM_VER}/parabam-${PARABAM_VER}.tar.gz

# build the tools in this repo, separate to reduce build time on errors
COPY . .
# git clone https://github.com/cancerit/telomerecat.git
# cd telomerecat
# git checkout feature/explore_to_py3
RUN python3 setup.py sdist
RUN pip3 install --install-option="--prefix=$CGP_OPT/python-lib" dist/$(ls -1 dist/)

FROM ubuntu:18.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="3.4.0" \
      description="telomerecat docker"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
ca-certificates \
time \
unattended-upgrades \
python3 \
zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV CGP_OPT /opt/wtsi-cgp
ENV PATH $CGP_OPT/bin:$CGP_OPT/python-lib/bin:$PATH
ENV PYTHONPATH $CGP_OPT/python-lib/lib/python3.6/site-packages
ENV LD_LIBRARY_PATH $OPT/lib
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN mkdir -p $CGP_OPT
COPY --from=builder $CGP_OPT $CGP_OPT

## USER CONFIGURATION
RUN adduser --disabled-password --gecos '' ubuntu && chsh -s /bin/bash && mkdir -p /home/ubuntu

USER    ubuntu
WORKDIR /home/ubuntu

CMD ["/bin/bash"]
