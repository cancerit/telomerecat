FROM ubuntu:20.04 as builder
USER root

# Version of tools that are going to be installed.
# can specify  "hotfix/X.X.X", "feature/fixstuff" or "3.4.1"
ARG BRANCH_OR_TAG_PARABAM="hotfix/2.3.2"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
locales \
g++ \
make \
gcc \
pkg-config \
python3 python3-dev python3-pip python3-setuptools \
zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev \
curl
# zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev are for building pysam

ENV CGP_OPT /opt/wtsi-cgp
ENV PYTHONPATH $CGP_OPT/lib/python3.8/site-packages
RUN mkdir -p $PYTHONPATH

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8
ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8

RUN pip3 install wheel cython

RUN curl -sSL https://github.com/cancerit/parabam/archive/${BRANCH_OR_TAG_PARABAM}.tar.gz | tar zx \
&& cd parabam* \
&& python3 setup.py sdist \
&& python3 setup.py install --prefix=$CGP_OPT

# build the tools in this repo, separate to reduce build time on errors
COPY setup.py .
COPY telomerecat telomerecat
RUN python3 setup.py sdist
RUN python3 setup.py install --prefix=$CGP_OPT

FROM ubuntu:20.04

LABEL maintainer="cgphelp@sanger.ac.uk" \
      uk.ac.sanger.cgp="Cancer, Ageing and Somatic Mutation, Wellcome Trust Sanger Institute" \
      version="3.4.1" \
      description="telomerecat docker"

RUN apt-get -yq update
RUN apt-get install -yq --no-install-recommends \
apt-transport-https \
locales \
ca-certificates \
time \
unattended-upgrades \
python3 python3-pkg-resources \
zlib1g-dev libbz2-dev liblzma-dev libcurl4-gnutls-dev && \
unattended-upgrade -d -v && \
apt-get remove -yq unattended-upgrades && \
apt-get autoremove -yq

RUN locale-gen en_US.UTF-8
RUN update-locale LANG=en_US.UTF-8

ENV CGP_OPT /opt/wtsi-cgp
ENV PATH $CGP_OPT/bin:$CGP_OPT/bin:$PATH
ENV PYTHONPATH $CGP_OPT/lib/python3.8/site-packages
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
