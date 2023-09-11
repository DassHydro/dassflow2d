FROM debian 

RUN apt-get update && DEBIAN_FRONTEND=noninteractive TZ=Europe/UTC apt-get install -y  \
     tzdata  \
     default-jdk  \
     git  \
     mpich  \
     python3  \
     python3-pip \
     python3-venv \
     wget \
     && rm -rf /var/lib/apt/lists/*

# Install f90wrap
RUN python3 -m venv /root/venv \
    && /root/venv/bin/pip3 install f90wrap
ENV PATH=/root/venv/bin:$PATH

COPY . /opt/dassflow/dassflow2d/

# Install tapenade
RUN wget https://tapenade.gitlabpages.inria.fr/tapenade/distrib/tapenade_3.16.tar \
    && tar xvfz tapenade_3.16.tar \
    && mv tapenade_3.16 /opt/dassflow/tapenade/ \
    && rm tapenade_3.16.tar

ENV PATH=$PATH:/opt/dassflow/tapenade/bin

# Compile project
RUN cd /opt/dassflow/dassflow2d/code/ \
    && rm -r ./bin_A/* && cp -r ../cases/tuto_case/1_lake-at-rest/bin_A/* ./bin_A \
    && make install

CMD ["/bin/bash"]

