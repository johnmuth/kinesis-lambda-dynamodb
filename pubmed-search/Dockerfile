FROM amazonlinux:latest
ENV LANG=en_GB.UTF-8
RUN   yum -y update && \
      yum -y install \
        git             \
        gcc             \
        gcc-c++         \
        make            \
        zlib-devel      \
        bzip2-devel     \
        openssl-devel   \
        readline-devel  \
        ncurses-devel   \
        sqlite-devel && \
      curl --remote-name --progress https://www.python.org/ftp/python/3.6.1/Python-3.6.1.tgz && \
      echo '2d0fc9f3a5940707590e07f03ecb08b9 Python-3.6.1.tgz' | md5sum -c - && \
      tar zxf Python-3.6.1.tgz && \
      cd Python-3.6.1 && \
        ./configure --prefix=/usr/local --enable-shared && \
        make && \
        make altinstall && \
      echo '/usr/local/lib' > /etc/ld.so.conf.d/local-libs.conf && \
      ldconfig && \
      cd .. && \
      rm -rf Python-3.6.1.tgz Python-3.6.1

RUN easy_install-3.6 pip
RUN pip3 install virtualenv
ADD . /src
WORKDIR /src
RUN virtualenv -p /usr/local/bin/python3.6 /app
RUN source /app/bin/activate && pip3 install -r requirements.txt
