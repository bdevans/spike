# https://hub.docker.com/_/gcc/
# https://github.com/bdevans/spike

# docker build -t spike .
# docker run -it --rm --name spike spike

FROM gcc:latest
LABEL maintainer="Ben Evans <ben.d.evans@gmail.com>"
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
                    git \
                    nano \
                    libgsl-dev && \
                    apt-get clean && \
    apt-get autoremove && \
    rm -rf /var/lib/apt/lists/*

COPY . /usr/src/spike
WORKDIR /usr/src/spike
#RUN dpkg -L libgsl-dev
#RUN pkg-config --libs gsl
#RUN gsl-config --prefix

ENV GSL_INCLUDE=/usr/include GSL_LIB=/usr/lib/x86_64-linux-gnu \
    LD_LIBRARY_PATH=$GSL_LIB:$LD_LIBRARY_PATH

WORKDIR /usr/src/spike/source
#RUN make -j "$(nproc)"
RUN gcc -fopenmp -lm -ldl -lgomp -lgsl -lgslcblas -I$GSL_INCLUDE -L$GSL_LIB \
    -O3 -Wall main.c spike.c read_parameters.c array_utils.c utils.c -o Spike
RUN mv Spike ../
WORKDIR /usr/src/spike

ENTRYPOINT ["./Spike"]
# Specify default arguments (in this case no arguments prints the help)
#CMD ["-f defaults.m"]
#ENTRYPOINT ["/bin/bash"]
