FROM julia:1.6.2
RUN cat /etc/resolv.conf
RUN  apt-get update -y && \
     apt-get upgrade -y && \
     apt-get dist-upgrade -y && \
     apt-get -y autoremove && \
     apt-get clean

RUN apt-get install -y unzip libghc-setenv-dev git build-essential libnlopt0 liblapack-dev libblas-dev openmpi-bin libopenmpi-dev gfortran cmake python
RUN apt install -y htop

RUN apt-get update && apt-get install -y openssh-server
RUN sed -i 's/PermitRootLogin prohibit-password/PermitRootLogin yes/' /etc/ssh/sshd_config

#RUN echo "sshuser:root" | 123456
ENTRYPOINT service ssh start && bash
COPY id_rsa.pub /home/root/.ssh/authorized_keys
RUN chown -R root:root /home/root/.ssh
RUN chmod 600 /home/root/.ssh/authorized_keys


RUN julia -e 'using Pkg; Pkg.update()'

EXPOSE 8888
WORKDIR /srv
CMD julia
