# FROM openmicroscopy/octave:latest
FROM schickling/octave:latest

RUN mkdir -p /home/octave/kihasa
WORKDIR /home/octave/kihasa

ADD . /home/octave/kihasa

CMD octave KIHASA_Main.m