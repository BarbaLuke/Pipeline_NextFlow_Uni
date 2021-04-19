FROM ubuntu:latest

RUN apt update -y
RUN apt install -y trimmomatic
RUN apt install -y fastqc

RUN apt upgrade -y
