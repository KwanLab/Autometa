#!/bin/bash

# Steps to update dockerhub images

rm -rf /tmp/dockerdev
mkdir /tmp/dockerdev
cd /tmp/dockerdev
git clone git@github.com:KwanLab/Autometa.git
git switch main
git pull
cd /tmp/dockerdev/Autometa/docker/modules

# One or both of of these may take a while to build

docker build -t jasonkwan/autometa-nf-modules-get_genomes_for_mock:main -f get_genomes_for_mock.Dockerfile .
docker build -t jasonkwan/autometa-nf-modules-mock_data_reporter:main -f mock_data_reporter.Dockerfile .

docker push jasonkwan/autometa-nf-modules-get_genomes_for_mock:main
docker push jasonkwan/autometa-nf-modules-mock_data_reporter:main
