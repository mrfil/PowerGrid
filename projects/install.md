---
layout: project
title: Quick Start
category: install
permalink: install/
excerpt: Getting Started
---

# Quick Start
{: .content-subhead }

To test out PowerGrid or to try modifying or creating reconstruction routines or objects, Vagrant provides an easy way to setup a virtual machine with PowerGrid and dependencies for CPU based reconstructions.

## Getting Started with Vagrant

Get the latest version of [Vagrant](http://vagrantup.com).

Clone the PowerGrid repository and Vagrant up

```shell
git clone https://github.com/mrfil/PowerGrid.git ./PG

cd PG

vagrant up

vagrant ssh
```

## Building PowerGrid in the Virtual machine

Run the following commands from the ssh session into the virtual machine.

```shell
cd /vagrant
cd /PG
mkdir build
cmake ../
make
```

## Running PowerGrid on a test dataset

Run the following commands from the ssh session into the virtual machine.

```shell
./PowerGridGnufft ~/PowerGridTestData/192_192_1_32coils/
```

## View images in MATLAB

Run the following command from a terminal on the host machine.

```shell
scp -P 2222 vagrant@127.0.0.1:/home/vagrant/PowerGridTestData/192_192_1_32coils/test_pwls.mat .
```

Open MATLAB on the host.

```MATLAB
load test_pwls.mat
temp = reshape(test_pwls,[192,192,1]);
im(temp) %Or another image viewer
```
