# -*- mode: ruby -*-
# vi: set ft=ruby :

# All Vagrant configuration is done below. The "2" in Vagrant.configure
# configures the configuration version (we support older styles for
# backwards compatibility). Please don't change it unless you know what
# you're doing.

Vagrant.configure(2) do |config|
    # The most common configuration options are documented and commented below.
    # For a complete reference, please see the online documentation at
    # https://docs.vagrantup.com.

    # Every Vagrant development environment requires a box. You can search for
    # boxes at https://atlas.hashicorp.com/search.

    config.vm.box = "boxcutter/ubuntu1604"


    # Enable provisioning with a shell script. Additional provisioners such as
    # Puppet, Chef, Ansible, Salt, and Docker are also available. Please see the
    # documentation for more information about their specific syntax and use.
    vagrant_root = File.dirname(__FILE__)
    config.vm.provision "shell" do |s|
      s.privileged = false
      s.inline = <<-SHELL
        sudo apt-get update
        sudo apt-get install -y software-properties-common python-software-properties
        sudo apt-get update
        sudo add-apt-repository universe
        sudo apt-get update
        sudo apt-get install -y sshfs autossh
        #  Installing ismrmrd
        sudo apt-get install -y libhdf5-serial-dev h5utils cmake cmake-curses-gui libboost-all-dev doxygen git libfftw3-dev g++
        git clone https://github.com/ismrmrd/ismrmrd
        cd ismrmrd/
        mkdir build
        cd build
        cmake ../
        make
        sudo make install
        #  Installing the PowerGrid Dependencies
        sudo apt-get install -y libmatio-dev libopenblas-dev libxerces-c-dev libarmadillo-dev xsdcxx
        cd ~
        git clone http://github.com/mrfil/PowerGridTestData.git
      SHELL

    end

end
