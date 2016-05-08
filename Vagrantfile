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

    config.vm.provider :aws do |aws, override|

        aws.access_key_id = "#{ENV['AWS_ACCESS_KEY']}"
        aws.secret_access_key = "#{ENV['AWS_SECRET_KEY']}"
        aws.keypair_name = "#{ENV['AWS_SSH_KEY']}"
        aws.ami = "emi-dc40a7ee"
        aws.instance_ready_timeout = 300
        aws.instance_type = "m1.medium"
        aws.tags = {
            "Name" => "VagrantPowerGrid",
        }
        aws.security_groups = ["Default"]
        aws.region = ""
        aws.endpoint = ""
        override.vm.box = "dummy"
        override.vm.box_url = "https://github.com/mitchellh/vagrant-aws/raw/master/dummy.box"
        override.ssh.username = "ubuntu"
        # override.ssh.private_key_path = "#{ENV['AWS_SSH_KEY']}"
        override.vm.synced_folder '.', '/vagrant', disabled: true
        # override.sync.host_folder = ""  #relative to the folder your Vagrantfile is in
        # override.sync.guest_folder = "/vagrant" # relative to the vagrant home folder -> /home/vagrant
        # override.sshfs.enabled = false
        # override.sshfs.use_ssh_key = true
        # override.sshfs.mount_on_guest = true
        # override.sshfs.paths = { "" => "/vagrant" }
        # override.sshfs.host_addr = "#{ENV['HOSTNAME']}"
        # override.ssh.private_key_path = "~/.ssh/id_rsa"
        override.ssh.forward_agent = true
        override.ssh.insert_key = false
        override.ssh.private_key_path = ["#{ENV['HOME']}/#{ENV['AWS_SSH_KEY']}.pem", '~/.ssh/id_rsa']

    end

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

    config.vm.provider :aws do |aws, override|
      override.vm.provision "shell" do |s|
        s.inline = <<-SHELL
          sudo mkdir -p /vagrant
	  sudo chmod a+rwx /vagrant
	  sudo chmod a+r /etc/fuse.conf
	  sudo echo 'user_allow_other' >> /etc/fuse.conf
          sshfs -o StrictHostKeyChecking=no -o allow_other \
           -o reconnect \
           -o ServerAliveInterval=45 \
           -o ServerAliveCountMax=2 \
           -o ssh_command='autossh -M 0' $1@$2:$3 /vagrant
          SHELL
          s.args = ["#{ENV['USER']}","#{ENV['HOSTNAME']}","#{vagrant_root}"]
        end
    end
end
