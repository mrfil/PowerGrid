"""
HPC Base image
Contents:
  CUDA version 9.0
  PGI compilers version 18.10
"""
# pylint: disable=invalid-name, undefined-variable, used-before-assignment

# The PGI End-User License Agreement (https://www.pgroup.com/doc/LICENSE)
# must be accepted.
pgi_eula=False
if USERARG.get('pgi_eula_accept', False):
  pgi_eula=True
else:
  raise RuntimeError('PGI EULA not accepted. To accept, use "--userarg pgi_eula_accept=yes"\nSee PGI EULA at https://www.pgroup.com/doc/LICENSE')

# Choose between either Ubuntu 16.04 (default) or CentOS 7
# Add '--userarg centos=true' to the command line to select CentOS
devel_image = 'nvidia/cuda:9.0-devel-ubuntu16.04'
runtime_image = 'nvidia/cuda:9.0-runtime-ubuntu16.04'
if USERARG.get('centos', False):
    devel_image = 'nvidia/cuda:9.0-devel-centos7'
    runtime_image = 'nvidia/cuda:9.0-runtime-centos7'

######
# Devel stage
######

Stage0 += comment(__doc__, reformat=False)

Stage0 += baseimage(image=devel_image, _as='devel')

# PGI compilers
compiler = pgi(eula=pgi_eula, version='18.10')
Stage0 += compiler
