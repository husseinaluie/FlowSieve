import subprocess, os

read_the_docs_build = os.environ.get('READTHEDOCS', None) == 'True'

# Extract version number from version file
with open('VERSION') as f:
    contents = f.readlines()
contents = [ x.strip() for x in contents ]
version_nums = [ x[14:] for x in contents ]
VERSION = version_nums[0] + '.' + version_nums[1] + '.' + version_nums[2]

if read_the_docs_build:
    subprocess.call('DOXY_VERSION={0} doxygen Doxyfile'.format(VERSION), shell=True)
