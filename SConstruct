import os
import socket
import getpass
import time

build_dir = time.strftime("%Y%m%d%H%M%S", time.gmtime()) + '_' + getpass.getuser() + '_' + socket.gethostname()

try:
 env = Environment(ENV = { 'PATH' : os.environ['PATH'], 'LD_LIBRARY_PATH' : os.environ['LD_LIBRARY_PATH']})
except:
 env = Environment(ENV = { 'PATH' : os.environ['PATH'],})

#modify cpu architecture before compiling
#modify library paths before compiling

env.Append(CCFLAGS=['-O3', '-std=c++11', '-march=native'])
env.Append(CPPPATH=['.', '/opt/eigen-3.2.3', '/usr/include/libxml2/', '/usr/include', ])
env.Append(LIBPATH=['/usr/lib/x86_64-linux-gnu',])

env.Append(LINKFLASGS=['-lz', '-lxml2', '-pthread', '-lm', '-lhdf5_cpp', '-lhdf5'])
env.Append(LIBS = ['libxml2','pthread','z','m', 'hdf5', 'hdf5_cpp'])

Export('env')

trj2hdf5 = env.SConscript("trj2hdf5/SConscript", variant_dir="built/trj2hdf5", duplicate=0)
frameout  = env.SConscript("frameout/SConscript", variant_dir="built/frameout", duplicate=0)
rmsd      = env.SConscript("rmsd/SConscript", variant_dir="built/rmsd", duplicate=0)
rmsf      = env.SConscript("rmsf/SConscript", variant_dir="built/rmsf", duplicate=0)
checkbox  = env.SConscript("checkbox/SConscript", variant_dir="built/checkbox", duplicate=0)
hbond  = env.SConscript("hbond/SConscript", variant_dir="built/hbond", duplicate=0)

#note: programs order by release data
#lower on this list means more recent
env.Install(dir='built/bin', source=trj2hdf5)
env.Install(dir='built/bin', source=frameout)
env.Install(dir='built/bin', source=rmsd)
env.Install(dir='built/bin', source=checkbox)
env.Install(dir='built/bin', source=rmsf)
env.Install(dir='built/bin', source=hbond)
