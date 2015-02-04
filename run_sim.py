#!/usr/bin/python3
import subprocess
import numpy as np

for ii in np.arange(0.2,2,0.2):
    cmd = './simulator  -k 5  -g 400 -n 4000 -t 5000 -d 0.2 -l'+' elastic%02u'%int(ii*10) +' -s' + ' %1.2f;' % ii
    print(cmd)
    p = subprocess.Popen(cmd, shell=True)
    p.wait()
