description = """a module that holds records of known PSDs"""
author = "Reed Essick"

#-------------------------------------------------

import os
import numpy as np

#-------------------------------------------------

def path2dict( path ):
    freqs, asd = np.loadtxt(path).transpose()
    return {'freqs':freqs, 'vals':asd**2}

#-------------------------------------------------

dirname = os.path.join(os.path.dirname(__file__), '../etc')

### LIGOs
aLIGO = path2dict(os.path.join(dirname, 'aLIGO.txt'))
aLIGO_O1 = path2dict(os.path.join(dirname, 'O1.txt'))
aLIGO_O2 = path2dict(os.path.join(dirname, 'O2_proj.txt'))
aLIGO_O3 = path2dict(os.path.join(dirname, 'O3_proj.txt'))
aLIGO_design = path2dict(os.path.join(dirname, 'aLIGO_design.txt'))

### Aplus
aPlus = path2dict(os.path.join(dirname, 'Aplus.txt'))
aPlus_sqzonly = path2dict(os.path.join(dirname, 'Aplus_sqzonly.txt'))

### Virgos
aVirgo = path2dict(os.path.join(dirname, 'AdVirgo.txt'))
aVirgo_sqz = path2dict(os.path.join(dirname, 'AdVirgo_sqz.txt'))
aVirgo_wb = path2dict(os.path.join(dirname, 'AdVirgo_wb.txt'))

### CE
CE = path2dict(os.path.join(dirname, 'CE.txt'))
CE_wb = path2dict(os.path.join(dirname, 'CE_wb.txt'))

### ET
ET = path2dict(os.path.join(dirname, 'ET_D.txt'))

### Voyager
Voyager = path2dict(os.path.join(dirname, 'Voyager.txt'))
