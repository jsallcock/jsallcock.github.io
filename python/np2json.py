import os
import json
import numpy as np

"""
testing a quick way to export numpy arrays into importable javascript. This is probably horrible
"""

# FPATH_OUT = os.path.join(os.path.dirname(os.path.realpath(__file__)), 'example.json')
FPATH_OUT = '/Users/jsallcock/repo/jsallcock.github.io/public/js/example.js'
try:
    os.remove(FPATH_OUT)
except FileNotFoundError:
    pass
arr1 = np.arange(36.).reshape(3, 12)
arr1[0, 0] = np.nan
data = {
    'var1': arr1.tolist(),
    'var2': np.arange(36).reshape(9, 4).tolist(),
}
with open(FPATH_OUT, 'w+') as f:
    f.write('export default\n' + json.dumps(data, indent=4))
