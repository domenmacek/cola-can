import time
import os
import shutil
from filecmp import dircmp
import pickle
from json_tricks import dump, load
import numpy as np

def save_to_json(a, outpath):

    # Save to file
    with open(outpath, "w") as f:
        dump(a, f, indent=4)

    return


def load_from_json(path):

    # Save to file
    with open(path, "r") as f:
        a = load(f)

    return a


def makeDirIfNotExists(filePath):
    '''if no directory for files create it'''
    if not os.path.exists(os.path.dirname(filePath)):
        try:
            os.makedirs(os.path.dirname(filePath))
        except OSError as exc: # Guard against race condition
            if exc.errno != errno.EEXIST:
                raise
