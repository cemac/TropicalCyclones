#!/usr/bin/env python
'''
Script to check for all Python code modified as part of this branch
and give pylint report for this code and a summary report of scores.
'''
import subprocess
import os
import sys


def prelude():
    '''Set path to toolbox Python code.'''
    lib_dir = os.path.normpath(os.path.join(os.path.dirname(
        os.path.realpath(__file__)), "..", "lib"))
    if sys.path[0:1] != [lib_dir]:
        if lib_dir in sys.path:
            sys.path.remove(lib_dir)
        sys.path.insert(0, lib_dir)


prelude()


def pylint():
    '''
    Run pylint on all code modified as part of this branch.
    '''
    summary = []

    # Set up command to determine which .py files have been modified in branch
    command = 'fcm diff -b | grep "^Index: " | cut -c8- | grep ".py"'
    p1 = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE)
    for pyfile in p1.stdout:
        pyfile = pyfile.strip()
        print '*******************************************'
        print 'Pylint report for {0:s}'.format(pyfile)
        command = 'pylint {0:s}'.format(pyfile)
        p2 = subprocess.Popen(
            command, shell=True, stdout=subprocess.PIPE)

        # Also get score for summary report
        for out in p2.stdout:
            print out
            if out[:27] == 'Your code has been rated at':
                summary.append(pyfile + ': ' + out.split()[6].split('/')[0])

    # Print summary
    print '*******************************************'
    print '*******************************************'
    print 'Summary:'
    for line in summary:
        print line


if __name__ == '__main__':
    pylint()
