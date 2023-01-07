#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  6 23:27:02 2023

@author: nawaz
"""

def readPath(paths_file, path_name):
    '''Read paths from a text file using a given path name'''
    with open(paths_file) as f: 
        lines = f.readlines()
        

        for line in lines: 
            if line.find(path_name) != -1: 
                line = line.strip()
                x = line.split(": ")
                return x[1]
        
        content = f.read()
        if path_name not in content: 
            print ('path does not exists')
                
# spath = readPath('paths.txt', 'path_val_signal:')
# print(spath)