#! /usr/bin/env python3
import os, glob, re

# Put this in the top directory of pyJose.
dir = "./optspecextr_source/"
ext = ".pro"

def readfunc(file):
    # Reads an IDL function and checks for all of the variables in that function.
    with open(file) as f:
        data = f.read().splitlines()
    
    for line in range(len(data)):
        # skip comments
        if data[line].__contains__(";"):
           continue
        
        # TODO regex match for only capitalized function definitions
        if re.match("function\\b", data[line]):
            while data[line].__contains__("$"):
                print(data[line])
                line += 1
            print(data[line], "\n")


for file in glob.glob(dir + "*" + ext):
    readfunc(file)