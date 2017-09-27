#!/usr/bin/env python

# Program to convert variables.cfg in the original BCVTB format (XML) to JSON format used by OBN-EnergyPlus
# Copyright (C) 2016 by Truong X. Nghiem (truong.nghiem@gmail.com)

from __future__ import print_function
import sys

# Print help
def print_help():
    print("Convert original BCVTB variable XML file to JSON file for OpenBuildNet.")
    print("{} <input file> <output file>".format(sys.argv[0]))

# Process input arguments
if len(sys.argv) < 3:
    print("ERROR: Not enough input arguments.")
    print_help()
    sys.exit(1)

    
import xml.parsers.expat
import json

thelist = []

# Where am I right now?
inlist = False
invariable = False
inep = False

# The current variable
cur_source = ""
cur_attrs = dict()


# Function to read BCVTB file in to a list of dictionary objects
def readBCVTB(f):
    try:
        infile = open(f, 'r')
    except IOError as e:
        print("ERROR: Could not open input file: ", e.strerror)
        sys.exit(2)
            
    def start_elem(name, attrs):
        # print("Start:", name, attrs)
        global inlist, invariable, inep, cur_source, cur_attrs
        if name.lower() == "variable" and inlist and not invariable:
            # Enter variable
            invariable = True
            if "source" in attrs:
                cur_source = attrs["source"]
            else:
                raise RuntimeError("Syntax error: A variable definition without source.")
        elif name.lower() == "energyplus" and invariable:
            # The definition of the variable
            inep = True
            cur_attrs = attrs
        elif name.lower() == "bcvtb-variables" and not inlist:
            # Enter list
            inlist = True
        else:
            # Syntax Error
            raise RuntimeError("Syntax error while processing input file.")

    def end_elem(name):
        # print("End:", name)
        global inlist, invariable, inep, cur_source, cur_attrs
        if name.lower() == "variable" and inlist and invariable and not inep:
            # Exit variable
            invariable = False
            cur_attrs["source"] = cur_source
            thelist.append(cur_attrs)  # Add the variable to the list
        elif name.lower() == "energyplus" and inep:
            # Exit variable definition
            inep = False
        elif name.lower() == "bcvtb-variables" and inlist and not invariable:
            # Exit list
            inlist = False
        else:
            # Syntax Error
            raise RuntimeError("Syntax error while processing input file.")

    p = xml.parsers.expat.ParserCreate()
    
    p.StartElementHandler = start_elem
    p.EndElementHandler = end_elem

    p.ParseFile(infile)
    infile.close()

    if inlist or invariable or inep:
        raise RuntimeError("Syntax error: XML file not closed properly.")

# Function to write JSON file from thelist
def writeJSON(f):
    # Go through the list and modify its attributes so that we can dump to JSON
    for idx, item in enumerate(thelist):
        if item["source"].lower() == "energyplus":
            # An output: name => object, type => name
            del item["source"]            
            item["object"] = item["name"]
            item["name"] = item["type"]
            item["type"] = "out"
        elif item["source"].lower() == "ptolemy":
            # An input: key name => key, key value => name
            del item["source"]
            keys = item.keys()
            if len(keys) != 1:
                raise RuntimeError("Syntax error: a Ptolemy variable must have exactly one key defined.")
            item["key"] = keys[0]
            item["name"] = item[keys[0]]
            item["type"] = "in"
            del item[keys[0]]
        else:
            raise RuntimeError("Syntax error: unrecognized variable source {}".format(item["source"]))

        
    try:
        outfile = open(f, 'w')
        json.dump(thelist, outfile, indent=2, separators=(',', ': '), sort_keys=True)
        outfile.close()
    except IOError as e:
        print("ERROR: Could not open or write to output file: ", e.strerror)
        sys.exit(3)


        
infilename = sys.argv[1]
outfilename = sys.argv[2]

print("Convert variable list from BCVTB input file {} to JSON output file {}...".format(infilename, outfilename))

readBCVTB(infilename)
writeJSON(outfilename)

print("Done!")
