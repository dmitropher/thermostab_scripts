#!/usr/bin/env python
import os,optparse
from rosetta import *
from toolbox import *
import rosetta.core.init

def makeSymm (pos_obj, symm_file, outname):
    core.pose.symmetry.make_symmetric_pose(pos_obj,symm_file)
    pos_obj.dump_pdb(outname)

def main (argv):
    
    parser = optparse.OptionParser(usage="make an oligo using a pdb and .symm file")\
    
    parser.add_option('--pdb', dest='pdb_path' , help = "the path to the pdb")
    
    parser.add_option('--symm',dest='symm',help = 'path to the symm file')

    (options,args) = parser.parse_args() 

    pdb_filename   = options.pdb_path
    symm           = options.symm
    pos_obj        = pose_from_pdb(pdb_filename)


    name           = pdb_filename.split('/')
    name           = name[-1]
    name           = name.split('.')
    name           = name[0]
    
    outname        =  name + "_oligo.pdb"
    
    makeSymm(pos_obj,symm, outname)

if __name__ == "__main__":
    rosetta.init ( extra_options='-ignore_zero_occupancy false @/home/ddz6/Rosetta/flags/subs_YC.flags -extra_res_fa /home/ddz6/Rosetta/params/TYZ.params /home/ddz6/Rosetta/params/CYZ.params  -run:preserve_header')
    main(sys.argv[1:])

