#!/usr/bin/env python
from rosetta import *
import math, csv, sys, optparse, os, ast, shutil
from toolbox import *
from transform import *
import rosetta.core.init

def main(argv):
    #command line options
    #
    
    parser = optparse.OptionParser(usage="\n\nCopy TYZ and CYX on all chains but A to chain A and symmeterize chain A. Requires an existing .symm file.")

    parser.add_option('--pdb-path', dest = 'pdb_path',
        help = 'The path to the pdbs listed in the outcsv \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('-o', dest = 'final_out_path',
        default = False,
        help = 'The path to where the output INPUT files should be stored. \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('--cst-path', dest = 'cst_path',
        default = './',
        help="Output cst files path: ( default = 'execute user's home dir' )")

    parser.add_option('--params-file-TYZ', dest = 'params_file_TYZ',
        default = '',
        help="Specify path to params file when mutating to ncAA.")
 
    parser.add_option('--params-file-CYX', dest = 'params_file_CYX',
        default = '',
        help="Specify path to params file when mutating to ncAA.")

    
    (options,args) = parser.parse_args()
    #Read pdb file from command line option
    pdb_filename = options.pdb_path
    pdb          = read_file(pdb_filename)
    pdb_obj      = Transform(pdb)
    res_to_check = "TYZ:CYZ" #This string contains the residues to
                             #+ search for in chain non-A

    
    #The lists to be populated with residues to be copied
    #+ over to chain A from the other chainsi
    #+ They are unified by index number, so bad input may break them
    pdb_resnum          = [] 
    resnum2pose         = []
    
    #this for loop uses res_to_check, a string which must contain all
    #+ residue types to check for in chains other than A
    #+ this loop also assumes all atoms belonging to a residue
    #+ will be sequential
    #+ it checks only if the previous addition was a duplicate before appending
    #+ not the entire list.


#BEGIN FOR LOOP
    for name in pdb_obj.xyz_names:  #Iterates through the entire pdb, and only saves
        split_name = name.split('_')
        if (split_name[2] in res_to_check and  #+ checks if a residue is CYZ or TYZ
            (len(pdb_resnum) == 0         or
             pdb_resnum[-1] != split_name[0])):
            pdb_resnum.append(split_name[0])
            print pdb_resnum
    pos_obj = pose_from_pdb(pdb_filename) 
    #ADD THE POS NUMBERING FOR TYZ AND CYZ AS A REMARK 7 [pos numbers]    remark = "REMARK 7"
    name    = pdb_filename.split('/')
    name    = name[-1]
    name    = name.split('.')
    name    = name[0]
    out     = options.final_out_path 
    if not out:
        out = ""
    outname = out +  name + ".pdb"
    pos_obj.dump_pdb(outname)
    remark = "REMARK 7"
    for i in range (len (pdb_resnum)):
        posenum=pos_obj.pdb_info().pdb2pose('A',int( pdb_resnum[i]))
 
        remark= remark + " " + str(posenum)
    with open(outname, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(remark.rstrip('\r\n') + '\n' + content)


if __name__ == "__main__":
    rosetta.init ( extra_options='-ignore_zero_occupancy false @/home/ddz6/Rosetta/flags/subs_YC.flags -extra_res_fa /home/ddz6/Rosetta/params/TYZ.params /home/ddz6/Rosetta/params/CYX.params  -run:preserve_header')
    main(sys.argv[1:])
