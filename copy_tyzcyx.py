#!/usr/bin/env python
from rosetta import *
import math, sys, optparse, os, ast, shutil
from toolbox import *
from transform import *
import rosetta.core.init

#A_pos_obj is a pos with only one chain
#name is the name to output the pdb to
#res_type is the string representing the residue type to be swapped in
#res_set is residue object built by pyrosetta residue factory
#resnum is the original pdb residue number on the other chain
def replace_chain_A_res (A_pos_obj,name,res_type,resnum,all_chi):
#    ncAA = rosetta.core.conformation.ResidueFactory.create_residue(res_set.name_map( res_type))
    posnum = A_pos_obj.pdb_info().pdb2pose( 'A' ,int( resnum ))
    print 1
    print posnum
    print resnum
    print res_type
    A_pos_obj.replace_residue( posnum, res_type, True)
    print 2
    A_pos_obj.residue(posnum).set_all_chi(all_chi)
    print 3
    A_pos_obj.dump_pdb(name)
    return
    
def makeSymm (pos_obj, symm_file, dump=False):
    if dump != False:
        core.pose.symmetry.make_symmetric_pose(pos_obj,symm_file)
        pos_obj.dump_pdb(dump)
    else:
        core.pose.symmetry.make_symmetric_pose(pos_obj, symm_file)
        owd=os.getcwd()
        if dump == owd:
            pos_obj.dump_pdb("name")
        else:
            os.chdir(dump)
            pos_obj.dump_pdb("name")
            os.chdir(owd)

def main(argv):
    #command line options
    #
    
    parser = optparse.OptionParser(usage="\n\nCopy TYZ and CYX on all chains but A to chain A and symmeterize chain A. Requires an existing .symm file.")

    parser.add_option('--pdb-path', dest = 'pdb_path',
        help = 'The path to the pdbs listed in the outcsv \n \
                Example: /home/user/my_pdbs/')

    parser.add_option ('--chain-a', dest = 'A_pdb',
        default='',
        help = "the file location of the chain a only pdb")

    parser.add_option('-o', dest = 'final_out_path',
        default = False,
        help = 'The path to where the output INPUT files should be stored. \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('--rot-lib', dest = 'rot_lib',
        default = '',
        help="Rotamer library dictionary file: ( default = '' )")

    parser.add_option('--cst-path', dest = 'cst_path',
        default = './',
        help="Output cst files path: ( default = 'execute user's home dir' )")

    parser.add_option('--flag-path', dest = 'flag_path',
        default = './',
        help="Output flag files path: ( default = 'execute user's home dir' )")

    parser.add_option('--params-file-TYZ', dest = 'params_file_TYZ',
        default = '',
        help="Specify path to params file when mutating to ncAA.")
 
    parser.add_option('--params-file-CYX', dest = 'params_file_CYX',
        default = '',
        help="Specify path to params file when mutating to ncAA.")
    parser.add_option('--symm-file', dest = 'symm_file',
        default = '',
        help="Specify path to the symm file to use when oligomerizing")

    
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
    non_A_resnum        = [] 
    non_A_restype       = []
    non_A_reschain      = []
    resnum2pose         = []
    
    #this for loop uses res_to_check, a string which must contain all
    #+ residue types to check for in chains other than A
    #+ this loop also assumes all atoms belonging to a residue
    #+ will be sequential
    #+ it checks only if the previous addition was a duplicate before appending
    #+ not the entire list.

    #Populate non_A residue lists

#BEGIN FOR LOOP
    for name in pdb_obj.xyz_names:  #Iterates through the entire pdb, and only saves
        split_name = name.split('_')
        if (split_name[1] == "A" and
            split_name[2] in res_to_check and
            (len(resnum2pose) == 0 or resnum2pose[-1] != split_name[0])) :
           resnum2pose.append(split_name[0])
        if (split_name[1] != "A" and         # if a residue is not already in chain A
            split_name[2] in res_to_check and #+ if a residue is
            (len(non_A_reschain) == 0 or      #+ actually to be copied
            non_A_reschain[-1]   != split_name[1] or    #+if the residue was
                                                        #+ added already
            non_A_resnum[-1]     != split_name[0] ) ):  #+if the residue 
                                                     #+was added already
           non_A_resnum.append(split_name[0])
           non_A_reschain.append(split_name[1])
           non_A_restype.append(split_name[2])
#END FOR LOOP
#    res_set = Vector1("'%s' '%s'" % (options.params_file_CYX,options.params_file_TYZ))
#    res_set = generate_nonstandard_residue_set(res_set)
#    res_set_TYZ = generate_nonstandard_residue_set( Vector1( [options.params_file_TYZ] ))

    resnum2pose.extend(non_A_resnum)
    print resnum2pose
    non_A_chi            = [None] * len(non_A_reschain)
    #This makes a pose object from the whole protein
    #+gets the chis of the residue to be replaced in chains B, C, D etc
    #Requires that options.RESIDUE_params_file points to the
    #+correct params file
    pos_obj = pose_from_pdb(pdb_filename)
    for i in range(0,len(non_A_reschain)):
        #Gets the pos number of the each residue not in A that should be
        posnum = pos_obj.pdb_info().pdb2pose(non_A_reschain[i], int(non_A_resnum[i]))
        #and copies the chis to a python list
        non_A_chi[i] = pos_obj.residue(posnum).chi()
        
    #takes A_pdb, a pdb generated only from chain A and replaces the target
    #+residues (non_A_resnums[i]) with the residue (non_A_restypes[i])
    #+setting the chi to (non_A_chi[i])
    #+writes the final pdb to SOURCEFILENAME_all2A_mono.pdb
    #Requires that options.RESIDUE_params_file points to the
    #+correct params files

    name    = pdb_filename.split('/')
    name    = name[-1]
    name    = name.split('.')
    name    = name[0]

    repname = name + '_allreplace_mono.pdb'
   
    A_pdb     = options.A_pdb
    A_pos_obj = pose_from_pdb(A_pdb)
    for i in range(0, len(non_A_resnum)):
        all_chi = non_A_chi[i]
        resnum   = non_A_resnum[i]
        res_type = pos_obj.residue(pos_obj.pdb_info().pdb2pose(non_A_reschain[i],int(non_A_resnum[i])))
        replace_chain_A_res (A_pos_obj,repname,res_type,resnum,all_chi) 
    #ADD THE POS NUMBERING FOR TYZ AND CYZ AS A REMARK 7 [pos numbers]
    remark = "REMARK 7"
    for i in range (len (resnum2pose)):
        posenum=pos_obj.pdb_info().pdb2pose('A',int( resnum2pose[i]))
        remark= remark + " " + str(posenum)
    with open(repname, 'r+') as f:
        content = f.read()
        f.seek(0, 0)
        f.write(remark.rstrip('\r\n') + '\n' + content)
    #Symmeterizes SOURCEFILENAME_allreplace_mono.pdb
    #+ writes the result to SOURCEFILENAME_allreplace_oligo.pdb
    #+ uses the symm file proved as symm_file
    out_path = options.final_out_path
    
    oligoname = name + '_allreplace_oligo.pdb'
    symm_file = options.symm_file
    out_path  = out_path + oligoname
    
    makeSymm(A_pos_obj, symm_file,out_path )

if __name__ == "__main__":
    rosetta.init ( extra_options='-ignore_zero_occupancy false @/home/ddz6/Rosetta/flags/subs_YC.flags -extra_res_fa /home/ddz6/Rosetta/params/TYZ.params /home/ddz6/Rosetta/params/CYZ.params  -run:preserve_header')
    main(sys.argv[1:])
