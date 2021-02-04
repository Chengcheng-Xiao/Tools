#!/usr/bin/env python
"""
A script to make generate PP for VASP

set env -> VASP_PP_PATH to your PP dir
LDA -> $VASP_PP_PATH/potpaw
PBE -> $VASP_PP_PATH/potpaw_PBE
"""

import argparse
import os
import sys
import datetime
import time

# Command line praser
#----------------------------
parser = argparse.ArgumentParser(description='A script to make supercell(transformed supercell) and measure atomic distances.')
parser.add_argument('-v', '--verbose', action="store_false", default=True, dest="prt",
                    help='print out info? (Default=False)')
parser.add_argument('-o','--overwrite', action="store_true", default=False, dest='over',
                    help='Overwrite POTCAR?. (Default=False)')
parser.add_argument('-c','--POSCAR', action="store", default="POSCAR", dest="POSCAR",
                    help='Input file name.')
parser.add_argument('-i', action="store", default=None, dest="atom_list",
                    help='atom list (Default=None)')
parser.add_argument('-s', action="store", default="recommended", dest="setup",
                    help='setup (minimal | recommended) | materialsproject | gw )')
parser.add_argument('-xc', action="store", default="PBE", dest="xc",
                    help='setup (PBE | LDA)')

prm = parser.parse_args()

# Starting
#----------------------------
if prm.prt == True:
    starttime = time.clock()
    print "Starting at",
    print time.strftime("%H:%M:%S on %a %d %b %Y")

# check input
if prm.atom_list == None and not os.path.isfile(prm.POSCAR):
    sys.stdout.write("\033[1;31m" ) # set color red
    print "\n** ERROR: Initial position file %s was not found." % prm.POSCAR
    sys.stdout.write("\033[0;0m") # reset color
    sys.exit(0)

if prm.atom_list == None:
    f = open(prm.POSCAR)
    a = f.readlines()
    atom_list = a[5].strip("\n").split()
    f.close()
elif prm.POSCAR != "POSAR" and prm.atom_list != None:
    sys.stdout.write("\033[1;31m" ) # set color red
    print "\n** WARNING: -i override -c"
    sys.stdout.write("\033[0;0m") # reset color
    atom_list = prm.atom_list.split()
else:
    sys.stdout.write("\033[1;31m" ) # set color red
    print "\n** ERROR: no atom species specified."
    sys.stdout.write("\033[0;0m") # reset color
    sys.exit(0)

# print out atom_list
if prm.prt == True:
    print "\ninput atom species: %s " % atom_list
    print "input xc:           %s " % prm.xc
    print "input setup:        %s " % prm.setup

# Dictionary containing all data
#----------------------------
if prm.xc == "PBE":
    PP_dir = os.environ.get("VASP_PP_PATH")+"/potpaw_PBE/"
    minima = {
        "Ac" : ["Ac"],
        "Ag" : ["Ag","Ag_pv"],
        "Al" : ["Al"],
        "Am" : ["Am"],
        "Ar" : ["Ar"],
        "As" : ["As","As_d"],
        "At" : ["At"],
        "Au" : ["Au"],
        "B"  : ["B","B_h","B_s"],
        "Ba" : ["Ba_sv"],
        "Be" : ["Be", "Be_sv"],
        "Bi" : ["Bi", "Bi_d"],
        "Br" : ["Br"],
        "C"  : ["C","C_h","C_s"],
        "Ca" : ["Ca_pv","Ca_sv"],
        "Cd" : ["Cd"],
        "Ce" : ["Ce","Ce_3","Ce_h"],
        "Cf" : ["Cf"],
        "Cl" : ["Cl","Cl_h"],
        "Cm" : ["Cm"],
        "Co" : ["Co","Co_pv","Co_sv"],
        "Cr" : ["Cr","Cr_pv","Cr_sv"],
        "Cs" : ["Cs_sv"],
        "Cu" : ["Cu","Cu_pv"],
        "Dy" : ["Dy","Dy_3"],
        "Er" : ["Er","Er_2","Er_3"],
        "Eu" : ["Eu","Eu_2","Eu_3"],
        "F"  : ["F","F_h","F_s"],
        "Fe" : ["Fe","Fe_pv","Fe_sv"],
        "Fr" : ["Fr_sv"],
        "Ga" : ["Ga","Ga_d","Ga_h"],
        "Gd" : ["Gd","Gd_3"],
        "Ge" : ["Ge","Ge_d","Ge_h"],
        "H"  : ["H","H.25","H.33","H.42","H.5","H.58","H.66","H.75","H1.25","H1.33","H1.5","H1.66","H1.75","H_AE","H_h","H_s"],
        "He" : ["He","He_AE"],
        "Hf" : ["Hf","Hf_pv","Hf_sv"],
        "Hg" : ["Hg"],
        "Ho" : ["Ho","Ho_3"],
        "I"  : ["I"],
        "In" : ["In","In_d"],
        "Ir" : ["Ir"],
        "K"  : ["K_pv","K_sv"],
        "Kr" : ["Kr"],
        "La" : ["La","La_s"],
        "Li" : ["Li","Li_sv"],
        "Lu" : ["Lu","Lu_3"],
        "Mg" : ["Mg","Mg_pv","Mg_sv"],
        "Mn" : ["Mn","Mn_pv","Mn_sv"],
        "Mo" : ["Mo","Mo_pv","Mo_sv"],
        "N"  : ["N","N_h","N_s"],
        "Na" : ["Na","Na_pv","Na_sv"],
        "Nb" : ["Nb_pv","Nb_sv"],
        "Nd" : ["Nd","Nd_3"],
        "Ne" : ["Ne"],
        "Ni" : ["Ni","Ni_pv"],
        "Np" : ["Np","Np_s"],
        "O"  : ["O","O_h","O_s"],
        "Os" : ["Os","Os_pv"],
        "P"  : ["P","P_h"],
        "Pa" : ["Pa","Pa_s"],
        "Pb" : ["Pb","Pb_d"],
        "Pd" : ["Pd","Pd_pv"],
        "Pm" : ["Pm","Pm_3"],
        "Po" : ["Po","Po_d"],
        "Pr" : ["Pr","Pr_3"],
        "Pt" : ["Pt","Pt_pv"],
        "Pu" : ["Pu","Pu_s"],
        "Ra" : ["Ra_sv"],
        "Rb" : ["Rb_pv","Rb_sv"],
        "Re" : ["Re","Re_pv"],
        "Rh" : ["Rh","Rh_pv"],
        "Rn" : ["Rn"],
        "Ru" : ["Ru","Ru_pv","Ru_sv"],
        "S"  : ["S","S_h"],
        "Sb" : ["Sb"],
        "Sc" : ["Sc","Sc_sv"],
        "Se" : ["Se"],
        "Si" : ["Si"],
        "Sm" : ["Sm","Sm_3"],
        "Sn" : ["Sn","Sn_d"],
        "Sr" : ["Sr_sv"],
        "Ta" : ["Ta","Ta_pv"],
        "Tb" : ["Tb","Tb_3"],
        "Tc" : ["Tc","Tc_pv","Tc_sv"],
        "Te" : ["Te"],
        "Th" : ["Th","Th_s"],
        "Ti" : ["Ti","Ti_pv","Ti_sv"],
        "Tl" : ["Tl","Tl_d"],
        "Tm" : ["Tm","Tm_3"],
        "U"  : ["U","U_s"],
        "V"  : ["V","V_pv","V_sv"],
        "W"  : ["W","W_pv","W_sv"],
        "Xe" : ["Xe"],
        "Y"  : ["Y_sv"],
        "Yb" : ["Yb","Yb_2","Yb_3"],
        "Zn" : ["Zn"],
        "Zr" : ["Zr_sv"]
    }

    #https://cms.mpi.univie.ac.at/vasp/vasp/Recommended_PAW_potentials_DFT_calculations_using_vasp_5_2.html
    recommended = {
        "Ac" : ["Ac"],
        "Ag" : ["Ag_pv","Ag"],
        "Al" : ["Al"],
        "Am" : ["Am"],
        "Ar" : ["Ar"],
        "As" : ["As","As_d"],
        "At" : ["At"],
        "Au" : ["Au"],
        "B"  : ["B","B_h","B_s"],
        "Ba" : ["Ba_sv"],
        "Be" : ["Be", "Be_sv"],
        "Bi" : ["Bi_d","Bi"],
        "Br" : ["Br"],
        "C"  : ["C","C_h","C_s"],
        "Ca" : ["Ca_sv","Ca_pv"],
        "Cd" : ["Cd"],
        "Ce" : ["Ce","Ce_3","Ce_h"],
        "Cf" : ["Cf"],
        "Cl" : ["Cl","Cl_h"],
        "Cm" : ["Cm"],
        "Co" : ["Co","Co_pv","Co_sv"],
        "Cr" : ["Cr_pv","Cr","Cr_sv"],
        "Cs" : ["Cs_sv"],
        "Cu" : ["Cu","Cu_pv"],
        "Dy" : ["Dy_3","Dy"],
        "Er" : ["Er_3","Er","Er_2"],
        "Eu" : ["Eu_2","Eu","Eu_3"],
        "F"  : ["F","F_h","F_s"],
        "Fe" : ["Fe","Fe_pv","Fe_sv"],
        "Fr" : ["Fr_sv"],
        "Ga" : ["Ga_d","Ga","Ga_h"],
        "Gd" : ["Gd_3","Gd"],
        "Ge" : ["Ge_d","Ge","Ge_h"],
        "H"  : ["H","H.25","H.33","H.42","H.5","H.58","H.66","H.75","H1.25","H1.33","H1.5","H1.66","H1.75","H_AE","H_h","H_s"],
        "He" : ["He","He_AE"],
        "Hf" : ["Hf_pv","Hf","Hf_sv"],
        "Hg" : ["Hg"],
        "Ho" : ["Ho_3","Ho"],
        "I"  : ["I"],
        "In" : ["In_d","In"],
        "Ir" : ["Ir"],
        "K"  : ["K_sv","K_pv"],
        "Kr" : ["Kr"],
        "La" : ["La","La_s"],
        "Li" : ["Li_sv","Li"],
        "Lu" : ["Lu_3","Lu"],
        "Mg" : ["Mg","Mg_pv","Mg_sv"],
        "Mn" : ["Mn_pv","Mn","Mn_sv"],
        "Mo" : ["Mo_sv","Mo","Mo_pv"],
        "N"  : ["N","N_h","N_s"],
        "Na" : ["Na_pv","Na","Na_sv"],
        "Nb" : ["Nb_sv","Nb_pv"],
        "Nd" : ["Nd_3","Nd"],
        "Ne" : ["Ne"],
        "Ni" : ["Ni","Ni_pv"],
        "Np" : ["Np","Np_s"],
        "O"  : ["O","O_h","O_s"],
        "Os" : ["Os","Os_pv"],
        "P"  : ["P","P_h"],
        "Pa" : ["Pa","Pa_s"],
        "Pb" : ["Pb_d","Pb"],
        "Pd" : ["Pd","Pd_pv"],
        "Pm" : ["Pm_3","Pm"],
        "Po" : ["Po_d","Po"],
        "Pr" : ["Pr_3","Pr"],
        "Pt" : ["Pt","Pt_pv"],
        "Pu" : ["Pu","Pu_s"],
        "Ra" : ["Ra_sv"],
        "Rb" : ["Rb_sv","Rb_pv"],
        "Re" : ["Re","Re_pv"],
        "Rh" : ["Rh_pv","Rh"],
        "Rn" : ["Rn"],
        "Ru" : ["Ru_pv","Ru","Ru_sv"],
        "S"  : ["S","S_h"],
        "Sb" : ["Sb"],
        "Sc" : ["Sc_sv","Sc"],
        "Se" : ["Se"],
        "Si" : ["Si"],
        "Sm" : ["Sm_3","Sm"],
        "Sn" : ["Sn_d","Sn"],
        "Sr" : ["Sr_sv"],
        "Ta" : ["Ta_pv","Ta"],
        "Tb" : ["Tb_3","Tb"],
        "Tc" : ["Tc_pv","Tc","Tc_sv"],
        "Te" : ["Te"],
        "Th" : ["Th","Th_s"],
        "Ti" : ["Ti_sv","Ti","Ti_pv"],
        "Tl" : ["Tl_d","Tl"],
        "Tm" : ["Tm_3","Tm"],
        "U"  : ["U","U_s"],
        "V"  : ["V_sv","V","V_pv"],
        "W"  : ["W_sv","W","W_pv"],
        "Xe" : ["Xe"],
        "Y"  : ["Y_sv"],
        "Yb" : ["Yb_2","Yb","Yb_3"],
        "Zn" : ["Zn"],
        "Zr" : ["Zr_sv"]
    }
    gw = {
        "Ag" : ["Ag_sv_GW","Ag_GW"],
        "Al" : ["Al_GW","Al_sv_GW"],
        "Ar" : ["Ar_GW"],
        "As" : ["As_GW","As_sv_GW"],
        "At" : ["At_d_GW","At_sv_GW"],
        "Au" : ["Au_sv_GW","Au_GW"],
        "B"  : ["B_GW"],
        "Ba" : ["Ba_sv_GW"],
        "Be" : ["Be_sv_GW","Be_GW"],
        "Bi" : ["Bi_d_GW","Bi_GW","Bi_sv_GW"],
        "Br" : ["Br_GW","Br_sv_GW"],
        "C"  : ["C_GW","C_GW_new","C_h_GW"],
        "Ca" : ["Ca_sv_GW"],
        "Cd" : ["Cd_sv_GW","Cd_GW"],
        "Ce" : ["Ce_GW"],
        "Cl" : ["Cl_GW"],
        "Co" : ["Co_sv_GW","Co_GW"],
        "Cr" : ["Cr_sv_GW"],
        "Cs" : ["Cs_sv_GW"],
        "Cu" : ["Cu_sv_GW","Cu_GW"],
        "F"  : ["F_GW","F_GW_new","F_h_GW"],
        "Fe" : ["Fe_sv_GW","Fe_GW",],
        "Ga" : ["Ga_d_GW","Ga_GW","Ga_sv_GW"],
        "Ge" : ["Ge_d_GW","Ge_GW","Ge_sv_GW"],
        "H"  : ["H_GW","H_h_GW"],
        "He" : ["He_GW"],
        "Hf" : ["Hf_sv_GW"],
        "Hg" : ["Hg_sv_GW"],
        "I"  : ["I_GW","I_sv_GW"],
        "In" : ["In_d_GW","In_sv_GW"],
        "Ir" : ["Ir_sv_GW"],
        "K"  : ["K_sv_GW"],
        "Kr" : ["Kr_GW"],
        "La" : ["La_GW"],
        "Li" : ["Li_sv_GW","Li_AE_GW","Li_GW"],
        "Mg" : ["Mg_sv_GW","Mg_GW","Mg_pv_GW"],
        "Mn" : ["Mn_sv_GW","Mn_GW"],
        "Mo" : ["Mo_sv_GW"],
        "N"  : ["N_GW","N_GW_new","N_h_GW","N_s_GW"],
        "Na" : ["Na_sv_GW"],
        "Nb" : ["Nb_sv_GW"],
        "Ne" : ["Ne_GW","Ne_s_GW"],
        "Ni" : ["Ni_sv_GW","Ni_GW"],
        "O"  : ["O_GW","O_GW_new","O_h_GW","O_s_GW"],
        "Os" : ["Os_sv_GW"],
        "P"  : ["P_GW"],
        "Pb" : ["Pb_sv_GW","Pb_d_GW","Pd_GW"],
        "Pb" : ["Pd_sv_GW"],
        "Po" : ["Po_d_GW","Po_sv_GW"],
        "Pt" : ["Pt_sv_GW","Pt_GW"],
        "Rb" : ["Rb_sv_GW"],
        "Re" : ["Re_sv_GW"],
        "Rh" : ["Rh_sv_GW","Rh_GW"],
        "Rn" : ["Rn_d_GW","Rn_sv_GW"],
        "Ru" : ["Ru_sv_GW"],
        "S"  : ["S_GW"],
        "Sb" : ["Sb_d_GW","Sb_GW","Sb_sv_GW"],
        "Sc" : ["Sc_sv_GW"],
        "Se" : ["Se_GW","Se_sv_GW"],
        "Si" : ["Si_GW","Si_sv_GW"],
        "Sn" : ["Sn_d_GW","Sn_sv_GW"],
        "Sr" : ["Sr_sv_GW"],
        "Ta" : ["Ta_sv_GW"],
        "Tc" : ["Tc_sv_GW"],
        "Te" : ["Te_GW","Te_sv_GW"],
        "Ti" : ["Ti_sv_GW"],
        "Tl" : ["Tl_d_GW","Tl_sv_GW"],
        "V"  : ["V_sv_GW"],
        "W"  : ["W_sv_GW"],
        "Xe" : ["Xe_GW","Xe_sv_GW"],
        "Y"  : ["Y_sv_GW"],
        "Zn" : ["Zn_sv_GW","Zn_GW"],
        "Zr" : ["Zr_sv_GW"]
    }
elif prm.xc == "LDA":
    PP_dir = os.environ.get("VASP_PP_PATH")+"/potpaw/"
    minima = {
        "Ac" : ["Ac"],
        "Ag" : ["Ag","Ag_pv"],
        "Al" : ["Al"],
        "Am" : ["Am"],
        "Ar" : ["Ar"],
        "As" : ["As","As_d"],
        "At" : ["At"],
        "Au" : ["Au"],
        "B"  : ["B","B_h","B_s"],
        "Ba" : ["Ba_sv"],
        "Be" : ["Be", "Be_sv"],
        "Bi" : ["Bi", "Bi_d"],
        "Br" : ["Br"],
        "C"  : ["C","C_h","C_s"],
        "Ca" : ["Ca_pv","Ca_sv"],
        "Cd" : ["Cd"],
        "Ce" : ["Ce","Ce_h"],
        "Cl" : ["Cl","Cl_h"],
        "Cm" : ["Cm"],
        "Co" : ["Co","Co_pv","Co_sv"],
        "Cr" : ["Cr","Cr_pv","Cr_sv"],
        "Cs" : ["Cs_sv"],
        "Cu" : ["Cu","Cu_pv"],
        "F"  : ["F","F_h","F_s"],
        "Fe" : ["Fe","Fe_pv","Fe_sv"],
        "Fr" : ["Fr_sv"],
        "Ga" : ["Ga","Ga_d","Ga_h"],
        "Ge" : ["Ge","Ge_d","Ge_h"],
        "H"  : ["H","H.25","H.33","H.42","H.5","H.58","H.66","H.75","H1.25","H1.33","H1.5","H1.66","H1.75","H_AE","H_h","H_s"],
        "He" : ["He"],
        "Hf" : ["Hf","Hf_pv","Hf_sv"],
        "Hg" : ["Hg"],
        "I"  : ["I"],
        "In" : ["In","In_d"],
        "Ir" : ["Ir"],
        "K"  : ["K_pv","K_sv"],
        "Kr" : ["Kr"],
        "La" : ["La","La_s"],
        "Li" : ["Li","Li_sv"],
        "Mg" : ["Mg","Mg_pv","Mg_sv"],
        "Mn" : ["Mn","Mn_pv","Mn_sv"],
        "Mo" : ["Mo","Mo_pv","Mo_sv"],
        "N"  : ["N","N_h","N_s"],
        "Na" : ["Na","Na_pv","Na_sv"],
        "Nb" : ["Nb_pv","Nb_sv"],
        "Ne" : ["Ne"],
        "Ni" : ["Ni","Ni_pv"],
        "Np" : ["Np","Np_s"],
        "O"  : ["O","O_h","O_s"],
        "Os" : ["Os","Os_pv"],
        "P"  : ["P","P_h"],
        "Pa" : ["Pa","Pa_s"],
        "Pb" : ["Pb","Pb_d"],
        "Pd" : ["Pd","Pd_pv"],
        "Po" : ["Po","Po_d"],
        "Pt" : ["Pt","Pt_pv"],
        "Pu" : ["Pu","Pu_s"],
        "Ra" : ["Ra_sv"],
        "Rb" : ["Rb_pv","Rb_sv"],
        "Re" : ["Re","Re_pv"],
        "Rh" : ["Rh","Rh_pv"],
        "Rn" : ["Rn"],
        "Ru" : ["Ru","Ru_pv","Ru_sv"],
        "S"  : ["S","S_h"],
        "Sb" : ["Sb"],
        "Sc" : ["Sc","Sc_sv"],
        "Se" : ["Se"],
        "Si" : ["Si"],
        "Sn" : ["Sn","Sn_d"],
        "Sr" : ["Sr_sv"],
        "Ta" : ["Ta","Ta_pv"],
        "Tc" : ["Tc","Tc_pv","Tc_sv"],
        "Te" : ["Te"],
        "Th" : ["Th","Th_s"],
        "Ti" : ["Ti","Ti_pv","Ti_sv"],
        "Tl" : ["Tl","Tl_d"],
        "U"  : ["U","U_s"],
        "V"  : ["V","V_pv","V_sv"],
        "W"  : ["W","W_pv","W_sv"],
        "Xe" : ["Xe"],
        "Y"  : ["Y_sv"],
        "Zn" : ["Zn"],
        "Zr" : ["Zr_sv"]
    }

    #https://cms.mpi.univie.ac.at/vasp/vasp/Recommended_PAW_potentials_DFT_calculations_using_vasp_5_2.html
    recommended = {
        "Ac" : ["Ac"],
        "Ag" : ["Ag_pv","Ag"],
        "Al" : ["Al"],
        "Am" : ["Am"],
        "Ar" : ["Ar"],
        "As" : ["As","As_d"],
        "At" : ["At"],
        "Au" : ["Au"],
        "B"  : ["B","B_h","B_s"],
        "Ba" : ["Ba_sv"],
        "Be" : ["Be", "Be_sv"],
        "Bi" : ["Bi_d","Bi"],
        "Br" : ["Br"],
        "C"  : ["C","C_h","C_s"],
        "Ca" : ["Ca_sv","Ca_pv"],
        "Cd" : ["Cd"],
        "Ce" : ["Ce","Ce_h"],
        "Cl" : ["Cl","Cl_h"],
        "Cm" : ["Cm"],
        "Co" : ["Co","Co_pv","Co_sv"],
        "Cr" : ["Cr_pv","Cr","Cr_sv"],
        "Cs" : ["Cs_sv"],
        "Cu" : ["Cu","Cu_pv"],
        "F"  : ["F","F_h","F_s"],
        "Fe" : ["Fe","Fe_pv","Fe_sv"],
        "Fr" : ["Fr_sv"],
        "Ga" : ["Ga_d","Ga","Ga_h"],
        "Ge" : ["Ge_d","Ge","Ge_h"],
        "H"  : ["H","H.25","H.33","H.42","H.5","H.58","H.66","H.75","H1.25","H1.33","H1.5","H1.66","H1.75","H_AE","H_h","H_s"],
        "He" : ["He"],
        "Hf" : ["Hf_pv","Hf","Hf_sv"],
        "Hg" : ["Hg"],
        "I"  : ["I"],
        "In" : ["In_d","In"],
        "Ir" : ["Ir"],
        "K"  : ["K_sv","K_pv"],
        "Kr" : ["Kr"],
        "La" : ["La","La_s"],
        "Li" : ["Li_sv","Li"],
        "Mg" : ["Mg","Mg_pv","Mg_sv"],
        "Mn" : ["Mn_pv","Mn","Mn_sv"],
        "Mo" : ["Mo_sv","Mo","Mo_pv"],
        "N"  : ["N","N_h","N_s"],
        "Na" : ["Na_pv","Na","Na_sv"],
        "Nb" : ["Nb_sv","Nb_pv"],
        "Ne" : ["Ne"],
        "Ni" : ["Ni","Ni_pv"],
        "Np" : ["Np","Np_s"],
        "O"  : ["O","O_h","O_s"],
        "Os" : ["Os","Os_pv"],
        "P"  : ["P","P_h"],
        "Pa" : ["Pa","Pa_s"],
        "Pb" : ["Pb_d","Pb"],
        "Pd" : ["Pd","Pd_pv"],
        "Po" : ["Po_d","Po"],
        "Pt" : ["Pt","Pt_pv"],
        "Pu" : ["Pu","Pu_s"],
        "Ra" : ["Ra_sv"],
        "Rb" : ["Rb_sv","Rb_pv"],
        "Re" : ["Re","Re_pv"],
        "Rh" : ["Rh_pv","Rh"],
        "Rn" : ["Rn"],
        "Ru" : ["Ru_pv","Ru","Ru_sv"],
        "S"  : ["S","S_h"],
        "Sb" : ["Sb"],
        "Sc" : ["Sc_sv","Sc"],
        "Se" : ["Se"],
        "Si" : ["Si"],
        "Sn" : ["Sn_d","Sn"],
        "Sr" : ["Sr_sv"],
        "Ta" : ["Ta_pv","Ta"],
        "Tc" : ["Tc_pv","Tc","Tc_sv"],
        "Te" : ["Te"],
        "Th" : ["Th","Th_s"],
        "Ti" : ["Ti_sv","Ti","Ti_pv"],
        "Tl" : ["Tl_d","Tl"],
        "U"  : ["U","U_s"],
        "V"  : ["V_sv","V","V_pv"],
        "W"  : ["W_sv","W","W_pv"],
        "Xe" : ["Xe"],
        "Y"  : ["Y_sv"],
        "Zn" : ["Zn"],
        "Zr" : ["Zr_sv"]
    }
    gw = {
        "Ag" : ["Ag_sv_GW","Ag_GW"],
        "Al" : ["Al_GW","Al_sv_GW"],
        "Ar" : ["Ar_GW"],
        "As" : ["As_GW","As_sv_GW"],
        "At" : ["At_d_GW","At_sv_GW"],
        "Au" : ["Au_sv_GW","Au_GW"],
        "B"  : ["B_GW"],
        "Ba" : ["Ba_sv_GW"],
        "Be" : ["Be_sv_GW","Be_GW"],
        "Bi" : ["Bi_d_GW","Bi_GW","Bi_sv_GW"],
        "Br" : ["Br_GW","Br_sv_GW"],
        "C"  : ["C_GW","C_GW_new","C_h_GW"],
        "Ca" : ["Ca_sv_GW"],
        "Cd" : ["Cd_sv_GW","Cd_GW"],
        "Ce" : ["Ce_GW"],
        "Cl" : ["Cl_GW"],
        "Co" : ["Co_sv_GW","Co_GW"],
        "Cr" : ["Cr_sv_GW"],
        "Cs" : ["Cs_sv_GW"],
        "Cu" : ["Cu_sv_GW","Cu_GW"],
        "F"  : ["F_GW","F_GW_new","F_h_GW"],
        "Fe" : ["Fe_sv_GW","Fe_GW",],
        "Ga" : ["Ga_d_GW","Ga_GW","Ga_sv_GW"],
        "Ge" : ["Ge_d_GW","Ge_GW","Ge_sv_GW"],
        "H"  : ["H_GW","H_h_GW"],
        "He" : ["He_GW"],
        "Hf" : ["Hf_sv_GW"],
        "Hg" : ["Hg_sv_GW"],
        "I"  : ["I_GW","I_sv_GW"],
        "In" : ["In_d_GW","In_sv_GW"],
        "Ir" : ["Ir_sv_GW"],
        "K"  : ["K_sv_GW"],
        "Kr" : ["Kr_GW"],
        "La" : ["La_GW"],
        "Li" : ["Li_sv_GW","Li_AE_GW","Li_GW"],
        "Mg" : ["Mg_sv_GW","Mg_GW","Mg_pv_GW"],
        "Mn" : ["Mn_sv_GW","Mn_GW"],
        "Mo" : ["Mo_sv_GW"],
        "N"  : ["N_GW","N_GW_new","N_h_GW","N_s_GW"],
        "Na" : ["Na_sv_GW"],
        "Nb" : ["Nb_sv_GW"],
        "Ne" : ["Ne_GW","Ne_s_GW"],
        "Ni" : ["Ni_sv_GW","Ni_GW"],
        "O"  : ["O_GW","O_GW_new","O_h_GW","O_s_GW"],
        "Os" : ["Os_sv_GW"],
        "P"  : ["P_GW"],
        "Pb" : ["Pb_sv_GW","Pb_d_GW","Pd_GW"],
        "Pb" : ["Pd_sv_GW"],
        "Po" : ["Po_d_GW","Po_sv_GW"],
        "Pt" : ["Pt_sv_GW","Pt_GW"],
        "Rb" : ["Rb_sv_GW"],
        "Re" : ["Re_sv_GW"],
        "Rh" : ["Rh_sv_GW","Rh_GW"],
        "Rn" : ["Rn_d_GW","Rn_sv_GW"],
        "Ru" : ["Ru_sv_GW"],
        "S"  : ["S_GW"],
        "Sb" : ["Sb_d_GW","Sb_GW","Sb_sv_GW"],
        "Sc" : ["Sc_sv_GW"],
        "Se" : ["Se_GW","Se_sv_GW"],
        "Si" : ["Si_GW","Si_sv_GW"],
        "Sn" : ["Sn_d_GW","Sn_sv_GW"],
        "Sr" : ["Sr_sv_GW"],
        "Ta" : ["Ta_sv_GW"],
        "Tc" : ["Tc_sv_GW"],
        "Te" : ["Te_GW","Te_sv_GW"],
        "Ti" : ["Ti_sv_GW"],
        "Tl" : ["Tl_d_GW","Tl_sv_GW"],
        "V"  : ["V_sv_GW"],
        "W"  : ["W_sv_GW"],
        "Xe" : ["Xe_GW","Xe_sv_GW"],
        "Y"  : ["Y_sv_GW"],
        "Zn" : ["Zn_sv_GW","Zn_GW"],
        "Zr" : ["Zr_sv_GW"]
    }

file_list = []
if prm.setup == "recommended":
    for i in atom_list:
        file_list.append(recommended[i][0])
elif prm.setup == "minimal":
    for i in atom_list:
        file_list.append(minimal[i][0])
elif prm.setup == "gw":
    for i in atom_list:
        file_list.append(gw[i][0])
elif prm.setup == "materialsproject":
    sys.stdout.write("\033[1;31m" ) # set color red
    print "\n** ERROR: Not implemented yet."
    sys.stdout.write("\033[0;0m") # reset color
    sys.exit(0)
elif prm.setup == "manual":
    for i in atom_list:
        print i+": ",recommended[i]
        inputstr = raw_input("--->>>:")
        file_list.append(inputstr)
else:
    sys.stdout.write("\033[1;31m" ) # set color red
    print "** ERROR: wrong setup."
    sys.stdout.write("\033[0;0m") # reset color
    sys.exit(0)

# Combine file_list
#----------------------------
if prm.over == True or prm.over == False and os.path.exists("POTCAR") == False:
    with open('POTCAR', 'w+') as outfile:
        for fname in file_list:
            with open(PP_dir+fname+"/POTCAR") as infile:
                for line in infile:
                    outfile.write(line)
elif prm.over == False and os.path.exists("POTCAR") == True:
    sys.stdout.write("\033[1;31m" ) # set color red
    print "\n** ERROR: POTCAR exist. use -o to override."
    sys.stdout.write("\033[0;0m") # reset color
    sys.exit(0)



if prm.prt == True:
    # Post process
    #----------------------------
    endtime = time.clock()
    runtime = endtime-starttime
    print "\nEnd."
    print "Program was running for %.2f seconds." % runtime
