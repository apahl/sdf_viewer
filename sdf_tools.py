#!/usr/bin/env python3
#-*- coding: utf-8 -*-
# sdf_tools.py
# version: 2015-05-28
# author:  Axel Pahl (APL)
# contact: firstnamelastname at gmx dot de
# license: BSD, see license.txt in this folder

#==============================================================================
# current development branch (05-Jun-2015): 
# * make the tools less dependent from the type in the property name 
#   (k_molid, n_LogP, ...)
# * add SAR table to the tools
# -----------------------------------------------------------------------------
# - sort_sdf: done (05-Jun-2015)
#==============================================================================

from __future__ import absolute_import, division, print_function # , unicode_literals

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
import rdkit.Chem.Descriptors as Desc

# imports for similarity search
from rdkit.Chem.Fingerprints import FingerprintMols

# imports for clustering
from rdkit import DataStructs
from rdkit.ML.Cluster import Butina

from PIL import Image, ImageChops

import os
import sys
import random
import math
import csv
import os.path as op

from time import sleep, strftime
from collections import Counter

if sys.version_info[0] > 2:
    import io
    PY3 = True
    file_type = io.IOBase
else:
    PY3 = False
    file_type = file

if hasattr(sys, "ps1"): # <True> when called from an interactive session
    print("  > interactive session, trying to import display_png from IPython...", end="")
    try:
        from IPython.core.display import display_png
        IPYTHON = True
        print(" success!")

    except:
        IPYTHON = False
        print()
        print("  * no inline display of molecule structures supported.")

try:
    import pylab
    PYLAB = True

except:
    PYLAB = False
    print("  * failed to import pylab,")
    print("    plotting will not work.")

try:
    from rdkit_ipynb_tools.tools import Mol_List
    
except ImportError:
    print("  * failed to load RDKit IPython tools")
    Mol_List = list


class NoFieldTypes(Exception):
    def __str__(self):
        return repr("FieldTypeError: field types could not be extracted from sdf list")


def set_sdf_report_folder():

    if "SDF_VIEWER_REPORTS" in os.environ:
        folder = os.environ["SDF_VIEWER_REPORTS"]
        print("  > found environment var", folder)
    else:
        if "HOME" in os.environ:
            folder = op.join(os.environ["HOME"], "sdf_reports")
        elif "HOMEPATH" in os.environ: # Windows
            folder = op.join(os.environ["HOMEPATH"], "sdf_reports")
        else:
            folder = None
        print("  > setting default folder", folder)

    if not op.exists(folder):
        print("  * folder does not exist, creating...")
        os.mkdir(folder)

    for directory in ["sdf", "reports", "html", "session"]:
        subfolder = op.join(folder, directory)
        if not op.exists(subfolder):
            print("    - creating subfolder", subfolder)
            os.mkdir(subfolder)

    return folder


MISSING_VAL = -999
POINTSIZE = 40
REPORT_FOLDER = set_sdf_report_folder()


def create_dir_if_not_exist(dir_name):
    if not op.exists(dir_name):
        print("  * target folder does not exist, creating {}...".format(dir_name))
        os.makedirs(dir_name)


def autocrop(im, bgcolor="white"):
    if im.mode != "RGB":
        im = im.convert("RGB")
    bg = Image.new("RGB", im.size, bgcolor)
    diff = ImageChops.difference(im, bg)
    bbox = diff.getbbox()
    if bbox:
        return im.crop(bbox)
    return None # no contents


def print_dict(d):
    max_len = max(map(len, d.keys()))
    for k in sorted(d.keys()):
        print("  {d_key:{width}s}: {num:7d}".format(d_key=k, width=max_len, num=d[k]))


def load_sdf(file_name_or_obj="testset.sdf", large_sdf=False):
    """load small sdf completely in memory as list; return large sdf as file object
    function accepts a string filename or a file object"""

    if isinstance(file_name_or_obj, str):
        if PY3:
            file_obj = open(file_name_or_obj, "rb")
        else:
            file_obj = open(file_name_or_obj)
    else:
        file_obj = file_name_or_obj

    reader = Chem.ForwardSDMolSupplier(file_obj)

    if large_sdf:
        if isinstance(file_name_or_obj, str):
            print("  > large sdf {} loaded as file object.".format(file_name_or_obj.split(".")[0]))
        else:
            print("  > large sdf loaded as file object.")

        return reader

    sdf_list = Mol_List()

    for mol in reader:
        if mol:
            sdf_list.append(mol)

    if isinstance(file_name_or_obj, str):
        print("  > sdf {} loaded with {} records.".format(file_name_or_obj.split(".")[0], len(sdf_list)))
    else:
        print("  > sdf loaded with {} records.".format(len(sdf_list)))

    return sdf_list


def write_sdf(sdf_list, fn, conf_id=-1):
    if not (isinstance(sdf_list, list) or isinstance(sdf_list, file_type)):
        sdf_list = [sdf_list]

    writer = Chem.SDWriter(fn)

    for mol in sdf_list:
        writer.write(mol, confId=conf_id)

    writer.close()


def write_pdb(sdf_list, fn):
    """sdf_list can be a single molecule or a list of molecules"""

    if not isinstance(sdf_list, list):
        sdf_list = list(sdf_list)

    writer = Chem.PDBWriter(fn)

    for mol in sdf_list:
        writer.write(mol)

    writer.close()


def smiles_supplier(fn, smiles_prop=""):
    """Process a text file containing Smiles or Smarts and other properties (separated by tab)
    and yield mol objects as generator.
    If <smiles_prop> is empty, the first column header containing "smiles" (case-insensitive) 
    is assumed to be the Smiles column.
    If <smiles_prop> contains "smarts" (case-insensitive),
    then the molecule is generated from the Smarts in this column."""
    
    with open(fn) as f:
        first_line = True
        smiles_idx = -1
        props_dict = {}
        for line in f:
            line = line.strip().split("\t")
            if first_line: # determine available properties
                first_line = False
                for idx, prop in enumerate(line):
                    if not smiles_prop:
                        if "smiles" in prop.lower():
                            props_dict[idx] = "Smiles"
                            smiles_idx = idx
                            
                    else:
                        if prop == smiles_prop:
                            props_dict[idx] = smiles_prop
                            smiles_idx = idx
                    
                    props_dict[idx] = prop
                
                if smiles_idx < 0:
                    raise ValueError("no Smiles column found in {}".format(fn))
                
            if "smarts" in smiles_prop.lower():
                mol = Chem.MolFromSmarts(line[smiles_idx])
            else:
                mol = Chem.MolFromSmiles(line[smiles_idx])
            
            if not mol: continue
            
            for idx, prop in enumerate(line):
                if idx != smiles_idx:
                    mol.SetProp(props_dict[idx], prop)
            
            yield mol


def smiles_writer(sdf_list, fn="smiles.csv", smiles_prop="Smiles"):
    """Process an iterator containing RDKit mols and write a tab-separated file with
    the Smiles and the mol properties. The properties have to be present in every mol."""

    with open(fn, "w") as f:
        first_line = True
        for mol in sdf_list:
            if not mol: continue
            if first_line: # determine and write header
                first_line = False
                props = mol.GetPropNames()
                line = [smiles_prop]
                line.extend(props)
                line_str = "\t".join(line) + "\n"
                f.write(line_str)

            line = [Chem.MolToSmiles(mol)]
            line.extend([mol.GetProp(prop) for prop in props])
            line_str = "\t".join(line) + "\n"
            f.write(line_str)


def prepare_for_viewer(sdf_list):
    """deprecated. All functions detect the field types"""
    
    if not isinstance(sdf_list, list):
        print("  * function prepare_for_viewer currently only handles lists.")
        return
    
    return


def iterate_over_reagents_file(fn="testset.sdf", supplier="__guess__",
                               max_num_recs=1000000, mw_low=100, mw_high=350, dryrun=False):

    if not supplier in ["aldrich", "chemspider", "__guess__"]:
        print("  * unknown supplier.")
        print("    aborting.")
        return

    reader = Chem.SDMolSupplier(fn)
    if not dryrun:
        writer = Chem.SDWriter(op.join(REPORT_FOLDER, "sdf", "output.sdf"))
    mol_counter_in = 0
    mol_counter_out = 0
    removals = Counter()

    for mol in reader:
        mol_counter_in += 1

        if mol_counter_in == 1 and supplier == "__guess__":
            if mol.HasProp("CSID"):
                supplier = "chemspider"
            elif mol.HasProp("CAS_NUMBER"):
                supplier = "aldrich"
            else:
                print("  * supplier could not be guessed.")
                print("    aborting.")
                return

            print("  > guessed supplier:", supplier)

        if not mol:
            removals["rejected_by_rdkit"] += 1
            continue

        if supplier == "aldrich":
            remove_props_from_mol(mol, ["ASSAY_NAME", "COMMON_NAME", "MOLECULAR_FORMULA", "MOLECULAR_WEIGHT",
                                      "MOLECULAR_WEIGHT", "BOILING_POINT", "FLASH_POINT", "PRODUCTS", "DENSITY"])

            rename_prop_in_mol(mol, "Similarity", "sim")
            rename_prop_in_mol(mol, "IUPAC_NAME", "name")
            rename_prop_in_mol(mol, "MDL_NUMBER", "mdl")
            rename_prop_in_mol(mol, "CAS_NUMBER", "cas")

        elif supplier == "chemspider":
            remove_props_from_mol(mol, ["MF", "MW", "SMILES", "InChI", "InChIKey", "Data Sources", "References", "PubMed",
                                        "RSC", "CSURL"])

            rename_prop_in_mol(mol, "CSID", "csid")
            # rename_prop_in_mol(mol, "CSURL", "s_url")

        # remove organometallics
        WRITE_TO_OUTPUT = True
        calc_props_in_mol(mol, include_date = False)
        formula = mol.GetProp("formula")
        for element in ["Hg", "Pd", "Pt", "Os", "Mg", "Mn", "Ti", "Zn"]:
            if element in formula:
                WRITE_TO_OUTPUT = False
                removals["organometallic"] += 1
                break

        # remove low or high molwt
        if WRITE_TO_OUTPUT:
            molwt = float(mol.GetProp("molwt"))
            if molwt > mw_high:
                WRITE_TO_OUTPUT = False
                removals["molwt_high"] += 1
            elif molwt < mw_low:
                WRITE_TO_OUTPUT = False
                removals["molwt_low"] += 1


        if WRITE_TO_OUTPUT:
            mol_counter_out += 1
            if not dryrun:
                writer.write(mol)

        if mol_counter_in >= max_num_recs:
            break

        if dryrun and mol_counter_in < 10:
            if not WRITE_TO_OUTPUT:
                print("  *** REMOVED ***")
            show_record(mol)

        if mol_counter_in % 500 == 0:
            print("  > processed: {:7d}   found: {:6d}\r".format(mol_counter_in, mol_counter_out), end="")
            sys.stdout.flush()

    print("  > processed: {:7d}   found: {:6d}".format(mol_counter_in, mol_counter_out))
    print("    done.")

    if not dryrun:
        writer.close()

    print("Molecules removed for the following reasons:")
    for reason in removals:
        print("{:20s}: {:4d}".format(reason, removals[reason]))


def iterate_over_sdf_file(fn="testset.sdf", max_num_recs=1000000, actives_only=False, dryrun=False):
    reader = Chem.SDMolSupplier(fn)
    if not dryrun:
        removals_sdf = []
        writer = Chem.SDWriter("output.sdf")
    mol_counter_in = 0
    mol_counter_out = 0
    removals = Counter()

    for mol in reader:
        mol_counter_in += 1
        if not mol:
            removals["rejected_by_rdkit"] += 1
            continue
        remove_props_from_mol(mol, ["ASSAY_NAME"])

        if mol.HasProp("% activity@ClpP"):
            rename_prop_in_mol(mol, "% activity@ClpP", "clpp_percact")
            try:
                old_value = float(mol.GetProp("clpp_percact"))

            except ValueError:
                removals["valueerror"] += 1
                mol.ClearProp("clpp_percact")
                if not dryrun:
                    mol.SetProp("s_reason", "valueerror")
                    removals_sdf.append(mol)
                continue

            new_value = 100 - old_value
            if actives_only and new_value < 0:
                removals["activity_low"] += 1
                continue
            mol.SetProp("clpp_percinh", str(new_value))

        else:
            removals["no_activity"] += 1
            if not dryrun:
                mol.SetProp("s_reason", "no activity")
                removals_sdf.append(mol)

            continue

        if mol.HasProp("pIC50@ClpP"):
            rename_prop_in_mol(mol, "pIC50@ClpP", "clpp_pic50")
            try:
                old_value = float(mol.GetProp("clpp_pic50"))
            except ValueError: # not a number
                mol.SetProp("clpp_pic50", "n.d.")

        rename_prop_in_mol(mol, "COMPOUND_ID", "molid")
        rename_prop_in_mol(mol, "SUPPLIER", "supplier")
        rename_prop_in_mol(mol, "BATCH_ID", "batchid")

        # remove organometallics
        WRITE_TO_OUTPUT = True
        calc_props_in_mol(mol, include_date = False)
        formula = mol.GetProp("s_formula")
        for element in ["Hg", "Pd", "Pt", "Os", "Mn", "Ti"]:
            if element in formula:
                WRITE_TO_OUTPUT = False
                removals["organometallic"] += 1
                if not dryrun:
                    mol.SetProp("reason", "organometallic")
                    removals_sdf.append(mol)
                break

        # remove low or high molwt
        if WRITE_TO_OUTPUT:
            molwt = float(mol.GetProp("molwt"))
            if molwt > 600:
                WRITE_TO_OUTPUT = False
                removals["molwt_high"] += 1
                if not dryrun:
                    mol.SetProp("s_reason", "molwt high")
                    removals_sdf.append(mol)

            if molwt < 200:
                WRITE_TO_OUTPUT = False
                removals["molwt_low"] += 1
                if not dryrun:
                    mol.SetProp("s_reason", "molwt low")
                    removals_sdf.append(mol)



        if WRITE_TO_OUTPUT:
            mol_counter_out += 1
            if not dryrun:
                writer.write(mol)

        if mol_counter_in >= max_num_recs:
            break

        if dryrun and mol_counter_in < 10:
            if not WRITE_TO_OUTPUT:
                print("  *** REMOVED ***")
            show_record(mol)

        if mol_counter_in % 500 == 0:
            print("  > processed: {:7d}   found: {:6d}\r".format(mol_counter_in, mol_counter_out), end="")
            sys.stdout.flush()

    print("  > processed: {:6d}   found: {:6d}".format(mol_counter_in, mol_counter_out))
    print("    done.")
    print("Molecules removed for the following reasons:", end="")

    if not dryrun:
        writer.close()

        if len(removals_sdf) > 0:
            write_sdf(removals_sdf, "removals.sdf")
            print(" (written to removals.sdf)")
        else:
            print(" (no molecules in removals_sdf)")
    else:
        print()

    print_dict(removals)


def enum_racemates(sdf_list_or_file, find_only=True, mol_id="molid"):
    """returns: result_sdf::list<mol>, racemic_molids::list<int>
    find_only==True:  return new sdf as list which contains all the racemates of the input sdf.
    find_only==False: return new sdf as list with ALL input structures, where the racemates are
                      replaced by their two enantiomers. The returned sdf is always
                      equal in size or larger as the input sdf.
    Multiple stereo centers are not yet handled.
    In the new sdf the molids are no longer unique and should be reassigned
    (remove molid and run calc_props(sdf))."""

    result_sdf = []
    racemic_molids = []

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return None

    for mol in sdf_list_or_file:
        first_undefined = False
        chiral_centers  = Chem.FindMolChiralCenters(mol, includeUnassigned=True)

        if chiral_centers:
            first_center    = chiral_centers[0][0]
            first_undefined = chiral_centers[0][1] == "?"

        if first_undefined:
            racemic_molids.append(int(mol.GetProp(mol_id)))
            if find_only:
                result_sdf.append(mol)
                continue

            else:
                mol1 = Chem.Mol(mol)
                mol2 = Chem.Mol(mol)
                mol1.GetAtomWithIdx(first_center).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CW)
                mol2.GetAtomWithIdx(first_center).SetChiralTag(Chem.rdchem.ChiralType.CHI_TETRAHEDRAL_CCW)
                result_sdf.append(mol1)
                result_sdf.append(mol2)

        else:
            if not find_only: # return ALL mols
                result_sdf.append(mol)

    return result_sdf, racemic_molids


def list_fields(sdf_list_or_file):
    field_list = []

    if isinstance(sdf_list_or_file, list):
        if len(sdf_list_or_file) > 100:
            sdf_sample = random.sample(sdf_list_or_file, len(sdf_list_or_file)//2)
        else:
            sdf_sample = sdf_list_or_file

        for mol in sdf_sample:
            field_list.extend(mol.GetPropNames())

    else: # sdf is file
        if sdf_list_or_file.atEnd():
            print("  * file object is at end, please reload.")
            return None
        index = 0
        while sdf_list_or_file and index < 500:
            try:
                mol = sdf_list_or_file.next()
                field_list.extend(mol.GetPropNames())
            except StopIteration:
                break

    return list(set(field_list))


def get_field_types(sdf_list_or_file):
    """detect all the property field types and return as dict"""
    
    print("  > detecting field types...")
    field_types = {}

    if isinstance(sdf_list_or_file, list):
        if len(sdf_list_or_file) > 100:
            sdf_sample = random.sample(sdf_list_or_file, len(sdf_list_or_file)//2)
        else:
            sdf_sample = sdf_list_or_file

    else: # sdf is file
        sdf_sample = []
        if sdf_list_or_file.atEnd():
            raise Exception("  * file object is at end, please reload.")

        index = 0
        while sdf_list_or_file and index < 500:
            try:
                mol = sdf_list_or_file.next()
                sdf_sample.append(mol)
            except StopIteration:
                break

    for mol in sdf_sample:
        prop_names = mol.GetPropNames()
        
        for prop in prop_names:
            prop_type = "number"
            prop_str = mol.GetProp(prop)
            
            try:
                val = float(prop_str)
                if prop.lower().endswith("id"):
                    prop_type = "key"
                    
            except ValueError:
                prop_type = "str"
                
            if prop in field_types:
                if field_types[prop] in ["number", "key"] and prop_type == "str":
                    # "str" overrides everything: if one string is among the values 
                    # of a property, all values become type "str"
                    field_types[prop] = prop_type 
            else:
                field_types[prop] = prop_type
                
    if not field_types:
        raise NoFieldTypes()
        
    return field_types


def logp_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if IPYTHON:
        display_png(mol)
    logp = Desc.MolLogP(mol)
    return logp


def show_record(mol):
    if IPYTHON:
        display_png(mol)

    for field in mol.GetPropNames():
        print("  {:13s}: {}".format(field, mol.GetProp(field)))

    print("_" * 75)


def show_sdf(sdf_list, force=False):
    l = len(sdf_list)
    if l > 20 and not force:
        print("  * will not show more than 20 records.")
        print("    to force, use force=True.")
    else:
        for mol in sdf_list:
            show_record(mol)


def merge_prop_from_file(sdf_list, fn, prop):
    lines_in_file = 0
    counter_sdf = 0

    firstline = True
    f_in = open(fn, "rb")
    dr = csv.reader(f_in, delimiter="\t")

    for row in dr:
        if firstline:
            headers = row

            firstline = False
            index_of = headers.index(prop)

        else:
            lines_in_file += 1
            for mol in sdf_list:
                if mol.HasProp(prop) and mol.GetProp(prop) == row[index_of]:
                    counter_sdf += 1
                    for index, new_prop in enumerate(row):
                        if index != index_of:
                            mol.SetProp(headers[index], row[index])

    print("  > {} lines from {} parsed. {} records modified in sdf.".format(lines_in_file-1, fn, counter_sdf))


def remove_props_from_mol(mol, prop_or_propslist):
    if not isinstance(prop_or_propslist, list):
        prop_or_propslist = [prop_or_propslist]
    for prop in prop_or_propslist:
        if prop in mol.GetPropNames():
            mol.ClearProp(prop)


def remove_props(mol_or_sdf_list, props):
    if isinstance(mol_or_sdf_list, file_type):
        print("  * operation not supported for file objects.")
        return

    if isinstance(mol_or_sdf_list, list):
        for mol in mol_or_sdf_list:
            if mol:
                remove_props_from_mol(mol, props)
    else:
        remove_props_from_mol(mol_or_sdf_list, props)


def rename_prop_in_mol(mol, old_prop, new_prop):
    if old_prop in mol.GetPropNames():
        value = mol.GetProp(old_prop)
        mol.SetProp(new_prop, value)
        mol.ClearProp(old_prop)


def rename_prop(mol_or_sdf_list, old_prop, new_prop):
    if isinstance(mol_or_sdf_list, list):
        for mol in mol_or_sdf_list:
            rename_prop_in_mol(mol, old_prop, new_prop)
    else:
        rename_prop_in_mol(mol_or_sdf_list, old_prop, new_prop)


def calc_props_in_mol(mol, dateprop="date", include_date=True, force2d=False):

    if force2d:
        mol.Compute2DCoords()
    else:
        try:
            mol.GetConformer()
        except ValueError: # no 2D coords... calculate them
            mol.Compute2DCoords()

    mol.SetProp("molwt", "{:.2f}".format(Desc.MolWt(mol)))
    mol.SetProp("formula", Chem.CalcMolFormula(mol))
    mol.SetProp("logp", "{:.2f}".format(Desc.MolLogP(mol)))
    mol.SetProp("hba", str(Desc.NOCount(mol)))
    mol.SetProp("hbd", str(Desc.NHOHCount(mol)))
    mol.SetProp("rotb", str(Desc.NumRotatableBonds(mol)))
    mol.SetProp("tpsa", str(int(Desc.TPSA(mol))))
    if include_date and not dateprop in mol.GetPropNames():
        mol.SetProp(dateprop, strftime("%Y%m%d"))


def get_highest_counter(mol_or_sdf, counterprop="molid"):
    # get highest counter in sdf
    molid_counter = 0

    for mol in mol_or_sdf:
        if counterprop in mol.GetPropNames():
            value = int(mol.GetProp(counterprop))
            if value > molid_counter:
                molid_counter = value

    return molid_counter


def calc_props(mol_or_sdf, counterprop="molid", dateprop="date",
               include_date=False, force2d=False):
    if not isinstance(mol_or_sdf, list):
        calc_props_in_mol(mol_or_sdf, dateprop, include_date, force2d)
        return

    molid_counter = get_highest_counter(mol_or_sdf, counterprop=counterprop) + 1

    for mol in mol_or_sdf:
        # continue counting if counterprop not present
        if not counterprop in mol.GetPropNames():
            mol.SetProp(counterprop, str(molid_counter))
            molid_counter += 1

        calc_props_in_mol(mol, dateprop, include_date, force2d)


def _key_get_prop(mol, field):
    try:
        val = float(mol.GetProp(field))
    except ValueError: # GetProp value could not be converted to float
        val = mol.GetProp(field)
    except KeyError:   # field is not present in the mol properties
        val = 10000000.0
    return val


def sort_sdf(sdf_list, field, reverse=True):
    sdf_list.sort(key=lambda x: _key_get_prop(x, field), reverse=reverse)


def activity_hist(sdf_list_or_file, activityprop):
    hist = Counter()
    act_oor = "OOR    (<0)"
    act_high = "high   (    >=50)"
    act_med = "medium (20 - <50)"
    act_low = "low    ( 0 - <20)"

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return

    for mol_counter_in, mol in enumerate(sdf_list_or_file):
        try:
            value = float(mol.GetProp(activityprop))
        except:
            hist["NaN"] += 1
            continue
        if value >= 50:
            hist[act_high] += 1
        elif value >= 20:
            hist[act_med] += 1
        elif value >= 0:
            hist[act_low] += 1
        else:
            hist[act_oor] += 1

        if mol_counter_in > 0 and mol_counter_in % 500 == 0:
            print("  > processed: {:7d}\r".format(mol_counter_in, end=""))
            sys.stdout.flush()

    mol_counter_in += 1
    print("  > processed: {:7d}".format(mol_counter_in))
    print("    done.")

    act_counter = 0
    for key in [act_high, act_med, act_low]:
        if hist.has_key(key):
            act_counter += hist[key]

    for key in [act_high, act_med, act_low]:
        if hist.has_key(key):
            perc_act = 100.0 * hist[key] / act_counter
            print("  {:<18s}: {:6d}  ({:4.1f}%)".format(key, hist[key], perc_act))

    for key in [act_oor, "NaN"]:
        if hist.has_key(key):
            print("  {:<18s}: {:6d}".format(key, hist[key]))


def factsearch(sdf_list_or_file, query, invert=False, max_hits=2000, count_only=False, sorted=True, reverse=True, field_types=None):
    result_list = []
    result_indexes_list = []
    mol_counter_out = 0
    
    if not field_types:
        field_types = get_field_types(sdf_list_or_file)
    
    if not field_types:
        print("  # no field type information available! -aborted.")
        return None

    field = None
    for el in query.split():
        if el in field_types:
            field = el
            break

    if not field:
        print("  # field could be extracted from query! -aborted.")
        return None
    
    print("  > field {} extracted from query: {}.".format(field, query))

    query_mod = query.replace(field, "val")

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return None

    print("  > searching...")
    for mol_counter_in, mol in enumerate(sdf_list_or_file):
        if not mol:
            continue
        hit = False
        if field in mol.GetPropNames():
            val = mol.GetProp(field).lower()
            if field_types[field] in ["number", "key"]:
                try:
                    val_float = float(val)

                except ValueError:
                    continue

                val_int = int(val_float)
                if val_int == val_float:
                    val = val_int
                else:
                    val = val_float

            if eval(query_mod):
                hit = True

            if invert:
                hit = not hit

            if hit:
                mol_counter_out += 1
                if not count_only:
                    result_list.append(mol)
                    result_indexes_list.append(mol_counter_in)

        if mol_counter_in > 0 and mol_counter_in % 500 == 0:
            print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in, mol_counter_out), end="")
            sys.stdout.flush()

        if not count_only and mol_counter_out >= max_hits:
            print()
            print("  * maximum number of hits reached.")
            break

    print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in+1, mol_counter_out))
    print("    done.")
    if count_only:
        print("  > option <count_only> was selected.")
        print("    no results were returned.")

    if sorted:
        sort_sdf(result_list, field, reverse=reverse)

    if 0 < mol_counter_out < 6:
        print("  > found {} matching records at {}.".format(mol_counter_out, str(result_indexes_list)))

    return result_list


def substruct_search(sdf_list_or_file, smarts, invert=False, max_hits=5000, count_only=False, add_h=False):
    result_list = []
    result_indexes_list = []

    mol_counter_out = 0
    query = Chem.MolFromSmarts(smarts)
    if not query:
        print("  * ERROR: could not generate query from SMARTS.")
        return None

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return None

    if not add_h and "[H]" in smarts:
        add_h = True
        print("  > explicit hydrogens turned on (add_h = True)")

    print("  > searching...")
    for mol_counter_in, mol in enumerate(sdf_list_or_file):
        if not mol: continue

        hit = False
        if add_h:
            mol_with_h = Chem.AddHs(mol)
            if mol_with_h.HasSubstructMatch(query):
                hit = True

        else:
            if mol.HasSubstructMatch(query):
                hit = True

        if invert:
            # reverse logic
            hit = not hit

        if hit:
            mol_counter_out += 1
            if not count_only:
                result_list.append(mol)
                result_indexes_list.append(mol_counter_in)

        if mol_counter_in > 0 and mol_counter_in % 500 == 0:
            print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in, mol_counter_out), end="")
            sys.stdout.flush()

        if not count_only and mol_counter_out >= max_hits:
            print()
            print("  * maximum number of hits reached.")
            break

    print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in+1, mol_counter_out))
    print("    done.")
    if count_only:
        print("  > option <count_only> was selected.")
        print("    no results were returned.")

    if 0 < mol_counter_out < 6:
        print("  > found {} matching records at {}.".format(mol_counter_out, str(result_indexes_list)))

    return Mol_List(result_list)


def similarity_search(sdf_list_or_file, smarts, similarity=0.8, max_hits=2000, count_only=False):
    result_list = []
    result_indexes_list = []

    mol_counter_out = 0
    query_mol = Chem.MolFromSmarts(smarts)
    if not query_mol:
        print("  * ERROR: could not generate query from SMARTS.")
        return None

    query_fp = FingerprintMols.FingerprintMol(query_mol)

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return None

    print("  > searching...")
    for mol_counter_in, mol in enumerate(sdf_list_or_file):

        mol_fp = FingerprintMols.FingerprintMol(mol)
        sim = DataStructs.FingerprintSimilarity(query_fp, mol_fp)
        if sim >= similarity:
            mol_counter_out += 1
            if not count_only:
                mol.SetProp("sim", "{:4.3f}".format(sim))
                result_list.append(mol)
                result_indexes_list.append(mol_counter_in)

        if mol_counter_in > 0 and mol_counter_in % 500 == 0:
            print( "\r    processed: {:7d}   found: {:6d}".format(mol_counter_in, mol_counter_out), end="")
            sys.stdout.flush()

        if not count_only and mol_counter_out >= max_hits:
            print()
            print("  * maximum number of hits reached.")
            break

    print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in+1, mol_counter_out))
    print("    done.")
    if count_only:
        print("  > option <count_only> was selected.")
        print("    no results were returned.")
    print

    if 0 < mol_counter_out < 6:
        print("  > found {} matching records at {}.".format(mol_counter_out, str(result_indexes_list)))

    return result_list



def similarity_hist(sdf_list_or_file, smarts, min_similarity=0.5):
    """use this to get a quick overview of the similarity distribution for a given compound (smarts)"""
    result_list = []

    mol_counter_out = 0
    query_mol = Chem.MolFromSmarts(smarts)
    if not query_mol:
        print("  * ERROR: could not generate query from SMARTS.")
        return None

    query_fp = FingerprintMols.FingerprintMol(query_mol)

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return None


    print("  > searching...")
    for mol_counter_in, mol in enumerate(sdf_list_or_file):

        mol_fp = FingerprintMols.FingerprintMol(mol)
        sim = DataStructs.FingerprintSimilarity(query_fp, mol_fp)
        if sim >= min_similarity:
            mol_counter_out += 1
            result_list.append(sim)

        if mol_counter_in > 0 and mol_counter_in % 500 == 0:
            print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in, mol_counter_out), end="")
            sys.stdout.flush()

    print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in+1, mol_counter_out))
    print("    done.")

    return result_list


def similarity_list(sdf_list_or_file, smarts):
    """return the molids and their similarity to the given compound (smarts)
    returns list([molid, similarity])"""
    result_list = []

    query_mol = Chem.MolFromSmarts(smarts)
    if not query_mol:
        print("  * ERROR: could not generate query from SMARTS.")
        return None

    query_fp = FingerprintMols.FingerprintMol(query_mol)

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return None

    print("  > processing...")
    for mol_counter_in, mol in enumerate(sdf_list_or_file):

        mol_fp = FingerprintMols.FingerprintMol(mol)
        sim = DataStructs.FingerprintSimilarity(query_fp, mol_fp)
        result_list.append([mol.GetProp("molid"), sim])

        if mol_counter_in > 0 and mol_counter_in % 500 == 0:
            print("\r    processed: {:7d}".format(mol_counter_in, ), end="")
            sys.stdout.flush()

    print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in+1))
    print("    done.")

    return result_list


def get_num_props_dict(sdf_list, fields=None, id_prop="molid", field_types=None):
    print("  > getting numerical fields...")
    DEBUG = False
    props_dict = {}
    molid_list = [] # list of molid for picking and structure display
    
    if not field_types:
        field_types = get_field_types(sdf_list)
        
    # only return the numerical database fields:
    num_props_list = [prop for prop in field_types if field_types[prop] == "number"]


    if fields == None or not isinstance(fields, list):
        props_list = num_props_list
    else:
        # make sure the <fields> are valid fields in the sdf
        props_list = list(set(fields).intersection(set(num_props_list)))

    # initialize dict
    for prop in props_list:
        props_dict[prop] = []

    if DEBUG:
        print("    ", props_list)
        DEBUG = False

    for mol in sdf_list:
        molid = int(mol.GetProp(id_prop))
        molid_list.append(molid)
        for prop in props_list:
            if mol.HasProp(prop):
                val = mol.GetProp(prop)

                try:
                    val = float(val)
                    props_dict[prop].append(val)

                except ValueError:
                    # prop is present in this record, but it is not a number
                    # insert MISSING_VAL
                    props_dict[prop].append(MISSING_VAL)

            else:
                # property is defined in sdf, but not in this record:
                # insert a MSSING_VAL to ensure identical length of the data lists
                props_dict[prop].append(MISSING_VAL)

    return props_dict, molid_list


def get_str_list(sdf_list, field):
    str_list = []

    for mol in sdf_list:
        if field in mol.GetPropNames():
            str_list.append(mol.GetProp(field))
        else:
            str_list.append("n.d.")

    return str_list


def remove_missing_values(l1, l2, l3, l4=None):
    # remove the missing values in one list and the corresponding data
                # at the same pos from the other list:
    while MISSING_VAL in l1:
        pos = l1.index(MISSING_VAL)
        l1.pop(pos)
        l2.pop(pos)
        l3.pop(pos)
        if l4:
            l4.pop(pos)

    # rinse and repeat for the other list
    while MISSING_VAL in l2:
        pos = l2.index(MISSING_VAL)
        l2.pop(pos)
        l1.pop(pos)
        l3.pop(pos)
        if l4:
            l4.pop(pos)

    if len(l1) != len(l2):
        print("  * ERROR: different length of value lists after removal of MISSING_VAL !")
    if len(l1) != len(l3):
        print("  * ERROR: different length of value list and molid list after removal of MISSING_VAL !")
    if l4 and len(l1) != len(l4):
        print("  * ERROR: different length of value list and color list after removal of MISSING_VAL !")


def show_hist(sdf_list, fields=None, show=True, savefig=True, id_prop="molid", field_types=None):
    """if fields==None, show histograms of all available fields in the sdf,
    otherwise use the supplied list of fields, e.g.
    fields=["pic50, "hba"]"""

    if not PYLAB:
        print("  * show_hist does not work because pylab could not be imported.")
        return

    if not isinstance(sdf_list, list):
        print("  * ERROR: plots are currently only implemented for sdf lists (no files).")
        return

    props_dict, molid_list = get_num_props_dict(sdf_list, fields, id_prop, field_types)

    num_of_plots = len(props_dict)
    num_rows = int(math.sqrt(num_of_plots))
    num_cols = num_of_plots // num_rows
    if num_of_plots % num_rows > 0:
        num_cols += 1

    if not savefig:
        pylab.figure(figsize=(3, 2))
    else:
        pylab.figure()

    pylab.subplots_adjust(hspace=0.3, wspace=0.2)
    for counter, prop in enumerate(props_dict):
        value_list = props_dict[prop]

        # remove the missing values:
        while MISSING_VAL in value_list:
            value_list.remove(MISSING_VAL)

        pylab.subplot(num_rows, num_cols, counter+1)
        pylab.hist(value_list)
        pylab.title(prop)

    if savefig:
        pylab.savefig("histogram.png")
    if show:
        pylab.show()


def show_scattermatrix(sdf_list, fields=None, colorby=None, mode="show", id_prop="molid", field_types=None):
    """if fields==None, show a scattermatrix of all available fields in the sdf,
    otherwise use the supplied list of fields, e.g.
    fields=["pic50, "hba"]"""

    if not PYLAB:
        print("  * show_scattermatrix does not work because pylab could not be imported.")
        return

    if not isinstance(sdf_list, list):
        print("  * ERROR: plots are currently only implemented for sdf lists (no files).")
        return

    props_dict, molid_list = get_num_props_dict(sdf_list, fields, id_prop=id_prop, field_types=field_types)

    props_dict_len = len(props_dict)
    num_of_plots = props_dict_len ** 2
    num_rows = props_dict_len
    num_cols = props_dict_len
    if not colorby:
        COLOR_IS_LIST = False
    else:
        avail_colors = ["blue", "green", "red", "cyan", "magenta", "yellow", "black"]
        str_list = get_str_list(sdf_list, colorby)
        str_set = list(set(str_list))
        if len(str_set) > len(avail_colors):
            print("  * not enough colors!")
            COLOR_IS_LIST = False
        else:
            # convert the string values to colors
            COLOR_IS_LIST = True
            color_list_all = [avail_colors[str_set.index(i)] for i in str_list]
            color_coding = ""
            for ind, item in enumerate(str_set):
                color_coding = color_coding +  "{}: {}  ".format(item, avail_colors[ind])
            print("  > color coding:", color_coding)

    props_list = props_dict.keys() # to have a defined order
    num_of_props = len(props_list)
    fig = pylab.figure()
    pylab.subplots_adjust(hspace=0, wspace=0)
    axes = []
    axes_molid_dict = {} # track the displayed molides for each plot (axe)
    for y_count, prop_y in enumerate(props_list):
        for x_count, prop_x in enumerate(props_list):
            if y_count == 0: # first row
                if x_count == 0 or x_count == 1: # top left, don't let the top left hist share the y axis
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1)
                else:
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1, sharey=axes[1])
            else:
                if x_count == 0:
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1, sharex=axes[x_count])
                elif x_count == y_count:
                    # don't let the histograms share the y axis
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1, sharex=axes[x_count])
                else:
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1, sharex=axes[x_count], sharey=axes[y_count*num_of_props])

            if prop_x == prop_y:
                value_list = props_dict[prop_x][:]

               # remove the missing values:
                while MISSING_VAL in value_list:
                    value_list.remove(MISSING_VAL)

                ax.hist(value_list)
            else:
                value_x_list = props_dict[prop_x][:]
                value_y_list = props_dict[prop_y][:]
                ax_molid_list = molid_list[:]

                if len(value_x_list) != len(value_y_list):
                    print("  * ERROR: different length of value lists !")
                    return

                if COLOR_IS_LIST:
                    color_list = color_list_all[:]
                    remove_missing_values(value_x_list, value_y_list, ax_molid_list, color_list)
                else:
                    color_list = "r"
                    remove_missing_values(value_x_list, value_y_list, ax_molid_list)

                if len(value_x_list) != len(value_y_list):
                    print(" * ERROR: different length of value lists after removal of MISSING_VAL !")
                    return
                if COLOR_IS_LIST and len(value_x_list) != len(color_list):
                    print("  * ERROR: different length of value list and color list after removal of MISSING_VAL !")
                    return

                axes_molid_dict[ax] = ax_molid_list
                ax.scatter(value_x_list, value_y_list, s=POINTSIZE, color=color_list, picker=5)

            if x_count % 2 == 0:
                if y_count == 0:
                    ax.xaxis.set_label_position("top")
                    pylab.xlabel(prop_x)
                ax.xaxis.set_ticks_position("bottom")
                pylab.xticks(fontsize="small", visible=(y_count==num_of_props-1))
            else:
                if y_count == num_of_props-1:
                    ax.xaxis.set_label_position("bottom")
                    pylab.xlabel(prop_x)
                ax.xaxis.set_ticks_position("top")
                pylab.xticks(fontsize="small", visible=(y_count==0))
            if y_count % 2 == 0:
                if x_count == 0:
                    ax.yaxis.set_label_position("left")
                    pylab.ylabel(prop_y)
                ax.yaxis.set_ticks_position("right")
                pylab.yticks(fontsize="small", visible=(x_count==num_of_props-1 and y_count != num_of_props-1))
            else:
                if x_count == num_of_props-1:
                    ax.yaxis.set_label_position("right")
                    pylab.ylabel(prop_y)
                ax.yaxis.set_ticks_position("left")
                pylab.yticks(fontsize="small", visible=(x_count==0))

            axes.append(ax)

    if COLOR_IS_LIST:
        pylab.suptitle(color_coding)
        # pylab.text(.5, 0, color_coding, horizontalalignment='center')
    pylab.savefig("scatter.png")
    if mode == "show":
        # fig.show()
        pylab.show()
    elif mode == "gui":
        return fig, axes_molid_dict


def show_scattermatrix2(sdf_list, sdf_base, fields=None, mode="show", id_prop="molid", field_types=None):
    """if fields==None, show a scattermatrix of all available fields in the sdf,
    otherwise use the supplied list of fields, e.g.
    fields=["pic50, "hba"]"""

    if not PYLAB:
        print("  * show_scattermatrix does not work because pylab could not be imported.")
        return

    if not isinstance(sdf_list, list):
        print("  * ERROR: plots are currently only implemented for sdf lists (no files).")
        return

    props_dict, molid_list = get_num_props_dict(sdf_list, fields, id_prop=id_prop, field_types=field_types)
    props_dict_base, molid_list_base = get_num_props_dict(sdf_base, fields, id_prop=id_prop, field_types=field_types)

    props_dict_len = len(props_dict)
    num_rows = props_dict_len
    num_cols = props_dict_len

    fig = pylab.figure()
    pylab.subplots_adjust(hspace=0, wspace=0)
    props_list = props_dict.keys() # to have a defined order
    num_of_props = len(props_list)
    axes = []
    axes_molid_dict = {} # track the displayed molides for each plot (axe)
    for y_count, prop_y in enumerate(props_list):
        for x_count, prop_x in enumerate(props_list):
            if y_count == 0: # first row
                if x_count == 0 or x_count == 1: # don't let the top left hist share the y axis
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1)
                else:
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1, sharey=axes[1])
            else:
                if x_count == 0:
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1, sharex=axes[x_count])
                elif x_count == y_count:
                    # don't let the histograms share the y axis
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1, sharex=axes[x_count])
                else:
                    ax = fig.add_subplot(num_rows, num_cols, y_count*num_of_props+x_count+1, sharex=axes[x_count], sharey=axes[y_count*num_of_props])


            if prop_x == prop_y:
                value_list = props_dict[prop_x][:]
                value_list_base = props_dict_base[prop_x][:]

               # remove the missing values:
                while MISSING_VAL in value_list:
                    value_list.remove(MISSING_VAL)

                while MISSING_VAL in value_list_base:
                    value_list_base.remove(MISSING_VAL)

                ax.hist(value_list_base)
                ax.hist(value_list, color="r")
                # pylab.yticks(visible=False)

            else:

                value_x_list = props_dict[prop_x][:]
                value_y_list = props_dict[prop_y][:]
                ax_molid_list = molid_list[:]
                value_x_list_base = props_dict_base[prop_x][:]
                value_y_list_base = props_dict_base[prop_y][:]
                ax_molid_list_base = molid_list_base[:]

                remove_missing_values(value_x_list, value_y_list, ax_molid_list)
                remove_missing_values(value_x_list_base, value_y_list_base, ax_molid_list_base)

                if len(value_x_list) != len(value_y_list):
                    print("  * ERROR: different length of value lists !")

                if len(value_x_list_base) != len(value_y_list_base):
                    print("  * ERROR: different length of base value lists !")

                color_list = ["b"] * len(value_x_list_base) + ["r"] * len(value_x_list)

                axes_molid_dict[ax] = ax_molid_list_base + ax_molid_list
                ax.scatter(value_x_list_base+value_x_list, value_y_list_base+value_y_list, s=POINTSIZE, color=color_list, picker=5)
                # ax.scatter(value_x_list, value_y_list, s=POINTSIZE, color="r", picker=5)

            if x_count % 2 == 0:
                if y_count == 0:
                    ax.xaxis.set_label_position("top")
                    pylab.xlabel(prop_x)
                ax.xaxis.set_ticks_position("bottom")
                pylab.xticks(fontsize="small", visible=(y_count==num_of_props-1))
            else:
                if y_count == num_of_props-1:
                    ax.xaxis.set_label_position("bottom")
                    pylab.xlabel(prop_x)
                ax.xaxis.set_ticks_position("top")
                pylab.xticks(fontsize="small", visible=(y_count==0))
            if y_count % 2 == 0:
                if x_count == 0:
                    ax.yaxis.set_label_position("left")
                    pylab.ylabel(prop_y)
                ax.yaxis.set_ticks_position("right")
                pylab.yticks(fontsize="small", visible=(x_count==num_of_props-1 and y_count != num_of_props-1))
            else:
                if x_count == num_of_props-1:
                    ax.yaxis.set_label_position("right")
                    pylab.ylabel(prop_y)
                ax.yaxis.set_ticks_position("left")
                pylab.yticks(fontsize="small", visible=(x_count==0))

            axes.append(ax)

            # pylab.title(prop[2:])

    pylab.savefig("scatter.png")
    if mode == "show":
        pylab.show()
    elif mode == "gui":
        return fig, axes_molid_dict


def cluster_from_sdf_list(sdf_list, cutoff=0.2):

    counter = Counter()

    # generate the fingerprints
    fp_list = [Chem.GetMorganFingerprintAsBitVect(mol, 3, 1024) for mol in sdf_list]

    # second generate the distance matrix:
    dists = []
    num_of_fps = len(fp_list)
    for i in range(1, num_of_fps):
        sims = DataStructs.BulkTanimotoSimilarity(fp_list[i],fp_list[:i])
        dists.extend([1-x for x in sims])

    # now cluster the data:
    cluster_idx_list = Butina.ClusterData(dists, num_of_fps, cutoff, isDistData=True)
    for cluster in cluster_idx_list:
        counter[len(cluster)] += 1
    print("    clustersize  num_of_clusters")
    print("    ===========  ===============")
    for length in sorted(counter.keys(), reverse=True):
        print("        {:4d}            {:3d}".format(length, counter[length]))

    # return sorted by cluster length
    return sorted(cluster_idx_list, key=len, reverse=True)


def get_sdf_from_index_list(orig_sdf, index_list):
    """generate sdf_lists after clustering"""
    cluster_sdf = [orig_sdf[x] for x in index_list]
    return cluster_sdf


def get_sdf_from_id_list(sdf_list_or_file, id_dict_or_list, calc_ex_mw=True):
    result_sdf = []
    result_id_list = []

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return None

    if isinstance(id_dict_or_list, dict):
        add_props = True
        id_list = id_dict_or_list.keys()
    else:
        add_props = False
        id_list = id_dict_or_list

    print("  > searching...")
    for mol_counter_in, mol in enumerate(sdf_list_or_file):
        molid = int(mol.GetProp("molid"))
        if molid in id_dict_or_list:
            if calc_ex_mw:
                exmw = Desc.ExactMolWt(mol)
                mol.SetProp("exmolwt", "{:.2f}".format(exmw))
            if add_props:
                mol.SetProp("s_pos", id_dict_or_list[molid])

            result_sdf.append(mol)
            result_id_list.append(molid)

        if mol_counter_in > 0 and mol_counter_in % 500 == 0:
            print("\r    processed: {:7d}   found: {:6d}".format(mol_counter_in, len(result_id_list)), end="")
            sys.stdout.flush()

    not_found_list = list(set(id_list).difference(result_id_list))

    print("\r  > result sdf generated with {} out of {} input records.".format(len(id_list), len(result_id_list)))

    return result_sdf, not_found_list


def get_max_act_in_cluster(orig_sdf, cluster, act_prop):
    max_act = -1000
    for x in cluster:
        mol = orig_sdf[x]
        try:
            value = float(mol.GetProp(act_prop))
        except:
            print("  * molecule at index {:6d} has not activity prop {}".format(x, act_prop))
            continue
        if value > max_act:
            max_act = value
    return max_act


def get_med_act_in_cluster(orig_sdf, cluster, act_prop):
    med_act = 0.0
    num_of_members = 0
    for x in cluster:
        mol = orig_sdf[x]
        try:
            value = float(mol.GetProp(act_prop))
        except:
            print("  * molecule at index {:6d} has not activity prop {}".format(x, act_prop))
            continue
        med_act += value
        num_of_members += 1

    if num_of_members == 0:
        num_of_members = 1

    return med_act / num_of_members


def analyze_cluster_idx_list(orig_sdf, cluster_idx_list, mode="remove_singletons", act_prop=None):
    """available modes:
       remove_singletons: remove clusters with only one member
       ind_activity:      sort the clusters by the member with the highest activity
       med_activity:      sort the cluster by their medium activity

       a NEW cluster_idx_list is returned, use get_sdf_from_index_list to generate a sdf from it"""

    if act_prop:
        print("  > sorting cluster members on {}...".format(act_prop), end="")
        cluster_idx_list_sorted = []
        for cluster in cluster_idx_list:
            cluster_dict = {}

            # map the position in the orig_sdf to the molid
            sdf_list = get_sdf_from_index_list(orig_sdf, cluster)
            for pos, idx in enumerate(cluster):
                mol = sdf_list[pos]
                cluster_dict[int(mol.GetProp("molid"))] = idx
            sort_sdf(sdf_list, act_prop, reverse=True)
            cluster_sorted = [cluster_dict[int(mol.GetProp("molid"))] for mol in sdf_list]
            cluster_idx_list_sorted.append(cluster_sorted)
        cluster_idx_list = cluster_idx_list_sorted
        print(" done.")

    if mode == "remove_singletons":
        new_idx_list = [x for x in cluster_idx_list if len(x) > 1]
        print("  > removed {} singletons from cluster list".format(len(cluster_idx_list) - len(new_idx_list)))
        print("    the resulting list has {} clusters".format(len(new_idx_list)))
        return new_idx_list

    if mode == "ind_activity":
        new_idx_list = sorted(cluster_idx_list, key=lambda x: get_max_act_in_cluster(orig_sdf, x, act_prop), reverse=True)
        print("  > max ind_activity and number of members in the first ten clusters:")
        print("    index  ind_act  #members")
        print("    =====  =======  ========")
        for cl_counter, cluster in enumerate(new_idx_list):
            print("      {:2d}   {:6.1f}      {:3d}".format(cl_counter, get_max_act_in_cluster(orig_sdf, cluster, act_prop), len(cluster)))
            if cl_counter >= 9:
                break
        return new_idx_list

    if mode == "med_activity":
        new_idx_list = sorted(cluster_idx_list, key=lambda x: get_med_act_in_cluster(orig_sdf, x, act_prop), reverse=True)
        print("  > max med_activity in the first ten clusters:")
        for cl_counter, cluster in enumerate(new_idx_list):
            print("{}: {}   ".format(cl_counter, get_med_act_in_cluster(orig_sdf, cluster, act_prop)), end="")
            if cl_counter >= 9:
                break
        print()
        return new_idx_list

    print("  * unsupported mode.")


def write_clusters_as_sdf(orig_sdf, cluster_idx_list, basename="cluster"):
    """write every cluster in the index list as individual sdf file"""
    basename = basename.split(".")[0]

    for counter, cluster in enumerate(cluster_idx_list):
        num_of_clusters = len(cluster_idx_list)
        sdf_list = get_sdf_from_index_list(orig_sdf, cluster)
        write_sdf(sdf_list, "{}_{:03d}_{:03d}.sdf".format(basename, counter, num_of_clusters-1))


def write_cluster_report(orig_sdf, cluster_idx_list, captionprop, reportname="cluster"):
    intro = "<html>\n<head>\n  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n  <link rel=\"stylesheet\" type=\"text/css\" href=\"css/style.css\" />\n  <title>ClusterViewer</title>\n</head>\n<body>\n<table width=\"\" cellspacing=\"1\" cellpadding=\"1\" border=\"1\" align=\"\" height=\"60\" summary=\"\">\n<tbody>"
    extro = "</tbody>\n</table>\n</body>\n</html>"
    filename = op.join(REPORT_FOLDER, "html", reportname + ".htm")
    f = open(filename, "w")
    f.write(intro)
    for counter, cluster in enumerate(cluster_idx_list):
        num_of_clusters = len(cluster_idx_list)
        sdf_list = get_sdf_from_index_list(orig_sdf, cluster)
        img = Draw.MolsToGridImage(sdf_list, molsPerRow=4, legends=[m.GetProp(captionprop) for m in sdf_list])
        img_filename = "img/{}_{:03d}.png".format(reportname, counter)
        img.save("html/"+img_filename)
        hist_filename = "img/{}_{:03d}_hist.png".format(reportname, counter)
        show_hist(sdf_list, fields=[captionprop], show=False, savefig=False)
        pylab.savefig("html/"+hist_filename)
        f.write("  <tr><td>\n<br><b>Cluster {:03d}:</b></td><td></td></tr>\n<tr><td><img src=\"{}\" alt=\"icon\" /></td><td><img src=\"{}\" alt=\"icon\" /></td></tr>\n".format(counter, img_filename, hist_filename))

    f.write(extro)
    f.close()


def neutralize_mol(mol, pattern=None, id_prop="molid", show=False):
    if not pattern:
        pattern = (
            # Imidazoles
            ('[n+;H]','n'),
            # Amines
            ('[N+;!H0]','N'),
            # Carboxylic acids and alcohols
            ('[$([O-]);!$([O-][#7])]','O'),
            # Thiols
            ('[S-;X1]','S'),
            # Sulfonamides
            ('[$([N-;X2]S(=O)=O)]','N'),
            # Enamines
            ('[$([N-;X2][C,N]=C)]','N'),
            # Tetrazoles
            ('[n-]','[nH]'),
            # Sulfoxides
            ('[$([S-]=O)]','S'),
            # Amides
            ('[$([N-]C=O)]','N'),
            )

        reactions = [(Chem.MolFromSmarts(x), Chem.MolFromSmiles(y, False)) for x, y in pattern]

        old_smiles = Chem.MolToSmiles(mol)
        replaced = False
        for reactant, product in reactions:
            while mol.HasSubstructMatch(reactant):
                replaced = True
                rms = Chem.ReplaceSubstructs(mol, reactant, product)
                mol = rms[0]

        if replaced:
            new_smiles = Chem.MolToSmiles(mol)
            mol_props = list(mol.GetPropNames())
            props_dict = {prop: mol.GetProp(prop) for prop in mol_props}
            mol = Chem.MolFromSmiles(new_smiles)
            for prop in props_dict:
                mol.SetProp(prop, props_dict[prop])

            calc_props_in_mol(mol, include_date=False)

            if show:
                if mol.HasProp(id_prop):
                    molid = mol.GetProp(id_prop)
                else:
                    molid = ""
                print("    {:5d}:  {:15s} --> {:15s}".format(molid, old_smiles, new_smiles))
                if IPYTHON:
                    display_png(mol)

        return mol, replaced


def neutralize_sdf(sdf_list_or_file, id_prop="molid", show=False):
    """returns:            neutral_sdf::list<mol>, neutralized_molids::list<int>
    neutral_sdf:         new sdf with all input mols, where the salts have been neutralized
    neutralized_molids: list with the neutralized molids"""
    neutral_sdf = []
    neutralized_molids = []
    counter_in = 0
    counter_out = 0

    if not isinstance(sdf_list_or_file, list) and sdf_list_or_file.atEnd(): # sdf is file
        print("  * file object is at end, please reload.")
        return None

    print("  > neutralizing...")
    for mol in sdf_list_or_file:
        counter_in += 1
        new_mol, replaced = neutralize_mol(mol, id_prop=id_prop, show=show)
        if replaced:
            counter_out += 1
            neutralized_molids.append(int(mol.GetProp(id_prop)))

        neutral_sdf.append(new_mol)

    print("    neutralizing finished:")
    print("    in:          {:5d}".format(counter_in))
    print("    neutralized: {:5d}".format(counter_out))

    return neutral_sdf, neutralized_molids


def scatter_props(sdf_list, props_list):
    """Plot a scatter plot for exactly two properties, x and y.
    Currently cmdline only, not yet implemented in the GUI"""
    if len(props_list) != 2:
        print("  * function scatter_props has to be used with a props_list of length 2")
        return
    sdf_list_sel = [mol for mol in sdf_list if mol.HasProp(props_list[0]) and mol.HasProp(props_list[1])]
    x = [float(mol.GetProp(props_list[0])) for mol in sdf_list_sel]
    y = [float(mol.GetProp(props_list[1])) for mol in sdf_list_sel]
    pylab.scatter(x, y, s=40, color="r")


def mol_grid(sdf_list, props, fn="img/grid.png", mols_per_row=5, sub_img_size=(200, 200)):
    """Draw a molecule grid from the input <sdf_list>. On IPython, an inline graphics will be returned
    in addition to writing the image to <fn>.
    The given sdf <props> (as a list) will be concatenated to the molecules' legends."""

    if not isinstance(props, list):
        props = [props]

    legends = []
    for mol in sdf_list:
        leg = [mol.GetProp(prop) for prop in props]
        leg_str = "_".join(leg)
        legends.append(leg_str)

    img = Draw.MolsToGridImage(sdf_list, molsPerRow=mols_per_row, subImgSize=sub_img_size, legends=legends)
    img.save(fn)

    if IPYTHON:
        return img

