#!/usr/bin/env python2
# -*- coding: utf-8 -*-

from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw
#from rdkit.Chem.Descriptors import MolWt
#from rdkit.Chem import Lipinski as Lip

import sdf_tools as sdft

import numpy as np
import re
import os.path as op


try:
    from IPython.core.display import display_png
    IPYTHON = True
except:
    IPYTHON = False
    print("  * no inline display of molecule structures supported.")


def get_res_pos(smiles):
    pat = re.compile('\[(.*?)\*\]')
    pos_str = re.findall(pat, smiles)[0]
    if pos_str:
        return int(pos_str)
    else:
        return 0


def combine_frags():
    # TEST Version
    mol = ha[0]
    core_mol = Chem.MolFromSmiles("CCC(O)C(C)C(N)=O")
    tmp = Chem.ReplaceCore(mol, core_mol, labelByIndex=True)
    frag_mols = list(Chem.GetMolFrags(tmp, asMols=True))
    core_frag = Chem.ReplaceSidechains(mol, core_mol)
    display_png(core_frag)
    display_png(frag_mols[0])
    display_png(frag_mols[1])
    nmol = Chem.CombineMols(Chem.CombineMols(frag_mols[0], frag_mols[1]), core_frag)
    emol = Chem.EditableMol(nmol)

    for at in nmol.GetAtoms():
        if at.GetAtomicNum() == 0:
            emol.ReplaceAtom(at.GetIdx(), Chem.Atom(6))

    gen_mol = emol.GetMol()

    Chem.SanitizeMol(gen_mol)
    if IPYTHON:
        print("  try to display molecule:")
        display_png(gen_mol)



def generate_sar_table(db_list, core, id_prop, act_prop, dir_name="sar_table", color_prop="logp"):
    """core: smiles string; act_prop: string"""

    act_xy = np.zeros([30, 30], dtype=np.float)    # coordinates for the activity
    color_xy = np.zeros([30, 30], dtype=np.float)
    molid_xy = np.zeros([30, 30], dtype=np.int)
    # molid_xy = np.arange(900, dtype=np.int).reshape(30, 30)  # coordinates for the molid
    rx_dict = {}  # axes for the residues
    ry_dict = {}
    max_x = -1  # keep track of the arraysize
    max_y = -1
    res_pos_x = -1
    res_pos_y = -1

    core_mol = Chem.MolFromSmiles(core)
    Draw.MolToFile(core_mol, "%s/img/core.png" % dir_name, [80, 80])

    for idx, mol in enumerate(db_list):
        act = float(mol.GetProp(act_prop))
        color = float(mol.GetProp(color_prop))
        molid = int(mol.GetProp(id_prop))
        tmp = Chem.ReplaceCore(mol, core_mol, labelByIndex=True)
        frag_mols = list(Chem.GetMolFrags(tmp, asMols=True))
        frag_smiles = [Chem.MolToSmiles(m, True) for m in frag_mols]
        if len(frag_mols) == 1:
            # one of the two residues is H:
            pos = get_res_pos(frag_smiles[0])
            if pos == res_pos_x:
                h_smiles = "[%d*]([H])" % res_pos_y
                frag_smiles.append(h_smiles)
                frag_mols.append(Chem.MolFromSmiles(h_smiles))
            else:
                h_smiles = "[%d*]([H])" % res_pos_x
                frag_smiles.insert(0, h_smiles)
                frag_mols.insert(0, Chem.MolFromSmiles(h_smiles))

            print(" adding H residue in pos {} to  mol #{} (molid: {})".format(pos, idx, mol.GetProp("k_molid")))

        elif len(frag_mols) > 2:
            print("*  incorrect number of fragments ({}) in mol #{} (molid: {})".format(len(frag_mols), idx, mol.GetProp("k_molid")))
            continue

        if res_pos_x == -1:
            # print frag_smiles[0], frag_smiles[1]
            res_pos_x = get_res_pos(frag_smiles[0])
            res_pos_y = get_res_pos(frag_smiles[1])
            # print "res_pos_x: {}     res_pos_y: {}".format(res_pos_x, res_pos_y)
        else:
            test_pos_x = get_res_pos(frag_smiles[0])
            if test_pos_x != res_pos_x: # switch residues
                frag_smiles = frag_smiles[::-1]
                frag_mols = frag_mols[::-1]
        if frag_smiles[0] in rx_dict:
            curr_x = rx_dict[frag_smiles[0]]
        else:
            max_x += 1
            rx_dict[frag_smiles[0]] = max_x
            curr_x = max_x
            Draw.MolToFile(frag_mols[0], "%s/img/frag_x_%02d.png" % (dir_name, max_x), [100, 100])
        if frag_smiles[1] in ry_dict:
            curr_y = ry_dict[frag_smiles[1]]
        else:
            max_y += 1
            ry_dict[frag_smiles[1]] = max_y
            curr_y = max_y
            Draw.MolToFile(frag_mols[1], "%s/img/frag_y_%02d.png" % (dir_name, max_y), [100, 100])

        act_xy[curr_x][curr_y] = act
        color_xy[curr_x][curr_y] = color
        molid_xy[curr_x][curr_y] = molid

    return act_xy, molid_xy, color_xy, max_x, max_y


def sar_table_report_html(act_xy, molid_xy, color_xy, max_x, max_y, color_by="logp", show_link=False, show_tooltip=True):
    # logp_colors = {2.7: "#5F84FF", 3.0: "#A4D8FF", 4.2: "#66FF66", 5.0: "#FFFF66", 1000.0: "#FF4E4E"}
    logp_colors = {2.7: "#98C0FF", 3.0: "#BDF1FF", 4.2: "#AAFF9B", 5.0: "#F3FFBF", 1000.0: "#FF9E9E"}

    # write horizontal residues
    line = ["<table width=\"\" cellspacing=\"1\" cellpadding=\"1\" border=\"1\" align=\"center\" height=\"60\" summary=\"\">\n<tbody>"]
    line.append("  <tr><td align=\"left\">Core:<br><img src=\"img/core.png\" alt=\"icon\" /></td>")
    for curr_x in range(max_x+1):
        line.append("<td><img src=\"img/frag_x_%02d.png\" alt=\"icon\" /></td>" % curr_x)

    line.append("</tr>\n")

    for curr_y in range(max_y+1):
        line.append("<tr><td><img src=\"img/frag_y_%02d.png\" alt=\"icon\" /></td>" % curr_y)
        for curr_x in range(max_x+1):
            molid = molid_xy[curr_x][curr_y]
            if molid > 0:
                link_in = ""
                link_out = ""
                bg_color = " "
                tooltip = " "
                if show_link:
                    link = "../reports/ind_stock_results.htm#cpd_%05d" % molid
                    link_in = "<a href=\"%s\">" % link
                    link_out = "</a>"
                if color_by == "n_logp":
                    logp = color_xy[curr_x][curr_y]
                    if show_tooltip:
                        tooltip = ' title="LogP: %.2f"' % logp
                    for limit in sorted(logp_colors):
                        if logp <= limit:
                            bg_color = ' bgcolor="%s"' % logp_colors[limit]
                            break
                else:
                    if show_tooltip:
                        tooltip = ' title="%.2f"' % color_xy[curr_x][curr_y]

                line.append("<td%s align=\"center\"%s><b>%.2f</b><br>(%s%d%s)</td>" % (tooltip, bg_color, act_xy[curr_x][curr_y], link_in, molid_xy[curr_x][curr_y], link_out))

            else: #empty value in numpy array
                line.append("<td></td>")


        line.append("</tr>\n")

    if color_by == "n_logp":
        line.append("</tbody>\n</table>\n<p><p>LogP color coding:</p>\n<table width=\"\" cellspacing=\"1\" cellpadding=\"1\" border=\"1\" align=\"left\" height=\"40\" summary=\"\">\n<tbody><tr>\n")
        for limit in sorted(logp_colors):
            line.append('<td align="center" bgcolor="%s">&le; %.2f</td>' % (logp_colors[limit], limit))
        line.append("</tr>")

    line.append("</tbody>\n</table>\n")
    html_table = "".join(line)

    return html_table


def write_html_page(html_content, dir_name="sar_table", page_name="sar_table", page_title="SAR Table"):
    intro = "<html>\n<head>\n  <meta http-equiv=\"Content-Type\" content=\"text/html; charset=utf-8\" />\n  <link rel=\"stylesheet\" type=\"text/css\" href=\"css/style.css\" />\n  <title>%s</title>\n</head>\n<body>\n" % page_title
    extro = "</body>\n</html>"

    sdft.create_dir_if_not_exist(op.join(sdft.REPORT_FOLDER, dir_name))
    sdft.create_dir_if_not_exist(op.join(sdft.REPORT_FOLDER, dir_name, "img"))

    filename = op.join(sdft.REPORT_FOLDER, dir_name, "%s.htm" % page_name)
    f = open(filename, "w")
    f.write(intro)

    f.write(html_content)

    f.write(extro)
    f.close()



if __name__ == "__main__":

    GEN_TABLE_FOR = "hydroxyamides"  # hydroxyamides, clpp, all

    if GEN_TABLE_FOR == "hydroxyamides" or GEN_TABLE_FOR == "all":
        print("  * generating SAR table for hydroxyamides...")
        dir_name = "sar_table_ha"
        db = sdft.load_db("sdf/aviru_db.sdf")
        ha_sss = sdft.factsearch(sdft.substruct_search(db, "[H]OC(C)C(C)C(=O)N([H])[H]"), 's_class == "hydroxyamide"')
        print("    %3d structures from substructure search" % len(ha_sss))

        #keep only most active stereoisomer:
        ha_dict = {}
        prop = "n_pEC50_hem"
        replaced_by_list = [] # [original, replaced_by]
        for mol in ha_sss:
            if not mol.HasProp(prop):
                continue
            smiles = Chem.MolToSmiles(mol)
            if smiles in ha_dict:
                result = float(mol.GetProp(prop))
                curr_molid = mol.GetProp("k_molid")
                if result > float(ha_dict[smiles].GetProp(prop)):
                    replaced = ha_dict[smiles].GetProp("k_molid")
                    ha_dict[smiles] = mol
                    replaced_by_list.append([replaced, curr_molid])
                    print("        %3s replaced by %3s" % (replaced, curr_molid))
                else:
                    dropped_for = ha_dict[smiles].GetProp("k_molid")
                    replaced_by_list.append([curr_molid, dropped_for])
                    print("        %3s dropped for %3s" % (curr_molid, dropped_for))
            else:
                ha_dict[smiles] = mol

        ha = ha_dict.values()
        print("    %3d structures from most active isomer" % len(ha))

        sdft.sort_db(ha, prop)
        act, molid, color_xy, max_x, max_y = generate_sar_table(ha, "CC(O)C(C)C(N)=O", "n_pEC50_hem", dir_name=dir_name)
        write_html_page(sar_table_report_html(act, molid, color_xy, max_x, max_y),
                        dir_name=dir_name, page_name="sar_table", page_title="SAR Table HA")
        molid_elim = []
        for mol in ha_sss:
            molid = mol.GetProp("k_molid")
            found = False
            for ha_mol in ha:
                if ha_mol.GetProp("k_molid") == molid:
                    found = True
                    break
            if not found:
                molid_elim.append(molid)

        print("    eliminated molid:", molid_elim)
        print("    done.")

    elif GEN_TABLE_FOR == "clpp" or GEN_TABLE_FOR == "all":
        print("  * generating SAR table for clpp_145...")
        dir_name="sar_table_clpp"
#    db = sdft.load_db("sdf/aviru_db.sdf")
#    clpp_sss = sdft.substruct_search(db, "Cc1cocn1")
#    clpp = sdft.factsearch(clpp_sss, "n_pIC50_clpp > 0")
        clpp = sdft.load_db("sdf/clpp_145_sar_tbl.sdf")
        print("    %3d structures from substructure search" % len(clpp))
        act, molid, color_xy, max_x, max_y = generate_sar_table(clpp, "Cc1cocn1", "n_pIC50_clpp", dir_name=dir_name)
        write_html_page(sar_table_report_html(act, molid, color_xy, max_x, max_y),
                        dir_name=dir_name, page_name="sar_table", page_title="SAR Table ClpP")
        print("    done.")
