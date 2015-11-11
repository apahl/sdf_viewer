#!/usr/bin/env python
#-*- coding: utf-8 -*-
# sdf_viewer.py
# version: 2015-05-28
# author:  Axel Pahl (APL)
# contact: firstnamelastname at gmx dot de
# license: BSD, see license.txt in this folder

from __future__ import absolute_import, division, print_function # , unicode_literals

import sys
from rdkit.Chem import AllChem as Chem
from rdkit.Chem import Draw

try:
    Draw.DrawingOptions.atomLabelFontSize = 18

except:
    print("  * drawing options are not available.")

from PyQt4 import QtGui, QtCore
from gui.app import Ui_MainWindow

import sdf_tools as sdft

import pickle
import webbrowser
import os.path as op
import os

if sys.version_info[0] > 2:
    import io
    from io import BytesIO as IO
    print("  > runnning on Python3")
    PY3 = True
    file_type = io.IOBase
else:
    from cStringIO import StringIO as IO
    PY3 = False
    file_type = file


# try to import local options file from ~/.sdf_viewer/ by adding the folder to sys.path:
folder = op.join(op.expanduser("~"), ".sdf_viewer")
if op.isdir(folder):
    if not folder in sys.path:
        sys.path.insert(0, folder)

import sdf_viewer_config as config

sys.stdout.flush() # flush the printing buffer...


class App(QtGui.QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super(App, self).__init__(parent)
        self.setupUi(self)

        self.FIRSTRUN = True
        self.SDFTYPE = ""
        self.le_sdf_name.setText(config.STARTUP_SDF)
        self.highlight_dict = config.HIGHLIGHT_DICT

        if not sys.platform.startswith("linux"): # to have some icons also on Windows
            icon = QtGui.QIcon()
            icon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_DirOpenIcon))
            self.btn_sdf_load.setIcon(icon)
            self.btn_session_load.setIcon(icon)

            icon = QtGui.QIcon()
            icon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_MediaPlay))
            self.btn_start.setIcon(icon)

            icon = QtGui.QIcon()
            icon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_ArrowLeft))
            self.btn_prev.setIcon(icon)
            self.btn_load_prev_cluster.setIcon(icon)

            icon = QtGui.QIcon()
            icon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_ArrowRight))
            self.btn_next.setIcon(icon)
            self.btn_load_next_cluster.setIcon(icon)

            icon = QtGui.QIcon()
            icon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_TrashIcon))
            self.btn_del_sdf.setIcon(icon)

            icon = QtGui.QIcon()
            icon.addPixmap(self.style().standardPixmap(QtGui.QStyle.SP_DialogSaveButton))
            self.btn_sdf_save.setIcon(icon)
            self.btn_session_save.setIcon(icon)

        self.btn_sdf_save.setEnabled(False)
        self.btn_del_sdf.setEnabled(False)
        self.btn_session_save.setEnabled(False)
        self.btn_prev.setEnabled(False)
        self.btn_next.setEnabled(False)
        self.btn_hist.setEnabled(False)
        self.btn_scatter.setEnabled(False)
        self.btn_scatter2.setEnabled(False)
        self.check_lipinski.setEnabled(False)
        self.check_rec_selected.setEnabled(False)
        self.le_recnumber.setReadOnly(False)
        self.label_molimage.setEnabled(True)
        self.btn_load_next_cluster.setEnabled(False)
        self.btn_load_prev_cluster.setEnabled(False)
        self.btn_browse.setEnabled(False)
        self.btn_sdf_from_selection.setEnabled(False)
        self.btn_copy_molids.setEnabled(False)
        self.btn_select_all.setEnabled(False)
        self.btn_select_none.setEnabled(False)
        self.btn_select_invert.setEnabled(False)
        self.label_molimage.customContextMenuRequested.connect(self.cntxt_menu)

        # can be omitted when @pyqtslot decorator is ***NOT*** used !!
        # QtCore.QObject.connect(self.table_query_view, QtCore.SIGNAL("cellClicked(int, int)"),
        #                          self.on_table_query_view_cellClicked)


    @QtCore.pyqtSlot()
    def on_le_sdf_name_editingFinished(self):
        if len(self.le_sdf_name.text()) > 0:
            self.btn_sdf_load.setEnabled(True)


    @QtCore.pyqtSlot()
    def on_btn_sdf_load_clicked(self):

        #TODO: check that the specified file does exist.

        self.sdf_dict = {}
        self.sdf_info_dict = {}
        self.sdf_name_order = []

        self.get_sdfname()

        self.statusbar.showMessage("loading sd file {}.".format(self.sdf_filename), 2000)

        self.sdf_dict[self.sdf_name] = sdft.load_sdf(self.sdf_filename)
        # [<origin>, "factsearch", "n_molwt < 500", <num_of_molecules>, <remark>]
        self.sdf_info_dict[self.sdf_name] = ["", "", "", str(len(self.sdf_dict[self.sdf_name])), "original sdf"]
        self.sdf_name_order.append(self.sdf_name)

        self.create_index()

        # for cluster enumerator:
        # match sdf names of the form cluster_001_009
        CLUSTERFILE = False
        split_name = self.sdf_name.split("_")
        if len(split_name) >= 3:
            if len(split_name[-2]) == 3 and len(split_name[-1]) == 3:
                try:
                    self.cluster_ind = int(split_name[-2])
                    self.cluster_max = int(split_name[-1])
                    CLUSTERFILE = True

                except ValueError:
                    CLUSTERFILE = False
                    print("  ** no cluster file")

        if not CLUSTERFILE:
            self.btn_load_next_cluster.setEnabled(False)
            self.btn_load_prev_cluster.setEnabled(False)
        else:
            self.cluster_name = "_".join(split_name[:-2])
            if self.cluster_ind > 0:
                self.btn_load_prev_cluster.setEnabled(True)
            else:
                self.btn_load_prev_cluster.setEnabled(False)
            if self.cluster_ind < self.cluster_max:
                self.btn_load_next_cluster.setEnabled(True)
            else:
                self.btn_load_next_cluster.setEnabled(False)

        # check, if the name of the sdf is in the config dict BROWSE_OPTIONS
        #   to enable linking to websites based on an identifier
        self.btn_browse.setEnabled(False)
        for name in config.BROWSE_OPTIONS:
            if name in self.sdf_name:
                self.btn_browse.setEnabled(True)
                self.browse_key = config.BROWSE_OPTIONS[name][0]
                self.browse_key_type = config.BROWSE_OPTIONS[name][1]
                self.url_templ = config.BROWSE_OPTIONS[name][2]
                break

        self.init_curr_sdf(self.sdf_name)

        self.statusbar.showMessage("{} molecules loaded into entry {} from file {}.".format(self.curr_sdf_num_of_mols, self.sdf_name, self.sdf_filename), 5000)


    @QtCore.pyqtSlot()
    def on_btn_load_prev_cluster_clicked(self):
        self.cluster_ind -= 1
        new_cluster_name = "_".join([self.cluster_name, "{:03d}".format(self.cluster_ind), "{:03d}".format(self.cluster_max)])
        self.le_sdf_name.setText(new_cluster_name)
        self.on_btn_sdf_load_clicked()


    @QtCore.pyqtSlot()
    def on_btn_load_next_cluster_clicked(self):
        self.cluster_ind += 1
        new_cluster_name = "_".join([self.cluster_name, "{:03d}".format(self.cluster_ind), "{:03d}".format(self.cluster_max)])
        self.le_sdf_name.setText(new_cluster_name)
        self.on_btn_sdf_load_clicked()


    @QtCore.pyqtSlot()
    def on_btn_sdf_save_clicked(self):
        self.statusbar.showMessage("saving current sdf...", 2000)
        sdft.write_sdf(self.curr_sdf, op.join(sdft.REPORT_FOLDER, "sdf", self.curr_sdf_name+".sdf"))
        self.statusbar.showMessage("sdf {} saved with {} records.".format(self.curr_sdf_name, len(self.curr_sdf)), 2000)


    @QtCore.pyqtSlot()
    def on_btn_del_sdf_clicked(self):
        #self.sdf_dict[sdf_name] = mol_list
        #self.sdf_info_dict[sdf_name] = [self.curr_sdf_name, querytype, query, str(len(mol_list)), ""]
        #self.sdf_name_order.append(sdf_name)
        #self.init_curr_sdf(sdf_name)

        self.sdf_dict.pop(self.curr_sdf_name)
        self.sdf_info_dict.pop(self.curr_sdf_name)
        index_of = self.sdf_name_order.index(self.curr_sdf_name)
        self.sdf_name_order.pop(index_of)
        l = len(self.sdf_name_order)
        if index_of > l-1:
            index_of = l-1
        sdf_name = self.sdf_name_order[index_of]

        self.init_curr_sdf(sdf_name)


    @QtCore.pyqtSlot()
    def on_btn_session_save_clicked(self):
        session_dir = op.join(sdft.REPORT_FOLDER, "session")
        if not op.isdir(session_dir):
            print("  > creating session folder {}...".format(session_dir))
            os.mkdir(session_dir)
        f = open(op.join(session_dir, self.sdf_name+".session"), "wb")
        self.statusbar.showMessage("saving session...", 2000)
        #session = [self.sdf_dict, self.sdf_info_dict, self.sdf_name_order, self.curr_sdf_name]
        #pickle.dump(session, f, protocol=2)
        session = [self.sdf_info_dict, self.sdf_name_order, self.curr_sdf_name]
        pickle.dump(session, f, protocol=2)
        f.close()
        for num, sdf_name in enumerate(self.sdf_name_order):
            sdf = self.sdf_dict[sdf_name]
            sdft.write_sdf(sdf, op.join(sdft.REPORT_FOLDER, "session", "session_{}_{:02d}".format(self.sdf_name, num)))

        self.statusbar.showMessage("session {} saved.".format(self.sdf_name), 2000)


    @QtCore.pyqtSlot()
    def on_btn_session_load_clicked(self):
        self.get_sdfname()

        f = open(op.join(sdft.REPORT_FOLDER, "session", self.sdf_name+".session"), "rb")
        self.statusbar.showMessage("loading session {}...".format(self.sdf_name), 2000)
        self.sdf_dict = {}
        [self.sdf_info_dict, self.sdf_name_order, self.curr_sdf_name] = pickle.load(f)
        f.close()
        for num, sdf_name in enumerate(self.sdf_name_order):
            sdf = sdft.load_sdf(op.join(sdft.REPORT_FOLDER, "session", "session_{}_{:02d}".format(self.sdf_name, num)))
            self.sdf_dict[sdf_name] = sdf

        self.create_index()
        self.init_curr_sdf(self.curr_sdf_name)

        self.statusbar.showMessage("session loaded.", 2000)


    @QtCore.pyqtSlot()
    def on_btn_start_clicked(self):
        querytype = str(self.combo_querytype.currentText())
        query = str(self.le_query.text())
        sdf_name = str(self.le_out.text())
        while sdf_name in self.sdf_dict:
            # make sure the name is not already in use
            sdf_name = sdf_name + "x"

        invert = "inverted" in querytype

        if "fact" in querytype:
            mol_list = sdft.factsearch(self.curr_sdf, query, invert=invert)
        else:
            mol_list = sdft.substruct_search(self.curr_sdf, query, invert=invert)
        self.statusbar.showMessage("found {} matching records.".format(len(mol_list)), 2000)
        if mol_list:
            self.sdf_dict[sdf_name] = mol_list
            self.sdf_info_dict[sdf_name] = [self.curr_sdf_name, querytype, query, str(len(mol_list)), ""]
            self.sdf_name_order.append(sdf_name)
            self.init_curr_sdf(sdf_name)


    @QtCore.pyqtSlot()
    def on_btn_prev_clicked(self):
        self.curr_sdf_mol_index -= 1
        self.display_mol()


    @QtCore.pyqtSlot()
    def on_btn_next_clicked(self):
        self.curr_sdf_mol_index += 1
        self.display_mol()


    @QtCore.pyqtSlot()
    def on_btn_browse_clicked(self):
        browse_id = self.browse_key_type(self.curr_sdf[self.curr_sdf_mol_index].GetProp(self.browse_key))
        url = self.url_templ.format(browse_id)
        webbrowser.open(url)
        #subprocess.call(["kde-open", url])


    @QtCore.pyqtSlot()
    def on_btn_select_all_clicked(self):
        self.selected_recs = range(len(self.curr_sdf))
        print("  > all records selected.")
        self.statusbar.showMessage("all records selected.", 2000)
        self.display_mol()


    @QtCore.pyqtSlot()
    def on_btn_select_none_clicked(self):
        self.selected_recs = []
        print("  > no records selected.")
        self.statusbar.showMessage("no records selected.", 2000)
        self.display_mol()


    @QtCore.pyqtSlot()
    def on_btn_select_invert_clicked(self):
        all_indexes = range(len(self.curr_sdf))
        self.selected_recs = list(set(all_indexes)-set(self.selected_recs))
        print("  > selection inverted.")
        self.statusbar.showMessage("selection inverted.", 2000)
        self.display_mol()


    @QtCore.pyqtSlot()
    def on_btn_copy_molids_clicked(self):
        molid_list = []
        for idx in self.selected_recs:
            molid_list.append(self.curr_sdf[idx].GetProp("k_molid"))
        QtGui.QApplication.clipboard().setText("\n".join(molid_list))


    @QtCore.pyqtSlot()
    def on_btn_sdf_from_selection_clicked(self):
        mol_list = []
        sdf_name = "selection"
        while sdf_name in self.sdf_dict:
            # make sure the name is not already in use
            sdf_name = sdf_name + "x"

        for idx in self.selected_recs:
            mol_list.append(self.curr_sdf[idx])

        self.sdf_dict[sdf_name] = mol_list
        self.sdf_info_dict[sdf_name] = [self.curr_sdf_name, "selection", "", str(len(mol_list)), ""]
        self.sdf_name_order.append(sdf_name)
        self.init_curr_sdf(sdf_name)


    @QtCore.pyqtSlot()
    def on_label_molimage_mousePressEvent(self, msevent):
        if self.check_rec_selected.isChecked():
            self.check_rec_selected.setChecked(False)
        else:
            self.check_rec_selected.setChecked(True)


    def on_check_rec_selected_stateChanged(self, state):
        if state == 2:  # isChecked()
            if not self.curr_sdf_mol_index in self.selected_recs:
                self.selected_recs.append(self.curr_sdf_mol_index)
                self.statusbar.showMessage("{} records selected.".format(len(self.selected_recs)), 1000)
                print("    {} was added to the list of selected records.".format(self.curr_sdf_mol_index))
        else:
            if self.curr_sdf_mol_index in self.selected_recs:
                self.selected_recs.remove(self.curr_sdf_mol_index)
                print("    {} was removed from the list of selected records.".format(self.curr_sdf_mol_index))

        if self.selected_recs:
            self.btn_sdf_from_selection.setEnabled(True)
            self.btn_copy_molids.setEnabled(True)
            self.btn_select_invert.setEnabled(True)
        else:
            self.btn_sdf_from_selection.setEnabled(False)
            self.btn_copy_molids.setEnabled(False)
            self.btn_select_invert.setEnabled(True)


    def on_le_recnumber_returnPressed(self):
        self.curr_sdf_mol_index = int(self.le_recnumber.text()) - 1
        self.display_mol()


    #@QtCore.pyqtSlot()
    def wheelEvent(self, event):
        if event.delta() < 0 and self.curr_sdf_mol_index < self.curr_sdf_num_of_mols - 1:
            self.curr_sdf_mol_index += 1
            self.display_mol()
            event.ignore()
            return
        if event.delta() > 0 and self.curr_sdf_mol_index > 0:
            self.curr_sdf_mol_index -= 1
            self.display_mol()
            event.ignore()
            return


    def cntxt_menu(self, pos):
        menu = QtGui.QMenu()
        copy_as_ctab = menu.addAction("Copy as CTab")
        copy_as_smiles = menu.addAction("Copy as Smiles")
        copy_as_image = menu.addAction("Copy as Image")
        action = menu.exec_(self.mapToGlobal(pos))
        if action == copy_as_ctab:
            self.do_copy_as_ctab()
        elif action == copy_as_smiles:
            self.do_copy_as_smiles()
        elif action == copy_as_image:
            self.do_copy_as_image()


    def do_copy_as_ctab(self):
        ctab = Chem.MolToMolBlock(self.curr_sdf[self.curr_sdf_mol_index])
        QtGui.QApplication.clipboard().setText(ctab)


    def do_copy_as_smiles(self):
        smiles = Chem.MolToSmiles(self.curr_sdf[self.curr_sdf_mol_index])
        QtGui.QApplication.clipboard().setText(smiles)


    def do_copy_as_image(self):
        QtGui.QApplication.clipboard().setPixmap(self.qpixmap)


    #@QtCore.pyqtSlot()
    def on_table_query_view_cellClicked(self, row, colmn):
        sdf_name = self.sdf_name_order[row]
        if sdf_name != self.curr_sdf_name:
            self.init_curr_sdf(sdf_name)


    def onpick(self, event):
        x = event.mouseevent.xdata
        y = event.mouseevent.ydata
        ind = event.ind
        ax = event.mouseevent.inaxes
        if ind:
            ind = ind[0]
            if self.curr_sdf_index != 0:
                # use root sdf to display structures
                self.init_curr_sdf(self.sdf_name_order[0])

            molid = self.axes_molindex_dict[ax][ind]
            molindex = self.sdf_mol_index[molid]

            if self.curr_sdf_mol_index != self.sdf_mol_index[molid]:
                self.curr_sdf_mol_index = molindex
                print("    x= {:6.2f}  y= {:6.2f}    index {:4d}: molid {}".format(x, y, molindex, molid))
                self.display_mol()

            else:
                # toggle the selection checkbox only on the second click on the data point:
                # self.curr_sdf_mol_index == self.sdf_mol_index[molid]
                # (this automagically triggers on_check_rec_selected_stateChanged)
                self.check_rec_selected.setChecked(not self.check_rec_selected.isChecked())


    @QtCore.pyqtSlot()
    def on_btn_hist_clicked(self):
        selected_fields = self.get_selected_fields()
        self.statusbar.showMessage("generating histogram...", 2000)
        sdft.show_hist(self.curr_sdf, selected_fields, show=True)
        self.statusbar.showMessage("histogram generated.", 2000)
        # self.label_plot_image.setPixmap(QtGui.QPixmap("histogram.png"))


    @QtCore.pyqtSlot()
    def on_btn_scatter_clicked(self):
        colorby = str(self.combo_colorby.currentText())
        print("  > colorby:", colorby)
        if colorby == "<none>":
            colorby = None
        selected_fields = self.get_selected_fields()
        self.statusbar.showMessage("generating scatter matrix...", 2000)
        self.fig, self.axes_molindex_dict = sdft.show_scattermatrix(self.curr_sdf, fields=selected_fields, colorby=colorby, mode="gui")
        self.statusbar.showMessage("scatter matrix generated.", 2000)
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.show()


    @QtCore.pyqtSlot()
    def on_btn_scatter2_clicked(self):
        compareto = str(self.combo_compareto.currentText())
        print("  > compareto:", compareto)
        selected_fields = self.get_selected_fields()
        self.statusbar.showMessage("generating scatter matrix...", 2000)
        self.statusbar.showMessage("scatter matrix generated.", 2000)
        self.fig, self.axes_molindex_dict = sdft.show_scattermatrix2(self.curr_sdf, self.sdf_dict[compareto], fields=selected_fields, mode="gui")
        self.fig.canvas.mpl_connect('pick_event', self.onpick)
        self.fig.show()


    #@QtCore.pyqtSlot()
    def on_check_lipinski_stateChanged(self, state):
        # lipinski = self.check_lipinski.isChecked()
        self.display_mol()


    def get_selected_fields(self):
        selected_fields = []
        for item in self.table_props.selectedIndexes():
            selected_fields.append(self.curr_sdf_fields[item.row()])
        if len(selected_fields) == 0:
            # if no cells are selected, plot all fields
            return None
        else:
            print("   ", selected_fields)
            return selected_fields


    def create_index(self):
        self.sdf_mol_index = {}
        sdf = self.sdf_dict[self.sdf_name_order[0]]
        fields = sdft.list_fields(sdf)

        prep_needed = False
        for field in fields:
            if field[:2] not in ["n_", "s_", "k_"]:
                prep_needed = True
                break

        if prep_needed:
            print("  > renaming existing fields for viewer:")
            sdft.prepare_for_viewer(sdf)

        if not "k_molid" in fields:
            print("  * no molid field found!")
            print("    calculating fields...",)
            sdft.calc_props(sdf, force2d=prep_needed)
            print("done.")

        for index, mol in enumerate(sdf):
            self.sdf_mol_index[int(mol.GetProp("k_molid"))] = index

        print("  > index created.")

        # show final list of available fields:
        fields = sdft.list_fields(sdf)
        print("  > available fields:", end="")

        for idx, field in enumerate(fields):
            if idx % 5 == 0:
                # new line
                print()
                print("    ", end="")

            print("{}   ".format(field), end="")

        print()



    def get_sdfname(self):
        self.sdf_filename = str(self.le_sdf_name.text())
        if ".sdf" in self.sdf_filename:
            self.sdf_name = self.sdf_filename.split(".sdf")[0]
            self.le_sdf_name.setText(self.sdf_name)
        else:
            self.sdf_name = self.sdf_filename
            self.sdf_filename = self.sdf_filename + ".sdf"
        self.sdf_filename = op.join(sdft.REPORT_FOLDER, "sdf", self.sdf_filename)


    #TODO: make sure the entry fields are filled before enabling the start button
    def edit_query(self):
        if len(self.le_sdf_name.text()) > 0:
            self.btn_sdf_load.setEnabled(True)


    def init_curr_sdf(self, sdf_name):

        self.SDF_CHANGED = True
        self.selected_recs = []
        self.btn_sdf_from_selection.setEnabled(False)
        self.btn_copy_molids.setEnabled(False)
        self.selected_fields = self.get_selected_fields()

        self.curr_sdf = self.sdf_dict[sdf_name][:] # make a copy
        self.curr_sdf_name = sdf_name
        self.curr_sdf_index = self.sdf_name_order.index(sdf_name)
        self.curr_sdf_num_of_mols = int(self.sdf_info_dict[sdf_name][3])
        self.curr_sdf_fields = sdft.list_fields(self.curr_sdf)
        self.curr_sdf_num_of_fields = len(self.curr_sdf_fields)
        self.le_curr_sdf.setText(sdf_name)

        self.table_props.setRowCount(self.curr_sdf_num_of_fields)

        self.curr_sdf_mol_index = 0
        self.display_mol()

        self.fill_query_table()

        self.btn_sdf_save.setEnabled(True)
        self.btn_select_all.setEnabled(True)
        self.btn_select_none.setEnabled(True)

        if len(self.sdf_dict) > 0:
            self.btn_session_save.setEnabled(True)
            self.btn_hist.setEnabled(True)
            self.btn_scatter.setEnabled(True)
            self.check_lipinski.setEnabled(True)
            self.check_rec_selected.setEnabled(True)

        if len(self.sdf_dict) > 1:
            # scatter2 can only be activated when two sdfs are avail. (after the 1. query)
            self.btn_scatter2.setEnabled(True)

        if self.curr_sdf_index == 0:
            self.btn_del_sdf.setEnabled(False) # the original sdf can not be deleted
        else:
            self.btn_del_sdf.setEnabled(True)

        self.combo_colorby.clear()
        text_fields = ["<none>"]
        text_fields.extend([item for item in self.curr_sdf_fields if item[:2] == "s_"])
        self.combo_colorby.addItems(text_fields)

        self.combo_compareto.clear()
        compare_sdfs = [sdf_entry for sdf_entry in self.sdf_dict if sdf_entry != self.curr_sdf_name]
        self.combo_compareto.addItems(compare_sdfs)


    def display_mol(self):
        highlight_lipinski = self.check_lipinski.isChecked()

        if self.curr_sdf_mol_index == 0:
            self.btn_prev.setEnabled(False)
        else:
            self.btn_prev.setEnabled(True)

        if self.curr_sdf_mol_index == self.curr_sdf_num_of_mols - 1:
            self.btn_next.setEnabled(False)
        else:
            self.btn_next.setEnabled(True)

        img_file = IO() # for structure depiction
        img = sdft.autocrop(Draw.MolToImage(self.curr_sdf[self.curr_sdf_mol_index]), "white")
        img.save(img_file, format='PNG')
        # qimg = QtGui.QImage.fromData(img_file.getvalue())
        self.qpixmap = QtGui.QPixmap()
        self.qpixmap.loadFromData(img_file.getvalue(), "PNG")
        self.label_molimage.setPixmap(self.qpixmap)
        self.le_recnumber.setText("{} of {}".format(self.curr_sdf_mol_index+1, self.curr_sdf_num_of_mols))

        if self.SDF_CHANGED:
            # self.SDF_CHANGED = False
            # self.selected_fields = self.get_selected_fields()
            self.table_props.clear()
            self.table_props.setHorizontalHeaderLabels(["prop", "value"])


        for row, prop in enumerate(self.curr_sdf_fields):
            tbl_item = QtGui.QTableWidgetItem(prop[2:])
            # tbl_item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table_props.setItem(row, 0, tbl_item)

            if self.SDF_CHANGED:
                if self.selected_fields and prop in self.selected_fields:
                    self.table_props.setItemSelected(tbl_item, True)

            if prop in self.curr_sdf[self.curr_sdf_mol_index].GetPropNames():
                value = self.curr_sdf[self.curr_sdf_mol_index].GetProp(prop)
                tbl_item = QtGui.QTableWidgetItem(value)
                # QtCore.Qt.ItemIsEditable is required to edit the cells
                # tbl_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable)
                tbl_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable)
                if highlight_lipinski and prop in self.highlight_dict:
                    if float(value) > self.highlight_dict[prop]:
                        tbl_item.setBackgroundColor(QtGui.QColor(255, 0, 0))
                    else:
                        tbl_item.setBackgroundColor(QtGui.QColor(0, 255, 0))
                self.table_props.setItem(row, 1, tbl_item)
            else:
                tbl_item = QtGui.QTableWidgetItem("n.d.")
                # see above for flag QtCore.Qt.ItemIsEditable
                tbl_item.setFlags(QtCore.Qt.ItemIsEnabled | QtCore.Qt.ItemIsEditable)
                self.table_props.setItem(row, 1, tbl_item)

            self.table_props.setRowHeight(row, 18)

        self.SDF_CHANGED = False

        if self.curr_sdf_mol_index in self.selected_recs:
            self.check_rec_selected.setChecked(True)
        else:
            self.check_rec_selected.setChecked(False)


    def fill_query_table(self):
        self.table_query_view.setRowCount(len(self.sdf_name_order))
        self.table_query_view.setHorizontalHeaderLabels(["sdf", "origin", "type", "query", "hits", "remark"])
        for row, sdf_name in enumerate(self.sdf_name_order):
            tbl_item = QtGui.QTableWidgetItem(sdf_name)
            tbl_item.setFlags(QtCore.Qt.ItemIsEnabled)
            self.table_query_view.setItem(row, 0, tbl_item)
            for colmn, item in enumerate(self.sdf_info_dict[sdf_name]):
                tbl_item = QtGui.QTableWidgetItem(item)
                if colmn != 4: # leave remark field editable
                    tbl_item.setFlags(QtCore.Qt.ItemIsEnabled)
                if colmn == 3: # num_of_mols
                    tbl_item.setTextAlignment(QtCore.Qt.AlignRight)
                self.table_query_view.setItem(row, colmn+1, tbl_item)
            self.table_query_view.setRowHeight(row, 18)

        if len(self.sdf_dict) > 1:
            self.btn_scatter2.setEnabled(True)


if __name__ == "__main__":

    app = QtGui.QApplication(sys.argv)
    app.setApplicationName("SDF Viewer")

    ui = App()
    ui.show()
    sys.exit(app.exec_())
