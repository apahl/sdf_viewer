# sdf_viewer.py - an interactive SDF viewer
### by Axel Pahl
### axelpahl at gmx dot de
### license: BSD, see [license.txt](license.txt) in this folder

## Usage

* For features and usage of the program, see also sdf_viewer.pdf, included in this folder (presented on RDKit UGM 2014).

* When performing fact searches, please remember that number field names are prefixed with "n_"
and text field names are prefixed with "s_". The program uses this to determine the plottable fields
and the fields that can be used for the "colorby" option.   
  * a search in the field that is displayed as "pIC50" in the table would therefore have to look 
  like this in the **query** entry field:  n_pIC50 >= 8
  * and a search in the text field "family" accordingly like this 
  (values are converted to lower-case in the text search):  "a.1" in s_family
  * field names beginning with "k_" (key, e.g. "k_molid") are used neither for plotting nor for coloring,
  but can be used in fact searching.
* When the opened SDF does not follow this field naming scheme, the program tries to guess the types
of the fields and renames them accordingly.
* Searches can be performed in the original sdf as well as in any other sublist derived from it by earlier searches.
Just click on the respective entry in the table that lists all the search results 
to make it the base for the next search. 

* The **open in browser** button is used to call a website for the displayed record.
The function of the **open in browser** button can be defined in the *sdf_viewer_config.py* file (2 examples are given)
  * when defined, a url template is formatted with an sdf property (e.g. "k_molid") 
  and the url is passed to webbrowser.open(url)


* A local version of the config file *sdf_viewer_config.py* can be stored in a *~/.sdf_viewer/* dir

* Saved sessions can be opened again under the name of the original sd file, e.g. when the original file was openend 
as *ugm2014_bzr* (with or without *.sdf*) then the saved session can also be loaded as *ugm2014_bzr*

* The functions used by the sdf_viewer are bundled in the module *sdf_tools.py* and can also be accessed
e.g. from an IPython session. A few examples are given in the included IPython notebook [nb_example.ipynb](http://nbviewer.ipython.org/github/apahl/sdf_viewer/blob/master/nb_example.ipynb)

### Planned Improvements (many based on discussions during the RDKit UGM 2014)

- [ ] combine search results with *or*
- [ ] show molecule grid of selected records
- [ ] include editable property ("s_remark")
- [ ] rename property
- [ ] select calculable properties
- [ ] create some unit tests
