#!/usr/bin/env python2
# -*- coding: utf-8 -*-
# sdf_viewer_config.py

from __future__ import absolute_import, division, print_function

# used for Lipinski highlighting, change as needed
HIGHLIGHT_DICT = {"n_molwt": 500.0, "n_hbd": 5, "n_hba": 10, "n_logp": 5.0, "n_rotb": 10, "n_tpsa": 130}

# {partial_sdf_name: [browse_key_ browse_key_type, url_templ]}
BROWSE_OPTIONS = {"aviru": ["k_molid", int, "file:///home/apl/aviru/db_reports/reports/ind_stock_results.htm#cpd_{:05d}"],
                  "cs":    ["k_csid",  str, "http://www.chemspider.com/Chemical-Structure.{}.html"]
                  }
