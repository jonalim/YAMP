#!/usr/bin/env python

import pandas as pd, sys, csv

old_data = pd.read_table(sys.argv[1], sep='\t', header=0, index_col=0, dtype=str, quotechar=None, quoting=3)
new_data = pd.read_table(sys.argv[2], sep='\t', header=0, index_col=0, dtype=str, quotechar=None, quoting=3)

old_data = old_data[[c for c in new_data.columns if c in old_data.columns]]
# Combine. New fields not in old table get assigned NaN for samples in old table
new_data = pd.concat([old_data, new_data], join='outer')
# new data overwrites old data for same sample
new_data = new_data[~new_data.index.duplicated(keep='last')] 
new_data.to_csv(sys.argv[2], sep='\t', quotechar='', quoting=csv.QUOTE_NONE, index=True)
