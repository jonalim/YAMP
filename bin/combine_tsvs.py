#!/usr/bin/env python

import pandas as pd, sys, csv

old_data = pd.read_table(sys.argv[1], sep='\t', header=0, dtype=str)
new_data = pd.read_table(sys.argv[2], sep='\t', header=0, dtype=str)

old_data = old_data[[c for c in new_data.columns if c in old_data.columns]]
new_data = pd.concat([old_data, new_data], join='outer')
new_data.to_csv(sys.argv[2], sep='\t', quotechar='', quoting=csv.QUOTE_NONE, index=False)
