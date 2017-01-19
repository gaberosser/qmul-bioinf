from load_data import microarray_data
import pandas as pd


if __name__ == "__main__":
    data, meta = microarray_data.load_annotated_microarray_sb_data(aggr_field='SYMBOL', aggr_method='median')
    print data.loc['Chd7']