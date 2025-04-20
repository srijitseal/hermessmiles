import pandas as pd
import os
from hermessmiles.processor import SMILESProcessor, find_overlaps

# prepare sample CSVs
CSV1 = "test1.csv"
CSV2 = "test2.csv"

def setup_module(module):
    data1 = pd.DataFrame({"smiles": ["CCOCCCN", "CCOCCCNC"]})
    data2 = pd.DataFrame({"smiles": ["CCOCCCNC", "CCOCCCNC"]})
    data1.to_csv(CSV1, index=False)
    data2.to_csv(CSV2, index=False)

def teardown_module(module):
    os.remove(CSV1)
    os.remove(CSV2)


def test_process_and_overlap(tmp_path):
    proc = SMILESProcessor()
    df1 = proc.process(CSV1)
    df2 = proc.process(CSV2)
    overlaps = find_overlaps(df1, df2)
    assert len(df1) == 2
    assert len(df2) == 2
    assert overlaps.iloc[0]['original_smiles_file1'] == 'CCOCCCNC'
    assert overlaps.iloc[0]['original_smiles_file2'] == 'CCOCCCNC'