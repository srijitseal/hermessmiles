import pandas as pd
import os
from hermessmiles.visualizer import Visualizer

def test_visualizer_creates_file(tmp_path):
    df = pd.DataFrame({
        'index_file1': [0],
        'index_file2': [0],
        'original_smiles_file1': ['CCO'],
        'original_smiles_file2': ['CCO'],
        'std_smiles_file1': ['CCO'],
        'std_smiles_file2': ['CCO'],
        'inchikey_prefix': ['LFQSCWFLJHTTH']
    })
    out = tmp_path / "out.html"
    viz = Visualizer(df, str(out))
    viz.render()
    assert out.exists()
    content = out.read_text()
    assert '<h2>Found 1 match(es)</h2>' in content