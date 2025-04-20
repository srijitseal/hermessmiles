import pandas as pd
from rdkit import Chem
from .athena_smiles_standardisation import SMILESStandardizer

def generate_inchikey(smiles: str) -> str:
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return None
        return Chem.MolToInchiKey(mol)
    except:
        return None

class SMILESProcessor:
    def __init__(self):
        self.standardizer = SMILESStandardizer()

    def process(self, file_path: str, smiles_col: str = None) -> pd.DataFrame:
        df = pd.read_csv(file_path, sep=None, engine='python')
        if smiles_col is None:
            for col in ['SMILES','smiles','Smiles','Canonical_SMILES']:
                if col in df.columns:
                    smiles_col = col
                    break
            else:
                smiles_col = df.columns[0]
        df = df.rename(columns={smiles_col: 'original_smiles'})
        df['std_smiles'] = df['original_smiles'].apply(self.standardizer.standardize)
        df = df.dropna(subset=['std_smiles'])
        df['inchikey'] = df['std_smiles'].apply(generate_inchikey)
        df['inchikey_prefix'] = df['inchikey'].str[:14]
        df = df.dropna(subset=['inchikey_prefix'])
        return df.reset_index(drop=True)


def find_overlaps(df1: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    overlaps = []
    for i, r1 in df1.iterrows():
        matches = df2[df2['inchikey_prefix'] == r1['inchikey_prefix']]
        for j, r2 in matches.iterrows():
            overlaps.append({
                'index_file1': i,
                'index_file2': j,
                'original_smiles_file1': r1['original_smiles'],
                'original_smiles_file2': r2['original_smiles'],
                'std_smiles_file1': r1['std_smiles'],
                'std_smiles_file2': r2['std_smiles'],
                'inchikey_prefix': r1['inchikey_prefix']
            })
    return pd.DataFrame(overlaps)