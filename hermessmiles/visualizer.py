import io
import base64
from datetime import datetime
from PIL import Image
from rdkit import Chem
from rdkit.Chem import Draw
import os

class Visualizer:
    def __init__(self, overlaps_df, output_path: str = None):
        self.overlaps = overlaps_df
        if output_path:
            self.output_html = output_path
        else:
            ts = datetime.now().strftime('%Y%m%d_%H%M%S')
            self.output_html = f"overlap_structures_{ts}.html"

    def render(self):
        # ensure output directory exists
        out_dir = os.path.dirname(self.output_html)
        if out_dir:
            os.makedirs(out_dir, exist_ok=True)
            
        # build the HTML (guaranteed to be a str)
        content = self._build_html()
        # write it out
        with open(self.output_html, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Saved visualization to {self.output_html}")

    def _build_html(self) -> str:
        df = self.overlaps
        # no matches → simple page
        if df is None or df.empty:
            return "<html><body><h2>No matching structures found.</h2></body></html>"

        # header + style
        header = (
            "<!DOCTYPE html><html><head><meta charset='utf-8'>"
            "<style>"
            "body{font-family:Arial,sans-serif;} "
            "table{border-collapse:collapse;width:100%;} "
            "th,td{border:1px solid #ddd;padding:8px;text-align:left;} "
            "th{background:#f2f2f2;} tr:nth-child(even){background:#f9f9f9;} "
            ".structure-img{width:200px;}"
            "</style>"
            f"<title>Matches: {len(df)}</title></head><body>"
            f"<h2>Found {len(df)} match(es)</h2>"
            "<table>"
            "<tr><th>#</th><th>Idx1</th><th>Structure1</th><th>SMILES1</th>"
            "<th>Idx2</th><th>Structure2</th><th>SMILES2</th><th>InChIKey Prefix</th></tr>"
        )

        rows_html = ""
        for idx, row in df.iterrows():
            m1 = Chem.MolFromSmiles(row['std_smiles_file1'])
            m2 = Chem.MolFromSmiles(row['std_smiles_file2'])
            img1 = Draw.MolToImage(m1, size=(300,200))
            img2 = Draw.MolToImage(m2, size=(300,200))
            b1 = base64.b64encode(self._to_bytes(img1)).decode('utf-8')
            b2 = base64.b64encode(self._to_bytes(img2)).decode('utf-8')
            rows_html += (
                "<tr>"
                f"<td>{idx+1}</td>"
                f"<td>{row['index_file1']}</td>"
                f"<td><img src='data:image/png;base64,{b1}' class='structure-img'/></td>"
                f"<td>{row['original_smiles_file1']}</td>"
                f"<td>{row['index_file2']}</td>"
                f"<td><img src='data:image/png;base64,{b2}' class='structure-img'/></td>"
                f"<td>{row['original_smiles_file2']}</td>"
                f"<td>{row['inchikey_prefix']}</td>"
                "</tr>"
            )

        footer = "</table></body></html>"
        return header + rows_html + footer

    @staticmethod
    def _to_bytes(img: Image.Image) -> bytes:
        buf = io.BytesIO()
        img.save(buf, format='PNG')
        return buf.getvalue()
