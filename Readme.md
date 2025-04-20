```markdown
# Hermessmiles

A tool to compare two CSV files of SMILES strings, find overlapping compounds by InChIKey prefix, and generate an HTML visualization of matching structures.

## Installation
```bash
pip install .
```

## Usage
```bash
hermessmiles --file1 path/to/df1.csv --file2 path/to/df2.csv [--smiles-col SMILES]
```

By default, the output HTML is saved as `overlap_structures_<timestamp>.html` in your working directory.
```