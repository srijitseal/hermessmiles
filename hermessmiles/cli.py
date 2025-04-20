import argparse
from .processor import SMILESProcessor, find_overlaps
from .visualizer import Visualizer


def main():
    parser = argparse.ArgumentParser(prog='hermessmiles', description='Compare SMILES files')
    parser.add_argument('--file1', required=True, help='Path to first CSV')
    parser.add_argument('--file2', required=True, help='Path to second CSV')
    parser.add_argument('--smiles-col', default=None, help='SMILES column name')
    parser.add_argument('--output', help='Output HTML path')
    args = parser.parse_args()

    proc = SMILESProcessor()
    df1 = proc.process(args.file1, args.smiles_col)
    df2 = proc.process(args.file2, args.smiles_col)
    overlaps = find_overlaps(df1, df2)
    viz = Visualizer(overlaps, args.output)
    viz.render()

if __name__ == '__main__':
    main()
