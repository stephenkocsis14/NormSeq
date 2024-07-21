import argparse
import pandas as pd
from norm_seq import cpm, tpm, fpkm, rpkm, log2cpm, zscore

def main():
    parser = argparse.ArgumentParser(description="NormSeq - RNA-seq Data Normalization Tool")
    parser.add_argument('--normMethod', type=str, required=True, help="Normalization method (e.g., cpm, tpm, fpkm, rpkm, log2cpm, zscore)")
    parser.add_argument('--input', type=str, required=True, help="Input CSV file with raw counts")
    parser.add_argument('--output', type=str, required=True, help="Output CSV file for normalized data")
    parser.add_argument('--lengths', type=str, help="CSV file with lengths for TPM, FPKM, and RPKM methods")

    args = parser.parse_args()

    data = pd.read_csv(args.input, index_col=0)
    
    if args.normMethod == 'cpm':
        normalized_data = cpm(data)
    elif args.normMethod == 'tpm':
        lengths = pd.read_csv(args.lengths, index_col=0).squeeze()
        normalized_data = tpm(data, lengths)
    elif args.normMethod == 'fpkm':
        lengths = pd.read_csv(args.lengths, index_col=0).squeeze()
        normalized_data = fpkm(data, lengths)
    elif args.normMethod == 'rpkm':
        lengths = pd.read_csv(args.lengths, index_col=0).squeeze()
        normalized_data = rpkm(data, lengths)
    elif args.normMethod == 'log2cpm':
        normalized_data = log2cpm(data)
    elif args.normMethod == 'zscore':
        normalized_data = zscore(data)
    else:
        raise ValueError("Unsupported normalization method")

    normalized_data.to_csv(args.output)

if __name__ == "__main__":
    main()
