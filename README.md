# NormSeq

NormSeq is a tool for performing normalization of RNA-seq read counts, written in C++ and callable from Python. It supports various normalization methods such as CPM, TPM, FPKM, RPKM, log2CPM, and Z-score normalization.

## Installation

To install NormSeq, clone the repository and build the C++ extension:

```bash
git clone https://github.com/your-username/NormSeq.git
cd NormSeq
python setup.py build_ext --inplace
