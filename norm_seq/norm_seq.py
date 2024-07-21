import pandas as pd
import normalization

class NormSeq:
    def __init__(self, raw_counts: pd.DataFrame, lengths: pd.Series = None):
        self.raw_counts = raw_counts
        self.lengths = lengths

    def cpm(self):
        data = self.raw_counts.values.astype(int)
        normalized_data = normalization.cpm(data)
        return pd.DataFrame(normalized_data, index=self.raw_counts.index, columns=self.raw_counts.columns)

    def tpm(self):
        data = self.raw_counts.values.astype(int)
        lengths = self.lengths.values
        normalized_data = normalization.tpm(data, lengths)
        return pd.DataFrame(normalized_data, index=self.raw_counts.index, columns=self.raw_counts.columns)

    def fpkm(self):
        data = self.raw_counts.values.astype(int)
        lengths = self.lengths.values
        normalized_data = normalization.fpkm(data, lengths)
        return pd.DataFrame(normalized_data, index=self.raw_counts.index, columns=self.raw_counts.columns)

    def rpkm(self):
        data = self.raw_counts.values.astype(int)
        lengths = self.lengths.values
        normalized_data = normalization.rpkm(data, lengths)
        return pd.DataFrame(normalized_data, index=self.raw_counts.index, columns=self.raw_counts.columns)

    def log2cpm(self):
        data = self.raw_counts.values.astype(int)
        normalized_data = normalization.log2cpm(data)
        return pd.DataFrame(normalized_data, index=self.raw_counts.index, columns=self.raw_counts.columns)

    def zscore(self):
        data = self.raw_counts.values.astype(int)
        normalized_data = normalization.zscore(data)
        return pd.DataFrame(normalized_data, index=self.raw_counts.index, columns=self.raw_counts.columns)

def cpm(dataframe):
    norm_seq = NormSeq(dataframe)
    return norm_seq.cpm()

def tpm(dataframe, lengths):
    norm_seq = NormSeq(dataframe, lengths)
    return norm_seq.tpm()

def fpkm(dataframe, lengths):
    norm_seq = NormSeq(dataframe, lengths)
    return norm_seq.fpkm()

def rpkm(dataframe, lengths):
    norm_seq = NormSeq(dataframe, lengths)
    return norm_seq.rpkm()

def log2cpm(dataframe):
    norm_seq = NormSeq(dataframe)
    return norm_seq.log2cpm()

def zscore(dataframe):
    norm_seq = NormSeq(dataframe)
    return norm_seq.zscore()
