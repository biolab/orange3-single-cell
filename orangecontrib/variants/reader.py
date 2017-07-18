import Orange
from Orange.data import StringVariable, DiscreteVariable
import vcf
import numpy as np


class VariantData:
    def __init__(self, filename):
        reader = vcf.Reader(open(filename, "r"))
        records = [r for r in reader]

        self.samples = np.array(reader.samples)

        self.gq = np.array([[s.data.GQ for s in r.samples] for r in records],
                           dtype="f")
        self.gq = np.nan_to_num(self.gq, 0)

        gt = np.array([[s.data.GT for s in r.samples] for r in records])
        self.gt = gt != "0/0"
        self.records = records

        self.variables = [
            DiscreteVariable("%s-%s" % (r.CHROM, r.POS), values=["0", "1"])
            for r in self.records
        ]

        for v, r in zip(self.variables, records):
            v.attributes["CHROM"] = str(r.CHROM)
            v.attributes["POS"] = str(r.POS)
            v.attributes["REF"] = str(r.REF)
            v.attributes["ALT"] = "".join(str(s) for s in r.ALT)

    def info(self):
        print("Samples: %d" % len(self.samples))
        unique_samples = set(s[:-2] for s in self.samples)
        print("Unique samples: %d" % len(unique_samples))
        print("Variants: %d" % len(self.records))

    def get_data(self, quality=None, frequency=None):
        """Orange data table with genotypes above quality and frequency
        threshold."""

        X = self.gt.astype(dtype="float", copy=True)
        if quality is not None:
            X[self.gq < quality] = np.nan
        selected = ~np.isnan(X).all(axis=1)
        if frequency is not None:
            selected &= np.nansum(X, axis=1) >= frequency
        X = X[selected].T

        variables = tuple(np.array(self.variables)[selected])
        metas = [StringVariable(s) for s in ["sample", "id"]]
        M = np.empty((len(X), 2), dtype="object")
        M[:, 0] = self.samples
        M[:, 1] = np.arange(1, len(X)+1)
        domain = Orange.data.Domain(variables, [], metas=metas)
        data = Orange.data.Table(domain, X, Y=None, metas=M)

        return data

if __name__ == "__main__":
    variants = VariantData("cells.vcf")
    data = variants.get_data(40, 10)
