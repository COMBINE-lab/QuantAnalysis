import pandas as pd

def readRSEMTruth(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t', skiprows=1,
                     names=['Name', 'Gene{}'.format(suffix),
                            'Length{}'.format(suffix),
                            'EffectiveLength{}'.format(suffix),
                            'NumReads{}'.format(suffix),
                            'TPM{}'.format(suffix),
                            'FPKM{}'.format(suffix),
                            'IsoPct{}'.format(suffix)], engine='c')
    df.set_index("Name", inplace=True)
    return df


def readStringTie(fn, suffix=""):
    """
    Not yet tested
    """
    df = pd.read_csv(fn, sep="\t", skiprows=1,
                     names=["tid{}".format(suffix),
                            "chr{}".format(suffix),
                            "strand{}".format(suffix),
                            "start{}".format(suffix),
                            "end{}".format(suffix),
                            "Name",
                            "num_exons{}".format(suffix),
                            "Length{}".format(suffix),
                            "gene_id{}".format(suffix),
                            "gene_name{}".format(suffix),
                            "cov{}".format(suffix),
                            "FPKM{}".format(suffix)])
    df.set_index('Name', inplace=True)
    df.convert_objects(convert_numeric=True)
    return df

def readExpress(fn, suffix=""):
    df = pd.read_csv(fn, sep="\t", skiprows=1,
                     names=["bundle_id{}".format(suffix),
                            "Name",
                            "Length{}".format(suffix),
                            "EffectiveLength{}".format(suffix),
                            "tot_counts{}".format(suffix),
                            "uniq_counts{}".format(suffix),
                            "NumReads{}".format(suffix),
                            "NumReadsEffective{}".format(suffix),
                            "ambig_distr_alpha{}".format(suffix),
                            "ambig_distr_beta{}".format(suffix),
                            "fpkm{}".format(suffix),
                            "fpkm_conf_low{}".format(suffix),
                            "fpkm_conf_high{}".format(suffix),
                            "solvable{}".format(suffix),
                            "TPM{}".format(suffix)])
    df.convert_objects(convert_numeric=True)
    df.set_index('Name', inplace=True)
    return df

def readSailfish(fn, suffix=""):
    df = pd.read_table(fn, engine='c')
    df.columns = [ "{}{}".format(cn, suffix) if cn != "Name" else cn for cn in df.columns.tolist()]
    df.dropna(how='all', inplace=True)
    df.convert_objects(convert_numeric=True)
    df.set_index('Name', inplace=True)
    return df

def readSalmon(fn, suffix=""):
    return readSailfish(fn, suffix)

def readSalmonDeprecated(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t', comment='#',
                     names=['Name',
                            'Length{}'.format(suffix),
                            'TPM{}'.format(suffix),
                            'NumReads{}'.format(suffix)])
    df.dropna(how='all', inplace=True)
    df.convert_objects(convert_numeric=True)
    df.set_index('Name', inplace=True)
    return df

def readSailfishDeprecated(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t', comment='#',
                     names=['Name',
                            'Length{}'.format(suffix),
                            'TPM{}'.format(suffix),
                            'RPKM{}'.format(suffix),
                            'KPKM{}'.format(suffix),
                            'NumKmers{}'.format(suffix),
                            'NumReads{}'.format(suffix)])
    df.dropna(how='all', inplace=True)
    df.convert_objects(convert_numeric=True)
    df.set_index('Name', inplace=True)
    return df

def readKallisto(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t', skiprows=1,
                     names=['Name',
                            'Length{}'.format(suffix),
                            'EffectiveLength{}'.format(suffix),
                            'NumReads{}'.format(suffix),
                            'TPM{}'.format(suffix)], engine='c')
    df.set_index('Name', inplace=True)
    df.convert_objects(convert_numeric=True)
    return df

def readProFile(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t',
                     names=['Locus{}'.format(suffix),
                            'Name',
                            'Coding{}'.format(suffix),
                            'Length{}'.format(suffix),
                            'ExpFrac{}'.format(suffix),
                            'ExpNum{}'.format(suffix),
                            'LibFrac{}'.format(suffix),
                            'LibNum{}'.format(suffix),
                            'SeqFrac{}'.format(suffix),
                            'SeqNum{}'.format(suffix),
                            'CovFrac{}'.format(suffix),
                            'ChiSquare{}'.format(suffix),
                            'CV{}'.format(suffix)])
    df.set_index('Name', inplace=True)
    df.convert_objects(convert_numeric=True)
    return df

def readResFile(fn, suffix=""):
    df = pd.read_csv(fn, sep='\t',
                     names=['Name',
                            'Length{}'.format(suffix),
                            'Abund{}'.format(suffix)])
    df.set_index('Name', inplace=True)
    df.convert_objects(convert_numeric=True)
    return df

