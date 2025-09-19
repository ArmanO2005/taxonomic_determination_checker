import pandas as pd
from BK_Tree import BKTree

def load_kew_data(path):
    """Load Kew data from a CSV file."""
    kew_data = pd.read_csv(path, sep='|')
    return kew_data[['taxonid', 'genus', 'specificepithet', 'infraspecificepithet', 'scientfiicname', 'scientfiicnameauthorship', 'taxonomicstatus', 'acceptednameusageid', 'parentnameusageid', 'originalnameusageid']]


