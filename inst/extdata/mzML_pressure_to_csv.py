import sys
import base64
import zlib
import re
import numpy as np
import pandas as pd
from lxml import etree

def convert_to_safe_filename(string):
    #replace all not alphanum chars with underscore
    return re.sub(r'[^a-zA-Z0-9]', '_', string)

def run(file):
    #parse mzML and read namespace
    mzml = etree.parse(file)
    namespace = mzml.getroot().nsmap
    
    #iterate through pressure chromatograms
    for chrom in mzml.findall('//chromatogram/cvParam[@accession="MS:1003019"]/..', namespaces=namespace):
        chrom_name = chrom.get('id')
        decoded_arrays = dict()
        #iterate through binary arrays
        for array in chrom.findall('binaryDataArrayList/binaryDataArray', namespaces=namespace):
            #decode binary data
            binary = base64.b64decode(array.find('binary', namespaces=namespace).text.encode('ascii'))
            
            #parsing array parameters (quick and dirty approach)
            #set all to None
            kind = None
            compression = None
            size = None
            
            for param in array.findall('cvParam', namespaces=namespace):
                name = param.get('name').split()
                if name[-1] == 'array':
                    kind = name[0]
                elif name[-1] == 'compression':
                    compression = name[0]
                elif name[-1] == 'float':
                    size = name[0]
            
            #try parsing binary format
            if size == '64-bit':
                datatype='float64'
            elif size == '32-bit':
                datatype='float32'
            else:
                raise ValueError(f'Unknown binary format {size}')
            
            #try parsing compression and decode the array
            if compression == 'zlib':
                decoded_array = np.frombuffer(zlib.decompress(binary), dtype=datatype)
            elif compression == 'no':
                decoded_array = np.frombuffer(binary, dtype=datatype)
            else:
                raise ValueError(f'Cannot parse compression from {compression}')
            
            if not kind is None:
                decoded_arrays[kind] = decoded_array
            else:
                raise ValueError(f'Kind of array is not set')
                
        csv_path = file.replace('.mzML', f'_{convert_to_safe_filename(chrom_name)}.csv')
        pd.DataFrame(decoded_arrays).to_csv(csv_path, index=False)

if __name__ == '__main__':
    for file in sys.argv[1:]:
        try:
            run(file)
        except Exception as e:
            print (f'Processing of {file} failed; raised {e}')