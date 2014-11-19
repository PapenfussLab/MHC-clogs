#!/usr/bin/env python

"""
picard.py

Supporting used functions of Picard suite.

Currently implemented:
    - SortSam
    - AddOrReplaceReadGroups
    - MarkDuplicate
    - MergeSamFiles
"""

import os
from massivetools.pipeline import Pipeline
from massivetools.toolkit import Toolkit
from massivetools.utilities import which
from massivetools.exceptions import MissingExecutableError

PICARD_PATH = "/usr/local/bioinfsoftware/picard-tools/current/jars/"


class Picard(Toolkit):
    def __init__(self, picard_jar_path=PICARD_PATH, name='Picard', 
                    pipeline=Pipeline.get_default_instance(), max_memory='4g', 
                    tmp_dir='.', verbosity=None, quiet=None, 
                    validation_stringency='SILENT', compression_level=None, 
                    max_records_in_ram=None, create_index=None, create_md5_file=None):
        Toolkit.__init__(self)
        
        #'''
        # Common Picard options
        #'''
        def default_value( allowed_list, value ):
            if value is None: return None
            
            value_str = str(value)
            if allowed_list is None:
                return value_str
            else:
                # case insensitive check
                if value_str.upper() in allowed_list: return value_str
            return None
        
        self.common_values = {}
        self.common_values['TMP_DIR']   = default_value( None, tmp_dir )
        self.common_values['VERBOSITY'] = default_value( ['ERROR','INFO','WARNING','DEBUG'], verbosity )
        self.common_values['QUIET']     = default_value( None, quiet )
        self.common_values['VALIDATION_STRINGENCY'] = default_value( ['STRICT','LENIENT', 'SILENT'], validation_stringency )
        self.common_values['COMPRESSION_LEVEL']     = default_value( None, compression_level )
        self.common_values['MAX_RECORDS_IN_RAM']    = default_value( None, max_records_in_ram )
        self.common_values['CREATE_INDEX']          = default_value( None, create_index )
        self.common_values['CREATE_MD5_FILE']       = default_value( None, create_md5_file )
        
        #'''
        # get Java runtime executable
        #'''
        java_exe = which('java')
        if not java_exe:
            raise MissingExecutableError("Unable to locate Java runtime environment")
        self.java_cmd = java_exe + ' -jar -Xmx%s' % max_memory
        
        # ensure path to Picard JAR files
        if not os.path.exists(picard_jar_path):
            raise MissingExecutableError("Picard JAR directory not found")
        self.picard = picard_jar_path
        
        # SortSam
        self.register_tool('SortSam', ' '.join( [self.java_cmd, 
            os.path.join(self.picard, 'SortSam.jar'),
            'INPUT=%(input)s', 'OUTPUT=%(output)s', 'SORT_ORDER=%(sort_order)s',
            'TMP_DIR=%(tmp_dir)s', self.common_options()] ),
            defaults={"sort_order": "coordinate", "tmp_dir": "/tmp"})
        
        # AddOrReplaceReadGroups
        self.register_tool('AddOrReplaceReadGroups', ' '.join( [self.java_cmd, 
            os.path.join(self.picard, 'AddOrReplaceReadGroups.jar'),
            'I=%(input)s', 'O=%(output)s', 'ID=%(id)s', 'LB=%(library)s', 
            'PL=%(platform)s', 'PU=%(unit)s', 'SM=%(sample)s', 
            'SO=%(sort_order)s', self.common_options()] ),
            defaults={'sort_order': 'coordinate'})
        
        # MarkDuplicates
        self.register_tool('MarkDuplicates', ' '.join( [self.java_cmd, os.path.join(self.picard, 'MarkDuplicates.jar'),
            'INPUT=%(input)s', 'OUTPUT=%(output)s', 
            'METRICS_FILE=%(metrics_file)s', 'REMOVE_DUPLICATES=%(remove_duplicates)s',
            'ASSUME_SORTED=%(assume_sorted)s', 'MAX_FILE_HANDLES=%(max_file_handles)d', 
            self.common_options()] ),
            defaults={'remove_duplicates': 'true', 'assume_sorted': 'true', 'max_file_handles': 8000})
        
        # Merge
        self.register_tool('MergeSamFiles', ' '.join( [self.java_cmd, os.path.join(self.picard, 'MergeSamFiles.jar'),
            '%(in_args)s', 'OUTPUT=%(output)s', 'USE_THREADING=%(use_threading)s', 'ASSUME_SORTED=%(assume_sorted)s'] ),
            defaults={'use_threading': 'true', 'assume_sorted': 'true'})
        
    def common_options(self):
        opt = ['%s=%s'%(k,v) for k,v in self.common_values.items() if v is not None]
        return ' '.join(opt)


def PicardToolkit(*args, **kw):
    import sys
    print sys.stderr, "Use Picard"
    return Picard(*args, **kw)


#from massivetools.utilities import add_suffix
#picard = Picard('/Users/hsu/sw/picard/current')
#
#data = '449B_anomalous.bam'
#data_sorted = add_suffix(data, '_s')
#data_rg     = add_suffix(data_sorted, '_rg')
#data_md     = add_suffix(data_rg, '_md')
#
#jid1 = picard.SortSam(input=data, output=data_sorted)
#jid2 = picard.AddOrReplaceReadGroups(input=data_sorted, output=data_rg, id='449B', library='600bp', platform='illumina', platform_unit='GA-IIx', sample='liposarcoma',
#            dependencies=[jid1])
#jid3 = picard.MarkDuplicates(input=data_rg, output=data_md, metrics_file='log.449B_anomalous_s_rg_md',
#            dependencies=[jid2])
#
#Pipeline.run()
#
