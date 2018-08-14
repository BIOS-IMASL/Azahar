from collections import namedtuple
Constants = namedtuple('Constants', [
                       'radians_to_degrees',
                       'degrees_to_radians',
                       'pi',
                       'peptide_bond_lenght'])

constants = Constants(57.29577951308232,
                      0.017453292519943295,
                      3.141592653589793,
                      1.364)

rotamers = {1:['R', 'K', 'Q', 'E', 'M', 'D', 'F', 'H', 'L','I', 'W','Y', 'N',
               'S', 'T', 'V', 'C'],
            2:['R', 'K', 'Q', 'E', 'M', 'D', 'F', 'H', 'L','I', 'W','Y', 'N'],
            3:['R', 'K', 'Q', 'E', 'M'],
            4:['R', 'K'],
            5:['R']}
            
chi_ijh = {2:(0, 1, ['H', 'HA', 'HB', 'HB2', 'HB3']),
           3:(1, 2, ['H', 'HA', 'HB', 'HB2', 'HB3', 'HG2', 'HG3']),
           4:(2, 3, ['H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3']),
           5:(3, 4, ['H', 'HA', 'HB2', 'HB3', 'HG2', 'HG3', 'HD2', 'HD3', 'HE'])}

# see scafold/gen_vdw_parameters.py for details on the generation of
# par_s_ij  and par_eps_ij
par_s_ij = {'CC': 3.816, 'HC': 3.395, 'OH': 3.1482, 'NH': 3.311, 'OO': 3.3224,
            'CN': 3.732, 'HN': 3.311, 'ON': 3.4852, 'NN': 3.648, 'SN': 3.824,
            'NO': 3.4852, 'HS': 3.487, 'NS': 3.824, 'OS': 3.6612, 'CH': 3.395,
            'SS': 4.0, 'SH': 3.487, 'CS': 3.908, 'SO': 3.6612, 'NC': 3.732,
            'HH': 2.974, 'OC': 3.5692, 'HO': 3.1482, 'SC': 3.908, 'CO': 3.5692}

par_eps_ij = {'SC': 0.14663, 'SH': 0.06265, 'SN': 0.20616, 'CC': 0.086,
              'CS': 0.14663, 'OC': 0.13439, 'HO': 0.05742, 'HH': 0.0157,
              'NS': 0.20616, 'HS': 0.06265, 'SO': 0.22913, 'SS': 0.25,
              'OO': 0.21, 'HC': 0.03675, 'OS': 0.22913, 'OH': 0.05742,
              'CH': 0.03675, 'NN': 0.17, 'HN': 0.05166, 'CN': 0.12091,
              'NO': 0.18894, 'NC': 0.12091, 'NH': 0.05166, 'ON': 0.18894,
'CO': 0.13439}