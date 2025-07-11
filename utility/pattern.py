import itertools
import os
from dataclasses import dataclass, field
from typing import Protocol,Tuple

@dataclass
class ExportDataColumns:
    """A class for defining the basic struture of the dataframes
    for the export of the simulation results."""
    _crackfront: list = field(init=False)
    nodemap: list = field(default_factory=lambda: ['ID',
                                                   'x', 'y', 'z',
                                                   'ux', 'uy', 'uz',
                                                   'epsx', 'epsy', 'epsxy', 'epseqv',
                                                   'sigx', 'sigy', 'sigxy', 'sigeqv',])
    nodemap_old: list = field(default_factory=lambda: ['ID',
                                                   'x', 'y', 'z',
                                                   'ux', 'uy', 'uz',
                                                   'epsx', 'epsy', 'epsxy', 'epseqv',])

    def crackfront(self, n_contours: int = 6) -> list:
        """A method for defining the column names for the crack front dataframes."""
        self._crackfront = ['ID',
                            'x', 'y', 'z',
                            'sigx', 'sigy', 'sigz', 'sigeqv',
                            'epelx', 'epely', 'epelz', 'epseqv',
                            *[f'SIFS_K{K + 1}_Kont{C + 1}' for K in range(3) for C in range(n_contours)],
                            'DLTN', 'DLTA', 'DLTK', 'R',
                            'K_I_SIFS', 'K_II_SIFS', 'K_III_SIFS',
                            'c_cum', 'rel_order', 'node_type',
                            ]
        return list(self._crackfront)