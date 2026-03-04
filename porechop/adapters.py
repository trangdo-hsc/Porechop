"""
Copyright 2017 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Porechop

This module contains the class and sequences for known adapters used in Oxford Nanopore library
preparation kits.

This file is part of Porechop. Porechop is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Porechop is distributed in
the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Porechop. If
not, see <http://www.gnu.org/licenses/>.
"""


class Adapter(object):

    def __init__(self, name, start_sequence=None, end_sequence=None, both_ends_sequence=None):
        self.name = name
        self.start_sequence = start_sequence if start_sequence else []
        self.end_sequence = end_sequence if end_sequence else []
        if both_ends_sequence:
            self.start_sequence = both_ends_sequence
            self.end_sequence = both_ends_sequence
        self.best_start_score, self.best_end_score = 0.0, 0.0

    def best_start_or_end_score(self):
        return max(self.best_start_score, self.best_end_score)

    def is_barcode(self):
        return self.name.startswith('Barcode ')

    def barcode_direction(self):
        if '_rev' in self.start_sequence[0]:
            return 'reverse'
        else:
            return 'forward'

    def get_barcode_name(self):
        """
        Gets the barcode name for the output files. We want a concise name, so it looks at all
        options and chooses the shortest.
        """
        possible_names = [self.name]
        if self.start_sequence:
            possible_names.append(self.start_sequence[0])
        if self.end_sequence:
            possible_names.append(self.end_sequence[0])
        barcode_name = sorted(possible_names, key=lambda x: len(x))[0]
        return barcode_name.replace(' ', '_')


# INSTRUCTIONS FOR ADDING CUSTOM ADAPTERS
# ---------------------------------------
# If you need Porechop to remove adapters that aren't included, you can add your own my modifying
# the ADAPTERS list below.
#
# Here is the format for a normal adapter:
#     Adapter('Adapter_set_name',
#             start_sequence=('Start_adapter_name', 'AAAACCCCGGGGTTTTAAAACCCCGGGGTTTT'),
#             end_sequence=('End_adapter_name', 'AACCGGTTAACCGGTTAACCGGTTAACCGGTT'))
#
# You can exclude start_sequence and end_sequence as appropriate.
#
# If you have custom Barcodes, make sure that the adapter set name starts with 'Barcode '. Also,
# remove the existing barcode sequences from this file to avoid conflicts:
#     Adapter('Barcode 1',
#             start_sequence=('Barcode_1_start', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT'),
#             end_sequence=('Barcode_1_end', 'AAAAAAAACCCCCCCCGGGGGGGGTTTTTTTT')),
#     Adapter('Barcode 2',
#             start_sequence=('Barcode_2_start', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'),
#             end_sequence=('Barcode_2_end', 'TTTTTTTTGGGGGGGGCCCCCCCCAAAAAAAA'))


ADAPTERS = [Adapter('SQK-NSK007',
                    start_sequence=('SQK-NSK007_Y_Top', 'AATGTACTTCGTTCAGTTACGTATTGCT'),
                    end_sequence=('SQK-NSK007_Y_Bottom', 'GCAATACGTAACTGAACGAAGT')),


            Adapter('Rapid',
                    start_sequence=('Rapid_adapter',
                                    'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA')),

            Adapter('RBK004_upstream',
                    start_sequence=('RBK004_upstream', 'AATGTACTTCGTTCAGTTACGGCTTGGGTGTTTAACC')),


            Adapter('SQK-MAP006',
                    start_sequence=('SQK-MAP006_Y_Top_SK63',    'GGTTGTTTCTGTTGGTGCTGATATTGCT'),
                    end_sequence=  ('SQK-MAP006_Y_Bottom_SK64', 'GCAATATCAGCACCAACAGAAA')),

            Adapter('SQK-MAP006 short',
                    start_sequence=('SQK-MAP006_Short_Y_Top_LI32',    'CGGCGTCTGCTTGGGTGTTTAACCT'),
                    end_sequence=  ('SQK-MAP006_Short_Y_Bottom_LI33', 'GGTTAAACACCCAAGCAGACGCCG')),


            # The PCR adapters are used both in PCR DNA kits and some cDNA kits.
            Adapter('PCR adapters 1',
                    start_sequence=('PCR_1_start', 'ACTTGCCTGTCGCTCTATCTTC'),
                    end_sequence=  ('PCR_1_end',   'GAAGATAGAGCGACAGGCAAGT')),

            Adapter('PCR adapters 2',
                    start_sequence=('PCR_2_start', 'TTTCTGTTGGTGCTGATATTGC'),
                    end_sequence=  ('PCR_2_end',   'GCAATATCAGCACCAACAGAAA')),

            Adapter('PCR adapters 3',
                    start_sequence=('PCR_3_start', 'TACTTGCCTGTCGCTCTATCTTC'),
                    end_sequence=  ('PCR_3_end',   'GAAGATAGAGCGACAGGCAAGTA')),


            # 1D^2 kit adapters are interesting. ONT provided the following sequences on their site:
            #   start: GGCGTCTGCTTGGGTGTTTAACCTTTTTGTCAGAGAGGTTCCAAGTCAGAGAGGTTCCT
            #   end:   GGAACCTCTCTCTGACTTGGAACCTCTCTGACAAAAAGGTTAAACACCCAAGCAGACGCCAGCAAT
            # But when looking at actual reads, I found two parts. The first corresponds to one end
            # of the provided sequences (through slightly different):
            Adapter('1D^2 part 1',
                    start_sequence=('1D2_part_1_start', 'GAGAGGTTCCAAGTCAGAGAGGTTCCT'),
                    end_sequence=  ('1D2_part_1_end',   'AGGAACCTCTCTGACTTGGAACCTCTC')),
            # and the second part corresponds to the other end, combined with a bit of standard 1D
            # adapter:
            Adapter('1D^2 part 2',
                    start_sequence=('1D2_part_2_start', 'CTTCGTTCAGTTACGTATTGCTGGCGTCTGCTT'),
                    end_sequence=  ('1D2_part_2_end',   'CACCCAAGCAGACGCCAGCAATACGTAACT')),
            # The middle part of the provided sequences is less common, so I've left it out of the
            # adapter sequences here.


            Adapter('cDNA SSP',
                    start_sequence=('cDNA_SSP',     'TTTCTGTTGGTGCTGATATTGCTGCCATTACGGCCGGG'),
                    end_sequence=  ('cDNA_SSP_rev', 'CCCGGCCGTAATGGCAGCAATATCAGCACCAACAGAAA')),


            # Some barcoding kits (like the native barcodes) use the rev comp barcode at the start
            # of the read and the forward barcode at the end of the read.
            Adapter('Barcode 1 (reverse)',
                    start_sequence=('BC01_rev', 'CACAAAGACACCGACAACTTTCTT'),
                    end_sequence=('BC01', 'AAGAAAGTTGTCGGTGTCTTTGTG')),
            Adapter('Barcode 2 (reverse)',
                    start_sequence=('BC02_rev', 'ACAGACGACTACAAACGGAATCGA'),
                    end_sequence=('BC02', 'TCGATTCCGTTTGTAGTCGTCTGT')),
            Adapter('Barcode 3 (reverse)',
                    start_sequence=('BC03_rev', 'CCTGGTAACTGGGACACAAGACTC'),
                    end_sequence=('BC03', 'GAGTCTTGTGTCCCAGTTACCAGG')),
            Adapter('Barcode 4 (reverse)',
                    start_sequence=('BC04_rev', 'TAGGGAAACACGATAGAATCCGAA'),
                    end_sequence=('BC04', 'TTCGGATTCTATCGTGTTTCCCTA')),
            Adapter('Barcode 5 (reverse)',
                    start_sequence=('BC05_rev', 'AAGGTTACACAAACCCTGGACAAG'),
                    end_sequence=('BC05', 'CTTGTCCAGGGTTTGTGTAACCTT')),
            Adapter('Barcode 6 (reverse)',
                    start_sequence=('BC06_rev', 'GACTACTTTCTGCCTTTGCGAGAA'),
                    end_sequence=('BC06', 'TTCTCGCAAAGGCAGAAAGTAGTC')),
            Adapter('Barcode 7 (reverse)',
                    start_sequence=('BC07_rev', 'AAGGATTCATTCCCACGGTAACAC'),
                    end_sequence=('BC07', 'GTGTTACCGTGGGAATGAATCCTT')),
            Adapter('Barcode 8 (reverse)',
                    start_sequence=('BC08_rev', 'ACGTAACTTGGTTTGTTCCCTGAA'),
                    end_sequence=('BC08', 'TTCAGGGAACAAACCAAGTTACGT')),
            Adapter('Barcode 9 (reverse)',
                    start_sequence=('BC09_rev', 'AACCAAGACTCGCTGTGCCTAGTT'),
                    end_sequence=('BC09', 'AACTAGGCACAGCGAGTCTTGGTT')),
            Adapter('Barcode 10 (reverse)',
                    start_sequence=('BC10_rev', 'GAGAGGACAAAGGTTTCAACGCTT'),
                    end_sequence=('BC10', 'AAGCGTTGAAACCTTTGTCCTCTC')),
            Adapter('Barcode 11 (reverse)',
                    start_sequence=('BC11_rev', 'TCCATTCCCTCCGATAGATGAAAC'),
                    end_sequence=('BC11', 'GTTTCATCTATCGGAGGGAATGGA')),
            Adapter('Barcode 12 (reverse)',
                    start_sequence=('BC12_rev', 'TCCGATTCTGCTTCTTTCTACCTG'),
                    end_sequence=('BC12', 'CAGGTAGAAAGAAGCAGAATCGGA')),

            # Other barcoding kits (like the PCR and rapid barcodes) use the forward barcode at the
            # start of the read and the rev comp barcode at the end of the read.
            Adapter('Barcode 1 (forward)',
                    start_sequence=('NB01', 'CACAAAGACACCGACAACTTTCTT'),
                    end_sequence=('NB01_rev', 'AAGAAAGTTGTCGGTGTCTTTGTG')),
            Adapter('Barcode 2 (forward)',
                    start_sequence=('NB02', 'ACAGACGACTACAAACGGAATCGA'),
                    end_sequence=('NB02_rev', 'TCGATTCCGTTTGTAGTCGTCTGT')),
            Adapter('Barcode 3 (forward)',
                    start_sequence=('NB03', 'CCTGGTAACTGGGACACAAGACTC'),
                    end_sequence=('NB03_rev', 'GAGTCTTGTGTCCCAGTTACCAGG')),
            Adapter('Barcode 4 (forward)',
                    start_sequence=('NB04', 'TAGGGAAACACGATAGAATCCGAA'),
                    end_sequence=('NB04_rev', 'TTCGGATTCTATCGTGTTTCCCTA')),
            Adapter('Barcode 5 (forward)',
                    start_sequence=('NB05', 'AAGGTTACACAAACCCTGGACAAG'),
                    end_sequence=('NB05_rev', 'CTTGTCCAGGGTTTGTGTAACCTT')),
            Adapter('Barcode 6 (forward)',
                    start_sequence=('NB06', 'GACTACTTTCTGCCTTTGCGAGAA'),
                    end_sequence=('NB06_rev', 'TTCTCGCAAAGGCAGAAAGTAGTC')),
            Adapter('Barcode 7 (forward)',
                    start_sequence=('NB07', 'AAGGATTCATTCCCACGGTAACAC'),
                    end_sequence=('NB07_rev', 'GTGTTACCGTGGGAATGAATCCTT')),
            Adapter('Barcode 8 (forward)',
                    start_sequence=('NB08', 'ACGTAACTTGGTTTGTTCCCTGAA'),
                    end_sequence=('NB08_rev', 'TTCAGGGAACAAACCAAGTTACGT')),
            Adapter('Barcode 9 (forward)',
                    start_sequence=('NB09', 'AACCAAGACTCGCTGTGCCTAGTT'),
                    end_sequence=('NB09_rev', 'AACTAGGCACAGCGAGTCTTGGTT')),
            Adapter('Barcode 10 (forward)',
                    start_sequence=('NB10', 'GAGAGGACAAAGGTTTCAACGCTT'),
                    end_sequence=('NB10_rev', 'AAGCGTTGAAACCTTTGTCCTCTC')),
            Adapter('Barcode 11 (forward)',
                    start_sequence=('NB11', 'TCCATTCCCTCCGATAGATGAAAC'),
                    end_sequence=('NB11_rev', 'GTTTCATCTATCGGAGGGAATGGA')),
            Adapter('Barcode 12 (forward)',
                    start_sequence=('NB12', 'TCCGATTCTGCTTCTTTCTACCTG'),
                    end_sequence=('NB12_rev', 'CAGGTAGAAAGAAGCAGAATCGGA')),
            Adapter('Barcode 13 (forward)',
                    start_sequence=('NB13', 'AGAACGACTTCCATACTCGTGTGA'),
                    end_sequence=('NB13_rev', 'TCACACGAGTATGGAAGTCGTTCT')),
            Adapter('Barcode 14 (forward)',
                    start_sequence=('NB14', 'AACGAGTCTCTTGGGACCCATAGA'),
                    end_sequence=('NB14_rev', 'TCTATGGGTCCCAAGAGACTCGTT')),
            Adapter('Barcode 15 (forward)',
                    start_sequence=('NB15', 'AGGTCTACCTCGCTAACACCACTG'),
                    end_sequence=('NB15_rev', 'CAGTGGTGTTAGCGAGGTAGACCT')),
            Adapter('Barcode 16 (forward)',
                    start_sequence=('NB16', 'CGTCAACTGACAGTGGTTCGTACT'),
                    end_sequence=('NB16_rev', 'AGTACGAACCACTGTCAGTTGACG')),
            Adapter('Barcode 17 (forward)',
                    start_sequence=('NB17', 'ACCCTCCAGGAAAGTACCTCTGAT'),
                    end_sequence=('NB17_rev', 'ATCAGAGGTACTTTCCTGGAGGGT')),
            Adapter('Barcode 18 (forward)',
                    start_sequence=('NB18', 'CCAAACCCAACAACCTAGATAGGC'),
                    end_sequence=('NB18_rev', 'GCCTATCTAGGTTGTTGGGTTTGG')),
            Adapter('Barcode 19 (forward)',
                    start_sequence=('NB19', 'GTTCCTCGTGCAGTGTCAAGAGAT'),
                    end_sequence=('NB19_rev', 'ATCTCTTGACACTGCACGAGGAAC')),
            Adapter('Barcode 20 (forward)',
                    start_sequence=('NB20', 'TTGCGTCCTGTTACGAGAACTCAT'),
                    end_sequence=('NB20_rev', 'ATGAGTTCTCGTAACAGGACGCAA')),
            Adapter('Barcode 21 (forward)',
                    start_sequence=('NB21', 'GAGCCTCTCATTGTCCGTTCTCTA'),
                    end_sequence=('NB21_rev', 'TAGAGAACGGACAATGAGAGGCTC')),
            Adapter('Barcode 22 (forward)',
                    start_sequence=('NB22', 'ACCACTGCCATGTATCAAAGTACG'),
                    end_sequence=('NB22_rev', 'CGTACTTTGATACATGGCAGTGGT')),
            Adapter('Barcode 23 (forward)',
                    start_sequence=('NB23', 'CTTACTACCCAGTGAACCTCCTCG'),
                    end_sequence=('NB23_rev', 'CGAGGAGGTTCACTGGGTAGTAAG')),
            Adapter('Barcode 24 (forward)',
                    start_sequence=('NB24', 'GCATAGTTCTGCATGATGGGTTAG'),
                    end_sequence=('NB24_rev', 'CTAACCCATCATGCAGAACTATGC')),
            Adapter('Barcode 25 (forward)',
                    start_sequence=('NB25', 'GTAAGTTGGGTATGCAACGCAATG'),
                    end_sequence=('NB25_rev', 'CATTGCGTTGCATACCCAACTTAC')),
            Adapter('Barcode 26 (forward)',
                    start_sequence=('NB26', 'CATACAGCGACTACGCATTCTCAT'),
                    end_sequence=('NB26_rev', 'ATGAGAATGCGTAGTCGCTGTATG')),
            Adapter('Barcode 27 (forward)',
                    start_sequence=('NB27', 'CGACGGTTAGATTCACCTCTTACA'),
                    end_sequence=('NB27_rev', 'TGTAAGAGGTGAATCTAACCGTCG')),
            Adapter('Barcode 28 (forward)',
                    start_sequence=('NB28', 'TGAAACCTAAGAAGGCACCGTATC'),
                    end_sequence=('NB28_rev', 'GATACGGTGCCTTCTTAGGTTTCA')),
            Adapter('Barcode 29 (forward)',
                    start_sequence=('NB29', 'CTAGACACCTTGGGTTGACAGACC'),
                    end_sequence=('NB29_rev', 'GGTCTGTCAACCCAAGGTGTCTAG')),
            Adapter('Barcode 30 (forward)',
                    start_sequence=('NB30', 'TCAGTGAGGATCTACTTCGACCCA'),
                    end_sequence=('NB30_rev', 'TGGGTCGAAGTAGATCCTCACTGA')),
            Adapter('Barcode 31 (forward)',
                    start_sequence=('NB31', 'TGCGTACAGCAATCAGTTACATTG'),
                    end_sequence=('NB31_rev', 'CAATGTAACTGATTGCTGTACGCA')),
            Adapter('Barcode 32 (forward)',
                    start_sequence=('NB32', 'CCAGTAGAAGTCCGACAACGTCAT'),
                    end_sequence=('NB32_rev', 'ATGACGTTGTCGGACTTCTACTGG')),
            Adapter('Barcode 33 (forward)',
                    start_sequence=('NB33', 'CAGACTTGGTACGGTTGGGTAACT'),
                    end_sequence=('NB33_rev', 'AGTTACCCAACCGTACCAAGTCTG')),
            Adapter('Barcode 34 (forward)',
                    start_sequence=('NB34', 'GGACGAAGAACTCAAGTCAAAGGC'),
                    end_sequence=('NB34_rev', 'GCCTTTGACTTGAGTTCTTCGTCC')),
            Adapter('Barcode 35 (forward)',
                    start_sequence=('NB35', 'CTACTTACGAAGCTGAGGGACTGC'),
                    end_sequence=('NB35_rev', 'GCAGTCCCTCAGCTTCGTAAGTAG')),
            Adapter('Barcode 36 (forward)',
                    start_sequence=('NB36', 'ATGTCCCAGTTAGAGGAGGAAACA'),
                    end_sequence=('NB36_rev', 'TGTTTCCTCCTCTAACTGGGACAT')),
            Adapter('Barcode 37 (forward)',
                    start_sequence=('NB37', 'GCTTGCGATTGATGCTTAGTATCA'),
                    end_sequence=('NB37_rev', 'TGATACTAAGCATCAATCGCAAGC')),
            Adapter('Barcode 38 (forward)',
                    start_sequence=('NB38', 'ACCACAGGAGGACGATACAGAGAA'),
                    end_sequence=('NB38_rev', 'TTCTCTGTATCGTCCTCCTGTGGT')),
            Adapter('Barcode 39 (forward)',
                    start_sequence=('NB39', 'CCACAGTGTCAACTAGAGCCTCTC'),
                    end_sequence=('NB39_rev', 'GAGAGGCTCTAGTTGACACTGTGG')),
            Adapter('Barcode 40 (forward)',
                    start_sequence=('NB40', 'TAGTTTGGATGACCAAGGATAGCC'),
                    end_sequence=('NB40_rev', 'GGCTATCCTTGGTCATCCAAACTA')),
            Adapter('Barcode 41 (forward)',
                    start_sequence=('NB41', 'GGAGTTCGTCCAGAGAAGTACACG'),
                    end_sequence=('NB41_rev', 'CGTGTACTTCTCTGGACGAACTCC')),
            Adapter('Barcode 42 (forward)',
                    start_sequence=('NB42', 'CTACGTGTAAGGCATACCTGCCAG'),
                    end_sequence=('NB42_rev', 'CTGGCAGGTATGCCTTACACGTAG')),
            Adapter('Barcode 43 (forward)',
                    start_sequence=('NB43', 'CTTTCGTTGTTGACTCGACGGTAG'),
                    end_sequence=('NB43_rev', 'CTACCGTCGAGTCAACAACGAAAG')),
            Adapter('Barcode 44 (forward)',
                    start_sequence=('NB44', 'AGTAGAAAGGGTTCCTTCCCACTC'),
                    end_sequence=('NB44_rev', 'GAGTGGGAAGGAACCCTTTCTACT')),
            Adapter('Barcode 45 (forward)',
                    start_sequence=('NB45', 'GATCCAACAGAGATGCCTTCAGTG'),
                    end_sequence=('NB45_rev', 'CACTGAAGGCATCTCTGTTGGATC')),
            Adapter('Barcode 46 (forward)',
                    start_sequence=('NB46', 'GCTGTGTTCCACTTCATTCTCCTG'),
                    end_sequence=('NB46_rev', 'CAGGAGAATGAAGTGGAACACAGC')),
            Adapter('Barcode 47 (forward)',
                    start_sequence=('NB47', 'GTGCAACTTTCCCACAGGTAGTTC'),
                    end_sequence=('NB47_rev', 'GAACTACCTGTGGGAAAGTTGCAC')),
            Adapter('Barcode 48 (forward)',
                    start_sequence=('NB48', 'CATCTGGAACGTGGTACACCTGTA'),
                    end_sequence=('NB48_rev', 'TACAGGTGTACCACGTTCCAGATG')),
            Adapter('Barcode 49 (forward)',
                    start_sequence=('NB49', 'ACTGGTGCAGCTTTGAACATCTAG'),
                    end_sequence=('NB49_rev', 'CTAGATGTTCAAAGCTGCACCAGT')),
            Adapter('Barcode 50 (forward)',
                    start_sequence=('NB50', 'ATGGACTTTGGTAACTTCCTGCGT'),
                    end_sequence=('NB50_rev', 'ACGCAGGAAGTTACCAAAGTCCAT')),
            Adapter('Barcode 51 (forward)',
                    start_sequence=('NB51', 'GTTGAATGAGCCTACTGGGTCCTC'),
                    end_sequence=('NB51_rev', 'GAGGACCCAGTAGGCTCATTCAAC')),
            Adapter('Barcode 52 (forward)',
                    start_sequence=('NB52', 'TGAGAGACAAGATTGTTCGTGGAC'),
                    end_sequence=('NB52_rev', 'GTCCACGAACAATCTTGTCTCTCA')),
            Adapter('Barcode 53 (forward)',
                    start_sequence=('NB53', 'AGATTCAGACCGTCTCATGCAAAG'),
                    end_sequence=('NB53_rev', 'CTTTGCATGAGACGGTCTGAATCT')),
            Adapter('Barcode 54 (forward)',
                    start_sequence=('NB54', 'CAAGAGCTTTGACTAAGGAGCATG'),
                    end_sequence=('NB54_rev', 'CATGCTCCTTAGTCAAAGCTCTTG')),
            Adapter('Barcode 55 (forward)',
                    start_sequence=('NB55', 'TGGAAGATGAGACCCTGATCTACG'),
                    end_sequence=('NB55_rev', 'CGTAGATCAGGGTCTCATCTTCCA')),
            Adapter('Barcode 56 (forward)',
                    start_sequence=('NB56', 'TCACTACTCAACAGGTGGCATGAA'),
                    end_sequence=('NB56_rev', 'TTCATGCCACCTGTTGAGTAGTGA')),
            Adapter('Barcode 57 (forward)',
                    start_sequence=('NB57', 'GCTAGGTCAATCTCCTTCGGAAGT'),
                    end_sequence=('NB57_rev', 'ACTTCCGAAGGAGATTGACCTAGC')),
            Adapter('Barcode 58 (forward)',
                    start_sequence=('NB58', 'CAGGTTACTCCTCCGTGAGTCTGA'),
                    end_sequence=('NB58_rev', 'TCAGACTCACGGAGGAGTAACCTG')),
            Adapter('Barcode 59 (forward)',
                    start_sequence=('NB59', 'TCAATCAAGAAGGGAAAGCAAGGT'),
                    end_sequence=('NB59_rev', 'ACCTTGCTTTCCCTTCTTGATTGA')),
            Adapter('Barcode 60 (forward)',
                    start_sequence=('NB60', 'CATGTTCAACCAAGGCTTCTATGG'),
                    end_sequence=('NB60_rev', 'CCATAGAAGCCTTGGTTGAACATG')),
            Adapter('Barcode 61 (forward)',
                    start_sequence=('NB61', 'AGAGGGTACTATGTGCCTCAGCAC'),
                    end_sequence=('NB61_rev', 'GTGCTGAGGCACATAGTACCCTCT')),
            Adapter('Barcode 62 (forward)',
                    start_sequence=('NB62', 'CACCCACACTTACTTCAGGACGTA'),
                    end_sequence=('NB62_rev', 'TACGTCCTGAAGTAAGTGTGGGTG')),
            Adapter('Barcode 63 (forward)',
                    start_sequence=('NB63', 'TTCTGAAGTTCCTGGGTCTTGAAC'),
                    end_sequence=('NB63_rev', 'GTTCAAGACCCAGGAACTTCAGAA')),
            Adapter('Barcode 64 (forward)',
                    start_sequence=('NB64', 'GACAGACACCGTTCATCGACTTTC'),
                    end_sequence=('NB64_rev', 'GAAAGTCGATGAACGGTGTCTGTC')),
            Adapter('Barcode 65 (forward)',
                    start_sequence=('NB65', 'TTCTCAGTCTTCCTCCAGACAAGG'),
                    end_sequence=('NB65_rev', 'CCTTGTCTGGAGGAAGACTGAGAA')),
            Adapter('Barcode 66 (forward)',
                    start_sequence=('NB66', 'CCGATCCTTGTGGCTTCTAACTTC'),
                    end_sequence=('NB66_rev', 'GAAGTTAGAAGCCACAAGGATCGG')),
            Adapter('Barcode 67 (forward)',
                    start_sequence=('NB67', 'GTTTGTCATACTCGTGTGCTCACC'),
                    end_sequence=('NB67_rev', 'GGTGAGCACACGAGTATGACAAAC')),
            Adapter('Barcode 68 (forward)',
                    start_sequence=('NB68', 'GAATCTAAGCAAACACGAAGGTGG'),
                    end_sequence=('NB68_rev', 'CCACCTTCGTGTTTGCTTAGATTC')),
            Adapter('Barcode 69 (forward)',
                    start_sequence=('NB69', 'TACAGTCCGAGCCTCATGTGATCT'),
                    end_sequence=('NB69_rev', 'AGATCACATGAGGCTCGGACTGTA')),
            Adapter('Barcode 70 (forward)',
                    start_sequence=('NB70', 'ACCGAGATCCTACGAATGGAGTGT'),
                    end_sequence=('NB70_rev', 'ACACTCCATTCGTAGGATCTCGGT')),
            Adapter('Barcode 71 (forward)',
                    start_sequence=('NB71', 'CCTGGGAGCATCAGGTAGTAACAG'),
                    end_sequence=('NB71_rev', 'CTGTTACTACCTGATGCTCCCAGG')),
            Adapter('Barcode 72 (forward)',
                    start_sequence=('NB72', 'TAGCTGACTGTCTTCCATACCGAC'),
                    end_sequence=('NB72_rev', 'GTCGGTATGGAAGACAGTCAGCTA')),
            Adapter('Barcode 73 (forward)',
                    start_sequence=('NB73', 'AAGAAACAGGATGACAGAACCCTC'),
                    end_sequence=('NB73_rev', 'GAGGGTTCTGTCATCCTGTTTCTT')),
            Adapter('Barcode 74 (forward)',
                    start_sequence=('NB74', 'TACAAGCATCCCAACACTTCCACT'),
                    end_sequence=('NB74_rev', 'AGTGGAAGTGTTGGGATGCTTGTA')),
            Adapter('Barcode 75 (forward)',
                    start_sequence=('NB75', 'GACCATTGTGATGAACCCTGTTGT'),
                    end_sequence=('NB75_rev', 'ACAACAGGGTTCATCACAATGGTC')),
            Adapter('Barcode 76 (forward)',
                    start_sequence=('NB76', 'ATGCTTGTTACATCAACCCTGGAC'),
                    end_sequence=('NB76_rev', 'GTCCAGGGTTGATGTAACAAGCAT')),
            Adapter('Barcode 77 (forward)',
                    start_sequence=('NB77', 'CGACCTGTTTCTCAGGGATACAAC'),
                    end_sequence=('NB77_rev', 'GTTGTATCCCTGAGAAACAGGTCG')),
            Adapter('Barcode 78 (forward)',
                    start_sequence=('NB78', 'AACAACCGAACCTTTGAATCAGAA'),
                    end_sequence=('NB78_rev', 'TTCTGATTCAAAGGTTCGGTTGTT')),
            Adapter('Barcode 79 (forward)',
                    start_sequence=('NB79', 'TCTCGGAGATAGTTCTCACTGCTG'),
                    end_sequence=('NB79_rev', 'CAGCAGTGAGAACTATCTCCGAGA')),
            Adapter('Barcode 80 (forward)',
                    start_sequence=('NB80', 'CGGATGAACATAGGATAGCGATTC'),
                    end_sequence=('NB80_rev', 'GAATCGCTATCCTATGTTCATCCG')),
            Adapter('Barcode 81 (forward)',
                    start_sequence=('NB81', 'CCTCATCTTGTGAAGTTGTTTCGG'),
                    end_sequence=('NB81_rev', 'CCGAAACAACTTCACAAGATGAGG')),
            Adapter('Barcode 82 (forward)',
                    start_sequence=('NB82', 'ACGGTATGTCGAGTTCCAGGACTA'),
                    end_sequence=('NB82_rev', 'TAGTCCTGGAACTCGACATACCGT')),
            Adapter('Barcode 83 (forward)',
                    start_sequence=('NB83', 'TGGCTTGATCTAGGTAAGGTCGAA'),
                    end_sequence=('NB83_rev', 'TTCGACCTTACCTAGATCAAGCCA')),
            Adapter('Barcode 84 (forward)',
                    start_sequence=('NB84', 'GTAGTGGACCTAGAACCTGTGCCA'),
                    end_sequence=('NB84_rev', 'TGGCACAGGTTCTAGGTCCACTAC')),
            Adapter('Barcode 85 (forward)',
                    start_sequence=('NB85', 'AACGGAGGAGTTAGTTGGATGATC'),
                    end_sequence=('NB85_rev', 'GATCATCCAACTAACTCCTCCGTT')),
            Adapter('Barcode 86 (forward)',
                    start_sequence=('NB86', 'AGGTGATCCCAACAAGCGTAAGTA'),
                    end_sequence=('NB86_rev', 'TACTTACGCTTGTTGGGATCACCT')),
            Adapter('Barcode 87 (forward)',
                    start_sequence=('NB87', 'TACATGCTCCTGTTGTTAGGGAGG'),
                    end_sequence=('NB87_rev', 'CCTCCCTAACAACAGGAGCATGTA')),
            Adapter('Barcode 88 (forward)',
                    start_sequence=('NB88', 'TCTTCTACTACCGATCCGAAGCAG'),
                    end_sequence=('NB88_rev', 'CTGCTTCGGATCGGTAGTAGAAGA')),
            Adapter('Barcode 89 (forward)',
                    start_sequence=('NB89', 'ACAGCATCAATGTTTGGCTAGTTG'),
                    end_sequence=('NB89_rev', 'CAACTAGCCAAACATTGATGCTGT')),
            Adapter('Barcode 90 (forward)',
                    start_sequence=('NB90', 'GATGTAGAGGGTACGGTTTGAGGC'),
                    end_sequence=('NB90_rev', 'GCCTCAAACCGTACCCTCTACATC')),
            Adapter('Barcode 91 (forward)',
                    start_sequence=('NB91', 'GGCTCCATAGGAACTCACGCTACT'),
                    end_sequence=('NB91_rev', 'AGTAGCGTGAGTTCCTATGGAGCC')),
            Adapter('Barcode 92 (forward)',
                    start_sequence=('NB92', 'TTGTGAGTGGAAAGATACAGGACC'),
                    end_sequence=('NB92_rev', 'GGTCCTGTATCTTTCCACTCACAA')),
            Adapter('Barcode 93 (forward)',
                    start_sequence=('NB93', 'AGTTTCCATCACTTCAGACTTGGG'),
                    end_sequence=('NB93_rev', 'CCCAAGTCTGAAGTGATGGAAACT')),
            Adapter('Barcode 94 (forward)',
                    start_sequence=('NB94', 'GATTGTCCTCAAACTGCCACCTAC'),
                    end_sequence=('NB94_rev', 'GTAGGTGGCAGTTTGAGGACAATC')),
            Adapter('Barcode 95 (forward)',
                    start_sequence=('NB95', 'CCTGTCTGGAAGAAGAATGGACTT'),
                    end_sequence=('NB95_rev', 'AAGTCCATTCTTCTTCCAGACAGG')),
            Adapter('Barcode 96 (forward)',
                    start_sequence=('NB96', 'CTGAACGGTCATAGAGTCCACCAT'),
                    end_sequence=('NB96_rev', 'ATGGTGGACTCTATGACCGTTCAG'))
        ]


def make_full_native_barcode_adapter(barcode_num):
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (reverse)'][0]
    start_barcode_seq = barcode.start_sequence[1]
    end_barcode_seq = barcode.end_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACGTATTGCTAAGGTTAA' + start_barcode_seq + 'CAGCACCT'
    end_full_seq = 'AGGTGCTG' + end_barcode_seq + 'TTAACCTTAGCAATACGTAACTGAACGAAGT'

    return Adapter('Native barcoding ' + str(barcode_num) + ' (full sequence)',
                   start_sequence=('NB' + '%02d' % barcode_num + '_start', start_full_seq),
                   end_sequence=('NB' + '%02d' % barcode_num + '_end', end_full_seq))


def make_old_full_rapid_barcode_adapter(barcode_num):  # applies to SQK-RBK001
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'TATTGCT' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence, old)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))


def make_new_full_rapid_barcode_adapter(barcode_num):  # applies to SQK-RBK004
    barcode = [x for x in ADAPTERS if x.name == 'Barcode ' + str(barcode_num) + ' (forward)'][0]
    start_barcode_seq = barcode.start_sequence[1]

    start_full_seq = 'AATGTACTTCGTTCAGTTACG' + 'GCTTGGGTGTTTAACC' + start_barcode_seq + \
                     'GTTTTCGCATTTATCGTGAAACGCTTTCGCGTTTTTCGTGCGCCGCTTCA'

    return Adapter('Rapid barcoding ' + str(barcode_num) + ' (full sequence, new)',
                   start_sequence=('RB' + '%02d' % barcode_num + '_full', start_full_seq))
