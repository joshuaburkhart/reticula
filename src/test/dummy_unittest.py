# coding=utf-8
import unittest

from hw3.src.lib.ScoringFunction import entropy, sum_of_pairs
from hw3.src.lib.ScoringMatrix import blosum_62

__author__ = 'burkhart'

S_FUNC_EPSILON = 0.01


class TestScoringFunction(unittest.TestCase):
    def setUp(self):
        pass

    def tearDown(self):
        pass

    def te_helper(self, bp_alignment, expect):
        actual = entropy(bp_alignment, ())
        self.assertLess(abs(actual - expect), S_FUNC_EPSILON,
                        "Incorrect entropy of {0}: {1} "
                        "(expected {2} to be within {3})".format(
                            bp_alignment,
                            actual,
                            expect,
                            S_FUNC_EPSILON))

    def sp_helper(self, bp_alignment, expect, *matrix_wrapper):
        actual = sum_of_pairs(bp_alignment, matrix_wrapper)
        self.assertLess(abs(actual - expect), S_FUNC_EPSILON,
                        "Incorrect sop of {0}: {1} "
                        "(expected {2} to be within {3})".format(
                            bp_alignment,
                            actual,
                            expect,
                            S_FUNC_EPSILON))

    def test_entropy_1(self):
        """Entropy Compeau, P, Pevzner, P,c Active
        Learning Publishers,Bioinformatics Algorithms:
        An Active Learning Approach,2nd Ed., Vol. I, p. 76,
        Example 1
        """
        self.te_helper('CCCTATCCAC', 1.371)

    def test_entropy_2(self):
        """Entropy Compeau, P, Pevzner, P,c Active
        Learning Publishers,Bioinformatics Algorithms:
        An Active Learning Approach,2nd Ed., Vol. I, p. 76,
        Example 2
        """
        self.te_helper('TCCTCCTTCC', 0.971)

    def test_entropy_3(self):
        """Entropy Compeau, P, Pevzner, P,c Active
        Learning Publishers,Bioinformatics Algorithms:
        An Active Learning Approach,2nd Ed., Vol. I, p. 76,
        Example 3
        """
        self.te_helper('GTGGGGGGGG', 0.467)

    def test_entropy_minimum_possible(self):
        """Entropy Compeau, P, Pevzner, P,c Active
        Learning Publishers,Bioinformatics Algorithms:
        An Active Learning Approach,2nd Ed., Vol. I, p. 76,
        Last Paragraph
        """
        self.te_helper('GGGGGGGGGG', 0.0)

    def test_entropy_maximum_possible(self):
        """Entropy Compeau, P, Pevzner, P,c Active
        Learning Publishers,Bioinformatics Algorithms:
        An Active Learning Approach,2nd Ed., Vol. I, p. 76,
        Last Paragraph
        """
        self.te_helper('AATTCCGG', 2.0)

    def test_sum_of_pairs_1(self):
        """http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        AR = -1
        AN = -2
        RN = 0
        """
        self.sp_helper('ARN', -3, blosum_62)

    def test_sum_of_pairs_2(self):
        """http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        DC = -3
        DQ = 0
        CQ = -3
        """
        self.sp_helper('DCQ', -6, blosum_62)

    def test_sum_of_pairs_3(self):
        """http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        EG = -2
        EH = 0
        GH = -2
        """
        self.sp_helper('EGH', -4, blosum_62)

    def test_sum_of_pairs_partial_1(self):
        """http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        IL = 2
        """
        self.sp_helper('IL-', 2, blosum_62)

    def test_sum_of_pairs_partial_2(self):
        """http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        FP = -4
        """
        self.sp_helper('-FP', -4, blosum_62)

    def test_sum_of_pairs_partial_3(self):
        """http://www.ncbi.nlm.nih.gov/Class/FieldGuide/BLOSUM62.txt
        SW = -3
        """
        self.sp_helper('S-W', -3, blosum_62)

    def test_sum_of_pairs_invalid_1(self):
        self.sp_helper('Y--', 0, blosum_62)

    def test_sum_of_pairs_invalid_2(self):
        self.sp_helper('-V-', 0, blosum_62)

    def test_sum_of_pairs_invalid_3(self):
        self.sp_helper('--B', 0, blosum_62)


if __name__ == '__main__':
    unittest.main()