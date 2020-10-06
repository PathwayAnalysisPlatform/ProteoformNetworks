import os
import unittest

from lib.data_read_write import get_graph


class TestPlotModulePair(unittest.TestCase):
    def test_get_graph_correct_edges_gene(self):
        print("WD: ", os.getcwd())

        G = get_graph("Cleft Lip", "genes", "../../../reports/All_modules/")
        self.assertEqual(5, len(G.edges))
        self.assertTrue(("FGF10", "FGFR1") in G.edges)
        self.assertTrue(("PTPN11", "SPRY2") in G.edges)
        self.assertTrue(("NTN1", "PTPN11") in G.edges)
        self.assertTrue(("FGF10", "PTPN11") in G.edges)
        self.assertTrue(("FGFR1", "PTPN11") in G.edges)

    def test_get_graph_correct_vertices_gene(self):
        G = get_graph("Cleft Lip", "genes", "../../../reports/All_modules/")
        self.assertEqual(5, len(G.nodes))
        self.assertTrue("FGF10" in G.nodes)
        self.assertTrue("FGFR1" in G.nodes)
        self.assertTrue("NTN1" in G.nodes)
        self.assertTrue("PTPN11" in G.nodes)
        self.assertTrue("SPRY2" in G.nodes)

    def test_get_graph_invalid_level_raises_error(self):
        with self.assertRaises(ValueError):
            get_graph("Cleft Lip", "something", "../../../reports/All_modules/")

    def test_get_graph_correct_protein_level_graph(self):
        G = get_graph("\"Carcinoma, Basal Cell\"", "proteins", "../../../reports/All_modules/")
        self.assertEqual(4, len(G.nodes))
        self.assertEqual(2, len(G.edges))

        self.assertTrue("P14679" in G.nodes)
        self.assertTrue("Q14790" in G.nodes)
        self.assertTrue("Q8IUC6" in G.nodes)
        self.assertTrue("Q9UMX9" in G.nodes)

        self.assertTrue(("P14679", "Q9UMX9") in G.edges)
        self.assertTrue(("Q14790", "Q8IUC6") in G.edges)

    def test_get_graph_correct_proteoform_level_graph(self):
        G = get_graph("\"Lung Diseases, Interstitial\"", "proteoforms", "../../../reports/All_modules/")
        self.assertEqual(8, len(G.nodes))
        self.assertEqual(28, len(G.edges))

        self.assertTrue("P98088;" in G.nodes)
        self.assertTrue("P98088;00163:null" in G.nodes)
        self.assertTrue("Q02817;" in G.nodes)
        self.assertTrue("Q02817;00163:null" in G.nodes)
        self.assertTrue("Q6W4X9;" in G.nodes)
        self.assertTrue("Q9HC84;00163:null" in G.nodes)

        self.assertTrue(("P98088;", "P98088;00163:null") in G.edges)
        self.assertTrue(("P98088;", "Q02817;") in G.edges)
        self.assertTrue(("P98088;", "Q02817;00163:null") in G.edges)
        self.assertTrue(("P98088;", "Q6W4X9;") in G.edges)
        self.assertTrue(("P98088;", "Q6W4X9;00163:null") in G.edges)
        self.assertTrue(("P98088;", "Q9HC84;") in G.edges)
        self.assertTrue(("P98088;", "Q9HC84;00163:null") in G.edges)

        self.assertTrue(("P98088;00163:null", "Q02817;") in G.edges)
        self.assertTrue(("P98088;00163:null", "Q02817;00163:null") in G.edges)
        self.assertTrue(("P98088;00163:null", "Q6W4X9;") in G.edges)
        self.assertTrue(("P98088;00163:null", "Q6W4X9;00163:null") in G.edges)
        self.assertTrue(("P98088;00163:null", "Q9HC84;") in G.edges)
        self.assertTrue(("P98088;00163:null", "Q9HC84;00163:null") in G.edges)

        self.assertTrue(("Q02817;", "Q02817;00163:null") in G.edges)
        self.assertTrue(("Q02817;", "Q6W4X9;") in G.edges)
        self.assertTrue(("Q02817;", "Q6W4X9;00163:null") in G.edges)
        self.assertTrue(("Q02817;", "Q9HC84;") in G.edges)
        self.assertTrue(("Q02817;", "Q9HC84;00163:null") in G.edges)

        self.assertTrue(("Q02817;00163:null", "Q6W4X9;") in G.edges)
        self.assertTrue(("Q02817;00163:null", "Q6W4X9;00163:null") in G.edges)
        self.assertTrue(("Q02817;00163:null", "Q9HC84;") in G.edges)
        self.assertTrue(("Q02817;00163:null", "Q9HC84;00163:null") in G.edges)

        self.assertTrue(("Q6W4X9;", "Q6W4X9;00163:null") in G.edges)
        self.assertTrue(("Q6W4X9;", "Q9HC84;") in G.edges)
        self.assertTrue(("Q6W4X9;", "Q9HC84;00163:null") in G.edges)

        self.assertTrue(("Q6W4X9;00163:null", "Q9HC84;") in G.edges)
        self.assertTrue(("Q6W4X9;00163:null", "Q9HC84;00163:null") in G.edges)

        self.assertTrue(("Q9HC84;", "Q9HC84;00163:null") in G.edges)

    if __name__ == '__main__':
        unittest.main()
