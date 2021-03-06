from complexes.utils.get_annotated_name import GetAnnotatedName
from unittest import TestCase


expected_result_one = [
    {
        "accession": "P43773",
        "stoichiometry": "12",
        "accession_stoichiometry": "P43773_12",
    },
    {
        "accession": "P43772",
        "stoichiometry": "12",
        "accession_stoichiometry": "P43772_12",
    },
]

expected_result_two = [
    {
        "accession": "P19483",
        "stoichiometry": "3",
        "accession_stoichiometry": "P19483_3",
    },
    {
        "accession": "P00829",
        "stoichiometry": "3",
        "accession_stoichiometry": "P00829_3",
    },
    {
        "accession": "P05631",
        "stoichiometry": "1",
        "accession_stoichiometry": "P05631_1",
    },
    {
        "accession": "P05630",
        "stoichiometry": "1",
        "accession_stoichiometry": "P05630_1",
    },
    {
        "accession": "P05632",
        "stoichiometry": "1",
        "accession_stoichiometry": "P05632_1",
    },
    {
        "accession": "P32876",
        "stoichiometry": "8",
        "accession_stoichiometry": "P32876_8",
    },
    {
        "accession": "P13621",
        "stoichiometry": "1",
        "accession_stoichiometry": "P13621_1",
    },
    {
        "accession": "P13619",
        "stoichiometry": "1",
        "accession_stoichiometry": "P13619_1",
    },
    {
        "accession": "P13620",
        "stoichiometry": "1",
        "accession_stoichiometry": "P13620_1",
    },
    {
        "accession": "P02721",
        "stoichiometry": "1",
        "accession_stoichiometry": "P02721_1",
    },
    {
        "accession": "P00847",
        "stoichiometry": "1",
        "accession_stoichiometry": "P00847_1",
    },
]


class TestGetAnnotatedName(TestCase):
    # TODO Convert all the other tests to this format

    def setUp(self) -> None:
        self.gan = GetAnnotatedName()
        self.gan.get_data()

    def test_read_molecule_names(self):
        # Test if method returns the correct name
        self.assertEqual("AAA+ proteases", self.gan.molecule_names["1"])
        self.assertEqual("ATP Synthase", self.gan.molecule_names["27"])

    def test_read_components(self):
        # Test if method returns the correct components for a complex
        self.assertEqual(expected_result_one, self.gan.molecule_components["1"])
        self.assertEqual(expected_result_two, self.gan.molecule_components["27"])
