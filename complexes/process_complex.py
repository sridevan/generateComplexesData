import argparse

import csv
import hashlib

# import re
from collections import OrderedDict

# import cx_Oracle
# import neo4j

from py2neo import Graph
import os
import time

# from pdbecommon.graph import neo4j
# import pandas as pd

# from get_complex_name import ProcessComplexName

# from utils.database import create_db_connection
# from orc.base.utils import close_db_connection, create_db_connection
# from utils.get_annotated_name import GetAnnotatedName
# from utils.get_data_from_graph_db import GetComplexData

# from datetime import datetime


# from orc.base.log import logger

MERGE_ACCESSION_QUERY = """
WITH $accession_params_list AS batch
UNWIND batch AS row
MATCH (u:UniProt {ACCESSION:row.accession})
MERGE (c:PDBComplex {COMPLEX_ID:row.complex_id})
MERGE (c)<-[:IS_PART_OF_PDB_COMPLEX {STOICHIOMETRY:row.stoichiometry}]-(u)
"""

MERGE_ENTITY_QUERY = """
WITH $entity_params_list AS batch
UNWIND batch AS row
MATCH (e:Entry {ID:row.entry_id})-[:HAS_ENTITY]->(en:Entity {ID:row.entity_id})
MERGE (c:PDBComplex {COMPLEX_ID:row.complex_id})
MERGE (c)<-[:IS_PART_OF_PDB_COMPLEX {STOICHIOMETRY:row.stoichiometry}]-(en)
"""

MERGE_ASSEMBLY_QUERY = """
WITH $assembly_params_list AS batch
UNWIND batch AS row
MATCH (e:Entry {ID:row.entry_id})-[:HAS_ENTITY]->(:Entity)-
    [:IS_PART_OF_ASSEMBLY]->(assembly:Assembly {UNIQID:row.assembly_id})
MERGE (c:PDBComplex {COMPLEX_ID:row.complex_id})
MERGE (c)<-[:IS_PART_OF_PDB_COMPLEX]-(assembly)
"""

MERGE_RFAM_QUERY = """
WITH $rfam_params_list AS batch
UNWIND batch AS row
MATCH (rfam:RfamFamily {RFAM_ACC:row.rfam_acc})
MERGE (c:PDBComplex {COMPLEX_ID:row.complex_id})
MERGE (c)<-[:IS_PART_OF_PDB_COMPLEX]-(rfam)
"""

COMMON_COMPLEX_QUERY = """
WITH $complex_params_list AS batch
UNWIND batch AS row
MATCH
    (p:PDBComplex {COMPLEX_ID:row.pdb_complex_id}),
    (c:Complex {COMPLEX_ID:row.complex_portal_id})
CREATE (p)-[:SAME_AS]->(c)
"""

MERGE_UNMAPPED_POLYMER_QUERY = """
WITH $unmapped_polymer_params_list AS batch
UNWIND batch AS row
MERGE (up:UnmappedPolymer {TYPE:row.polymer_type})
MERGE (c:PDBComplex {COMPLEX_ID:row.complex_id})
MERGE (c)<-[:IS_PART_OF_PDB_COMPLEX]-(up)
"""

COMPLEX_PORTAL_DATA_QUERY = """
MATCH
    (complex:Complex)<-[rel:IS_PART_OF_COMPLEX]-(unp:UniProt)-[:HAS_TAXONOMY]->(tax:Taxonomy)
OPTIONAL MATCH
    (complex)<-[:IS_PART_OF_COMPLEX]-(entry:Entry)
WITH
    complex.COMPLEX_ID AS complex_id,
    unp.ACCESSION +'_' + rel.STOICHIOMETRY +'_' +tax.TAX_ID AS uniq_accessions,
    COLLECT(entry.ID) AS entries ORDER BY uniq_accessions
WITH
    complex_id AS complex_id,
    COLLECT(DISTINCT uniq_accessions) AS uniq_accessions,
    entries
WITH
    complex_id AS complex_id,
    REDUCE(s = HEAD(uniq_accessions),
    n in TAIL(uniq_accessions) | s +',' +n) AS uniq_accessions,
    entries
RETURN
    complex_id,
    uniq_accessions,
    REDUCE(s = HEAD(entries), n in TAIL(entries) | s +',' +n) AS entries_str
LIMIT 50
"""

PDB_ASSEMBLY_DATA_QUERY = """
MATCH (assembly:Assembly {PREFERED: 'True'})<-[rel:IS_PART_OF_ASSEMBLY]
-(entity:Entity {TYPE:'p'})
WITH assembly, rel, entity
OPTIONAL MATCH (entity)-[:HAS_UNIPROT {BEST_MAPPING:'1'}]->(uniprot:UniProt)
-[:HAS_TAXONOMY]->(tax:Taxonomy)
OPTIONAL MATCH (entity)-[:HAS_RFAM]->(rfam:RfamFamily)
WITH assembly.UNIQID AS assembly_id,
CASE uniprot
    WHEN null
        THEN
            CASE rfam
                WHEN null
                    THEN
                        CASE entity.POLYMER_TYPE
                            WHEN 'R'
                                THEN 'RNA' +':UNMAPPED'
                            WHEN 'D'
                                THEN 'DNA' +':UNMAPPED'
                            WHEN 'D/R'
                                THEN 'DNA/RNA' +':UNMAPPED'
                            WHEN 'P'
                                THEN 'NA_' +entity.UNIQID +'_' +rel.NUMBER_OF_CHAINS
                        END
                ELSE
                    rfam.RFAM_ACC
            END
    ELSE uniprot.ACCESSION +'_' +rel.NUMBER_OF_CHAINS +'_' +tax.TAX_ID
END AS accession ORDER BY accession
WITH assembly_id AS assembly_id, COLLECT (DISTINCT accession) AS accessions
WITH assembly_id AS assembly_id, REDUCE(s = HEAD(accessions),
n in TAIL(accessions) | s +',' +n) AS accessions
WITH accessions, COLLECT(DISTINCT assembly_id) AS assemblies
WITH accessions AS accessions, REDUCE(s = HEAD(assemblies),
n in TAIL(assemblies) | s +',' +n) AS assemblies
RETURN accessions, assemblies
"""

DROP_PDB_COMPLEX_NODES_QUERY = """
MATCH (p:PDBComplex) DETACH DELETE p
"""

DROP_SUBCOMPLEX_RELATION_QUERY = """
MATCH (:PDBComplex)-[r:IS_SUB_COMPLEX_OF]->(:PDBComplex) DELETE r
"""

CREATE_SUBCOMPLEX_RELATION_QUERY = """
MATCH
(src_complex:PDBComplex)<-[rel1:IS_PART_OF_PDB_COMPLEX]-()-
[rel2:IS_PART_OF_PDB_COMPLEX]->(dest_complex:PDBComplex)
WHERE rel1.STOICHIOMETRY=rel2.STOICHIOMETRY
    WITH DISTINCT src_complex, dest_complex, rel1
    WITH src_complex, startNode(rel1) AS relRelations, dest_complex
    WITH src_complex, COUNT(relRelations) AS relRelationsAmount, dest_complex
    MATCH (src_complex)<-[allRelations:IS_PART_OF_PDB_COMPLEX]-()
        WITH src_complex, relRelationsAmount, count(allRelations) AS allRelationsAmount,
            dest_complex
                WHERE relRelationsAmount = allRelationsAmount
CREATE (dest_complex)<-[:IS_SUB_COMPLEX_OF]-(src_complex)
"""


class Neo4JProcessComplex:
    def __init__(self, bolt_uri, username, password, csv_path):
        self.graph = None
        self.bolt_uri = bolt_uri
        self.username = username
        self.password = password
        self.csv_path = csv_path
        # self._driver = neo4j.GraphDatabase.driver(bolt_uri, auth=(username, password))
        # self.session = self._driver.session()
        # self.complex_subcomplex_outcsv = outcsv
        self.dict_complex_portal_id = {}
        self.dict_complex_portal_entries = {}
        self.dict_pdb_complex = {}
        self.common_complexes = []
        self.reference_mapping = OrderedDict()
        # self.db_conn_str = db_connection_str
        self.accession_params_list = []
        self.entity_params_list = []
        self.assembly_params_list = []
        self.rfam_params_list = []
        self.complex_params_list = []
        self.unmapped_polymer_params_list = []

    def close(self):
        self._driver.close()

    def use_persistent_identifier(
        self, hash_str, accession, complex_portal_id, entries
    ):
        """
        This method generates a new PDB complex ID if the complex composition
        is new and not present in reference_mapping. On the other hand, if it's
        an existing complex composition then the existing ID is returned

        Args:
            hash_str (hash): md5 hash obj of complex composition
            accession (string): complex composition string
            complex_portal_id (string): Complex Portal identifier
            entries (string): assemblies

        Returns:
            string: PDB Complex ID
        """
        basic_PDB_complex_str = "PDB-CPX-"
        initial_num = 100001
        if hash_str in self.reference_mapping:
            pdb_complex_id = self.reference_mapping.get(hash_str).get("pdb_complex_id")
        # when the dict is empty
        elif len(self.reference_mapping) == 0:
            # identifier = len(self.reference_mapping) + 1
            pdb_complex_id = basic_PDB_complex_str + str(initial_num)
            # print("This is new complex id {pdb_complex_id}")
            # print(complex_portal_id)
            self.reference_mapping[hash_str] = {
                "pdb_complex_id": pdb_complex_id,
                "complex_portal_id": complex_portal_id,
                "accession": accession,
                "entries": entries,
            }
        # when the dict has one or more elems
        elif len(self.reference_mapping) >= 1:
            # last_dict_key = list(self.reference_mapping.keys())[-1]
            last_dict_key = next(reversed(self.reference_mapping))
            last_pdb_complex_id = self.reference_mapping[last_dict_key][
                "pdb_complex_id"
            ]
            _, _, last_pdb_complex_id_num = last_pdb_complex_id.split("-")
            current_num = int(last_pdb_complex_id_num) + 1
            pdb_complex_id = basic_PDB_complex_str + str(current_num)
            self.reference_mapping[hash_str] = {
                "pdb_complex_id": pdb_complex_id,
                "complex_portal_id": complex_portal_id,
                "accession": accession,
                "entries": entries,
            }

        return pdb_complex_id

    # def check_reference_data(self):
    #     conn = create_db_connection(self.db_conn_str)

    #     with conn as connection:
    #         query = connection.execute("SELECT * FROM COMPLEX_REFERENCE_MAPPING")

    #         if query:
    #             for row in query:
    #                 self.reference_mapping[row[0]] = {"pdb_complex_id": row[2],
    #                                             "accession": row[1]}

    #             return self.reference_mapping

    # def export_data(self):
    #     data = [(k, v['accession'], v['pdb_complex_id']) for k,
    #              v in self.reference_mapping.items()]
    #     conn = create_db_connection(self.db_conn_str)

    #     with conn as connection:
    #         empty_table_command = "TRUNCATE TABLE COMPLEX_REFERENCE_MAPPING"
    #         insert_data_command = "INSERT INTO COMPLEX_REFERENCE_MAPPING VALUES (:1, :2, :3)"
    #         connection.execute(empty_table_command)
    #         for row in data:
    #             connection.execute(insert_data_command, row)

    def export_csv(self):
        """
        Generate a csv file that contains complexes related information such as
        complex_composition, pdb_complex_id, complex_portal_id, assemblies and
        md5_obj representing the complex_composition
        """
        headers = [  # noqa: F841
            "pdb_complex_id",
            "complex_portal_id",
            "accession",
            "entries",
        ]
        # base_path = Path.cwd()
        base_path = self.csv_path
        filename = "complexes_mapping.csv"
        # csv_filename = (
        #     base_path.joinpath(filename)
        # )  # noqa: F841
        complete_path = os.path.join(base_path, filename)
        with open(complete_path, "w", newline="") as reference_file:
            file_csv = csv.writer(reference_file)
            file_csv.writerow(["md5_obj", *headers])
            for key, val in self.reference_mapping.items():
                file_csv.writerow([key] + [val.get(i, "") for i in headers])
        print("Complexes_mapping file has been produced")

    def get_graph(self):
        """
        Create an instance of Neo4j with the given connection
        parameters
        """
        self.graph = Graph(self.bolt_uri, user=self.username, password=self.password)

    def run_query(self, query, param=None):
        """
        Run Neo4j query

        Args:
            query (string): Neo4j query

        Returns:
            object: Neo4j query result
        """
        if not self.graph:
            self.get_graph()
        if param is not None:
            return self.graph.run(query, parameters=param)
        else:
            return self.graph.run(query)

    def get_complex_portal_data(self):
        """
        Gets and processes Complex Portal data - Complex Portal ID,
        complex composition str and entries
        """
        print("Querying Complex Portal data")
        mappings = self.run_query(query=COMPLEX_PORTAL_DATA_QUERY)
        # rdict = [rec.data() for rec in mappings]
        # print(rdict)
        for row in mappings:
            accessions = row.get("uniq_accessions")
            complex_id = row.get("complex_id")
            entries = row.get("entries_str")
            self.dict_complex_portal_id[accessions] = complex_id
            self.dict_complex_portal_entries[complex_id] = entries

    def create_nodes_relationship(
        self, query_name, n1_name, n2_name, param_name, param_val
    ):
        """
        Runs Neo4j query to create a relationship between a given pair of nodes

        Args:
            query_name (string): Neo4j query
            n1_name (string): first node name
            n2_name (string): second node name
            param_name (string): parameter name
            param_val (string): parameter value
        """
        print(f"Creating relationship between {n1_name} and {n2_name} nodes - START")
        self.run_query(
            query_name,
            param={param_name: param_val},
        )
        print(f"Creating relationship between {n1_name} and {n2_name} nodes - DONE")

    def create_subcomplex_relationships(self):
        """
        Create subcomplex relationships in the graph db
        """
        return self.run_query(CREATE_SUBCOMPLEX_RELATION_QUERY)

    def drop_PDBComplex_nodes(self):
        """
        Drop any existing PDB complex nodes in the graph db
        """
        return self.run_query(DROP_PDB_COMPLEX_NODES_QUERY)

    def drop_existing_subcomplex_relationships(self):
        """
        Drop any existing subcomplex relationships in the graph db
        """
        return self.run_query(DROP_SUBCOMPLEX_RELATION_QUERY)

    def process_assembly_data(self):
        """
        Aggregate unique complex compositions from PDB data, compares them to
        Complex Portal data and processes them for use later.
        """
        # Had to include this to run in local machine. We don\t need this now
        # cx_Oracle.init_oracle_client(lib_dir="/Users/sria/Downloads/instantclient_19_8")

        # Update self.reference_mapping if reference is available
        # self.reference_mapping = self.check_reference_data()
        print("Querying PDB Assembly data")
        mappings = self.run_query(PDB_ASSEMBLY_DATA_QUERY)
        for row in mappings:
            # count += 1
            # (uniq_accessions, assemblies) = row

            uniq_accessions = row.get("accessions")
            assemblies = row.get("assemblies")

            # remove all occurences of NA_ from the unique complex combination
            tmp_uniq_accessions = uniq_accessions.replace("NA_", "")

            # pdb_complex_id = basic_complex_string + str(uniq_id)
            accession_hash = hashlib.md5(
                tmp_uniq_accessions.encode("utf-8")
            ).hexdigest()
            complex_portal_id = self.dict_complex_portal_id.get(tmp_uniq_accessions)

            # print(f"Printing Complex Portal ID: {complex_portal_id}")
            pdb_complex_id = self.use_persistent_identifier(
                accession_hash, tmp_uniq_accessions, complex_portal_id, assemblies
            )

            # common complex; delete from dictionary else will be processed again
            if complex_portal_id is not None:
                del self.dict_complex_portal_id[tmp_uniq_accessions]
                self.common_complexes.append((pdb_complex_id, complex_portal_id))

            # keep data for each PDB complex in dict_pdb_complex to be used later
            self.dict_pdb_complex[pdb_complex_id] = (
                tmp_uniq_accessions,
                assemblies,
            )

            for uniq_accession in uniq_accessions.split(","):
                tokens = uniq_accession.split("_")

                # handle cases of PDB entity
                if len(tokens) == 4:
                    [_, entry_id, entity_id, stoichiometry] = tokens
                    self.entity_params_list.append(
                        {
                            "complex_id": str(pdb_complex_id),
                            "entry_id": str(entry_id),
                            "entity_id": str(entity_id),
                            "stoichiometry": str(stoichiometry),
                        }
                    )

                # handle cases of UniProt
                elif len(tokens) == 3:
                    [accession, stoichiometry, tax_id] = tokens
                    self.accession_params_list.append(
                        {
                            "complex_id": str(pdb_complex_id),
                            "accession": str(accession),
                            "stoichiometry": str(stoichiometry),
                        }
                    )

                # handle unmapped polymers and Rfam accessions
                elif len(tokens) == 1:
                    token = tokens[0]

                    # check for unmapped polymers (:UNMAPPED string)
                    if ":UNMAPPED" in token:
                        polymer_type = token.replace(":UNMAPPED", "")
                        self.unmapped_polymer_params_list.append(
                            {
                                "complex_id": str(pdb_complex_id),
                                "polymer_type": str(polymer_type),
                            }
                        )

                    # handle Rfam
                    else:
                        self.rfam_params_list.append(
                            {
                                "complex_id": str(pdb_complex_id),
                                "rfam_acc": str(token),
                            }
                        )

            for uniq_assembly in assemblies.split(","):
                # [entry, assembly_id] = uniq_assembly.split("_")
                [entry, _] = uniq_assembly.split("_")
                self.assembly_params_list.append(
                    {
                        "complex_id": str(pdb_complex_id),
                        "assembly_id": str(uniq_assembly),
                        "entry_id": str(entry),
                    }
                )

                # uniq_id += 1
            # print(f"Number of records: {count}")

        # create list of common Complex and PDB_Complex nodes
        for common_complex in self.common_complexes:
            (pdb_complex_id, complex_portal_id) = common_complex
            self.complex_params_list.append(
                {
                    "pdb_complex_id": str(pdb_complex_id),
                    "complex_portal_id": str(complex_portal_id),
                }
            )

        print("Done querying PDB Assembly data")
        # print(self.csv_path)

        # for accessions in self.dict_complex_portal_id.keys():
        #     accession_hash = hashlib.md5(accessions.encode("utf-8")).hexdigest()
        #     complex_portal_id = self.dict_complex_portal_id[accessions]
        #     entries = self.dict_complex_portal_entries.get(complex_portal_id)
        #     pdb_complex_id = self.use_persistent_identifier(
        #         accession_hash, accessions, complex_portal_id, entries
        #     )

        #     # keep data for each PDB complex in dict_pdb_complex to be used later
        #     self.dict_pdb_complex[pdb_complex_id] = (accessions, entries)

        #     for item in accessions.split(","):
        #         [accession, stoichiometry, tax_id] = item.split("_")

        #         # this is the data from complex portal, there won't be any PDB entity
        #         # as a participant
        #         accession_params_list.append(
        #             {
        #                 "complex_id": str(pdb_complex_id),
        #                 "accession": str(accession),
        #                 "stoichiometry": str(stoichiometry),
        #             }
        #         )

    def create_graph_relationships(self):
        """
        Create relationships between complexes related nodes in the
        graph db - Uniprot, PDBComplex, Entity, Rfam, Unmapped
        Polymer, Assembly and Complex
        """
        print("Start creating relationships betweeen nodes")
        # Create relationship between Uniprot and PDBComplex nodes
        self.create_nodes_relationship(
            query_name=MERGE_ACCESSION_QUERY,
            n1_name="Uniprot",
            n2_name="PDBComplex",
            param_name="accession_params_list",
            param_val=self.accession_params_list,
        )
        # Create relationship between Entity and PDBComplex nodes
        self.create_nodes_relationship(
            query_name=MERGE_ENTITY_QUERY,
            n1_name="Entity",
            n2_name="PDBComplex",
            param_name="entity_params_list",
            param_val=self.entity_params_list,
        )
        # Create relationship between UnmappedPolymer and PDBComplex nodes
        self.create_nodes_relationship(
            query_name=MERGE_UNMAPPED_POLYMER_QUERY,
            n1_name="UnmappedPolymer",
            n2_name="PDBComplex",
            param_name="unmapped_polymer_params_list",
            param_val=self.unmapped_polymer_params_list,
        )
        # Create relationship between Rfam and PDBComplex nodes
        self.create_nodes_relationship(
            query_name=MERGE_RFAM_QUERY,
            n1_name="Rfam",
            n2_name="PDBComplex",
            param_name="rfam_params_list",
            param_val=self.rfam_params_list,
        )
        # Create relationship between Assembly and PDBComplex nodes
        self.create_nodes_relationship(
            query_name=MERGE_ASSEMBLY_QUERY,
            n1_name="Assembly",
            n2_name="PDBComplex",
            param_name="assembly_params_list",
            param_val=self.assembly_params_list,
        )
        # Create relationship between PDBComplex and Complex nodes
        self.create_nodes_relationship(
            query_name=COMMON_COMPLEX_QUERY,
            n1_name="PDBComplex",
            n2_name="Complex",
            param_name="complex_params_list",
            param_val=self.complex_params_list,
        )
        print("Done creating relationships between nodes")

    # def generate_complexes_names(self):
    #     complex_data = GetComplexData(self.bolt_uri, self.username, self.password)
    #     complex_data.get_pdb_complex_data()

    #     complex_names = ProcessComplexName(complex_data.pdb_complexes)
    #     complex_names.process_complex_name()
    #     print("Complexes name file has been produced")

    # def create_and_merge_csv_files(self):
    #     self.export_csv()
    #     self.merge_csv_files()

    # def process_subcomplex_data(self):

    #     print(
    #         f"Processing Complex-Subcomplex relationships - Started at {datetime.now()}"
    #     )

    #     query = """
    #     MATCH (complex:PDBComplex)<-[:IS_SUB_COMPLEX_OF]-(sub_complex:PDBComplex)
    #     WITH complex AS complex, sub_complex AS sub_complex
    #     MATCH
    #         (complex)<-[rel1:IS_PART_OF_PDB_COMPLEX]-(u1),
    #         (sub_complex)<-[rel2:IS_PART_OF_PDB_COMPLEX]-(u2)
    #     WITH complex.COMPLEX_ID AS complex_id,
    #     CASE u1.ACCESSION
    #         WHEN null
    #         THEN u1.UNIQID +'_' +rel1.STOICHIOMETRY +'_' +u1.POLYMER_TYPE
    #         ELSE u1.ACCESSION +'_' +rel1.STOICHIOMETRY
    #     END AS unique_complex,
    #     sub_complex.COMPLEX_ID AS sub_complex_id,
    #     CASE u2.ACCESSION
    #         WHEN null
    #         THEN u2.UNIQID +'_' +rel2.STOICHIOMETRY +'_' +u2.POLYMER_TYPE
    #         ELSE u2.ACCESSION +'_' +rel2.STOICHIOMETRY
    #     END AS unique_sub_complex
    #     ORDER BY sub_complex, unique_sub_complex
    #     WITH
    #         complex_id AS complex_id,
    #         COLLECT(DISTINCT unique_complex) AS unique_complex,
    #         sub_complex_id AS sub_complex_id,
    #         COLLECT(DISTINCT unique_sub_complex) AS unique_sub_complex
    #     RETURN
    #         complex_id,
    #         REDUCE(s = HEAD(unique_complex), n in TAIL(unique_complex) | s +',' +n)
    #             AS unique_complex,
    #         sub_complex_id,
    #         REDUCE(s = HEAD(unique_sub_complex), n in TAIL(unique_sub_complex) | s +',' +n)
    #             AS unique_sub_complex
    #     """

    #     with self._driver.session() as session:
    #         mappings = session.run(query)

    #         with open(self.complex_subcomplex_outcsv, "w") as complex_subcomplex_file:
    #             complex_subcomplex_file_csv = csv.writer(
    #                 complex_subcomplex_file, dialect="excel"
    #             )
    #             complex_subcomplex_file_csv.writerow(
    #                 (
    #                     "PDB_COMPLEX",
    #                     "PDB_COMPLEX_PARTICIPANTS",
    #                     "PDB_SUBCOMPLEX",
    #                     "PDB_SUBCOMPLEX_PARTICIPANTS",
    #                     "PDB_ENTRIES",
    #                 )
    #             )

    #             for row in mappings:
    #                 (
    #                     complex_id,
    #                     unique_complex,
    #                     sub_complex_id,
    #                     unique_sub_complex,
    #                 ) = row

    #                 (participants, assemblies) = self.dict_pdb_complex.get(complex_id)

    #                 complex_subcomplex_file_csv.writerow(
    #                     (
    #                         complex_id,
    #                         unique_complex,
    #                         sub_complex_id,
    #                         unique_sub_complex,
    #                         assemblies,
    #                     )
    #                 )

    #     print(
    #         f"Processing Complex-Subcomplex relationships - Ended at {datetime.now()}"
    #     )

    # def find_subcomplexes(self):

    #     print(
    #         f"Checking for sub complexes and making relationships - Started now"
    #     )

    #     # drop existing IS_SUB_COMPLEX_OF relationship, if any
    #     print("Dropping IS_SUB_COMPLEX_OF relationships, if any - START")

    #     with self._driver.session() as session:
    #         session.run(
    #             "MATCH (:PDBComplex)-[r:IS_SUB_COMPLEX_OF]->(:PDBComplex) DELETE r"
    #         )

    #     print("Dropping IS_SUB_COMPLEX_OF relationships, if any - DONE")

    #     query = """
    #     MATCH
    #         (src_complex:PDBComplex)<-[rel1:IS_PART_OF_PDB_COMPLEX]-()-
    #         [rel2:IS_PART_OF_PDB_COMPLEX]->(dest_complex:PDBComplex)
    #         WHERE rel1.STOICHIOMETRY=rel2.STOICHIOMETRY
    #     WITH DISTINCT src_complex, dest_complex, rel1
    #     WITH src_complex, startNode(rel1) AS relRelations, dest_complex
    #     WITH src_complex, COUNT(relRelations) AS relRelationsAmount, dest_complex
    #     MATCH (src_complex)<-[allRelations:IS_PART_OF_PDB_COMPLEX]-()
    #     WITH src_complex, relRelationsAmount, count(allRelations) AS allRelationsAmount,
    #     dest_complex
    #         WHERE relRelationsAmount = allRelationsAmount
    #     CREATE (dest_complex)<-[:IS_SUB_COMPLEX_OF]-(src_complex)
    #     """  # noqa: B950

    #     print(
    #         f"Creating IS_SUB_COMPLEX_OF relationships - Started at {datetime.now()}"
    #     )
    #     with self._driver.session() as session:
    #         session.run(query)
    #     print(
    #         f"Creating IS_SUB_COMPLEX_OF relationships - Ended at {datetime.now()}"
    #     )

    #     print(
    #         f"Checking for sub complexes and making relationships - Ended at {datetime.now()}"
    #     )


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "-b",
        "--bolt-url",
        required=True,
        help="BOLT url",
    )

    parser.add_argument(
        "-u",
        "--username",
        required=True,
        help="DB username",
    )

    parser.add_argument(
        "-p",
        "--password",
        required=True,
        help="DB password",
    )

    # parser.add_argument(
    #     "-d",
    #     "--connection",
    #     required=False,
    #     help="DB connection string",
    # )

    parser.add_argument(
        "-o",
        "--csv-path",
        required=True,
        help="Path to CSV file containing complexes information",
    )

    args = parser.parse_args()

    complex = Neo4JProcessComplex(
        bolt_uri=args.bolt_url,
        username=args.username,
        password=args.password,
        csv_path=args.csv_path,
    )

    complex.get_complex_portal_data()
    complex.drop_PDBComplex_nodes()
    complex.process_assembly_data()
    complex.create_graph_relationships()
    complex.drop_existing_subcomplex_relationships()
    complex.create_subcomplex_relationships()
    complex.export_csv()


if __name__ == "__main__":
    start_time = time.time()
    main()
    print("--- %s seconds ---" % (time.time() - start_time))
