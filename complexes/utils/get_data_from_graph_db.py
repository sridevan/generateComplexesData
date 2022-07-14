#!/usr/bin/env python3

import logging

from py2neo import Graph

# from pdbecommon.graph import neo4j

# from pdbe.complex_analysis.settings import Settings

# logger = logging.getLogger()
# FORMAT = "%(funcName)s [%(levelname)s] - %(message)s"
# logging.basicConfig(format=FORMAT)


class GetComplexData:
    def __init__(self, bolt_uri, username, password):
        self.graph = None
        self.molecule_names = {}
        self.pdb_complexes = {}
        # self._driver = neo4j.GraphDatabase.driver(bolt_uri, auth=(username, password))
        # self.cfg = Settings().dict()
        self.bolt_host = bolt_uri
        self.neo4j_username = username
        self.neo4j_password = password

        self.pdb_complex_query = """
        MATCH (pdb_complex:PDBComplex)<-[rel:IS_PART_OF_PDB_COMPLEX]-(component)
        OPTIONAL MATCH (component)-[:HAS_TAXONOMY]->(tax:Taxonomy)
        OPTIONAL MATCH (component)-[:HAS_ANTIBODY_MAPPING]->(ab:Antibody)
        WITH pdb_complex.COMPLEX_ID AS complex_id,
        labels(component) AS component_db,
        rel.STOICHIOMETRY AS stoichiometry,
        component.TYPE AS component_type,
        component.ACCESSION AS accession,
        component.RFAM_ACC AS rfam_accession,
        component.POLYMER_TYPE AS polymer_type,
        component.NUMBER_POLY_SEQ AS entity_length,
        tax.TAX_ID AS taxonomy,
        CASE labels(component)
            WHEN ['Assembly']
            THEN component.UNIQID
        END AS entry_assembly,
        CASE ab.ID
            WHEN '1'
                THEN True
            ELSE
                False
            END AS antibody,
        component.UNIQID AS entity
        RETURN complex_id, component_db, component_type, stoichiometry, accession,
        rfam_accession, polymer_type, taxonomy, entry_assembly, entity, antibody,
        entity_length
        """

        self.entity_query = """
        MATCH (e:Entity)
        WITH e.UNIQID AS entity_uniqid,
            e.DESCRIPTION AS description
        RETURN entity_uniqid, description
        """

        self.uniprot_query = """
        MATCH (u:UniProt)
        WITH u.ACCESSION AS accession,
            u.DESCR AS description
        RETURN accession, description
        """

        self.rfam_query = """
        MATCH (u:RfamFamily)
        WITH u.RFAM_ACC AS accession,
            u.DESCRIPTION AS description
        RETURN accession, description
        """

    def populate_molecule_names_from_entity(self):
        """
        Stores the molecules description obtained from
        entity node into a dictionary
        """
        query = self.entity_query
        mappings = self.run_query(query=query)

        for row in mappings:
            entity_uniqid = row.get("entity_uniqid")
            description = row.get("description")
            self.molecule_names[entity_uniqid] = description

    def populate_molecule_names_from_uniprot_or_rfam(self, db_type):
        """
        Stores the molecules description obtained from
        uniprot/rfam node into a dictionary


        Args:
            db_type (str): either uniprot or rfam
        """
        if db_type == "uniprot":
            query = self.uniprot_query
        else:
            query = self.rfam_query

        mappings = self.run_query(query=query)

        for row in mappings:
            accession = row.get("accession")
            name = row.get("description")
            self.molecule_names[accession] = name

    # def populate_molecule_names(self):
    #     """
    #     Get molecules description from Uniprot, RFAM
    #     and entity nodes from the graph database
    #     """
    #     logging.info("Populate molecule names - START")

    #     query = """
    #     MATCH (e:Entity)
    #     WITH e.UNIQID AS entity_uniqid,
    #         e.DESCRIPTION AS description
    #     RETURN entity_uniqid, description
    #     """

    #     mappings = self.run_query(query=query)

    #     for row in mappings:
    #         entity_uniqid = row.get('entity_uniqid')
    #         description = row.get('description')
    #         self.molecule_names[entity_uniqid] = description

    #     query = """
    #     MATCH (u:UniProt)
    #     WITH u.ACCESSION AS accession,
    #         u.DESCR AS description
    #     RETURN accession, description
    #     """

    #     mappings = self.run_query(query=query)

    #     for row in mappings:
    #         accession = row.get('accession')
    #         name = row.get('description')
    #         self.molecule_names[accession] = name

    #     query = """
    #             MATCH (u:RfamFamily)
    #             WITH u.RFAM_ACC AS accession,
    #                 u.DESCRIPTION AS description
    #             RETURN accession, description
    #             """

    #     mappings = self.run_query(query=query)

    #     for row in mappings:
    #         accession = row.get('accession')
    #         name = row.get('description')
    #         self.molecule_names[accession] = name

    #     logging.info('{} names obtained'.format(len(self.molecule_names)))
    #     logging.info("Populate molecule names - END")

    def get_graph(self):
        """
        Create an instance of Neo4j with the connection
        parameters
        """
        self.graph = Graph(
            self.bolt_host, user=self.neo4j_username, password=self.neo4j_password
        )

    def run_query(self, query):
        """
        Runs Neo4J database query and returns its result

        Args:
            query (str): Neo4j query

        Returns:
            obj: result of Neo4j query
        """
        if not self.graph:
            self.get_graph()
        return self.graph.run(query)

    def get_pdb_complex_data(self):
        """
        Returns information related to complexes that are stored
        in the graph database

        Returns:
            dict: The key of the dict is the pdb_complex_id and the
            value is a nested dictionary containing information
            related to complexes
        """
        self.populate_molecule_names_from_uniprot_or_rfam("uniprot")
        self.populate_molecule_names_from_uniprot_or_rfam("rfam")
        self.populate_molecule_names_from_entity()

        logging.info("Get PDB Complex Data - START")
        mappings = self.run_query(self.pdb_complex_query)
        for row in mappings:

            pdb_complex_id = row.get("complex_id")
            component_db = row.get("component_db", [])
            component_type = row.get("component_type")
            stoichiometry = row.get("stoichiometry", 0)
            accession = row.get("accession", "")
            rfam_accession = row.get("rfam_accession", "")
            polymer_type = row.get("polymer_type")
            tax_id = row.get("taxonomy", "")
            entry_assembly = row.get("entry_assembly", "")
            entity = row.get("entity")
            antibody = row.get("antibody", False)
            entity_length = row.get("entity_length", 0)

            self.pdb_complexes.setdefault(pdb_complex_id, {})[
                "pdb_complex_id"
            ] = pdb_complex_id

            for db in component_db:
                stoichiometry = stoichiometry if stoichiometry else 0
                tax_id = tax_id if tax_id else ""
                component = {"stoichiometry": stoichiometry, "tax_id": tax_id}

                if db == "Assembly":
                    if entry_assembly:
                        entry_info = entry_assembly.split("_")
                        pdbid = entry_info[0]
                        # assembly_id = entry_info[1]

                        self.pdb_complexes.setdefault(pdb_complex_id, {}).setdefault(
                            "pdb_entries", []
                        ).append(pdbid)
                        self.pdb_complexes.setdefault(pdb_complex_id, {}).setdefault(
                            "pdb_entries_with_assemblies", []
                        ).append(entry_assembly)
                    continue

                elif db == "UniProt":
                    accession_without_isoform = accession.split("-")[0]
                    component.update(
                        {
                            "accession": accession_without_isoform,
                            "accession_with_isoform": accession,
                            "database": "UNP",
                            "polymer_type": "PROTEIN",
                            "molecule_name": self.molecule_names.get(accession),
                        }
                    )

                elif db == "RfamFamily":
                    component.update(
                        {
                            "accession": rfam_accession,
                            "database": "Rfam",
                            "polymer_type": "RNA",
                            "molecule_name": self.molecule_names.get(rfam_accession),
                        }
                    )

                elif db == "UnmappedPolymer":
                    component.update(
                        {
                            "polymer_type": component_type,
                            "molecule_name": component_type,
                        }
                    )

                elif db == "Entity":
                    polymer_type_dict = {"P": "PROTEIN"}
                    antibody_map = {"True": True, "False": False}
                    entity_length = int(entity_length) if entity_length else 0
                    if entity_length > 19:
                        molecule_name = self.molecule_names.get(entity)
                    elif entity_length > 9:
                        if tax_id in (None, ""):
                            molecule_name = "peptide (unknown source)"
                        elif tax_id == "32630":
                            molecule_name = "synthetic peptide"
                        else:
                            molecule_name = "peptide"
                    else:
                        if tax_id in (None, ""):
                            molecule_name = "short peptide (unknown source)"
                        elif tax_id == "32630":
                            molecule_name = "synthetic short peptide"
                        else:
                            molecule_name = "short peptide"
                    # molecule_name = self.molecule_names.get(entity) if
                    # entity_length > 19 else 'peptide'
                    component.update(
                        {
                            "polymer_type": polymer_type_dict.get(
                                polymer_type, polymer_type
                            ),
                            "molecule_name": molecule_name,
                            "is_antibody": antibody_map.get(antibody, antibody),
                            "pdb_entity": entity,
                            "entity_length": entity_length,
                        }
                    )

                else:
                    logging.info("unhandled db type")
                    logging.info(db)
                    continue

                self.pdb_complexes.setdefault(pdb_complex_id, {}).setdefault(
                    "components", []
                ).append(component)

        logging.info("Get PDB Complex Data - END")
        return self.pdb_complexes

    def return_pdb_data(self):
        # Return the complexes data dictionary
        return self.pdb_complexes
