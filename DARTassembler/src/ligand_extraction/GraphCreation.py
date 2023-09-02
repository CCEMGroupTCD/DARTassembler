"""
A collection of all graph creating methods
"""
import networkx as nx
from ase import Atoms, neighborlist
from ase.data import covalent_radii  # THE basic covalent radii data
import warnings
from DARTassembler.src.constants.Paths import csd_graph_path
from DARTassembler.src.constants.Periodic_Table import DART_Element
from DARTassembler.src.ligand_extraction.io_custom import load_json
from DARTassembler.src.ligand_extraction.utilities_graph import graph_from_graph_dict


class GraphCreation:

    def __init__(self,
                 graph_creating_strategy: str,
                 molecule: Atoms,
                 atomic_props: dict,
                 **kwargs
                 ):
        """
        :param kwargs: potential arguments for other graph creation methods, see the specific method for more information
        """
        self.G = None

        if graph_creating_strategy == "smiles":
            if "csd_code" not in kwargs:
                warnings.warn("Missing Information, no graph could be created")
            else:
                self.smiles_graph(identifier=kwargs["csd_code"])

        elif graph_creating_strategy == "CSD":
            if "csd_code" not in kwargs:
                warnings.warn("Missing Information, no graph could be created")
            else:
                self.CSD_graphs(identifier=kwargs["csd_code"])

        elif graph_creating_strategy in ["default", "ase_cutoff"]:
            self.ase_cutoff_graph(mol=molecule, **kwargs)

        ##########  deprecated graph creation methods
        # elif graph_creating_strategy == "pymatgen_NN":
        #     if atomic_props is None or atomic_props == {}:
        #         warnings.warn("Empty atomic props, no graph could be created")
        #     else:
        #         self.pymat_JmolNN_graph(mol=molecule, atomic_props=atomic_props)
        #
        # elif graph_creating_strategy == "molsimplifyGraphs":
        #     if atomic_props is None or atomic_props == {}:
        #         warnings.warn("Empty atomic props, no graph could be created")
        #     else:
        #         self.molsimiply_Graph(mol=molecule, atomic_props=atomic_props)

        if self.G is None:
            # then all graph strategies have failed so far, we use the default creation
            # without any additional parameters
            # Note that this here defines the default, so if we want to change that we have to change the function
            # here as well
            warnings.warn(f'Graph creating strategy {graph_creating_strategy} not found. Fallback to ASE cutoff graphs.')
            self.ase_cutoff_graph(mol=molecule)

    #
    #
    #
    # 1. ase based methods
    @staticmethod
    def convert_ase_cutoff_corrections(dict_):
        """
        We need to bring the cutoff directions from a readible and tunable format
        as I would like to have it
        into the format ase needs to actually overwrite the build-in cutoffs
        """

        return {
            DART_Element(an).symbol: covalent_radii[an] + correction
            for an, correction in dict_.items()
        }

    def ase_cutoff_graph(self,
                         mol: Atoms,
                         skin_: float = 0.2,
                         cutoff_corrections_for_metals: dict = None,
                         **kwargs
                         ):
        """
        Standard graph creation method using the cutoffs from THE covalent radii paper.
        :param cutoff_corrections_for_metals: potential corrections for metal cutoffs
                                                format: {atom_number: difference to covalent_radii}
                                                example: {22: -0.05} if we want to reduce the titanium cutoff by 0.05 A
        :param kwargs: is just to dump potential not required arguments which may be passed into the function
        """
        if cutoff_corrections_for_metals is None:
            cutoff_corrections_for_metals = {}

        cutOff = neighborlist.natural_cutoffs(mol,
                                              **self.convert_ase_cutoff_corrections(cutoff_corrections_for_metals)
                                              )
        neighborList = neighborlist.NeighborList(cutOff, skin=skin_, self_interaction=False, bothways=True)
        neighborList.update(mol)

        A = neighborList.get_connectivity_matrix(sparse=False)

        labels = {i: el for i, el in enumerate(mol.get_chemical_symbols())}

        self.G = nx.Graph(A)

        for node in self.G.nodes():
            self.G.nodes[node]["node_label"] = labels[node]

    #
    #
    #
    # 2. pymatgen based methods
    # @staticmethod
    # def atomic_props_to_pymatmol(atomic_props):
    #     return PyMatMol(species=atomic_props["atoms"],
    #                     coords=[[atomic_props[key_][i] for key_ in ["x", "y", "z"]] for i, _ in
    #                             enumerate(atomic_props["x"])])
    #
    # def pymat_JmolNN_graph(self, mol: Atoms, atomic_props: dict):
    #
    #     pymat_mol = self.atomic_props_to_pymatmol(atomic_props)
    #
    #     pymat_graph = MoleculeGraph.with_local_env_strategy(molecule=pymat_mol, strategy=JmolNN())
    #
    #     self.G = nx.Graph(pymat_graph.graph)  # is a networkx object
    #
    #     labels = {i: el for i, el in enumerate(mol.get_chemical_symbols())}
    #     for node in self.G.nodes():
    #         self.G.nodes[node]["node_label"] = labels[node]

    #
    #
    #
    # 3. molsimplify based methods
    @staticmethod
    def get_xyz_file_format_string_from_atomic_props(atomic_props):
        """
        returns a string that can be written into an .xyz file
        is based on the method get_xyz_file_format_string
        from an RCA_Molecule
        """
        str_ = f"{len(atomic_props['x'])}\n \n"
        for i, _ in enumerate(atomic_props['x']):
            str_ += f"{atomic_props['atoms'][i]}  {atomic_props['x'][i]}  {atomic_props['y'][i]}  {atomic_props['z'][i]} \n"

        return str_

    # def molsimiply_Graph(self, mol: Atoms, atomic_props: dict):
    #     """
    #     Here we use basically copy the method HK and her group use to create graphs
    #     """
    #     mol_mol = mol3D()
    #     mol_mol.readfromstring(xyzstring=self.get_xyz_file_format_string_from_atomic_props(atomic_props))
    #
    #     # now we create the mol.graph object
    #     # as we dont know if a complex is octahedral in general or not, we default this prop to false
    #     mol_mol.createMolecularGraph(oct=False)
    #
    #     self.G = nx.Graph(mol_mol.graph)  # is a networkx object
    #
    #     labels = {i: el for i, el in enumerate(mol.get_chemical_symbols())}
    #     for node in self.G.nodes():
    #         self.G.nodes[node]["node_label"] = labels[node]

    #
    #
    # 4. SmilesString based methods
    def smiles_graph(self, identifier):
        """
        simple graph creating method using pysmiles to convert the extracted smiles (from the CSD)
        to graphs
        """
        warnings.warn("Method not feasible as the order in the graph creation is crucial and smiles string "
                      "dont preserve that order")
        """
        try:
            smiles = smiles_df.set_index("CSD_code").loc[identifier, "smiles"]
            self.G = nx.Graph(read_smiles(smiles, explicit_hydrogen=True))
            # the nodes need to be renamed, as the default is "element" rather than "node_label"
            for node in self.G.nodes:
                self.G.nodes[node]["node_label"] = self.G.nodes[node]["element"]
        except Exception as e:
            print(f"Smiles based Graph creation not possible, {e}")
            return
        """
        return

    #
    # 5. CSD Graphs
    def CSD_graphs(self, identifier):
        """
        Here we shall use graphs created from the .mol2 files if possible
        if not, we shall use the default graphs, by setting self.G = None,
        because then the main programm will take care of it
        """
        try:
            self.G = graph_from_graph_dict(load_json(f"{csd_graph_path}/{identifier}_g.json"))

        except FileNotFoundError:
            warnings.warn("Graph directory not found, standard graphs are getting created", UserWarning)
            return