from src01.Molecule import RCA_Ligand


def get_reactant():
    """
    In this method we return a list which only contains the reactant.
    Need the list format to proceed.
    Will decode the reactant by denticity 0
    In our case it is only OH
    """
    reactants = []
    atomic_props = {"atoms": ["O", "H"], "x": [0, 0.2096], "y": [0, -0.5615], "z": [1.4361, 2.1227]}
    Hydroxi = RCA_Ligand(atomic_props=atomic_props,
                         denticity=0,
                         ligand_to_metal=[0],
                         name="ActiveSite",
                         unique_name="ActiveSite",
                         graph=None
                         )
    reactants.append(Hydroxi)

    return reactants


def get_monodentate_list():
    monodentate_ligands = []

    #
    #
    # 1.
    atomic_props = {"atoms": ["O", "H"], "x": [0, 0.2096], "y": [0, -0.5615], "z": [1.4361, 2.1227]}
    Hydroxi = RCA_Ligand(atomic_props=atomic_props,
                         denticity=1,
                         ligand_to_metal=[0],
                         name="Hydroxi",
                         unique_name="Hydroxi",
                         graph=None
                         )
    monodentate_ligands.append(Hydroxi)

    #
    #
    # 2.
    atomic_props = {"atoms": ["C", "O"], "x": [3.18739, 4.26068], "y": [-0.87756, -0.53207], "z": [0.64611, 0.60959]}
    CO = RCA_Ligand(atomic_props=atomic_props,
                    denticity=1,
                    ligand_to_metal=[0],
                    name="CO",
                    unique_name="CO",
                    graph=None
                    )
    monodentate_ligands.append(CO)

    #
    #
    # 3.
    atomic_props = {"atoms": ["C", "N"], "x": [-0.18655, -0.02880], "y": [0.84136, 0.97282], "z": [1.69712, 2.83591]}
    CN = RCA_Ligand(atomic_props=atomic_props,
                    denticity=1,
                    ligand_to_metal=[1],
                    name="CN",
                    unique_name="CN",
                    graph=None
                    )
    monodentate_ligands.append(CN)

    #
    #
    # 4.
    atomic_props = {"atoms": ["N", "H", "H", "H"],
                    "x": [1.64690, 0.69907, 1.61955, 1.83798],
                    "y": [0.44248, 0.35414, 1.23713, -0.43173],
                    "z": [2.25092, 1.82123, 2.92813, 2.78948]
                    }

    Ammonia = RCA_Ligand(atomic_props=atomic_props,
                         denticity=1,
                         ligand_to_metal=[0],
                         name="NH3",
                         unique_name="NH3"
                         )
    monodentate_ligands.append(Ammonia)

    return monodentate_ligands

