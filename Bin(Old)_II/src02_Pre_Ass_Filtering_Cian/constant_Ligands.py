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
                         unique_name="NH3",
                         graph = None
                         )
    monodentate_ligands.append(Ammonia)
    #
    #
    # 5.
    atomic_props = {"atoms": ["Cl"],
                    "x": [1.0],
                    "y": [1.0],
                    "z": [1.0],
                    }
    Chlorine = RCA_Ligand(atomic_props=atomic_props,
                          denticity=1,
                          ligand_to_metal=[0],
                          name="Cl",
                          unique_name="Cl",
                          graph=None
                          )
    monodentate_ligands.append(Chlorine)
    #
    #
    # 6.
    atomic_props = {"atoms": ["F"],
                    "x": [0],
                    "y": [0],
                    "z": [0],
                    }

    Fluorine = RCA_Ligand(atomic_props=atomic_props,
                          denticity=1,
                          ligand_to_metal=[0],
                          name="F",
                          unique_name="F",
                          graph=None
                          )
    monodentate_ligands.append(Fluorine)
    #
    #
    # 7.
    atomic_props = {"atoms": ["Br"],
                    "x": [0],
                    "y": [0],
                    "z": [0],
                    }

    Bromine = RCA_Ligand(atomic_props=atomic_props,
                         denticity=1,
                         ligand_to_metal=[0],
                         name="Br",
                         unique_name="Br",
                         graph=None
                         )
    monodentate_ligands.append(Bromine)
    #
    #
    # 8.
    atomic_props = {"atoms": ["I"],
                    "x": [1.0, ],
                    "y": [1.0, ],
                    "z": [1.0, ],
                    }

    Iodide = RCA_Ligand(atomic_props=atomic_props,
                        denticity=1,
                        ligand_to_metal=[0],
                        name="I",
                        unique_name="I",
                        graph=None
                        )
    monodentate_ligands.append(Iodide)
    #
    #
    # 9.
    atomic_props = {"atoms": ["N", "C", "C", "C", "H", "H", "H", "H", "H", "H", "H", "H", "H"],
                    "x": [-7.01591, -5.83377, -6.93084, -8.23178, - 5.73266, -4.92471, - 5.87628, - 6.04180, - 6.88690, - 7.79583, - 8.31804, - 9.12049, - 8.25563],
                    "y": [-0.65827, -1.38182, 0.75351, -1.25868, - 1.35195, - 0.95903, - 2.43084, 1.21526, 0.89466, 1.30534, - 2.30545, - 0.74356, - 1.22239],
                    "z": [-0.50304, -0.03640, -0.13104, 0.04461, 1.05466, - 0.47857, - 0.35010, - 0.57495, 0.95509, - 0.51569, - 0.26761, - 0.33684, 1.13989]
                    }
    Tetramethylamine = RCA_Ligand(atomic_props=atomic_props,
                                  denticity=1,
                                  ligand_to_metal=[0],
                                  name="NH3",
                                  unique_name="NH3",
                                  graph=None
                                  )
    monodentate_ligands.append(Tetramethylamine)

    return monodentate_ligands
