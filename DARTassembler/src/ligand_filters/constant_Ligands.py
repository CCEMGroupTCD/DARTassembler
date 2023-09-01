from DARTassembler.src.ligand_extraction.Molecule import RCA_Ligand


def get_reactant():
    """
    In this method we return a list which only contains the reactant.
    Need the list format to proceed.
    Will decode the reactant by denticity 0
    In our case it is only OH
    """
    reactants = []
    atomic_props = {"atoms": ["O", "H"], "x": [0, 0.2096], "y": [0, -0.5615], "z": [1.4361, 2.1227]}
    Hydroxide = RCA_Ligand(atomic_props=atomic_props,
                           denticity=0,
                           ligand_to_metal=[0],
                           name="ActiveSite_OH",
                           unique_name="ActiveSite_OH",
                           graph=None,
                           original_metal_position = [0, 0, 0]
                           )
    Hydroxide.global_props.update({"LCS_pred_charge": -1})
    #reactants.append(Hydroxide)

    return reactants


def get_monodentate_list():
    monodentate_ligands = []

    #
    #
    # 1.
    atomic_props = {"atoms": ["O", "H"], "x": [0, 0.2096], "y": [0, -0.5615], "z": [1.4361, 2.1227]}
    Hydroxide = RCA_Ligand(atomic_props=atomic_props,
                           denticity=1,
                           ligand_to_metal=[0],
                           name="OH",
                           unique_name="OH",
                           graph=None,
                           original_metal_position = [0, 0, 0]

                           )
    Hydroxide.global_props.update({"LCS_pred_charge": -1})
    #monodentate_ligands.append(Hydroxide)

    #
    #
    # 2.
    atomic_props = {"atoms": ["C", "O"], "x": [3.18739, 3.18739], "y": [-0.87756, -0.87756], "z": [0.64611, 1.77411]}
    CO = RCA_Ligand(atomic_props=atomic_props,
                    denticity=1,
                    ligand_to_metal=[0],
                    name="CO",
                    unique_name="CO",
                    graph=None,
                    original_metal_position = [0, 0, 0]
                    )
    CO.global_props.update({"LCS_pred_charge": 0})
    # monodentate_ligands.append(CO)

    #
    #
    # 3.
    atomic_props = {"atoms": ["C", "N"], "x": [-0.18655, -0.02880], "y": [0.84136, 0.97282], "z": [1.69712, 2.83591]}
    CN = RCA_Ligand(atomic_props=atomic_props,
                    denticity=1,
                    ligand_to_metal=[1],
                    name="CN",
                    unique_name="CN",
                    graph=None,
                    original_metal_position = [0, 0, 0]
                    )
    CN.global_props.update({"LCS_pred_charge": -1})
    # monodentate_ligands.append(CN)

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
                         graph=None,
                         original_metal_position = [0, 0, 0]
                         )
    Ammonia.global_props.update({"LCS_pred_charge": 0})
    # monodentate_ligands.append(Ammonia)
    #
    #
    # 5.
    atomic_props = {"atoms": ["Cl"],
                    "x": [1.0],
                    "y": [1.0],
                    "z": [1.0],
                    }
    Chloride = RCA_Ligand(atomic_props=atomic_props,
                          denticity=1,
                          ligand_to_metal=[0],
                          name="Cl",
                          unique_name="Cl",
                          graph=None,
                          original_metal_position = [0, 0, 0]
                          )
    Chloride.global_props.update({"LCS_pred_charge": -1})
    # monodentate_ligands.append(Chloride)
    #
    #
    # 6.
    atomic_props = {"atoms": ["F"],
                    "x": [0],
                    "y": [0],
                    "z": [0],
                    }

    Fluoride = RCA_Ligand(atomic_props=atomic_props,
                          denticity=1,
                          ligand_to_metal=[0],
                          name="F",
                          unique_name="F",
                          graph=None,
                          original_metal_position = [0, 0, 0]
                          )
    Fluoride.global_props.update({"LCS_pred_charge": -1})
    # monodentate_ligands.append(Fluoride)
    #
    #
    # 7.
    atomic_props = {"atoms": ["Br"],
                    "x": [1.0],
                    "y": [1.0],
                    "z": [1.0],
                    }

    Bromide = RCA_Ligand(atomic_props=atomic_props,
                         denticity=1,
                         ligand_to_metal=[0],
                         name="Br",
                         unique_name="Br",
                         graph=None,
                         original_metal_position = [0, 0, 0]
                         )
    Bromide.global_props.update({"LCS_pred_charge": -1})
    # monodentate_ligands.append(Bromide)
    #
    #
    # 8.
    atomic_props = {"atoms": ["I"],
                    "x": [1.0],
                    "y": [1.0],
                    "z": [1.0],
                    }

    Iodide = RCA_Ligand(atomic_props=atomic_props,
                        denticity=1,
                        ligand_to_metal=[0],
                        name="I",
                        unique_name="I",
                        graph=None,
                        original_metal_position = [0, 0, 0]
                        )
    Iodide.global_props.update({"LCS_pred_charge": -1})
    # monodentate_ligands.append(Iodide)
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
                                  graph=None,
                                  original_metal_position = [0, 0, 0]
                                  )
    Tetramethylamine.global_props.update({"LCS_pred_charge": 0})
    # monodentate_ligands.append(Tetramethylamine)

    #
    #
    # 10.
    atomic_props = {"atoms": ["C", "C", "N", "C", "C", "C", "H", "H", "H", "H", "H"],
                    "x": [-2.78072, -2.64593, -1.45555, -1.62507, -0.38426, -0.35119, -3.76261, -3.51801, -1.69090, 0.53485, 0.59281],
                    "y": [-1.73705, -0.35782, 0.28087, -2.51117, -1.88256, -0.49715, -2.19633, 0.28963, -3.59530, -2.45726, 0.04004],
                    "z": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    }
    Pyridine = RCA_Ligand(atomic_props=atomic_props,
                          denticity=1,
                          ligand_to_metal=[2],
                          name="Pyridine",
                          unique_name="Pyridine",
                          graph=None,
                          original_metal_position = [0, 0, 0]
                          )
    Pyridine.global_props.update({"LCS_pred_charge": 0})
    # monodentate_ligands.append(Pyridine)
                                        #N
    atomic_props = {"atoms": ["C", "C", "C", "C", "C", "C", "F", "F", "F", "F", "F"],
                    "x": [-2.78072, -2.64593, -1.45555, -1.62507, -0.38426, -0.35119, -3.76261, -3.51801, -1.69090, 0.53485, 0.59281],
                    "y": [-1.73705, -0.35782, 0.28087, -2.51117, -1.88256, -0.49715, -2.19633, 0.28963, -3.59530, -2.45726, 0.04004],
                    "z": [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
                    }
    pentafluorophenyl = RCA_Ligand(atomic_props=atomic_props,
                                   denticity=1,
                                   ligand_to_metal=[2],
                                   name="pentafluorophenyl",
                                   unique_name="pentafluorophenyl",
                                   graph=None,
                                   original_metal_position=[0, 0, 0]
                                   )
    pentafluorophenyl.global_props.update({"LCS_pred_charge": -1})
    monodentate_ligands.append(pentafluorophenyl)

    return monodentate_ligands
