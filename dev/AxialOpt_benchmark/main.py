from DARTassembler.src.assembler.isomer import AxialOptModifier, DuplicateIsomerFilter
from DARTassembler.src.assembler.isomer import get_atomic_props_from_ase_atoms
from scipy.optimize import differential_evolution
import logging
import pickle
import numpy as np
import time
from tqdm import tqdm
from ase import Atom
import copy


class AxialOptModifierCustom(AxialOptModifier):
    """
    Custom AxialOptModifier that allows additional parameters for differential_evolution.
    """

    def modify(self, target_vectors_list, ligand_origins_list, maxiter=1000, popsize=15, workers=None):
        """
        Optimize each isomer independently, with its own target_vectors and ligand_origins.
        """
        if not self.opt_command:
            return self.input_isomers

        # Clear output isomers each run
        self.output_isomers = []

        # Sanity check lengths
        if len(self.input_isomers) != len(target_vectors_list) or len(self.input_isomers) != len(ligand_origins_list):
            raise ValueError("Each isomer must have its own set of target_vectors and ligand_origins.")

        # Optimize each isomer independently
        for isomer, target_vectors, ligand_origins in zip(self.input_isomers, target_vectors_list, ligand_origins_list):
            atoms = isomer.atoms.copy()

            # Each ligand rotation angle gets its own bound
            bounds = [[0, 360] for _ in target_vectors]

            geometries = [ligand.geometry for ligand in isomer.ligands]

            # Run the optimizer
            result = differential_evolution(
                self.objective_function,
                bounds=bounds,
                args=(target_vectors, ligand_origins, atoms.copy(), isomer.ligand_idc, geometries),
                seed=42,
                maxiter=maxiter,
                popsize=popsize,
                workers=workers,
                polish=True,
                vectorized=True,
                updating='deferred'
            )

            best_ligand_angles = list(result.x)

            # Correctly apply rotations to this isomer's ligands
            for angle, axis, origin, idc, ligand in zip(best_ligand_angles, target_vectors, ligand_origins, isomer.ligand_idc, isomer.ligands):
                if ligand.geometry not in ['1_monodentate', '2_trans']:
                    continue
                self.rotate(atoms=atoms, vector=np.asarray(axis).squeeze(), origin=origin, idc=idc, angle=angle)

            # Copy isomer before modification to avoid unintended side effects
            new_isomer = deepcopy(isomer)
            new_isomer.atoms = atoms
            new_isomer.atomic_props = get_atomic_props_from_ase_atoms(atoms)

            for ligand, idc in zip(new_isomer.ligands, new_isomer.ligand_idc):
                ligand.atoms = deepcopy(atoms[idc])
                ligand.atomic_props = get_atomic_props_from_ase_atoms(ligand.atoms)

            self.output_isomers.append(new_isomer)

        logging.debug(f"Optimized {len(self.output_isomers)} complexes correctly.")
        return self.output_isomers


def plot_color_map(time_dict):
    """
    Plot a 3D scatter plot of timing and similarity results using Plotly.

    :param time_dict: Dictionary containing maxiter, popsize, time, and diff_sum data.
    """
    import pandas as pd
    import plotly.express as px

    # Convert the timing dictionary to a DataFrame
    df = pd.DataFrame(time_dict)

    # Create 3D scatter plot with color = diff_sum
    fig = px.scatter_3d(
        df,
        x='popsize',             # X-axis: Population Size
        y='num_workers',             # Y-axis: Max Iterations
        z='time',                # Z-axis: Timing Results
        color='max_diff',        # Color by similarity (diff_sum)
        color_continuous_scale='RdBu_r',  # Red-Blue reversed: blue=low, red=high
        size_max=10,
        title='Timing & Similarity Results for Axial Optimization (3D Scatter)',
    )

    # Update layout for better visuals
    fig.update_layout(
        scene=dict(
            xaxis_title='Population Size',
            yaxis_title='Number of worker threads',
            zaxis_title='Time (s)'
        ),
        coloraxis_colorbar=dict(
            title='Max of Differences (max_diff)',
            tickformat=".2f"
        ),
        margin=dict(l=0, r=0, b=0, t=40),
    )

    fig.show()


if __name__ == "__main__":
    # inputs
    ligand_origins = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    metal_centers = [[Atom('Ir', [0.0, 0.0, 0.0])],
                     [Atom('Ir', [0.0, 0.0, 0.0])],
                     [Atom('Ir', [0.0, 0.0, 0.0])],
                     [Atom('Ir', [0.0, 0.0, 0.0])]]

    print("Loading assembled isomers from file...")
    with open("assembled_isomers_list.pkl", "rb") as file:
        assembled_isomers_list = pickle.load(file)

    with open("target_vectors_lists.pkl", "rb") as file:
        target_vectors = pickle.load(file)

    # scan range
    maxiter_list = np.arange(10, 50, 10)
    popsize_list = np.arange(1, 35, 1)
    workers_list = [1, 2, 4, 8, -1]



    timing_dict = {
        'maxiter': [],
        'popsize': [],
        'time': [],
        'max_diff': [],
        'num_workers': []
    }

    from copy import deepcopy
    modifier = AxialOptModifierCustom(isomers=deepcopy(assembled_isomers_list))
    rotated_isomers_bench = modifier.modify(target_vectors_list=target_vectors,
                                            ligand_origins_list=[ligand_origins] * len(assembled_isomers_list),
                                            maxiter=1000, popsize=15,
                                            workers=1  # Use single worker for benchmarking
                                            )

    dupe_bench = DuplicateIsomerFilter(isomers=rotated_isomers_bench, metal_centers=metal_centers)
    dupe_bench.get_duplicate_groups()
    df_bench_matrix = dupe_bench.diff_matrix

    for workers in tqdm(workers_list, desc="num_workers Scan"):
        for popsize_val in tqdm(popsize_list, desc=f"Pop_size Scan"):
            logging.info(f"Running optimization with workers={workers}, popsize={popsize_val}")

            start = time.perf_counter()
            rotated_isomers = modifier.modify(target_vectors_list=target_vectors,
                                              ligand_origins_list=[ligand_origins] * len(assembled_isomers_list),
                                              maxiter=1000,
                                              popsize=popsize_val,
                                              workers=workers)

            dupe = DuplicateIsomerFilter(isomers=rotated_isomers, metal_centers=metal_centers)
            dupe.get_duplicate_groups()
            dupe_matrix = dupe.diff_matrix

            end = time.perf_counter()

            elapsed_time = end - start

            abs_diff = np.abs(df_bench_matrix - dupe_matrix)
            max_diff = np.max(abs_diff)

            timing_dict['maxiter'].append(1000)
            timing_dict['popsize'].append(popsize_val)
            timing_dict['time'].append(elapsed_time)
            timing_dict['num_workers'].append(workers)
            timing_dict['max_diff'].append(max_diff)

            logging.info(f"Number of isomers optimized: {len(rotated_isomers)}")

    plot_color_map(timing_dict)
