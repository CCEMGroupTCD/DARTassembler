import functools
import json
import sys
from copy import deepcopy
from tqdm import tqdm
import pandas as pd
import numpy as np
import warnings
from scipy.special import comb
from typing import Union
from pathlib import Path
import jsonlines
import itertools
import ase
from collections import defaultdict
from DARTassembler.src.constants.chem import Element
from DARTassembler.src.metalig.utils_molecule import get_all_ligands_by_graph_hashes
from DARTassembler.src.misc.io import save_json, NumpyEncoder, load_unique_ligand_db
from DARTassembler.src.misc.io import save_to_xyz


def update_old_ligand_db_to_new(path: str, csv=False) -> None:
    """
    Updates a ligand database in the old format before DARTassembler version 1.1.0 to the new format from version 1.1.0 on.
    :param path: Path to the ligand database
    :param csv: If True, saves new ligand database also as csv.
    :return: None
    """
    path = Path(path)
    # Ignore warnings about deprecated formats in the ligand database because now we are updating the database to the new format.
    warnings.filterwarnings("ignore", message='The provided ligand database has a format that was deprecated.*')
    db = LigandDB.from_json(path)
    db.save_to_file(path)
    if csv:
        db.save_to_csv(outpath=path.with_suffix('.csv'))

    return

class BaseDB(object):

    def __init__(self, db: dict):
        self.db = db

    def to_dict(self, desc: str= 'Convert DB to dict') -> dict:
        """
        Converts all molecules in the database to a dictionary format.
        :param desc: str: Description for the progress bar.
        :return: A dictionary where keys are molecule identifiers and values are dictionaries representing the molecules.
        """
        json_dict = {}
        for key, mol in tqdm(self.db.items(), desc):
            json_dict[key] = mol.to_dict()
        return json_dict

    def _to_json(self, path, desc: str= 'Save DB to json', json_lines: bool=False) -> None:
        """
        Saves the database to a .json or .jsonlines file.
        :param path: Path to the output file where the database will be saved.
        :param desc: str: Description for the progress bar.
        :param json_lines: bool: If True, saves the database as jsonlines. If False, saves the database as a single json file.
        :return: None
        """
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)

        if json_lines:
            with jsonlines.open(path, mode='w', dumps=functools.partial(json.dumps, cls=NumpyEncoder)) as writer:
                for key, mol in tqdm(self.db.items(), desc, file=sys.stdout):
                    data = {'key': key, 'value': mol.to_dict()}
                    writer.write(data)
        else:
            d = self.to_dict(desc=desc)
            save_json(d, path=path, indent=4)

        return

    def save_to_file(self, path: Union[str, Path], desc: Union[str, None]=None, json_lines=True) -> None:
        """
        Saves the database to a .json or .jsonlines file.
        :param path: Path to the output file where the database will be saved.
        :param desc: Description for the progress bar. If None, no progress bar will be shown.
        :param json_lines: If True, saves the database as jsonlines. If False, saves the database as a single json file.
        :return: None
        """
        self._to_json(path=path, desc=desc, json_lines=json_lines)
        return


class LigandDB(BaseDB):

    def __init__(self, db: dict):
        super().__init__(db=db)

    @classmethod
    def from_json(cls,
                  path: Union[str, Path] = 'metalig',
                  n_max: Union[int,None] = None,
                  show_progress: bool = True
                  ) -> 'LigandDB':
        """
        Loads a ligand database from a .jsonlines file.
        :param path: Path to the .jsonlines file of the ligand database. Alternatively, the string 'metalig' can be used to load the default ligand database.
        :param n_max: Maximum number of ligands to load. If None, all ligands will be loaded.
        :param show_progress: If True, a progress bar will be shown.
        :return: A LigandDB object
        """
        return cls(load_unique_ligand_db(path=path, n_max=n_max, show_progress=show_progress, molecule='class'))

    def get_df(self, max_entries: int=5) -> pd.DataFrame:
        """
        Returns a DataFrame with important ligand information for all ligands in the database, such as charge, stoichiometry, archetype, and more.
        :param max_entries: Maximum number of entries for long lists in the DataFrame, such as the list of CSD codes in which this ligand is present.
        :return: A DataFrame with ligand information, with the index being `unique_name`.
        """
        ligands = {uname: ligand.get_csv_info(max_entries=max_entries) for uname, ligand in self.db.items()}
        df_ligand_info = pd.DataFrame.from_dict(ligands, orient='index')
        assert df_ligand_info.index.tolist() == df_ligand_info['unique_name'].tolist()
        df_ligand_info = df_ligand_info.reset_index(drop=True).set_index('unique_name')

        return df_ligand_info

    def save_to_csv(self, outpath: Union[str, Path], max_entries: int=5) -> pd.DataFrame:
        """
        Saves a csv file with important ligand information for all ligands in the database, such as charge, stoichiometry, archetype, and more.
        :param outpath: Path to the output csv file.
        :param max_entries: Maximum number of entries for long lists in the DataFrame, such as the list of CSD codes in which this ligand is present.
        :return: A DataFrame with ligand information.
        """
        df = self.get_df(max_entries=max_entries)
        df.to_csv(outpath, index=True)  # the index is the unique_names of the ligands

        return df

    def save_to_concat_xyz(self, outpath: Union[str, Path], with_metal: bool=True, comment: str=None) -> str:
        """
        Save a concatenated xyz file with all ligands in the database.
        :param outpath: Path to the output xyz file.
        :param with_metal: If True, the output structure of each ligand will contain a metal center at the original position for better visualization. If False, only the ligand structure itself will be saved.
        :param comment: A comment to be added to each xyz structure. If None, a default comment will be used.
        :return: A concatenated xyz string with all ligands in the database.
        """
        outpath = Path(outpath)
        outpath.parent.mkdir(parents=True, exist_ok=True)
        xyz_string = self.get_concat_xyz_string(with_metal=with_metal, comment=comment)
        with open(outpath, 'w') as f:
            f.write(xyz_string)

        return xyz_string

    def get_concat_xyz_string(self, with_metal: bool=True, comment: str=None) -> str:
        """
        Get a concatenated xyz string with all ligands in the database.
        :param with_metal: If True, the output structure of each ligand will contain a metal center at the original position for better visualization. If False, only the ligand structure itself will be saved.
        :param comment: A comment to be added to each xyz structure. If None, a default comment will be used.
        :return: A concatenated xyz string with all ligands in the database.
        """
        return '\n'.join([lig.get_xyz_string(comment=comment, with_metal=with_metal) for lig in self.db.values()])

    def _get_ligand_archetypes(self, sort_by_rssd: bool=False) -> dict:
        """
        Assigns archetypes to all ligands in the database and returns a dictionary with archetypes as keys and a list of tuples as values.
        :param sort_by_rssd:  If True, sorts the ligands by the weight necessary for change (rssd) within each archetype. This is useful to have the ligands with the lowest rssd first.
        :return: A dictionary with archetypes as keys and a list of tuples as values. Each tuple contains the ligand name, rssd, weight necessary for change, second archetype, isomers and isomer_idc.
        """
        data = defaultdict(list)
        for name, ligand in tqdm(self.db.items(), desc='Assigning ligand archetypes'):
            archetype, isomers, isomer_idc, rssd, second_archetype, weight_necessary_for_change = ligand.get_ligand_archetype_and_isomers()
            data[archetype].append((name, rssd, weight_necessary_for_change, second_archetype, isomers, isomer_idc))

        # Sort data by archetype name to have archetypes of the same n_eff_donors together
        data = dict(sorted(data.items(), key=lambda x: x[0]))
        if sort_by_rssd:
            # Sort each archetype by weight necessary for change
            for archetype, names_rssd in data.items():
                names_rssd.sort(key=lambda x: x[2])

        return data

    def _save_ligand_archetype_concat_xyz_files(self,
                                               outdir: Union[str, Path],
                                               output_all_isomers: bool=True,
                                               sort_by_rssd: bool=False
                                               ) -> dict:
        """
        Assigns archetypes to all ligands in the database and saves the structures with all isomers for each archetype in a different concatenated xyz file.
        :param outdir: Path to the output directory where the xyz files will be saved.
        :param output_all_isomers: If True, all isomers will be saved in the xyz file. If False, only the first isomer will be saved.
        :param sort_by_rssd: If True, sorts the ligands by the weight necessary for change (rssd) within each archetype. This is useful to have the ligands with the lowest rssd first.
        :return: A dictionary with archetypes as keys and a list of tuples as values. Each tuple contains the ligand name, rssd, weight necessary for change, second archetype, isomers and isomer_idc.
        """
        data = self._get_ligand_archetypes(sort_by_rssd=sort_by_rssd)

        # Save structures with all isomers for each archetype in a different concatenated xyz file
        for archetype, info in data.items():
            atoms, comments, weights, second_archetypes = [], [], [], []
            for name, rssd, weight, second_archetype, isomers, isomer_idc in info:
                for isomer_idx, (isomer, idc) in enumerate(zip(isomers, isomer_idc)):
                    if not output_all_isomers and isomer_idx > 0:
                        continue
                    assert not 'Cu' in isomer.get_chemical_symbols(), 'There should be no Cu atoms in the isomers!'
                    isomer.append(ase.Atom('Cu', position=(0, 0, 0)))
                    atoms.append(isomer)
                    weights.append(weight)
                    comments.append(f'{name}-{isomer_idx} rssd={rssd:.3f} change:{weight:.3f}->{second_archetype} idc={idc}')
                    second_archetypes.append(second_archetype)
            outpath = Path(outdir, f'concat_{archetype}.xyz')
            n_isomers = np.unique([len(isomers) for _, _, _, _, isomers, _ in info])
            n_isomers = n_isomers[0] if len(n_isomers) == 1 else n_isomers
            print(f'{archetype}: {len(atoms)} structures, {n_isomers} isomer{"s" if n_isomers > 1 else ""}')
            save_to_xyz(outpath=outpath, structures=atoms, comments=comments)

        return data

    def _get_lig_db_in_old_format(self):
        """A little helper function to get the ligand database in the old format, where ligands are grouped by their number of donors."""
        lig_dict_old_format = {}
        for lig in self.db.values():
            if lig.n_donors not in lig_dict_old_format:
                lig_dict_old_format[lig.n_donors] = [lig]
            else:
                lig_dict_old_format[lig.n_donors].append(lig)

        return lig_dict_old_format

    @staticmethod
    def _check_property_and_print_if_not_same_for_all_same_ligands(check_props, unique_ligand, ligand):
        for prop in check_props:
            if getattr(unique_ligand, prop) != getattr(ligand, prop):
                print(
                    f'WARNING: Different {prop} for unique ligand {unique_ligand.name} ({getattr(unique_ligand, prop)}) and ligand {ligand.name} ({getattr(ligand, prop)}).')

        return

    def _get_unique_ligands_and_set_unique_ligand_name(self, grouped_unique_ligands):
        unique_ligands = []
        for same_ligands in grouped_unique_ligands:

            unique_ligand = same_ligands[0]
            unique_ligand_name = 'unq_' + unique_ligand.name

            for ligand in tqdm(same_ligands, desc="Filter Duplicates"):
                ligand.unique_name = unique_ligand_name
                ligand.n_total_unique_ligands = len(same_ligands)

                check_props = ['n_donors', 'graph_hash', 'hash', 'unique_name']
                self._check_property_and_print_if_not_same_for_all_same_ligands(check_props, unique_ligand, ligand)

            unique_ligands.append(deepcopy(unique_ligand))

        return unique_ligands

    def _filter_duplicates(self) -> dict:
        """
        Filters out duplicate ligands in the database based on their graph hashes.
        :return: A dictionary with unique ligands, where the keys are the unique ligand names and the values are the unique ligand objects.
        """

        print('Start filtering duplicates.')

        ligands_by_hash = get_all_ligands_by_graph_hashes(list(self.db.values()))
        grouped_unique_ligands = [ligand_list for ligand_list in ligands_by_hash.values()]

        unique_ligands = self._get_unique_ligands_and_set_unique_ligand_name(grouped_unique_ligands)
        unique_ligand_dict = {lig.name: lig for lig in unique_ligands}

        print(f'Number of unique ligands: {len(unique_ligands)}.')

        return unique_ligand_dict

    def _calc_number_of_possible_complexes(self, metals: list[str] = None) -> pd.DataFrame:
        """
        Calculates the number of possible complexes for each metal in the list of metals.
        :param metals: 
        :return:
        """
        if metals is None:
            metals = ['Cr', 'Mn', 'Fe', 'Ru', 'Co', 'Ni']

        df = []
        for metal in metals:
            df_n_metal_combs = self._calc_number_of_possible_complexes_for_metal(metal)
            df.append(df_n_metal_combs)
        df = pd.concat(df, ignore_index=True)

        return df

    def _calc_number_of_possible_complexes_for_metal(self, metal: str, archetypes: dict = None) -> pd.DataFrame:
        metal_oxi_states = Element(metal).common_oxidation_states
        results = []

        # possible archetypes for octahedral and square-planar complexes. This list needs to be expanded when adding new archetypes.
        if archetypes is None:
            archetypes = {
                'octahedral': [(3, 2, 1), (4, 1, 1), (5, 1)],
                'square_planar': [(2, 2), (2, 1, 1)]
            }

        for oxi_state in metal_oxi_states:
            target_charge = -oxi_state
            for archetype_name, archetype_list in archetypes.items():
                for archetype in archetype_list:
                    count = self._calc_number_of_combinations_of_ligands_for_topology(target_charge=target_charge, archetype=archetype)
                    results.append(
                        {'metal': metal, 'oxi_state': oxi_state, 'archetype': archetype_name, 'denticities': archetype,
                         'count': count})

        return pd.DataFrame(results)

    def _calc_number_of_combinations_of_ligands_for_topology(self, target_charge: int, archetype: tuple) -> int:
        """
        Calculates the number of possible ligand combinations for a given target charge and archetype.
        @target_charge: The targeted sum of charges of the ligands.
        @archetype: The topology of the complex, e.g. (3, 2, 1) for a octahedral complex with 3 bidentate, 2 monodentate and 1 tridentate ligand.
        """
        n_ligands = len(archetype)
        archetype = sorted(archetype)

        df = self.get_reduced_df().query('n_donors in @archetype and not charge.isnull()')[['n_donors', 'charge']].astype(int)

        df = df.groupby(['n_donors', 'charge']).size().reset_index().rename(columns={0: 'count'})
        count = 0
        for ligs in itertools.combinations_with_replacement(list(df.itertuples()), n_ligands):
            correct_denticities = sorted(lig.n_donors for lig in ligs) == archetype
            correct_charges = sum(lig.charge for lig in ligs) == target_charge
            if correct_charges and correct_denticities:
                # Group ligands which have the same charge and n_donors.
                groups = pd.DataFrame(ligs).groupby(['n_donors', 'charge'])['count']
                # Calculate the number of possible combinations for this combination of ligands. If there are multiple ligands with the same charge and n_donors, we need to pay attention that we don't count the same combination twice, e.g (lig1, lig2) and (lig2, lig1). That is because we define a complex here just in terms of its set of ligands, without caring about the order of the ligands.
                comb_count = 1
                for _, group in groups:
                    n_same_ligands = len(group)
                    lig_count = group.values[0]
                    if n_same_ligands == 1:
                        comb_count *= lig_count
                    else:
                        # Multiple ligands with same charge and n_donors: Avoid double counting.
                        comb_count *= comb(lig_count + n_same_ligands - 1, n_same_ligands, exact=True)
                count += comb_count

        return count