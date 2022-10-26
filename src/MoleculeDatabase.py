"""
Takes as input the tmqm downloaded data and makes them easier accessible.
"""
import pandas as pd
from pathlib import Path
import glob
from copy import deepcopy
import io
import numpy as np
from tqdm import tqdm
from Molecule import RCA_Molecule, RCA_Ligand
from Extracted_Molecule import Extracted_Molecule
import json
from src.read_database import read_local_tmqm_db


ATOMIC_PROPERTIES_COLUMN_SEPARATOR = '  ===  '
felix_molecules = ['DUCVIG', 'XUMHUJ', 'VEVTOH', 'AMAWAM', 'BALTUG', 'JIYZOG', 'LIRBUK', 'YAXYUU', 'ZAVZOO', 'COJRIF', 'VIFRIL', 'AFOYUP', 'KAYBOB', 'RATLOQ', 'IZICOJ', 'VIQMUE', 'DOGGIQ', 'PIPXIV', 'GAWQAW', 'ETATIC', 'QAFFAG', 'BOXJUU', 'TIDQAY', 'FIVXUE', 'PIBDEK', 'ABESUV', 'BACLUO', 'PURROJ', 'OSUCIP', 'WUXLUW', 'QEBSUL', 'AKIHIL', 'JUKRIQ', 'RAKMOH', 'UQAZAP', 'EBUFOX', 'MOYPIZ', 'MIKVOS', 'YEXFAL', 'NEKQID', 'PIWNIT', 'DINVIH', 'FOFZEF', 'PADQIU', 'UROPOH', 'WEKMAB', 'NEJMEU', 'ADUHIR', 'XACJET', 'WIBHOH', 'GIXMEF', 'QUBHEC', 'REPQUB', 'YIJTER', 'DEJHIJ', 'WOFYOG', 'OZUQIK', 'UGIFAR', 'IZUDOX', 'YUFBIL', 'KEYVUF', 'XESDEE', 'ADOXOH', 'VIXLOE', 'POBJUL', 'BATCUV', 'CEJZIC', 'KACHED', 'POVLIW', 'VEHLEZ', 'BABVOR', 'LANQUM', 'YEVVIF', 'TIBXEH', 'ZUKGOB', 'ZAGHEX', 'NAZBEV', 'ADEZIS', 'ROBQOQ', 'SAKDEN', 'JEBJEF', 'RISZUR', 'CEVVII', 'XIKRIS', 'PIMVAI', 'UJULAN', 'NUXPAX', 'XUSYOZ', 'DAXSOK', 'XIZQON', 'HAHSIU', 'MIJCAK', 'YEDPII', 'IZOVUP', 'AHUYOT', 'OPALOH', 'LUNBOM', 'WOPHAK', 'NOVVAV', 'URONAR', 'VACXUS', 'WOLQAP', 'WIWGEO', 'QAJCEL', 'YIWWOR', 'FUFMEY', 'WUGPAR', 'MIJKEV', 'XABBEI', 'CPCBMO', 'JOQKOP', 'RAGMIW', 'PIRJOP', 'RAFXAX', 'FOGMIX', 'QUGFAB', 'NAJWAY', 'QIJDUI', 'DOXQMO', 'MADREP', 'VOXBAL', 'FOPQAD', 'SUSXEJ', 'OKOWAM', 'DENPES', 'VETRAP', 'EFALOM', 'WAKNOL', 'XULGAO', 'JOPFOM', 'DITGET', 'FICHEH', 'EVIRUX', 'VIFTIN', 'CIJWUN', 'ZAZXAZ', 'FUTSUJ', 'YETBII', 'HAJNIP', 'NUFWUI', 'QASDIX', 'XEMPIP', 'HARKIU', 'YARBAV', 'ZOWPAC', 'CANGUT', 'SOHZIZ', 'CUGSUS', 'BAHYEO', 'MOMGAW', 'XAMKUT', 'WEDXIM', 'FAJGOO', 'MAZVAN', 'MAFXUN', 'LUFMAA', 'APUQAD', 'OYULEA', 'UMIGAY', 'XILXUL', 'SIXCIO', 'DERGAI', 'GENXAZ', 'IPUYUN', 'BAKVOY', 'CIHVAT', 'HECZOF', 'VAQFAU', 'EADJIC', 'WOLPIW', 'VALJEX', 'CIQSAY', 'YOXWUD', 'YEJFOJ', 'BUKJOH', 'OBORUS', 'RITKEN', 'YOKVOM', 'XELMIK', 'ROLDEE', 'XOVPUT', 'MEPKAU', 'SAWMIP', 'ELODUG', 'JUMMIN', 'WAGYIM', 'IHIMIV', 'VEQLAE', 'VOBMEH', 'SASYUI', 'IBUMUP', 'LUGBOE', 'ACOJEK', 'EKOFAM', 'LICPET', 'KANNES', 'ONEZOY', 'WEKREJ', 'EJUXOX', 'OKOVEP', 'GUPQEN', 'DAJWIX', 'KUCXIQ', 'GUDQEB', 'ACEBIV', 'MEPKOI', 'PITQUE', 'BACDER', 'GUPSAL', 'TACVAV', 'QEFCUC', 'GUGBUH', 'TUCSOA', 'PUFQEO', 'ECOSOH', 'WEDGAO', 'ZIFLIJ', 'TEBVAZ', 'COLXAD', 'PPTCFE', 'ZOJHOX', 'YAVNAL', 'KENVAA', 'GAVQEB', 'XOSPIG', 'IBUVUW', 'SUHMAK', 'KEPLAT', 'OGELOE', 'QOZVUY', 'ODACON', 'BUJSIK', 'XALSAE', 'KAQDIQ', 'YOGBAX', 'XIKXEX', 'QOVVON', 'EGETAL', 'BIDHOL', 'DAQQET', 'CEYRIH', 'DOVNUZ', 'VIBYEL', 'SIWVAV', 'SECHEN', 'VICLAU', 'AGIGEC', 'TUJHEM', 'WOKWOI', 'UFABIN', 'ZEBSAA', 'QEYTOD', 'XUXKUX', 'BOGCEH', 'VIRHEK', 'QULZIG', 'ADAQUU', 'QISPEO', 'ATOJUO', 'BOCWOI', 'ZOMMUM', 'FIZTIT', 'XERJIQ', 'KULPEO', 'ALIPIU', 'FEYRUZ', 'MESGEW', 'RAZHOS', 'MISQEM', 'DUNVAL', 'UWEPUJ', 'KOFKUO', 'LEXLOQ', 'HEVJIC', 'QUDDEZ', 'IHEQAO', 'MOWVOK', 'GAXZEN', 'FIHWUO', 'RESYUK', 'NOBBEN', 'PABZAV', 'WOYDOD', 'GIZZIA', 'IGESUI', 'VIRFEI', 'YAKXOZ', 'YAJSEH', 'WORNAS', 'KAQFAJ', 'DOPDET', 'ICETAN', 'QEYBEC', 'ZUNWEK', 'WUVPOT', 'TLSCRU', 'QUGMOW', 'VAMDUJ', 'ONEPAA', 'XAXHUC', 'MIGLOF', 'HEVKUP', 'LINTAF', 'CUHVAC', 'IHITAW', 'CPYFEM', 'SEVHOT', 'SOWVUW', 'CAYCAH', 'EZISUD', 'IBINUC', 'YIKZUN', 'SINXEV', 'KARRED', 'QOYFOB', 'TUMKIW', 'HUBVOR', 'CIKZOM', 'YUTSIQ', 'QOKKIL', 'ODECIJ', 'PAFFIN', 'QAHVOK', 'PIRDID', 'ZOWPOS', 'TAYMIS', 'QEXMEL', 'GAKLAI', 'LIRFAU', 'VAHBIR', 'AHUSOL', 'YENMUA', 'HUPMEL', 'SIVYIG', 'BASLAI', 'HEBPOT', 'SOTPEX', 'WOWPED', 'QAWXIW', 'QEZVIB', 'WAPDID', 'ZUZLIP', 'VUCCEC', 'CATSAU', 'BAZMEX', 'CASPAO', 'YUCMIU', 'FUNKUV', 'VIJGEB', 'GUKCOG', 'BISXIK', 'FIBPIR', 'CIYJUP', 'HIDKIP', 'RABYIF', 'KISZET', 'CUDKUI', 'NESCIX', 'UJELOL', 'AGABEQ', 'EVUYIE', 'KIKNAV', 'UHOWOE', 'DOVXUH', 'GAPMIT', 'BUCKOA', 'KIJDEN', 'MIHLAT', 'TAPDOE', 'PIYRIY', 'RURWAD', 'UTEGEG', 'VICDIV', 'NAPQEA', 'HIBVIA', 'XEZWAB', 'IPAPUK', 'EFOCOR', 'QADROC', 'UTOSED', 'KUMFAA', 'XACQAV', 'HAGTUH', 'EYEVIN', 'VECLOH', 'CENMEP', 'QEXSUH', 'BEJFIH', 'AWOQIN', 'GEGZEZ', 'WEYSEZ', 'NOGTIM', 'ZULROP', 'DICQUE', 'MAFSIY', 'MYNBCO', 'VAHHER', 'XOFZIB', 'ASELAL', 'LEHDOT', 'MIKHAQ', 'VIKNIO', 'MEGMES', 'ARIHIU', 'FOPWAI', 'FUQQAL', 'SIZWOP', 'UTOJEU', 'UKEDAR', 'KABZUJ', 'HITJAW', 'UMEGOK', 'HITDET', 'AVAQUK', 'ZAPQIT', 'CEMSEU', 'KIYTER', 'YEJCUL', 'TILBUO', 'REGTAY', 'UQIJOV', 'XANQUZ', 'PUCVUF', 'NAQDUH', 'CIYLAZ', 'ZIGLUW', 'VUGBUW', 'ZITFEN', 'RAXPEL', 'VULVON', 'YODMAH', 'ATOQOQ', 'FACFIA', 'POXQAV', 'METVIT', 'JOMJAW', 'DISDOA', 'LEDVUM', 'UCUQAK', 'VUXFAX', 'MAZVIV', 'GEBQEL', 'UQIJEL', 'TILBAU', 'HIPXOV', 'ZIKDUS', 'SEDYEI', 'FOJJUM', 'GAHRUC', 'SULBOQ', 'BAFCIW', 'SUQNIB', 'CBYPRH', 'NUYBEQ', 'FOKGER', 'LIWTOA', 'UCUMEM', 'LIPCOE', 'KUGTUD', 'HETNUQ', 'GIKJOB', 'TOLKEN', 'JOSYUL', 'JOVFOS', 'JUSQUJ', 'MEFTIB', 'MAPTUT', 'ABEVAH', 'MIDCOU', 'FEHNEL', 'OVETUF', 'QIFRAY', 'ZORNUR', 'LOLKON', 'YOZDUN', 'ZUFHUD', 'YUWNUA', 'CEFFOK', 'KUDKEB', 'NUXWEK', 'SIZTEC', 'LERKAX', 'PEZTIA', 'TUCWIZ', 'CESFUD', 'LERSIK', 'PIYBIL', 'KISKAZ', 'UQEPEM', 'BEJKOQ', 'HAYHEV', 'TUXGEY', 'ETESED', 'QISHIK', 'CIXHEY', 'QUQVED', 'LIYWUM', 'VUJCIM', 'LAVQOO', 'UFIYUG', 'LOFYEN', 'VAQJEE', 'JUXZOR', 'PUXYEO', 'GUWSOG', 'USOQUP', 'REDQIC', 'NEXBOH', 'HUHKAW', 'HESDOB', 'UQAKII', 'MEHNAO', 'ZOPJAP', 'BAZLAS', 'UPAKIH', 'OZUJAV', 'UFUPAO', 'IQCTNI', 'BAVLIT', 'HOLTAD', 'PIJNAY', 'QAZMIQ', 'ODIHUG', 'XOMXEE', 'FIJDOT', 'VARSIR', 'YOCFAZ', 'FEVROQ', 'WUTYOA', 'CUMWUC', 'TERCUQ', 'ROMRAO', 'XIQVOL', 'CAYWIJ', 'QATFIC', 'WERZAX', 'TOJHEF', 'LIQNEG', 'QAHCUY', 'XUSFOG', 'SURZIP', 'NIJCAN', 'GUKKII', 'KADGOM', 'CAKVOC', 'MALTID', 'CEFZAR', 'JUXVIJ', 'TOFMPD', 'SILCEV', 'VOQQAW', 'HUMPUC', 'NETTOX', 'GIYHEE', 'KOLTIQ', 'OYOPEY', 'UDIGAQ', 'ADUMAP', 'NIXANT', 'HOMWOX', 'MUFJUT', 'HAFMEH', 'BEXXOS', 'SOVBEL', 'TUYPUY', 'RICDOY', 'NILNUS', 'MAXGIC', 'LUXZOT', 'EQONII', 'NOXJAN', 'IDUFIY', 'AZUVIA', 'EFAGEX', 'WUMJEU', 'VAWKOT', 'SADRUN', 'QIKQAF', 'HEWDAO', 'QUMLEP', 'BADCOA', 'TIPWUK', 'PIJWUB', 'SESXUL', 'KABDAT', 'WEKZEU', 'BOLNEV', 'CIKBOQ', 'HIMFEP', 'ODABED', 'AHIHUU', 'BOZPEM', 'LEMBOX', 'NOXMEU', 'TEGZEM', 'LEQPOM', 'LUMSOC', 'MEHHAL', 'QEKCER', 'BUVNAI', 'BODRUI', 'WEXHOZ', 'PIBPIB', 'GAXDOY', 'WEPXAS', 'QOYMIC', 'UCECIQ', 'QEXRET', 'FAYVUZ', 'CAJNOQ', 'YEFLAW', 'AHABOB', 'YOSVOR', 'VUBXEX', 'CAWTEA', 'POXRAW', 'INICAJ', 'LIYXAT', 'OLIQAC', 'NUCDOG', 'CUCNAR', 'IMOSAG', 'DOFVEB', 'ZIWFUI', 'GADKUT', 'FAYDEO', 'XANHEB', 'JONCOE', 'RAZLAF', 'REFPOI', 'JIKBEN', 'QUCCOI', 'HOFBIN', 'QAMHOD', 'AQEZUS', 'RASYER', 'SIVBOQ', 'ACARUT', 'JIJQOL', 'KIHBUB', 'AZIDIY', 'BAHFOI', 'KEQSUU', 'TEHBEO', 'LINZUF', 'IXENEG', 'ZIBLEB', 'HUKSEL', 'UTIBOP', 'CUHXUZ', 'NUCPAC', 'XELMOR', 'QEJWUA', 'FITDOD', 'ZEYHIX', 'WEPYEY', 'OLOWUI', 'CAWLIV', 'DABZIS', 'MEZWAP', 'KEFCEF', 'QERTAL', 'XESPOB', 'NORKAG', 'RITHUZ', 'SUYXOB', 'PUPNUL', 'XEDCUG', 'DUDRIF', 'NEVXET', 'JOCVIJ', 'DARTIC', 'AVUSUF', 'REBQIY', 'IBESAJ', 'LIJTUW', 'MAFGIL', 'TULZOQ', 'QEGFOA', 'IBUGIX', 'FIRWEJ', 'MEBNEM', 'AWENEW', 'XENCUR', 'TEBDOV', 'MIMBUH', 'PIWJEL', 'ZODZID', 'OWEYOG', 'LOYXEC', 'NETFOJ', 'HUBYIM', 'NEYPUE', 'DUGWUX', 'QOZVEH', 'XABMIY', 'WULFOZ', 'BIKCUW', 'KEFJUC', 'YERYEC', 'COBPIT', 'TIYWIK', 'SUTPED', 'FOKFAM', 'TIZYUX', 'CMTPCU', 'FAYGAQ', 'ZURLON', 'TAHVEG', 'RISPIT', 'AMOKAQ', 'QAFYUU', 'LAQLUL', 'GEMLAM', 'DEFJON', 'UYOBER', 'UCIWIN', 'ROMFII', 'YUWWUI', 'FIRXEK', 'QAQVOU', 'PIJWIQ', 'ZAQHOR', 'DEDGAU', 'PATQAF', 'EPADOQ', 'ACIQEK', 'MILJOG', 'YEHJOL', 'MITKOQ', 'SUSBOZ', 'FIQTIL', 'OMAXUV', 'XUYVIX', 'ZARYIC', 'JAGWIY', 'EKALUY', 'FOTNOR', 'QURGAN', 'HOCLOA', 'TUQRIG', 'MIQVAL', 'DPIMZN', 'LAZGAU', 'PAZCUQ', 'GEQXEG', 'PADTIZ', 'IXIGIG', 'PEDKAK', 'NIPXUF', 'ICONOE', 'ESIRII', 'IMULOS', 'POXRIE', 'CODRIX', 'FOMQII', 'UGOKUW', 'XEYCOW', 'LIKXIN', 'SUNNAR', 'ZIYXEK', 'AGAMEB', 'PODVIO', 'SEZQAR', 'YUTNUY', 'TIQXAS', 'MIPLOM', 'BEZKOH', 'JONXIV', 'DIVLUQ', 'XUQCOB', 'SOKQUG', 'CONZEL', 'QODDAP', 'NOHZAM', 'NABMAF', 'MEHZUU', 'TAYQOA', 'DUVHIN', 'RALCUE', 'XAMFUN', 'HACFEY', 'IRETOP', 'JUBSUU', 'DUBPAT', 'BESXAB', 'SETMOW', 'SEMPUV', 'KOPYOD', 'GUWYED', 'YIHNEK', 'CEKWEX', 'RUQSOL', 'FOLRUW', 'TAMQIH', 'TAWJIM', 'GODCOU', 'NEKWIJ', 'QIFFES', 'XIRMER', 'KOFGET', 'UCASOG', 'XOPFUD', 'EBUSOL', 'YOHHIM', 'TUTCEQ', 'KEXNOR', 'MOLMEH', 'BERZII', 'QEZVUP', 'XIYGIW', 'CUCSUQ', 'PEYHGP', 'HIWCUN', 'HOLHOH', 'NUJSAN', 'VEFQUU', 'HOXYUQ', 'ISOSIU', 'MUHWUH', 'WIWGEQ']


class MoleculeDatabase:
    
    def __init__(self, data_path: str, id_col: str, TestSize=False):
        self.data_path = data_path
        
        self.global_props_path = Path(self.data_path, 'global_mol_properties.csv')
        self.atomic_props_dir = Path(self.data_path, 'atomic_properties')
        
        self.output_data_path = Path('..', 'data')
        self.default_atomic_properties_json_path = Path(self.atomic_props_dir, 'atomic_properties.json')
        
        self.id_col = id_col
        self.testing = TestSize
        
        self.has_df = False
        self.has_atomic_props = False
        self.extraction_error_count = 0
        self.extraction_errors = {}
        
        self.all_molecule_ids = self.get_all_molecule_ids()
        
        self.xyz_props = ['x', 'y', 'z']
        self.other_atomic_props = ['partial_charge']
        
    def get_global_properties(self) -> pd.DataFrame:
        """
        Returns a deepcopy of the internal global molecular properties dataframe.
        :return:
        """
        self.df = pd.read_csv(self.global_props_path)
        self.has_df = True
        
        return deepcopy(self.df)
    
    def get_all_molecule_ids(self) -> np.array:
        """
        Returns all molecular ids. The id is a unique string for each molecule, for example the CSD code.
        :return: molecular ids (np.array)
        """
        self.ensure_global_properties_present()
            
        ids = self.df[self.id_col]
        assert ids.is_unique, 'Molecular IDs are not unique.'
        ids = ids.values
        
        if self.testing:
            test_molecule_ids = self.get_test_molecule_ids(mol_ids=ids)
            ids = np.intersect1d(ids, test_molecule_ids)
        
        return ids
    
    def read_atomic_properties_file_header(self, filestr: str):
        lines = filestr.split('\n')
        num_atoms = int(lines[0])
    
        columns = lines[1].split(ATOMIC_PROPERTIES_COLUMN_SEPARATOR)
        mol_id = columns[0]
        comment = columns[-1]
        col_names = columns[1:-1]
        
        return num_atoms, mol_id, col_names, comment
        
    def read_single_atomic_properties_filestring(self, filestr: str,):
        num_atoms, mol_id, col_names, comment = self.read_atomic_properties_file_header(filestr)
        
        file = io.StringIO(filestr)
        txt_array = np.genfromtxt(file, skip_header=2, dtype=str)
        atoms, value_array = txt_array[:,0], txt_array[:,1:].astype(float)
        
        atoms = [str(a) for a in atoms]     # numpy to regular python list of strings for json compatibility
        col_names.pop(0)
        
        values = {name: column.tolist() for name, column in zip(col_names, value_array.T)}
        
        assert num_atoms == len(value_array) and num_atoms == len(atoms), 'Number of atoms specified in atomic properties file is not equal to the real number of atoms included.'
        return mol_id, atoms, values, comment
    
    def atomic_props_long_to_wide_format(self, atomic_props_long: dict, prop_names: list):
        if isinstance(prop_names, str):
            prop_names = [prop_names]
        assert all([name in atomic_props_long for name in prop_names]), f'At least one name of {prop_names} not in atomic property dict.'
        
        all_atoms = atomic_props_long['atoms']
        
        atomic_props_wide = {}
        for i, atom in enumerate(all_atoms):
            prop_values = [atomic_props_long[name][i] for name in prop_names]
            atomic_props_wide[i] = [atom, prop_values]
        
        return atomic_props_wide
        
    def read_full_atomic_properties_file(self, path, sep='\n\n'):
        with open(path, 'r') as file:
            full_file = file.read()
        files = full_file.split(sep)
        
        atomic_props = {}
        for filestr in tqdm(files, desc='Reading atomic properties file'):
            mol_id, atoms, values, comment = self.read_single_atomic_properties_filestring(filestr)
            assert not mol_id in atomic_props, f'Molecular id {mol_id} not unique in file{path}.'
            values['atoms'] = atoms
            values['comment'] = comment
            atomic_props[mol_id] = values
            
        return atomic_props
        
    def get_atomic_properties(self) -> dict:
        """
        Returns a dictionary of all atomic properties of all molecules. The form is {mol_id: {'prop1': list}}.
        :return: atomic properties (dict)
        """
        all_atomic_props = {}
        
        pattern = Path(self.atomic_props_dir, '*.xyz')
        self.all_atomic_props_paths = sorted(glob.glob(str(pattern)))
        print(f'Found {len(self.all_atomic_props_paths)} atomic property files. Start reading in.')
        
        for atm_path in self.all_atomic_props_paths:
            atomic_props = self.read_full_atomic_properties_file(path=atm_path)
            
            assert not any([mol_id in all_atomic_props for mol_id in atomic_props.keys()]), 'Molecular ids of atomic property files not unique.'
            all_atomic_props.update(atomic_props)
        
        self.atomic_props = all_atomic_props
        self.has_atomic_props = True
        return all_atomic_props
    
    def save_atomic_properties(self, outpath=None):
        print('Start saving atomic properties to json.')
        
        self.ensure_atomic_properties_present()
        
        self.atomic_properties_json = outpath or self.default_atomic_properties_json_path
        with open(self.atomic_properties_json, 'w') as file:
            json.dump(deepcopy(self.atomic_props), file)
            
        print(f'Saved atomic properties to {self.atomic_properties_json}.')
        return
    
    def load_atomic_properties(self, json_path=None):
        print('Start loading atomic properties from json file.')
        
        json_path = json_path or self.default_atomic_properties_json_path
        with open(json_path, 'r') as file:
            self.atomic_props = json.load(file)
        self.has_atomic_props = True
        
        print('Loaded atomic properties from json file.')
        return self.atomic_props
    
    def ensure_global_properties_present(self):
        if not self.has_df:
            self.get_global_properties()
    
    def ensure_atomic_properties_present(self):
        if not self.has_atomic_props:
            self.get_atomic_properties()
            
    def get_Extracted_Molecule(self, mol_id: str, atomic_props_long: dict, add_atomic_props: list) -> Extracted_Molecule:
        self.ensure_global_properties_present()
        
        mol_atomic_props = {name: atomic_props_long[name] for name in add_atomic_props}
        
        coordinates = self.atomic_props_long_to_wide_format(atomic_props_long=atomic_props_long, prop_names=['x', 'y', 'z'])
        
        extr_mol = Extracted_Molecule(
                                        coordinates=coordinates,
                                        csd_code=mol_id,
                                        atomic_props=mol_atomic_props,
                                        # global_props=mol_global_props
                                        )

        return extr_mol
    
    def get_test_molecule_ids(self, mol_ids):
        test_molecule_ids = felix_molecules
        
        return test_molecule_ids

    def get_all_coordinates(self):
        return self.get_all_atomic_props_in_wide_format(names=self.xyz_props)
    
    def get_atomic_props_individual_in_wide_format(self, names):
        names = [names] if isinstance(names, str) else names
        
        all_atomic_props = {}
        for name in names:
            name_props = self.get_all_atomic_props_in_wide_format(names=[name])
            all_atomic_props[name] = name_props
        
        return all_atomic_props
        
    def get_all_atomic_props_in_wide_format(self, names: list):
        self.ensure_atomic_properties_present()
        
        atomic_props = {
                        mol_id: self.atomic_props_long_to_wide_format(
                                                                        atomic_props_long=self.atomic_props[mol_id],
                                                                        prop_names=names)
                        for mol_id in self.all_molecule_ids
                        }
        
        return atomic_props
        
        
    def get_global_properties_of_mol(self, mol_id: str):
        self.ensure_global_properties_present()
        
        row = self.df[self.df[self.id_col] == mol_id].squeeze(axis=0)
        row = row.drop(self.id_col)
        row = row.rename(mol_id)
        
        return row
        
        
    def get_all_Extracted_Molecules(self):
        all_coords = self.get_all_atomic_props_in_wide_format(names=self.xyz_props)
        other_atomic_props = self.get_atomic_props_individual_in_wide_format(names=self.other_atomic_props)
        
        self.all_Extracted_Molecules = {}
        for mol_id, coordinates in tqdm(all_coords.items(), desc='Extracting molecules'):
            try:
                other_atomic_props_of_mol = {name: other_atomic_props[name][mol_id] for name in other_atomic_props}
                global_props = self.get_global_properties_of_mol(mol_id).to_dict()
                molecule = Extracted_Molecule(
                                                coordinates=coordinates,
                                                csd_code=mol_id,
                                                atomic_props=other_atomic_props_of_mol,
                                                global_props=global_props
                                                )

            except Exception as ex:
                print(f"An Error occured: {ex}")
                self.extraction_error_count += 1
                self.extraction_errors[mol_id] = str(ex)
                pass
        
            self.all_Extracted_Molecules[mol_id] = molecule
        
        return self.all_Extracted_Molecules
    
    def get_Extracted_Molecules_global_information(self):
        if not hasattr(self, 'all_Extracted_Molecules'):
            raise ValueError('No attribute all_Extracted_Molecules')
        
        infos = []
        possible_attrs = ['denticity_dict']
        
        for csd_code, mol in self.all_Extracted_Molecules.items():
            mol_infos = {
                            'CSD_code': csd_code,
                            'n_ligands': len(mol.ligands)
                        }
            
            for attr in possible_attrs:
                if hasattr(mol, attr):
                    mol_infos.update({'denticity_dict': getattr(mol, attr)})
            infos.append(mol_infos)
        
        if not self.has_df:
            self.get_global_properties()
            
        self.Extracted_Molecules_global_information = pd.DataFrame(infos)
        self.df = self.df.merge(self.Extracted_Molecules_global_information, on=self.id_col)
        
        return self.Extracted_Molecules_global_information
        
            
            
        

if __name__ == '__main__':

    id_col = 'CSD_code'

    tmqm = MoleculeDatabase(data_path='../database/tmQM/data', id_col=id_col)
    df = tmqm.get_global_properties()
    row = tmqm.get_global_properties_of_mol(mol_id='WIXKOE')
    tmqm.load_atomic_properties_from_json()
    # atm = tmqm.get_atomic_properties()
    # tmqm.save_atomic_properties()
    # ids = tmqm.get_all_molecule_ids()
    mol = tmqm.get_Extracted_Molecule('WIXKOE')
    print('Done!')

        
        
        

        
        
        

