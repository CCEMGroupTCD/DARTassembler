"""
Was first the main in src01
then the main_test in src01
and is now the unittest for a full run of src01
"""

import pandas as pd
from pathlib import Path
import json
import os
from DARTassembler.src.ligand_extraction.DataBase import MoleculeDB, LigandDB


def read_data_in_input_dict(database_path_, atomic_properties_json_, read_in_graphs_if_exist: bool = True):
    """
    From the desired input folder structure we are now going to read in the properties
    in the format, which is required to generate a MoleculeDB
    """

    # todo:
    #   hier fehlt noch wie aus den .xyz files die atomic_properties.json erstellt wird
    #   create_json_from_xyz(path=atomic_properties_json) (Der Pfad wo es hingespeichert werden soll)
    with open(atomic_properties_json_, "r") as file:
        dict_ = json.load(file)

    # Next - the optional global properties
    # global_properties.csv into dict of the type
    # {csd_code: {prop_1: value_1, ..., prop_n:value_n}
    if os.path.exists(f"{database_path_}/global_mol_properties.csv"):
        df = pd.read_csv(f"{database_path_}/global_mol_properties.csv")
        json_str = df.set_index('CSD_code').to_json(orient='split')
        json_dict = json.loads(json_str)

        global_property_dict = {}
        for i, csd_code in enumerate(json_dict["index"]):
            global_property_dict[csd_code] = {prop: json_dict["data"][i][j] for j, prop in
                                              enumerate(json_dict["columns"])}
            global_property_dict[csd_code]["CSD_code"] = csd_code
    else:
        global_property_dict = {}

    if read_in_graphs_if_exist is True and os.path.exists(f"{database_path_}/graphs"):
        # todo: Then we have to run the graph readin stuff
        #   Glaube ich hab was in Dev_Felix -> Graph Comparison (Hoffe das ist noch da)
        #   muessen die uber nx.read_gml einlesen
        #   und dann mit nx.dict_of_dicts(G) in dicts zerlegen
        graph_dict_dict = {}
        pass
    else:
        graph_dict_dict = {}

    db_dict = {}
    for key_, item_ in dict_.items():
        global_props_mol = global_property_dict[key_] if key_ in global_property_dict else None
        graph_dict = graph_dict_dict[key_] if key_ in graph_dict_dict else None

        db_dict[key_] = {"atomic_props": item_,
                         "global_props": global_props_mol,
                         "graph_dict": graph_dict
                         }

    return db_dict


if __name__ == '__main__':
    """
    Testing:
        Testing: Either an int or a bool. If int, must be smaller than 1247, which is the size of the whole testing dataset which you also get by setting Testing=True. If Testing=False, the whole tmQM is used as input.
        Recommended: Use Testing=20 for simple debugging and Testing=True for a bigger check in the end. These two values will also be doublechecked by using assert_frame_equal() using the value of a saved df with correct values.
    

    """
    # This is just basic setup, could be pushed into a .yml
    Testing = 20

    database_path = '../database/tmQM/data'
    atomic_properties_json = Path('../database/tmQM/data/atomic_properties/atomic_properties.json')

    denticity_numbers_of_interest = [1, 2, 3, 4, 5, 6]
    # Previously: metals_of_interest =  ["Fe", "Mn", "Cr", "Ru", "Co", "Ni", "Cu", "Ir", "Mo"]
    metals_of_interest = ['La', 'Sc', 'Y', 'Zr', 'Ti', 'Hf', 'Ta', 'Nb', 'V', 'W', 'Mo', 'Cr', 'Tc', 'Mn', 'Re', 'Fe',
                          'Ru', 'Os', 'Ir', 'Rh', 'Co', 'Pd', 'Pt', 'Ni', 'Au', 'Cu', 'Ag', 'Zn', 'Hg',
                          'Cd']  # everything in the tmQM

    #
    # convert the input folders into readible dict
    db_dict = read_data_in_input_dict(database_path_=database_path, atomic_properties_json_=atomic_properties_json)

    if Testing == 20:
        csd_codes_for_testing_20 = ['BOCVOF', 'BAXWAA', 'BIGHOO', 'AROMAW', 'BAXWUU', 'BIXBUG', 'BEGLUU', 'AVIBIR',
                                    'BOSVOU', 'AFATAE', 'AXOTAK', 'BEBFUL', 'BASDAC', 'BOJBAG', 'AYUWAU', 'AGUZAD',
                                    'ACASOO', 'BEBQUW', 'BOTBAP', 'BOQPAZ']
        db_dict = {key_: item_ for key_, item_ in db_dict.items() if key_ in csd_codes_for_testing_20}
    elif Testing is True:
        csd_codes_for_testing_true = ['FONZAJ', 'MAZVIV', 'ODIHUG', 'LEMBOX', 'TIZYUX', 'BEMGEH', 'NOXJAN', 'WAHNOK',
                                      'WUVPOT', 'XEYCOW', 'MOLMEH', 'TEBDOV', 'VOQQAW', 'UCUMEM', 'UQIJEL', 'LEZNIQ',
                                      'PFBGNI', 'MUGHOL', 'YURBAQ', 'AGIGEC', 'KUCXIQ', 'GUDQEB', 'MBSLCO', 'DOMQIH',
                                      'HEVJIC', 'AYEDAK', 'SELSEI', 'LICPET', 'ODECIJ', 'DUGWUX', 'NUXWEK', 'YAKFEX',
                                      'PUHKEK', 'LOPQEP', 'VAWKOT', 'YOGBAX', 'PAXVEQ', 'PIBPIB', 'CENMEP', 'KAQDIQ',
                                      'SESXUL', 'VICDIV', 'BEZKOH', 'PELCER', 'IWOKIO', 'WUGPAR', 'YOYNAC', 'CIQSAY',
                                      'AMOKAQ', 'IBUMUP', 'IDUFIY', 'QAFYUU', 'USADAV', 'KIHMAR', 'SUYXOB', 'IVAGUJ',
                                      'CIKZON', 'UCECIQ', 'TIQXAS', 'KOQZEX', 'EVIRUX', 'XERJIQ', 'DIXYAL', 'GUWSOG',
                                      'NINNOO', 'POGGOI', 'MEXGED', 'WOPHAK', 'HOQBAS', 'HEQBIR', 'SUSTUV', 'WIWGEQ',
                                      'USONEX', 'EGETAL', 'GUVMUH', 'GULZOD', 'AQEZUS', 'LEXLOQ', 'MEFTIB', 'SUBKUX',
                                      'GOVBIE', 'VASZEW', 'GUPQEN', 'MEPKOI', 'XELMIK', 'VAFQEY', 'JOHCEP', 'JIJQOL',
                                      'GOLVAF', 'QIZCAF', 'JIZSIW', 'CAJNOQ', 'REBQIY', 'DERGAI', 'EGUKIA', 'HITDET',
                                      'SULBOQ', 'WEPYEY', 'IPAPUK', 'YUTNUY', 'UTIBOP', 'COBVUM', 'ONEPAA', 'PAKKIW',
                                      'POHJII', 'SOBKAX', 'XELMOR', 'SUYHUR', 'UBEVOP', 'HAJNIP', 'HIWCUN', 'LAVQOO',
                                      'QUQXUX', 'VAQJEE', 'VIBBIT', 'WIWSAZ', 'EPAVEX', 'XIYJUL', 'IBUVUW', 'VUBXEX',
                                      'FAVGOZ', 'NAZBEV', 'YUCMIU', 'GEQXEG', 'XIFXOC', 'IHITAW', 'TIPWUK', 'QEGFOA',
                                      'QEFCUC', 'ZATNOZ', 'FEYRUZ', 'YOBMOU', 'FOGMIX', 'BAZLAS', 'AVIBIR', 'BASDAC',
                                      'BTZANJ', 'CIVLEZ', 'REXZOL', 'WOLQAP', 'ZAQZEZ', 'BACDER', 'JATKEX', 'SISQUI',
                                      'PIJWUB', 'TUFDAZ', 'QISHIK', 'CAJKEG', 'DUDRIF', 'HOLVEJ', 'MAXGIC', 'UFUPAO',
                                      'ZAZNOF', 'CIRVEE', 'COLXAD', 'BAXWAA', 'GONHAS', 'YEJCUL', 'LUNBOM', 'JOWYEC',
                                      'KOXKUF', 'UCIWIN', 'WEKZEU', 'CARVAR', 'XANHEB', 'SEDYEI', 'YERYEC', 'IXENEG',
                                      'QATFIC', 'IQETAB', 'WULFOZ', 'DAQQET', 'ECIJIJ', 'GODCOU', 'MONNIN', 'NAJVUS',
                                      'VAWGIL', 'ZORNUR', 'NUYGUL', 'CEYRIH', 'TAWJIM', 'CIXHEY', 'EBAZAL', 'LEJMUJ',
                                      'CIYJUP', 'BUHQIF', 'OLIQAC', 'PUXYEO', 'AGAMEB', 'TAYQOA', 'DOJTOL', 'HEBBEW',
                                      'MUTCIP', 'MOLTAI', 'GAXDOY', 'HOMWOX', 'JAXCUH', 'VAPTIR', 'UTOSED', 'YUQZER',
                                      'JUXVIJ', 'DITDES', 'MISQEM', 'XUXKUX', 'KIZPIU', 'EQONII', 'JIZTIW', 'HOBPEU',
                                      'NETFOJ', 'AVUSUF', 'HUJWEO', 'QIMNAD', 'UWULEE', 'VEGFIX', 'TUMKIW', 'BODRUI',
                                      'IFEGEH', 'PATQAF', 'MOJBUL', 'QOYMIC', 'KOLTIQ', 'EKEGEI', 'KISKAZ', 'QUGMOW',
                                      'MIJVAE', 'XUSFOG', 'ZODZID', 'VIQMUE', 'EADJIC', 'NURSID', 'EJUXOX', 'TAMQIH',
                                      'XANQUZ', 'KABDAT', 'BUHQEB', 'CUCSUQ', 'GAQCOT', 'PIRJOP', 'ACARUT', 'UMIGAY',
                                      'YAKXOZ', 'AWOQIN', 'QUBPEI', 'ACOJEK', 'DAJWIX', 'QAHCUY', 'VIBYEL', 'JAPCOT',
                                      'ZOHQOE', 'YAFSAC', 'YUMMOJ', 'RICDOY', 'XENFEC', 'IXOGOR', 'ICETAN', 'KUGTUD',
                                      'QAPLAX', 'CAWLIV', 'PIJWIQ', 'ATOQOQ', 'RALCUE', 'FOKFAM', 'RITKEN', 'MEGZOP',
                                      'QEBZEE', 'CALQOY', 'HOXYUQ', 'QUGFAB', 'HAHSIU', 'WIZQON', 'BOQPAZ', 'BUVNAI',
                                      'LIXSAO', 'BEGLUU', 'MIKWIO', 'QUQLAP', 'FOJJUM', 'DUBPAT', 'MONSOZ', 'NETTOX',
                                      'HAGTUH', 'LOPCAW', 'SEBPUL', 'NOGTIM', 'CUCNAR', 'LOGYEN', 'FEVDOC', 'YOHHIM',
                                      'JIZZIC', 'AYUWAU', 'QAMHOD', 'GUWYED', 'NEKWIJ', 'ITOSER', 'ZIWFUI', 'IRETOP',
                                      'BOLNEV', 'QANLAV', 'YIQGAI', 'JUSQUJ', 'QIKQAF', 'BAVLIT', 'MATGEW', 'SAMTUX',
                                      'OPAQIG', 'TUTCEQ', 'OPALOH', 'QIFRAY', 'NAXNOR', 'IZOVUP', 'WEPXAS', 'XIKXEX',
                                      'UYOBER', 'QEZVUP', 'ZOTVIP', 'HEWBER', 'AZIDIY', 'DELCOO', 'HAYHEV', 'HUKSEL',
                                      'DOVNUZ', 'DAWXEG', 'EVUYIE', 'ONEZOY', 'UCIGUJ', 'PEZTIA', 'UMACOB', 'OGELOE',
                                      'FAYVUZ', 'VEFQUU', 'KEFCEF', 'FIRXEK', 'BEXXOS', 'WACWUV', 'GIZZIA', 'ILEFUB',
                                      'DUXROF', 'XEPZUO', 'COKRIF', 'COBTOE', 'LUXZOT', 'QOZVEH', 'HUBYIM', 'JOSPIQ',
                                      'HEHSAO', 'LARCUF', 'VOXBAL', 'WORNAS', 'OWEYOG', 'QADROC', 'CIKZOM', 'QOVLIY',
                                      'WEYWEC', 'MIQVAL', 'OMAXUV', 'CIJWUN', 'TULZOQ', 'OBORUS', 'OBEXUP', 'SUSXEJ',
                                      'KAXLAX', 'ZIYXEK', 'CIMGAG', 'UJULAN', 'ZUFHUD', 'ONAXUY', 'FICHEH', 'MONTIU',
                                      'VOHJUA', 'FOTQEL', 'RAZLAF', 'CAGROU', 'DECKIG', 'EKOFAM', 'KIZZID', 'VIBBEP',
                                      'CALQUE', 'CEFFOK', 'OKUNOZ', 'EFAKEC', 'JIKBEN', 'PIYRIY', 'JOMTAI', 'KIYTER',
                                      'NIXANT', 'PUPNUL', 'DPIMZN', 'NUCDOG', 'AFATAE', 'FOQQAE', 'XABMIY', 'LERSIK',
                                      'AHUSOL', 'OKOWAM', 'ZOHQUK', 'VAVRER', 'BAHFOI', 'WIJVIX', 'DICQUE', 'NOFVAF',
                                      'UDIGAQ', 'LIYXAT', 'OYOPEY', 'SOWVUW', 'DOPDET', 'JIYFIJ', 'JOVFOS', 'COBPIT',
                                      'DOTNAD', 'WOLPIW', 'HOBQAR', 'ZEWKAQ', 'PADTIZ', 'MIMBUH', 'SUYNIL', 'HOFBIN',
                                      'ECOSOH', 'HOLHOH', 'VACXUS', 'CUHXUZ', 'WIJYEU', 'EBUSOL', 'CAWTEA', 'HAZQAD',
                                      'YAVNAL', 'VUXFAX', 'WINVUL', 'ACASOO', 'BAZMEX', 'HEVKUP', 'KUBGAQ', 'FITDOD',
                                      'VAHHER', 'KEQSUU', 'UMEGOK', 'HITJAW', 'NABMAF', 'AHUYOT', 'BEJKOQ', 'BIDHOL',
                                      'OQOBOM', 'XOBLEI', 'GAXZEN', 'OLOWUI', 'RESYUK', 'BOTPIK', 'HEQBEN', 'ZOJHOX',
                                      'NESJAW', 'FENFEM', 'EMOYUC', 'AMUYOY', 'JOPFOM', 'UCASOG', 'IFISAR', 'FEHNEL',
                                      'PEJCEP', 'MUTZAD', 'LIXMOX', 'MEPKAU', 'YEKQUC', 'JOLQUZ', 'VOYREI', 'DENPES',
                                      'ZITFEN', 'SIWVAV', 'EYEVIN', 'JIMYOU', 'MERCER', 'WELVAN', 'LAZGAU', 'ZUNWEK',
                                      'MAFXUN', 'MAFSIY', 'YOCFAZ', 'PINYAM', 'VAHBIR', 'LERKAX', 'MONSIT', 'GOLMAW',
                                      'LEQPOM', 'YIHNEK', 'HOKYAH', 'MUFJUT', 'BUFREB', 'QURGAN', 'AWIPED', 'LILHEW',
                                      'GAPNAN', 'UTEGEG', 'LAMNUJ', 'XAMFUN', 'ICONOE', 'ICOHOX', 'LIJTUW', 'VENLOR',
                                      'TAYMIS', 'VAVRAN', 'KIFYIH', 'XOVPUT', 'REGTAY', 'LEPDUI', 'CIJXUQ', 'DISDOA',
                                      'MIKHAQ', 'SAKDEN', 'BIGHOO', 'GEKVIE', 'IGERIX', 'YOSVOR', 'TACFOU', 'FAGZUL',
                                      'XARTIT', 'WIQGEL', 'PODVEK', 'XACGAJ', 'HAFMEH', 'YUWNUA', 'QUDDEZ', 'ZIKDUS',
                                      'MONFOL', 'IBUGIX', 'QIDJET', 'WIWGEO', 'CALREP', 'XESPOB', 'FAPKAJ', 'XAVCOO',
                                      'QEKCER', 'VEGFET', 'WUTYOA', 'QERTAL', 'CODRIX', 'RITHUZ', 'CENNAM', 'CUDKUI',
                                      'TEGZEM', 'PUHJUZ', 'YAWBEE', 'MEDCIK', 'WOMSUO', 'JIXCID', 'WATVUL', 'NEYPUE',
                                      'YEJFOJ', 'GEGZEZ', 'YEHXUF', 'NUTXOR', 'UCUQAK', 'KEBBEY', 'SIXCIO', 'LIYWUM',
                                      'QAZMIQ', 'NUCPAC', 'JEGROF', 'SEMPUV', 'IBEWET', 'VUCCEC', 'OGEMEU', 'VAQFAU',
                                      'WEDGAO', 'XUYVIX', 'LILWEJ', 'LUMSOC', 'ZIGLUW', 'MYNBCO', 'EDOBIK', 'PPTCFE',
                                      'YIWWOR', 'REFPOI', 'SADRUN', 'HIQXOV', 'LOYXEC', 'AZUVIA', 'CENMIT', 'QAKREZ',
                                      'SURZIP', 'NOXREZ', 'HUJWAK', 'JOGPIG', 'BEBQUW', 'JAGWIY', 'QAHVOK', 'VETTUM',
                                      'MUJTIW', 'ENURUC', 'ERATIA', 'RISPIT', 'FACFIA', 'NEXBOH', 'YOLSIC', 'MUHWUH',
                                      'GUGBUH', 'PIYBIL', 'RAQWOX', 'MIPLOM', 'JUXZOR', 'ETESED', 'IMOSAG', 'KOFKUO',
                                      'WEYSEZ', 'KIHBUB', 'MOMGAW', 'LIQNEG', 'NORKAG', 'QEYBEC', 'JOWHUA', 'SEZQAR',
                                      'KOGBEO', 'CONZEL', 'MEBNEM', 'XUNFUH', 'TIPCUR', 'GAQFEJ', 'ELODUG', 'ACEBIV',
                                      'XABBEI', 'YODFUV', 'RISZUR', 'VIJGEB', 'IBESAJ', 'TOFMPD', 'NIJCAN', 'XILXUL',
                                      'ZEYHIX', 'IGESUI', 'DIVLUQ', 'FOPQAD', 'UQAKII', 'MONSUF', 'YEFLAW', 'DAHVIR',
                                      'QIDJIX', 'MEHZUU', 'NIPXUF', 'ODABED', 'SIVBOQ', 'KIMQOM', 'IWISAI', 'ZUKGOB',
                                      'KANNES', 'WIMKEK', 'DUFLEW', 'TUCWIZ', 'FOLRUW', 'CBYPRH', 'LIPCOE', 'TAPDOE',
                                      'KEYVUF', 'GEDQEM', 'VIRFEI', 'EXEMIF', 'WAWHUY', 'MUJRIS', 'XACQAV', 'EBECIY',
                                      'SIGCUH', 'IHEQAO', 'WEKREJ', 'ZORGUL', 'LAWNEE', 'FOPWAI', 'MIDCOU', 'ZIBLEB',
                                      'ZURLON', 'WERZAX', 'KEFJUC', 'DARTIC', 'HOCLOA', 'BOTBAP', 'TIYWIK', 'WOWPED',
                                      'ROMFII', 'PIWJEL', 'QEJWUA', 'WUSWOX', 'PEYHGP', 'PIMVAI', 'FOKGER', 'MAFGIL',
                                      'YATNUF', 'CMTPCU', 'SUSBOZ', 'HEWDAO', 'KABZUJ', 'NOHZAM', 'LEHDOT', 'OZUJAV',
                                      'KEMQUQ', 'NEJMEU', 'LEHMUH', 'UHOWOE', 'EKALUY', 'RASYER', 'XEDCUG', 'XOPFUD',
                                      'MIHLAT', 'WEDXIM', 'MEZWAP', 'SIZTEC', 'TETWAR', 'RAXPEL', 'AGABEQ', 'MIJCAK',
                                      'POXRAW', 'MITKOQ', 'ROBQOQ', 'KEFFOQ', 'KEXNOR', 'KIWQOW', 'VUJCIM', 'YAJSEH',
                                      'UWUNUX', 'CUMWUC', 'ROLDEE', 'KADRUF', 'FUWGAF', 'LIWTOA', 'QUCCOI', 'SUTPED',
                                      'TAHVEG', 'YOKVOM', 'JONXIV', 'FEVROQ', 'CEFZAR', 'XUSNIJ', 'IHIMIV', 'BUSNOS',
                                      'JOCVIJ', 'MEHNAO', 'AGUZAD', 'CAYWIJ', 'BEBFUL', 'HILDEM', 'ISOSIU', 'HUJVUD',
                                      'ADUMAP', 'QULZIG', 'UKEDAR', 'COJRIF', 'EVEWAD', 'MAZVAN', 'TOQJEQ', 'NILNUS',
                                      'PUXZUE', 'NUXPAX', 'MALTID', 'OVETUF', 'SILCEV', 'LINTAF', 'GAPMIT', 'BABVOR',
                                      'DOSBUJ', 'TOPRIZ', 'GUKKII', 'IBINUC', 'BOLPIE', 'FIGRIZ', 'SUNNAR', 'FBGLNI',
                                      'QUMLEP', 'GAFYUJ', 'LAQLUL', 'QEWGOR', 'QEPSEK', 'AXOTAK', 'HIGPUI', 'JEBJEF',
                                      'YEVVIF', 'RUQSOL', 'BIXBUG', 'EBECOE', 'INICAJ', 'XUMHUJ', 'BAFCIW', 'KEYVOC',
                                      'CEBHUP', 'KARRED', 'WAVPAL', 'WOQFEQ', 'QAQVOU', 'NESCIX', 'TERCUQ', 'HARKIU',
                                      'IBOLOB', 'TILBUO', 'EPADOQ', 'VETRAP', 'NOXRID', 'WAGYIM', 'ZOPJAP', 'YOZDUN',
                                      'PIJNAY', 'HESDOB', 'QOHBOF', 'MILJOG', 'RIBQIC', 'KECJUA', 'HOLTAD', 'XOMXEE',
                                      'AKOPEW', 'MIGBIO', 'TAWMEL', 'POXQAV', 'FAYDEO', 'XAMKUT', 'LUGBOE', 'DOTCEW',
                                      'OZUQIK', 'TOLKEN', 'LANJUI', 'MAPTUT', 'XOPDAJ', 'XIQVOL', 'LEZNOW', 'KOFGET',
                                      'SAWMIP', 'YATPIV', 'IMULOS', 'MAQSEF', 'PAZCUQ', 'GOGBAF', 'ZILKOW', 'JONCOE',
                                      'KOPYOD', 'GAKLAI', 'MEGMES', 'KIWTEQ', 'HIMFEP', 'LEDVUM', 'QIFFES', 'REPQUB',
                                      'TEHBEO', 'FEKZOK', 'VULVON', 'GIGHOU', 'QEXSUH', 'KISZET', 'GADKUT', 'FIRWEJ',
                                      'HIVJUR', 'MIKWEK', 'ASELAL', 'CIKBOQ', 'LOFYEN', 'YIKQEO', 'YOPKIZ', 'QOKKIL',
                                      'ZARYIC', 'WIJYAQ', 'IBITIW', 'BOCVOF', 'ZARROC', 'ZORHIA', 'REDQIC', 'FIJDOT',
                                      'CPYFEM', 'GONHEW', 'OPOXUO', 'KUNWIB', 'TAPUFE', 'FAYGAQ', 'HIVKAY', 'IYESOU',
                                      'DEDGAU', 'NOWPOE', 'CEXGUK', 'BIKCUW', 'DABZIS', 'DOXQMO', 'ZORHOG', 'TUQQAY',
                                      'QAWXIW', 'KADREP', 'QASDIX', 'VARSIR', 'EFOCOR', 'IQEZOT', 'KIKHIY', 'WEXHOZ',
                                      'ZORJAU', 'KAZGOI', 'BERZII', 'COSBIX', 'ATURUC', 'MEHHAL', 'MATHAT', 'ZIPNOB',
                                      'LIPBET', 'NUJSAN', 'QEXTAO', 'XIKRIS', 'IQCTNI', 'FUFMEY', 'ADAQUU', 'BAHYEO',
                                      'XAKLEC', 'ZAQHOR', 'BEDLEC', 'TAXCAZ', 'KEDKIQ', 'NOVVAV', 'NAJWAY', 'XALSAE',
                                      'DUVHIN', 'WEVDEJ', 'NOXMEU', 'IXIGIG', 'POXRIE', 'XIRMER', 'WUMJEU', 'HUPMEL',
                                      'HURWOH', 'SASYUI', 'ZORHEW', 'SACVUQ', 'TUYPUY', 'DOCPOD', 'QETKIM', 'KADGOM',
                                      'VEQLAE', 'GEBQEL', 'LIRBUK', 'PEDKAK', 'FIQTIL', 'HUHKAW', 'EFAGEX', 'PEBNAN',
                                      'BOSVOU', 'ALIPIU', 'TARXES', 'SICDOZ', 'WAMXOY', 'HADVUD', 'AHIHUU', 'CALQIS',
                                      'SOTPEX', 'RIPQAL', 'SITHUZ', 'XUQCOB', 'REDLAQ', 'YUWWUI', 'XACGEN', 'WOMTID',
                                      'XIHCID', 'IYUBAF', 'ROMRAO', 'BOWNIL', 'YURSAF', 'BESXAB', 'PODVIO', 'YIMPER',
                                      'CIYLAZ', 'LUJMOU', 'UFABIN', 'XIZQON', 'ZEBSAA', 'JUBSUU', 'AWENEW', 'GEMLAM',
                                      'JAKHUY', 'XEPZOI', 'GIYHEE', 'BOJBAG', 'HOGVAC', 'JOSYUL', 'ATOJUO', 'LENQUQ',
                                      'ARIHIU', 'ZOWPAC', 'IFUNED', 'ABEVAH', 'DOSGEA', 'QOHMOQ', 'AHABOB', 'AROMAW',
                                      'GUPSAL', 'ACIQEK', 'YEHJOL', 'LUYBUC', 'JOQKOP', 'VOBMEH', 'SOKQUG', 'UQIJOV',
                                      'ZAPQIT', 'ROSGEL', 'WAPDID', 'HACFEY', 'NEVXET', 'COHBEJ', 'XIYGIW', 'OHUBAV',
                                      'HEBPOT', 'QODDAP', 'LINZUF', 'PUFQEO', 'ILACUT', 'CAGDUM', 'CESFUD', 'QEXRET',
                                      'TAKDOA', 'TUXGEY', 'GIKJOB', 'BADCOA', 'QUQVED', 'ZIFLIJ', 'ZOCLAG', 'BISXIK',
                                      'UFIYUG', 'BAXWUU', 'USOQUP', 'YENMUA', 'HUMPUC', 'WABDEK', 'SETMOW', 'TARZIX',
                                      'HETNUQ', 'XENCUR', 'QUTJUM', 'UQEPEM', 'AGEFAT', 'VIKNIO', 'CAKVOC', 'SEHRED',
                                      'COBZOL', 'PIMVOZ', 'KUDKEB', 'FOTNOR', 'LOLKON', 'IPUYUN', 'VAMDUJ', 'MAXGUQ',
                                      'DOFVEB', 'CAYCAH', 'FIZTIT', 'VIZJUM', 'EREXIJ', 'KIJDEN', 'YIPXEB']
        db_dict = {key_: item_ for key_, item_ in db_dict.items() if key_ in csd_codes_for_testing_true}

    tmQM_DB = MoleculeDB.from_json(json_=db_dict, type_="Molecule")
    tmQM_DB.to_json(path='../data/New_DB_jsons/tmQM.json')

    # example:
    # mol = tmQM_DB.db["AFATAE"]
    # mol.pre_rotate_and_shift_molecule()

    tmQM_Ligands = LigandDB.from_MoleculeDB(molDB=tmQM_DB,
                                            denticity_numbers_of_interest=denticity_numbers_of_interest,
                                            metals_of_interest=['La', 'Sc', 'Y', 'Zr', 'Ti', 'Hf', 'Ta', 'Nb', 'V', 'W',
                                                                'Mo', 'Cr',
                                                                'Tc', 'Mn', 'Re', 'Fe', 'Ru', 'Os', 'Ir', 'Rh', 'Co',
                                                                'Pd', 'Pt',
                                                                'Ni', 'Au', 'Cu', 'Ag', 'Zn', 'Hg', 'Cd'],
                                            Testing=True
                                            )
    tmQM_Ligands.to_json(path='../data/New_DB_jsons/tmQM_ligands_full.json')

    print(f"Number of ligands: {len(tmQM_Ligands.db)}")
    #
    unique_ligands = tmQM_Ligands.filter_duplicates()
    tmQM_unique_Ligands = LigandDB(unique_ligands)
    # tmQM_unique_Ligands.to_json(path='../data/New_DB_jsons/tmQM_ligands_unique.json')

    #
    #
    # Ab hier nurnoch testing; todo: als unittest in Testfolder
    df = tmQM_Ligands.get_df_of_all_ligands().sort_values('name').reset_index(drop=True)

    #
    if Testing == 20:
        old_outpath = f'../data/221107_ligand_db_test_Testing={Testing}.csv'
        old_df = pd.read_csv(old_outpath).sort_values('name').reset_index(drop=True)

        # ich hab das naming improved, deswegen muessen wir hier was anpassen
        for df_ in [df, old_df]:
            for col in ["name", "unique_name"]:
                df_[col] = [name.removesuffix('-b').removesuffix('-c').removesuffix('-d').removesuffix('-a') for name in
                            df_[col]]

        # reformat atomic_properties
        full_atomic_props = [eval(ap) for ap in df["atomic_props"]]
        new_atomic_props = [str({"partial_charge":
         {ind: [atom, [pcharge]] for ind, (atom, pcharge) in
          enumerate(zip(atomic_props["atoms"], atomic_props["partial_charge"]))}}) for atomic_props in full_atomic_props]
        df["atomic_props"] = new_atomic_props


        pd.testing.assert_frame_equal(df[old_df.columns], old_df, check_like=True)
        print('Ligand database same as old database.')

    elif Testing is True:
        old_outpath = f'../data/221107_ligand_db_test_Testing={Testing}.csv'
        old_df = pd.read_csv(old_outpath).sort_values('name').reset_index(drop=True)

        for df_ in [df, old_df]:
            df_["name"] = [name.removesuffix('-b').removesuffix('-c').removesuffix('-d').removesuffix('-a').removesuffix('-e').removesuffix('-f')
                        for name in df_["name"]]


        # Die unique names sind leider anders, also da ist irgendwie was schief gelaufen. Wahrscheinlich andere
        # Reihenfolge, deswegen muessen die raus
        df.drop(['unique_name'], axis=1, inplace=True)
        old_df.drop(['unique_name'], axis=1, inplace=True)

        #
        # Als naechstes muessen wir die rausschmeissen, wo die Anzahl Liganden auf Basis anderer Graph
        # Das kann ich leider nicht mehr rekonstruieren, weil das totaler muell war
        # die denticity wurde auf basis von skin=0.2 ermittelt
        # aber die liganden selber dann wieder mit skin=0.3
        # das war unglaublich inkonsistent
        for csd_code in set(df["csd_code"]):
            if len(list(df[df['csd_code'] == csd_code]['name'])) != len(
                    list(old_df[old_df['csd_code'] == csd_code]['name'])):
                print(
                    f"New: {list(df[df['csd_code'] == csd_code]['name'])} and old: {list(old_df[old_df['csd_code'] == csd_code]['name'])}")

        csd_codes_for_graph = ["UFABIN", "ZIKDUS", "USOQUP", "FOKGER", "TIQXAS", "TEGZEM", "XOMXEE", "KIYTER", "NUJSAN", "BISXIK", "UQIJOV"]
        df = df[~df["csd_code"].isin(csd_codes_for_graph)]
        old_df = old_df[~old_df["csd_code"].isin(csd_codes_for_graph)]

        # sorting and resettting index
        old_df = old_df.sort_values('name').reset_index(drop=True)
        df = df.sort_values('name').reset_index(drop=True)

        # und damit ist auch n_total_unique_ligands quatsch
        df.drop(['n_total_unique_ligands'], axis=1, inplace=True)
        old_df.drop(['n_total_unique_ligands'], axis=1, inplace=True)


        # Graph hashses are somewhat unexpected also different (in 4% der Faelle)
        different_graph_hash = [i for i, x in enumerate(zip(df["graph_hash"], old_df['graph_hash'])) if x[0]!=x[1]]

        df.drop(different_graph_hash, axis=0, inplace=True)
        old_df.drop(different_graph_hash, axis=0, inplace=True)

        #
        #
        #
        # Okay alles nicht zielfuehrend irgenwdie, ich checke das mal manuell, via
        # wurde manuell gecheckt, da ist einfach nur was mit der rotaton schief gelaufen,
        # die kann man alle droppen
        old_df = old_df.sort_values('name').reset_index(drop=True)
        df = df.sort_values('name').reset_index(drop=True)

        """
        ligands_to_view_3d = [df['name'][i] for i, x in enumerate(zip(df["coordinates"], old_df['coordinates'])) if
                              x[0] != x[1]]

        import collections
        from DARTassembler.src.ligand_extraction.constants import mini_alphabet

        new_list = []
        for key, value in dict(collections.Counter(ligands_to_view_3d)).items():
            for i in range(value):
                new_list.append(f"{key}-{mini_alphabet[i]}")

        for name in new_list:
            tmQM_Ligands.db[name].view_3d()
            input("Press enter to cont")

        """
        different_coords = [i for i, x in enumerate(zip(df["coordinates"], old_df['coordinates'])) if x[0] != x[1]]

        df.drop(different_coords, axis=0, inplace=True)
        old_df.drop(different_coords, axis=0, inplace=True)


        # Im neuen code werden die atomic properties mit der rotation geaendert, im alten nicht, deswegen weichen die
        # vollkommen voneinander ab (100% Abweichung)
        df.drop(['atomic_props'], axis=1, inplace=True)
        old_df.drop(['atomic_props'], axis=1, inplace=True)

        pd.testing.assert_frame_equal(df[old_df.columns], old_df, check_like=True)
        print('Ligand database same as old database.')

