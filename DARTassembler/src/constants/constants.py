# some global coordinates
#
from pathlib import Path

from DARTassembler.src.constants.Periodic_Table import DART_Element

project_path = Path(__file__).parent.resolve().parent

size_full_ligand_db = 41018

#
ATOMIC_PROPERTIES_COLUMN_SEPARATOR = '  ===  '

original_test_complexes = ['DUCVIG', 'XUMHUJ', 'VEVTOH', 'AMAWAM', 'BALTUG', 'JIYZOG', 'LIRBUK', 'YAXYUU', 'ZAVZOO',
                           'COJRIF', 'VIFRIL', 'AFOYUP', 'KAYBOB', 'RATLOQ', 'IZICOJ', 'VIQMUE', 'DOGGIQ', 'PIPXIV',
                           'GAWQAW', 'ETATIC', 'QAFFAG', 'BOXJUU', 'TIDQAY', 'FIVXUE', 'PIBDEK', 'ABESUV', 'BACLUO',
                           'PURROJ', 'OSUCIP', 'WUXLUW', 'QEBSUL', 'AKIHIL', 'JUKRIQ', 'RAKMOH', 'UQAZAP', 'EBUFOX',
                           'MOYPIZ', 'MIKVOS', 'YEXFAL', 'NEKQID', 'PIWNIT', 'DINVIH', 'FOFZEF', 'PADQIU', 'UROPOH',
                           'WEKMAB', 'NEJMEU', 'ADUHIR', 'XACJET', 'WIBHOH', 'GIXMEF', 'QUBHEC', 'REPQUB', 'YIJTER',
                           'DEJHIJ', 'WOFYOG', 'OZUQIK', 'UGIFAR', 'IZUDOX', 'YUFBIL', 'KEYVUF', 'XESDEE', 'ADOXOH',
                           'VIXLOE', 'POBJUL', 'BATCUV', 'CEJZIC', 'KACHED', 'POVLIW', 'VEHLEZ', 'BABVOR', 'LANQUM',
                           'YEVVIF', 'TIBXEH', 'ZUKGOB', 'ZAGHEX', 'NAZBEV', 'ADEZIS', 'ROBQOQ', 'SAKDEN', 'JEBJEF',
                           'RISZUR', 'CEVVII', 'XIKRIS', 'PIMVAI', 'UJULAN', 'NUXPAX', 'XUSYOZ', 'DAXSOK', 'XIZQON',
                           'HAHSIU', 'MIJCAK', 'YEDPII', 'IZOVUP', 'AHUYOT', 'OPALOH', 'LUNBOM', 'WOPHAK', 'NOVVAV',
                           'URONAR', 'VACXUS', 'WOLQAP', 'WIWGEO', 'QAJCEL', 'YIWWOR', 'FUFMEY', 'WUGPAR', 'MIJKEV',
                           'XABBEI', 'CPCBMO', 'JOQKOP', 'RAGMIW', 'PIRJOP', 'RAFXAX', 'FOGMIX', 'QUGFAB', 'NAJWAY',
                           'QIJDUI', 'DOXQMO', 'MADREP', 'VOXBAL', 'FOPQAD', 'SUSXEJ', 'OKOWAM', 'DENPES', 'VETRAP',
                           'EFALOM', 'WAKNOL', 'XULGAO', 'JOPFOM', 'DITGET', 'FICHEH', 'EVIRUX', 'VIFTIN', 'CIJWUN',
                           'ZAZXAZ', 'FUTSUJ', 'YETBII', 'HAJNIP', 'NUFWUI', 'QASDIX', 'XEMPIP', 'HARKIU', 'YARBAV',
                           'ZOWPAC', 'CANGUT', 'SOHZIZ', 'CUGSUS', 'BAHYEO', 'MOMGAW', 'XAMKUT', 'WEDXIM', 'FAJGOO',
                           'MAZVAN', 'MAFXUN', 'LUFMAA', 'APUQAD', 'OYULEA', 'UMIGAY', 'XILXUL', 'SIXCIO', 'DERGAI',
                           'GENXAZ', 'IPUYUN', 'BAKVOY', 'CIHVAT', 'HECZOF', 'VAQFAU', 'EADJIC', 'WOLPIW', 'VALJEX',
                           'CIQSAY', 'YOXWUD', 'YEJFOJ', 'BUKJOH', 'OBORUS', 'RITKEN', 'YOKVOM', 'XELMIK', 'ROLDEE',
                           'XOVPUT', 'MEPKAU', 'SAWMIP', 'ELODUG', 'JUMMIN', 'WAGYIM', 'IHIMIV', 'VEQLAE', 'VOBMEH',
                           'SASYUI', 'IBUMUP', 'LUGBOE', 'ACOJEK', 'EKOFAM', 'LICPET', 'KANNES', 'ONEZOY', 'WEKREJ',
                           'EJUXOX', 'OKOVEP', 'GUPQEN', 'DAJWIX', 'KUCXIQ', 'GUDQEB', 'ACEBIV', 'MEPKOI', 'PITQUE',
                           'BACDER', 'GUPSAL', 'TACVAV', 'QEFCUC', 'GUGBUH', 'TUCSOA', 'PUFQEO', 'ECOSOH', 'WEDGAO',
                           'ZIFLIJ', 'TEBVAZ', 'COLXAD', 'PPTCFE', 'ZOJHOX', 'YAVNAL', 'KENVAA', 'GAVQEB', 'XOSPIG',
                           'IBUVUW', 'SUHMAK', 'KEPLAT', 'OGELOE', 'QOZVUY', 'ODACON', 'BUJSIK', 'XALSAE', 'KAQDIQ',
                           'YOGBAX', 'XIKXEX', 'QOVVON', 'EGETAL', 'BIDHOL', 'DAQQET', 'CEYRIH', 'DOVNUZ', 'VIBYEL',
                           'SIWVAV', 'SECHEN', 'VICLAU', 'AGIGEC', 'TUJHEM', 'WOKWOI', 'UFABIN', 'ZEBSAA', 'QEYTOD',
                           'XUXKUX', 'BOGCEH', 'VIRHEK', 'QULZIG', 'ADAQUU', 'QISPEO', 'ATOJUO', 'BOCWOI', 'ZOMMUM',
                           'FIZTIT', 'XERJIQ', 'KULPEO', 'ALIPIU', 'FEYRUZ', 'MESGEW', 'RAZHOS', 'MISQEM', 'DUNVAL',
                           'UWEPUJ', 'KOFKUO', 'LEXLOQ', 'HEVJIC', 'QUDDEZ', 'IHEQAO', 'MOWVOK', 'GAXZEN', 'FIHWUO',
                           'RESYUK', 'NOBBEN', 'PABZAV', 'WOYDOD', 'GIZZIA', 'IGESUI', 'VIRFEI', 'YAKXOZ', 'YAJSEH',
                           'WORNAS', 'KAQFAJ', 'DOPDET', 'ICETAN', 'QEYBEC', 'ZUNWEK', 'WUVPOT', 'TLSCRU', 'QUGMOW',
                           'VAMDUJ', 'ONEPAA', 'XAXHUC', 'MIGLOF', 'HEVKUP', 'LINTAF', 'CUHVAC', 'IHITAW', 'CPYFEM',
                           'SEVHOT', 'SOWVUW', 'CAYCAH', 'EZISUD', 'IBINUC', 'YIKZUN', 'SINXEV', 'KARRED', 'QOYFOB',
                           'TUMKIW', 'HUBVOR', 'CIKZOM', 'YUTSIQ', 'QOKKIL', 'ODECIJ', 'PAFFIN', 'QAHVOK', 'PIRDID',
                           'ZOWPOS', 'TAYMIS', 'QEXMEL', 'GAKLAI', 'LIRFAU', 'VAHBIR', 'AHUSOL', 'YENMUA', 'HUPMEL',
                           'SIVYIG', 'BASLAI', 'HEBPOT', 'SOTPEX', 'WOWPED', 'QAWXIW', 'QEZVIB', 'WAPDID', 'ZUZLIP',
                           'VUCCEC', 'CATSAU', 'BAZMEX', 'CASPAO', 'YUCMIU', 'FUNKUV', 'VIJGEB', 'GUKCOG', 'BISXIK',
                           'FIBPIR', 'CIYJUP', 'HIDKIP', 'RABYIF', 'KISZET', 'CUDKUI', 'NESCIX', 'UJELOL', 'AGABEQ',
                           'EVUYIE', 'KIKNAV', 'UHOWOE', 'DOVXUH', 'GAPMIT', 'BUCKOA', 'KIJDEN', 'MIHLAT', 'TAPDOE',
                           'PIYRIY', 'RURWAD', 'UTEGEG', 'VICDIV', 'NAPQEA', 'HIBVIA', 'XEZWAB', 'IPAPUK', 'EFOCOR',
                           'QADROC', 'UTOSED', 'KUMFAA', 'XACQAV', 'HAGTUH', 'EYEVIN', 'VECLOH', 'CENMEP', 'QEXSUH',
                           'BEJFIH', 'AWOQIN', 'GEGZEZ', 'WEYSEZ', 'NOGTIM', 'ZULROP', 'DICQUE', 'MAFSIY', 'MYNBCO',
                           'VAHHER', 'XOFZIB', 'ASELAL', 'LEHDOT', 'MIKHAQ', 'VIKNIO', 'MEGMES', 'ARIHIU', 'FOPWAI',
                           'FUQQAL', 'SIZWOP', 'UTOJEU', 'UKEDAR', 'KABZUJ', 'HITJAW', 'UMEGOK', 'HITDET', 'AVAQUK',
                           'ZAPQIT', 'CEMSEU', 'KIYTER', 'YEJCUL', 'TILBUO', 'REGTAY', 'UQIJOV', 'XANQUZ', 'PUCVUF',
                           'NAQDUH', 'CIYLAZ', 'ZIGLUW', 'VUGBUW', 'ZITFEN', 'RAXPEL', 'VULVON', 'YODMAH', 'ATOQOQ',
                           'FACFIA', 'POXQAV', 'METVIT', 'JOMJAW', 'DISDOA', 'LEDVUM', 'UCUQAK', 'VUXFAX', 'MAZVIV',
                           'GEBQEL', 'UQIJEL', 'TILBAU', 'HIPXOV', 'ZIKDUS', 'SEDYEI', 'FOJJUM', 'GAHRUC', 'SULBOQ',
                           'BAFCIW', 'SUQNIB', 'CBYPRH', 'NUYBEQ', 'FOKGER', 'LIWTOA', 'UCUMEM', 'LIPCOE', 'KUGTUD',
                           'HETNUQ', 'GIKJOB', 'TOLKEN', 'JOSYUL', 'JOVFOS', 'JUSQUJ', 'MEFTIB', 'MAPTUT', 'ABEVAH',
                           'MIDCOU', 'FEHNEL', 'OVETUF', 'QIFRAY', 'ZORNUR', 'LOLKON', 'YOZDUN', 'ZUFHUD', 'YUWNUA',
                           'CEFFOK', 'KUDKEB', 'NUXWEK', 'SIZTEC', 'LERKAX', 'PEZTIA', 'TUCWIZ', 'CESFUD', 'LERSIK',
                           'PIYBIL', 'KISKAZ', 'UQEPEM', 'BEJKOQ', 'HAYHEV', 'TUXGEY', 'ETESED', 'QISHIK', 'CIXHEY',
                           'QUQVED', 'LIYWUM', 'VUJCIM', 'LAVQOO', 'UFIYUG', 'LOFYEN', 'VAQJEE', 'JUXZOR', 'PUXYEO',
                           'GUWSOG', 'USOQUP', 'REDQIC', 'NEXBOH', 'HUHKAW', 'HESDOB', 'UQAKII', 'MEHNAO', 'ZOPJAP',
                           'BAZLAS', 'UPAKIH', 'OZUJAV', 'UFUPAO', 'IQCTNI', 'BAVLIT', 'HOLTAD', 'PIJNAY', 'QAZMIQ',
                           'ODIHUG', 'XOMXEE', 'FIJDOT', 'VARSIR', 'YOCFAZ', 'FEVROQ', 'WUTYOA', 'CUMWUC', 'TERCUQ',
                           'ROMRAO', 'XIQVOL', 'CAYWIJ', 'QATFIC', 'WERZAX', 'TOJHEF', 'LIQNEG', 'QAHCUY', 'XUSFOG',
                           'SURZIP', 'NIJCAN', 'GUKKII', 'KADGOM', 'CAKVOC', 'MALTID', 'CEFZAR', 'JUXVIJ', 'TOFMPD',
                           'SILCEV', 'VOQQAW', 'HUMPUC', 'NETTOX', 'GIYHEE', 'KOLTIQ', 'OYOPEY', 'UDIGAQ', 'ADUMAP',
                           'NIXANT', 'HOMWOX', 'MUFJUT', 'HAFMEH', 'BEXXOS', 'SOVBEL', 'TUYPUY', 'RICDOY', 'NILNUS',
                           'MAXGIC', 'LUXZOT', 'EQONII', 'NOXJAN', 'IDUFIY', 'AZUVIA', 'EFAGEX', 'WUMJEU', 'VAWKOT',
                           'SADRUN', 'QIKQAF', 'HEWDAO', 'QUMLEP', 'BADCOA', 'TIPWUK', 'PIJWUB', 'SESXUL', 'KABDAT',
                           'WEKZEU', 'BOLNEV', 'CIKBOQ', 'HIMFEP', 'ODABED', 'AHIHUU', 'BOZPEM', 'LEMBOX', 'NOXMEU',
                           'TEGZEM', 'LEQPOM', 'LUMSOC', 'MEHHAL', 'QEKCER', 'BUVNAI', 'BODRUI', 'WEXHOZ', 'PIBPIB',
                           'GAXDOY', 'WEPXAS', 'QOYMIC', 'UCECIQ', 'QEXRET', 'FAYVUZ', 'CAJNOQ', 'YEFLAW', 'AHABOB',
                           'YOSVOR', 'VUBXEX', 'CAWTEA', 'POXRAW', 'INICAJ', 'LIYXAT', 'OLIQAC', 'NUCDOG', 'CUCNAR',
                           'IMOSAG', 'DOFVEB', 'ZIWFUI', 'GADKUT', 'FAYDEO', 'XANHEB', 'JONCOE', 'RAZLAF', 'REFPOI',
                           'JIKBEN', 'QUCCOI', 'HOFBIN', 'QAMHOD', 'AQEZUS', 'RASYER', 'SIVBOQ', 'ACARUT', 'JIJQOL',
                           'KIHBUB', 'AZIDIY', 'BAHFOI', 'KEQSUU', 'TEHBEO', 'LINZUF', 'IXENEG', 'ZIBLEB', 'HUKSEL',
                           'UTIBOP', 'CUHXUZ', 'NUCPAC', 'XELMOR', 'QEJWUA', 'FITDOD', 'ZEYHIX', 'WEPYEY', 'OLOWUI',
                           'CAWLIV', 'DABZIS', 'MEZWAP', 'KEFCEF', 'QERTAL', 'XESPOB', 'NORKAG', 'RITHUZ', 'SUYXOB',
                           'PUPNUL', 'XEDCUG', 'DUDRIF', 'NEVXET', 'JOCVIJ', 'DARTIC', 'AVUSUF', 'REBQIY', 'IBESAJ',
                           'LIJTUW', 'MAFGIL', 'TULZOQ', 'QEGFOA', 'IBUGIX', 'FIRWEJ', 'MEBNEM', 'AWENEW', 'XENCUR',
                           'TEBDOV', 'MIMBUH', 'PIWJEL', 'ZODZID', 'OWEYOG', 'LOYXEC', 'NETFOJ', 'HUBYIM', 'NEYPUE',
                           'DUGWUX', 'QOZVEH', 'XABMIY', 'WULFOZ', 'BIKCUW', 'KEFJUC', 'YERYEC', 'COBPIT', 'TIYWIK',
                           'SUTPED', 'FOKFAM', 'TIZYUX', 'CMTPCU', 'FAYGAQ', 'ZURLON', 'TAHVEG', 'RISPIT', 'AMOKAQ',
                           'QAFYUU', 'LAQLUL', 'GEMLAM', 'DEFJON', 'UYOBER', 'UCIWIN', 'ROMFII', 'YUWWUI', 'FIRXEK',
                           'QAQVOU', 'PIJWIQ', 'ZAQHOR', 'DEDGAU', 'PATQAF', 'EPADOQ', 'ACIQEK', 'MILJOG', 'YEHJOL',
                           'MITKOQ', 'SUSBOZ', 'FIQTIL', 'OMAXUV', 'XUYVIX', 'ZARYIC', 'JAGWIY', 'EKALUY', 'FOTNOR',
                           'QURGAN', 'HOCLOA', 'TUQRIG', 'MIQVAL', 'DPIMZN', 'LAZGAU', 'PAZCUQ', 'GEQXEG', 'PADTIZ',
                           'IXIGIG', 'PEDKAK', 'NIPXUF', 'ICONOE', 'ESIRII', 'IMULOS', 'POXRIE', 'CODRIX', 'FOMQII',
                           'UGOKUW', 'XEYCOW', 'LIKXIN', 'SUNNAR', 'ZIYXEK', 'AGAMEB', 'PODVIO', 'SEZQAR', 'YUTNUY',
                           'TIQXAS', 'MIPLOM', 'BEZKOH', 'JONXIV', 'DIVLUQ', 'XUQCOB', 'SOKQUG', 'CONZEL', 'QODDAP',
                           'NOHZAM', 'NABMAF', 'MEHZUU', 'TAYQOA', 'DUVHIN', 'RALCUE', 'XAMFUN', 'HACFEY', 'IRETOP',
                           'JUBSUU', 'DUBPAT', 'BESXAB', 'SETMOW', 'SEMPUV', 'KOPYOD', 'GUWYED', 'YIHNEK', 'CEKWEX',
                           'RUQSOL', 'FOLRUW', 'TAMQIH', 'TAWJIM', 'GODCOU', 'NEKWIJ', 'QIFFES', 'XIRMER', 'KOFGET',
                           'UCASOG', 'XOPFUD', 'EBUSOL', 'YOHHIM', 'TUTCEQ', 'KEXNOR', 'MOLMEH', 'BERZII', 'QEZVUP',
                           'XIYGIW', 'CUCSUQ', 'PEYHGP', 'HIWCUN', 'HOLHOH', 'NUJSAN', 'VEFQUU', 'HOXYUQ', 'ISOSIU',
                           'MUHWUH', 'WIWGEQ']

benchmark_complexes = ['ACASOO', 'AFATAE', 'AGUZAD', 'AROMAW', 'AVIBIR', 'AXOTAK', 'AYUWAU', 'BASDAC', 'BAXWAA', 'BAXWUU', 'BEBFUL', 'BEBQUW', 'BEGLUU', 'BIGHOO', 'BIXBUG', 'BOCVOF', 'BOJBAG', 'BOQPAZ', 'BOSVOU', 'BOTBAP', 'BOWNIL', 'BTZANJ', 'BUFREB', 'BUSNOS', 'CAGDUM', 'CAJKEG', 'CALQIS', 'CALQOY', 'CALQUE', 'CALREP', 'CIJXUQ', 'CIKZON', 'CIMGAG', 'COBTOE', 'COBVUM', 'COBZOL', 'COSBIX', 'CUMWUC', 'DAHVIR', 'DELCOO', 'DITDES', 'DOJTOL', 'DOSBUJ', 'DOSGEA', 'DOTCEW', 'DOTNAD', 'DUFLEW', 'DUXROF', 'EBECIY', 'EBECOE', 'EDOBIK', 'EFAKEC', 'EKEGEI', 'EPAVEX', 'ERATIA', 'EREXIJ', 'EVEWAD', 'FAVGOZ', 'FBGLNI', 'FEKZOK', 'FIGRIZ', 'FONZAJ', 'FOQQAE', 'FUWGAF', 'GAPNAN', 'GAQCOT', 'GEDQEM', 'GOLMAW', 'GONHAS', 'GONHEW', 'GULZOD', 'HADVUD', 'HEHSAO', 'HEWBER', 'HIQXOV', 'HIVJUR', 'HIVKAY', 'HOKYAH', 'HOLVEJ', 'HUJVUD', 'HUJWAK', 'HUJWEO', 'HURWOH', 'IBITIW', 'IBOLOB', 'ICOHOX', 'ILEFUB', 'IQETAB', 'IQEZOT', 'ITOSER', 'IVAGUJ', 'IWISAI', 'IWOKIO', 'IYESOU', 'JAPCOT', 'JEGROF', 'JIMYOU', 'JIXCID', 'JIYFIJ', 'JOMTAI', 'JOSPIQ', 'JOWHUA', 'KADREP', 'KADRUF', 'KEBBEY', 'KECJUA', 'KEDKIQ', 'KEFFOQ', 'KEMQUQ', 'KEYVOC', 'KIHMAR', 'KIKHIY', 'KIMQOM', 'KIWTEQ', 'KOGBEO', 'KUBGAQ', 'KUNWIB', 'LAMNUJ', 'LILHEW', 'LOGYEN', 'LOPCAW', 'MATGEW', 'MATHAT', 'MAXGUQ', 'MBSLCO', 'MEDCIK', 'MIGBIO', 'MIKWEK', 'MIKWIO', 'MUJRIS', 'MUJTIW', 'MUTCIP', 'NESJAW', 'NOWPOE', 'NOXREZ', 'NOXRID', 'NUYGUL', 'OBEXUP', 'OHUBAV', 'ONAXUY', 'OPAQIG', 'OPOXUO', 'PAKKIW', 'PEBNAN', 'PFBGNI', 'PIMVOZ', 'POGGOI', 'POHJII', 'QAPLAX', 'QEXTAO', 'QUBPEI', 'QUQXUX', 'SAMTUX', 'SEBPUL', 'SEHRED', 'SUYNIL', 'TACFOU', 'TAXCAZ', 'TETWAR', 'TOPRIZ', 'TOQJEQ', 'TUFDAZ', 'UBEVOP', 'USONEX', 'UWULEE', 'VAFQEY', 'VASZEW', 'VAVRAN', 'VAVRER', 'VEGFET', 'VEGFIX', 'VENLOR', 'WABDEK', 'WATVUL', 'WEVDEJ', 'WIMKEK', 'WIQGEL', 'WIZQON', 'XACGAJ', 'XACGEN', 'XAKLEC', 'XARTIT', 'XAVCOO', 'XEPZOI', 'XEPZUO', 'XIFXOC', 'XIYJUL', 'XOBLEI', 'YAFSAC', 'YEHXUF', 'YEKQUC', 'YIKQEO', 'YIQGAI', 'YOBMOU', 'YOLSIC', 'YUMMOJ', 'ZAQZEZ', 'ZAZNOF', 'ZEWKAQ', 'ZIPNOB', 'ZOHQOE', 'ZOHQUK', 'AGEFAT', 'AKOPEW', 'AMUYOY', 'ATURUC', 'AWIPED', 'AYEDAK', 'BEDLEC', 'BEMGEH', 'BOLPIE', 'BOTPIK', 'BUHQEB', 'BUHQIF', 'CAGROU', 'CARVAR', 'CEBHUP', 'CENMIT', 'CENNAM', 'CEXGUK', 'CIRVEE', 'CIVLEZ', 'COHBEJ', 'COKRIF', 'DAWXEG', 'DECKIG', 'DIXYAL', 'DOCPOD', 'DOMQIH', 'EBAZAL', 'ECIJIJ', 'EGUKIA', 'EMOYUC', 'ENURUC', 'EXEMIF', 'FAGZUL', 'FAPKAJ', 'FENFEM', 'FEVDOC', 'FOTQEL', 'GAFYUJ', 'GAQFEJ', 'GEKVIE', 'GIGHOU', 'GOGBAF', 'GOLVAF', 'GOVBIE', 'GUVMUH', 'HAZQAD', 'HEBBEW', 'HEQBEN', 'HEQBIR', 'HIGPUI', 'HILDEM', 'HOBPEU', 'HOBQAR', 'HOGVAC', 'HOQBAS', 'IBEWET', 'IFEGEH', 'IFISAR', 'IFUNED', 'IGERIX', 'ILACUT', 'IXENEG', 'IXOGOR', 'IYUBAF', 'JAKHUY', 'JATKEX', 'JAXCUH', 'JIZSIW', 'JIZTIW', 'JIZZIC', 'JOGPIG', 'JOHCEP', 'JOLQUZ', 'JOWYEC', 'KAXLAX', 'KAZGOI', 'KIFYIH', 'KIWQOW', 'KIZPIU', 'KIZZID', 'KOQZEX', 'KOXKUF', 'LANJUI', 'LARCUF', 'LAWNEE', 'LEHMUH', 'LEJMUJ', 'LENQUQ', 'LEPDUI', 'LEZNIQ', 'LEZNOW', 'LILWEJ', 'LIPBET', 'LIXMOX', 'LIXSAO', 'LOPQEP', 'LUJMOU', 'LUYBUC', 'MAQSEF', 'MEGZOP', 'MERCER', 'MEXGED', 'MIJVAE', 'MOJBUL', 'MOLTAI', 'MONFOL', 'MONNIN', 'MONSIT', 'MONSOZ', 'MONSUF', 'MONTIU', 'MUGHOL', 'MUTZAD', 'NAJVUS', 'NAXNOR', 'NINNOO', 'NOFVAF', 'NURSID', 'NUTXOR', 'OGEMEU', 'OKUNOZ', 'OQOBOM', 'PAXVEQ', 'PEJCEP', 'PELCER', 'PINYAM', 'PODVEK', 'PUHJUZ', 'PUHKEK', 'PUXZUE', 'QAKREZ', 'QANLAV', 'QEBZEE', 'QEPSEK', 'QETKIM', 'QEWGOR', 'QIDJET', 'QIDJIX', 'QIMNAD', 'QIZCAF', 'QOHBOF', 'QOHMOQ', 'QOVLIY', 'QUQLAP', 'QUTJUM', 'RAQWOX', 'REDLAQ', 'REXZOL', 'RIBQIC', 'RIPQAL', 'ROSGEL', 'SACVUQ', 'SELSEI', 'SICDOZ', 'SIGCUH', 'SISQUI', 'SITHUZ', 'SOBKAX', 'SUBKUX', 'SUSTUV', 'SUYHUR', 'TAKDOA', 'TAPUFE', 'TARXES', 'TARZIX', 'TAWMEL', 'TIPCUR', 'TUQQAY', 'UCIGUJ', 'UMACOB', 'USADAV', 'UWUNUX', 'VAPTIR', 'VAWGIL', 'VETTUM', 'VIBBEP', 'VIBBIT', 'VIZJUM', 'VOHJUA', 'VOYREI', 'WACWUV', 'WAHNOK', 'WAMXOY', 'WAVPAL', 'WAWHUY', 'WELVAN', 'WEYWEC', 'WIJVIX', 'WIJYAQ', 'WIJYEU', 'WINVUL', 'WIWSAZ', 'WOMSUO', 'WOMTID', 'WOQFEQ', 'WUSWOX', 'XENFEC', 'XIHCID', 'XOPDAJ', 'XUNFUH', 'XUSNIJ', 'YAKFEX', 'YATNUF', 'YATPIV', 'YAWBEE', 'YIMPER', 'YIPXEB', 'YODFUV', 'YOPKIZ', 'YOYNAC', 'YUQZER', 'YURBAQ', 'YURSAF', 'ZARROC', 'ZATNOZ', 'ZILKOW', 'ZOCLAG', 'ZORGUL', 'ZORHEW', 'ZORHIA', 'ZORHOG', 'ZORJAU', 'ZOTVIP']

specified_test_molecules = benchmark_complexes + original_test_complexes

# list of metal numbers in pse
metals_in_pse = [el for a in [[21, 31], [39, 49], [57, 81], [89, 113]] for el in range(a[0], a[1])]
transition_metal_symbols = [DART_Element(i).symbol for i in metals_in_pse]

# Extend mini-alphabet in excel style because some complexes in the CSD have more than 24 ligands.
mini_alphabet_0 = ["a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u",
                 "v", "w", "x", "y", "z"]
mini_alphabet = []
for i in range(1, 1001):
    mini_alphabet.extend([letter*i for letter in mini_alphabet_0])

# warnings
unconfident_charge_warning = 'Charge assignment not confident.'
similar_molecule_with_diff_n_hydrogens_warning = 'Similar molecule found with different number of hydrogens.'
odd_n_electrons_warning = 'Odd electron count.'