%chk=NOCOQOGE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
C  2.7706927117040063  1.135713848674156  -1.1910118415942672 
H  2.3213247064028555  1.1239483918710282  -2.0621726647709506 
C  3.32459010488357  -0.28364356847200123  -1.0361367621880881 
H  4.0260359651711966  -0.4007348066240439  -1.6951997387186986 
H  3.7408821570953954  -0.3442524737243293  -0.16185615627864292 
C  2.3614755579772675  -1.4510299079422693  -1.1729389696106343 
H  2.903059070044928  -2.26672230190529  -1.1795053908084492 
C  3.8819584238797917  2.149897960394884  -1.3587925451202294 
H  4.320222506616566  2.286107535025887  -0.517104666175937 
H  3.512381844676769  2.9801592545998616  -1.6662133405288393 
H  4.51837975055206  1.8243224506722582  -2.0011250543594845 
C  2.0814750791442247  -1.8315727687941459  1.3369250750953845 
H  2.6114315814780102  -1.0447646646535522  1.5813622917642733 
C  3.027312449600963  -3.031134259652937  1.2683609569567653 
H  3.347053833612211  -3.236167028251573  2.1503439109594766 
H  3.770311502160686  -2.8197324198401583  0.6994351533416366 
H  2.5574728042277286  -3.788975379217665  0.9124702642821196 
C  1.0234515244924718  -2.022570067173828  2.4014655251552894 
H  0.35694300439297644  -2.6376791068284553  2.0846521542727507 
H  0.6110367052504432  -1.1783703832528933  2.59987689211132 
H  1.429208947939967  -2.374527144917631  3.1967787708285136 
C  0.8531133811997773  3.1264288047924547  -0.673103659857415 
C  1.0294143558891165  4.352164451703374  -0.024437250818480867 
H  1.4576055277743727  4.377610415578981  0.8002889360108164 
C  0.5819820753653369  5.514090002432332  -0.5866201916024429 
H  0.7098809681664551  6.319908464740846  -0.14033433139287194 
C  -0.05557418703468575  5.507182386619196  -1.8137041402865914 
H  -0.3688567117491901  6.299653149340039  -2.1852723842757174 
C  -0.22077660806248045  4.328041562091777  -2.4649714402636063 
H  -0.6244489335831378  4.3185042794797805  -3.303517541470435 
C  0.207823163441484  3.132180368912319  -1.8942025356146712 
H  0.059869744068590025  2.3297007313840736  -2.3403977698433587 
C  1.582895953499062  -1.4326527284952413  -2.4676399732273593 
H  1.0625652060403143  -0.6271109547821435  -2.5174009423917907 
H  0.9988206217999234  -2.1934269304088407  -2.502205137819366 
H  2.1935512170461533  -1.4642907015291429  -3.2087535257714688 
C  2.085181140268685  1.9178007876969865  1.6505433127601916 
C  3.3452376668966934  2.434030491653965  1.9067233519018425 
H  3.9479898381697716  2.5481782318529866  1.2076978849601507 
C  3.712905103534382  2.779407360385805  3.1873001651314605 
H  4.568156488016234  3.1052109911274584  3.350723277136699 
C  2.8206630167254336  2.643687007324639  4.2262278403292095 
H  3.059304709059726  2.9112616147801025  5.084275713192313 
C  1.5950926368007134  2.1224122418930844  3.993417718608792 
H  1.000437828140589  2.0034278841048736  4.698924686071983 
C  1.2243634673894066  1.7657464366773374  2.7121808254709445 
H  0.37603281699705726  1.4150181591212199  2.564548705228018 
N  1.439777577440349  -1.5553586491545923  1.523585403891381e-16 
P  1.439777577440349  1.5553586491545923  -1.3496534789007444e-16 
H  0.9684673606864068  -2.271792701587245  -0.00844678035017385 
Br  -1.654629867976521  1.654629867976521  0.0 
C  -1.352040782177733  -1.3199206589279684  -0.4819701001757076 
C  -1.6380700204479186  -1.599217458033982  -1.8139382277462257 
H  -1.1579680913617836  -1.159339363685303  -2.503515938957912 
C  -2.63068216563725  -2.52632649780156  -2.142289956238523 
H  -2.8368222554424674  -2.694399528648736  -3.0538898079700485 
C  -3.3101544432602585  -3.1961573187854975  -1.1619874135498987 
H  -3.9960123560788676  -3.811990145619169  -1.3868663427944008 
C  -2.9974710358883443  -2.9701089012107795  0.12423560353005347 
H  -3.438659168087169  -3.468581175155742  0.8017084178296768 
C  -2.0389854250767776  -2.015940632791841  0.4912562833419403 
H  -1.8619130127668426  -1.8492165108145995  1.4101568590026907 

-Pd 0
lanl2dz
F 1 1.0
1.472 1.0
****
-H 0
6-31g(d,p)
****
-P 0
6-31+g(d)
****
-N 0
6-31+g(d)
****
-C 0
6-31g(d)
****
-Br 0
6-31+g(d)
****

-Pd 0
lanl2dz

