%chk=EBADAVIF_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
C  2.681928659072696  -0.32301282612470494  -0.027599480467513027 
C  3.644343831787773  -2.308472604503092  -0.0341235719352443 
H  4.2957169861083315  -3.0007420385814925  -0.045099202726564015 
C  2.2933267863450633  -2.4729961354273087  -0.005370125263106279 
H  1.8332399544140285  -3.304016018546009  0.008152853145897294 
C  5.21523000812663  -0.3555509704444768  0.01560380697607862 
C  5.827168971204836  -0.21546947941303946  1.2626244448715525 
C  7.060345134768337  0.4322638202558468  1.3098903111906384 
H  7.504255466261439  0.5443454430229584  2.141882701439247 
C  7.641858563486409  0.9133020671900387  0.14793170727778493 
H  8.458932672132047  1.3955698813267365  0.19710361418394579 
C  7.045814798519288  0.6987869941317778  -1.0798627884235903 
H  7.4758173368543375  1.0122579337001765  -1.8665401410094586 
C  5.824688703707608  0.03247801625824461  -1.1895172756788068 
C  5.229588278052019  -0.7997217680465024  2.535628285556118 
H  4.252126667512426  -0.9452969222419562  2.38025472151641 
C  5.868623201554005  -2.1533393259323916  2.828055632991086 
H  5.70508381469522  -2.760805940014832  2.0761182474854403 
H  5.476963352490026  -2.529967873476815  3.6429892191229682 
H  6.833544396971256  -2.0389803203613277  2.952292561797982 
C  5.383535163937372  0.1280968025408633  3.7322585402294055 
H  5.055773570103199  1.0211216924462585  3.497753713804143 
H  6.329726441448571  0.1831082010501795  3.9832042363386164 
H  4.865944321179935  -0.22235558740826356  4.486912912689164 
C  5.237599174195919  -0.3300947586316608  -2.5353122739034375 
H  4.241172864905149  -0.3241673002114183  -2.448328336363523 
C  5.619335596253947  0.653501116746551  -3.6385056014817576 
H  5.510480828296955  1.5701548264447438  -3.3101398206982355 
H  5.040519149343713  0.5130292736227443  -4.416266954388713 
H  6.55373611117907  0.5082200517456206  -3.896720240609522 
C  5.66424058462989  -1.7368607178178528  -2.9453638108812696 
H  5.4626103687023235  -2.3634000851263575  -2.2203518611606596 
H  6.626471093949348  -1.7458853725463723  -3.1271315267767013 
H  5.176659773476626  -2.003547586952959  -3.7533212889864505 
C  1.816212826496467  2.0656304780767107  -1.6685786066359956 
C  3.159136371374214  2.748333567170052  -1.9129833012624513 
H  3.2321663768705875  2.9961318110445756  -2.858077502191608 
H  3.885180754804197  2.133477627773097  -1.6784115061163485 
H  3.222032904799362  3.5539045045122415  -1.3577836425875436 
C  1.5541441342885371  0.9761302305021974  -2.709956646110389 
H  1.6472551716511101  1.354559945475514  -3.6097850661697173 
H  0.645101938125286  0.6265875732705483  -2.5956365850092777 
H  2.201543897486993  0.24966525472277037  -2.5929434606294772 
C  0.6644247090955899  3.088991694619461  -1.768493535902325 
H  0.5876603580401993  3.404892627920657  -2.6938688515156803 
H  0.8502021841586777  3.848677165610406  -1.1779371500047084 
H  -0.17533089502451626  2.6625610374351987  -1.4987711787231508 
C  2.21372011077362  2.2553215086471963  1.4787868522626006 
C  1.9720683949237863  1.3532975804162237  2.694861278901708 
H  2.170457989626726  1.8512433903674848  3.51517604718618 
H  2.5560667980767504  0.5683035906752704  2.6399274897972047 
H  1.036547814903916  1.0635119013095562  2.7061286896064605 
C  1.2960929467394215  3.471796586494025  1.5691421111252346 
H  1.3822833539598987  3.8815011422652637  2.4553450679491657 
H  0.36739330285646243  3.191649199826148  1.4275745188068567 
H  1.5491648758031715  4.123255203921236  0.8821977476577084 
C  3.6749983072354744  2.7103043938258113  1.433569745546237 
H  3.9245489954527724  3.0916763298973455  2.3010574344429338 
H  3.784539919916742  3.3896012827087927  0.7348593750456232 
H  4.25063379566981  1.9424584269874623  1.2356956298437376 
N  1.7178607852500734  -1.2413920905580156  2.015979405059486e-16 
N  3.893158989873213  -0.9474726288704609  -0.04532461808704213 
P  1.7178607852500734  1.2413920905580156  -2.0753783625012933e-16 
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
-P 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-N 0
6-31+g(d)
****
-Br 0
6-31+g(d)
****
-C 0
6-31g(d)
****

-Pd 0
lanl2dz
