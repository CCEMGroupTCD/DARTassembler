%chk=ANORALEB_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.6229630664312726  1.3631180744895142  -9.375600042134922e-17 
N  1.6229630664312726  -1.3631180744895137  1.1102230246251565e-16 
N  3.8440382105234967  -1.7002164804160018  0.7395996748360963 
N  2.966331967732391  0.44970840125465283  0.4942839506377417 
C  2.801832518414205  -0.9280171789034042  0.41067691174652604 
C  1.4404878345575234  -2.832773233122955  -0.04663672860237089 
H  0.7677650208551858  -3.0522235548509737  -0.7110913788366093 
H  1.1248841471938182  -3.145511243848585  0.8160109742113664 
C  2.7336954733978933  -3.536479953730584  -0.39309481295977683 
H  2.598589489777943  -4.496852715635496  -0.3857851357290193 
H  3.0234540962951018  -3.2767305221689798  -1.2819749979454507 
C  3.7750029810917836  -3.162711481626064  0.6076080518712732 
H  3.5592948640915214  -3.558445600405186  1.4669967564132538 
H  4.636975172725829  -3.5061525185848925  0.3248012538937951 
C  5.066031599214656  -1.1985520921302757  1.3952596625744413 
H  5.803989517867389  -1.2382267182775757  0.7659487264404479 
H  5.283912826991891  -1.7725965858862853  2.147122089997638 
C  4.907584970346528  0.2074114196643464  1.878958628153729 
H  5.775725053352601  0.5769904931343128  2.106698555990975 
H  4.354789478473095  0.21978699835830318  2.6757845120307686 
C  4.266173509422917  1.0392241813171679  0.7992355722975604 
H  4.155056187338361  1.9531988899766302  1.102911110765203 
H  4.824954115227315  1.0444695156237196  0.005922783389609831 
C  2.181558683790535  2.2779776624727655  -1.4616569332583444 
C  2.8649453128021047  3.4914020612501315  -1.4161906409089329 
H  3.0263337615549872  3.899215507791204  -0.5964073698979656 
C  3.3129363630899027  4.102803458236311  -2.593787129790189 
H  3.752854163453952  4.921236417844597  -2.558380780787727 
C  3.093805418603514  3.47980718895272  -3.8103250567070726 
H  3.3942090996660346  3.875745248245288  -4.596204851337102 
C  2.430953853875454  2.271870946440027  -3.856292623119703 
H  2.2974866722115  1.852930697838406  -4.674461236255351 
C  1.96474122846916  1.6753099623121215  -2.700350379697417 
H  1.5042762051242171  0.8682819318078852  -2.748257897441449 
C  1.4170038777093303  2.5379197112440703  1.3631660452270173 
C  1.0022199163355068  3.857203989440936  1.1536658096910686 
H  0.8490883732087168  4.162614326073665  0.28829801691304513 
C  0.8166644454492645  4.713868699456842  2.23169140856968 
H  0.5736850601347652  5.599839762860338  2.089913366814724 
C  0.9961809186050945  4.238843720243288  3.519910863860175 
H  0.8680049835339193  4.805487728899154  4.246399278855923 
C  1.3664272263172013  2.9179588910444  3.729335980829866 
H  1.4601517995693942  2.5977749112644757  4.595808942589782 
C  1.5962387299616734  2.0786500252820157  2.6662364743652343 
H  1.8719052396305935  1.2038759266064036  2.8158802571736996 
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


