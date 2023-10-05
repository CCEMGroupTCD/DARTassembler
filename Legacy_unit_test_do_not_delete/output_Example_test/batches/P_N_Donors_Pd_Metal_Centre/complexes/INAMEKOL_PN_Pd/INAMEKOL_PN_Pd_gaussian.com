%chk=INAMEKOL_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
N  1.5214811706031723  -1.4755321235066357  -2.74016299367768e-15 
N  1.3582936450664143  -2.646783984084804  0.6664343933211915 
H  0.6270009572296169  -2.7378648781290558  1.1792448704533847 
C  2.4933905277476462  -3.364242813435215  0.6318743714098207 
H  2.6317733087396897  -4.211348416957767  1.0390907787482393 
C  3.4054367677502357  -2.6638941774840497  -0.08536384858302717 
H  4.294359842033541  -2.926332547229416  -0.28944108908792243 
C  2.772000381322141  -1.4754446330359487  -0.46494944452947906 
C  3.3485708331348985  -0.357992089148593  -1.2014903980793217 
C  4.4588431240789  -0.6283054314647933  -2.018860814064535 
H  4.734303062992589  -1.52907114427174  -2.141999283622437 
C  5.156887605314941  0.3762490823857042  -2.6461532665287066 
H  5.903564691183478  0.16508347979650445  -3.191631549984624 
C  4.770821948105667  1.6927601695069228  -2.4807259235333854 
H  5.2655199463162585  2.3914230949538067  -2.89362725214889 
C  3.649437193479595  1.9868057874151859  -1.7039892650182598 
H  3.375168313582119  2.891480839034959  -1.6062348150094108 
C  2.9248409500963057  0.9785153196588533  -1.0693220871159304 
P  1.5214811706031728  1.4755321235066352  4.134403571588152e-17 
C  1.0484374631665383  3.0911292304911244  -0.6724965864284912 
C  0.5357555292133169  3.1067527398900343  -1.9651589555535602 
H  0.39946984270497765  2.291209482085393  -2.4322390089730397 
C  0.22125795994149633  4.320378943936517  -2.5759888956955455 
H  -0.1305624077864731  4.3323549658193885  -3.4588270313419627 
C  0.4219500642790672  5.506760176925056  -1.891999415316526 
H  0.2219225380669625  6.334739094182773  -2.3127950442555685 
C  0.9078502821749143  5.4922620411954455  -0.6078671821297942 
H  1.0328626227534854  6.310561053981777  -0.1421609867920192 
C  1.220315469733586  4.285335067780764  0.01934231031113836 
H  1.5474467153801243  4.279936164650774  0.9110488595563235 
C  2.2555857981406318  1.7061325046026223  1.6422094568192414 
C  3.3063988255479386  2.607663029630205  1.842119821419435 
H  3.6396986565513654  3.1217128750496417  1.115990329470959 
C  3.8587190069201878  2.7468551805349954  3.108701150184189 
H  4.568715293373058  3.3613699346056287  3.2467200266201766 
C  3.3840771739892244  2.001176667358895  4.170698977592261 
H  3.7604262266166257  2.11574340047634  5.0362358192398355 
C  2.362791902130323  1.0873502629594123  3.97512600884043 
H  2.0456804924629712  0.5644055789994159  4.702420344977328 
C  1.8019655123728122  0.9355001836351374  2.710838425273769 
H  1.1046088117675281  0.3026863023751657  2.576955103662711 
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
