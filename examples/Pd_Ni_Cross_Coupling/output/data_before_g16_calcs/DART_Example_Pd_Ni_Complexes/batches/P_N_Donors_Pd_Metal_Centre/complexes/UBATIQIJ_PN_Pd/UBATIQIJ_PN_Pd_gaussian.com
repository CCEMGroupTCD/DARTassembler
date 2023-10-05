%chk=UBATIQIJ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.5978050538160153  1.392522534826635  0.0 
N  1.5978050538160153  -1.392522534826635  0.0 
C  2.8128329248012416  -1.0064872308687547  0.06661851715014132 
C  3.1050520400363695  0.39486579660159915  0.24425334225735823 
C  4.293556377527445  0.8973306418740017  0.6299750152682327 
C  5.565742059950288  0.21068785076180063  0.7925922248494587 
C  6.041499610941235  -0.7517503909768517  -0.09102418712699303 
C  7.286679500205183  -1.3177095336761946  0.10230644202323566 
C  8.050391773081483  -0.9627639743488303  1.176372470799872 
C  7.60832473115747  -0.02358411384997039  2.0497625784138345 
C  6.384872917950409  0.5758302648221241  1.8604329362508232 
C  1.6611082481527026  2.79649786815465  1.129863035501775 
C  1.8318183232688374  2.5488150469621442  2.4836679218646536 
C  1.8508354085610774  3.598723068868602  3.3854740795301757 
C  1.6516548781272722  4.888670272648415  2.9336434721597504 
C  1.4707793483762615  5.134377283804682  1.6142493356690617 
C  1.4796258730851017  4.099562611897489  0.6930669032452768 
C  1.698321883993131  2.015500560893553  -1.70236349699145 
C  2.8890135706492  2.0734494399830004  -2.3922215148041377 
C  2.904960982630037  2.499341707553259  -3.7063304853522405 
C  1.761721187893051  2.8753695109849566  -4.324724819718307 
C  0.5777790056460435  2.8383057804461718  -3.6551146353761417 
C  0.5274796107116562  2.39442985417475  -2.342860196272983 
H  1.4436497976170601  -2.2275207132545045  -0.13299262764670647 
H  3.504843906601389  -1.623839200764901  0.0003795079987968242 
H  4.302545988253353  1.8092442578056058  0.8165859312874858 
H  5.520428797385162  -1.014241827099648  -0.8152706805408042 
H  7.608246086258497  -1.9436181909522574  -0.5060263069562317 
H  8.87744415144635  -1.3667493002093112  1.310473738245138 
H  8.13517046556055  0.21573170450516077  2.7787123717647133 
H  6.099406388525885  1.232652755751415  2.4532834114869786 
H  1.9343940886177677  1.6737514917775922  2.7870239944111828 
H  1.9964544213416966  3.4377916552682475  4.289030487970401 
H  1.6408639850528917  5.593097608987446  3.54115421646388 
H  1.3386803398007965  6.007474587146609  1.3230539359278668 
H  1.3657877041309028  4.278009971839809  -0.2128705951610917 
H  3.6822655120720404  1.825151540773028  -1.9744030366622336 
H  3.7097389318914873  2.530489935977008  -4.170570391272178 
H  1.7847252588861857  3.159160659865142  -5.210252332862319 
H  -0.20201850648262276  3.1121710935468165  -4.079917072231472 
H  -0.2863102749577293  2.351008273203691  -1.8953211831264194 
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


