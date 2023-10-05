%chk=HOLELEPU_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.6137414035402298  -1.3740228100362826  1.3660995767286241e-15 
N  1.6137414035402322  1.3740228100362821  -1.6826926362663845e-16 
N  1.5743851169583107  2.681785422763614  0.33439324374015494 
C  2.797061358535898  3.1122010284436863  0.6493214651951289 
C  3.656876227895375  2.0690428075389313  0.536929932652986 
C  2.8886107527479075  0.9864159296338109  0.13154386612647062 
C  3.1970312609427136  -0.45177024570909996  -0.16322023477014774 
C  1.751406515261673  -2.893107045577382  -0.9586592420381995 
C  1.1549073583005538  -4.055754749335529  -0.49578055092076223 
C  1.2242160477386421  -5.213019038097965  -1.2350347436877493 
C  1.8570944112886738  -5.217650646069474  -2.4415492339272253 
C  2.4321414175339795  -4.076431261975992  -2.9169180149568725 
C  2.3825861458334416  -2.9020300524666762  -2.1787493796507205 
C  1.4979296310579808  -1.881272813390736  1.72407140149458 
C  0.39691170173480206  -1.5063211122713083  2.48744922481443 
C  0.32946305530693554  -1.8604686755036477  3.8211142215627425 
C  1.3336768744445562  -2.603320946619181  4.394363075684176 
C  2.4214087562146167  -2.9762900605645806  3.635674632761931 
C  2.506817518288458  -2.619606062001951  2.320009934394446 
H  0.8644815871576262  3.1666082729624474  0.34264686728481486 
H  3.021142564257852  3.97830586040683  0.9018548754447353 
H  4.572280968055736  2.0778575811449684  0.6988009660660659 
H  3.852798362296557  -0.7934436091866711  0.46498431907508575 
H  3.5509183732555645  -0.5439770664944518  -1.0615246118826807 
H  0.7066597334774191  -4.054254439453523  0.31877615423284217 
H  0.8374566378781586  -5.993904735309695  -0.9106912424818059 
H  1.897664140761165  -6.000860744285811  -2.9413095755936025 
H  2.861007464031169  -4.084586082928131  -3.7421054258261397 
H  2.7737089664312897  -2.1256754117661822  -2.508220275285518 
H  -0.2934761270247068  -1.0171118028696322  2.100258280323902 
H  -0.39853514849750327  -1.593575897362419  4.334009922688818 
H  1.2786482133162729  -2.8530817529852266  5.2882123281075755 
H  3.103397044365353  -3.474325784572634  4.023070815547746 
H  3.2492731906606585  -2.8743480265773114  1.8204363701953008 
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
-C 0
6-31g(d)
****
-N 0
6-31+g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****
-P 0
6-31+g(d)
****

-Pd 0
lanl2dz


