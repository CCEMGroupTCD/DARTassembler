%chk=XOCITOLE_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5244288619020576  -1.47248655172127  -1.7566371353045757e-15 
N  1.524428861902057  1.47248655172127  -2.1787290385605793e-15 
H  1.3060007498804436  1.7692414385412607  0.8874016548480946 
H  1.568051734607827  2.1240450320880284  -0.5600179565327051 
C  2.87596534402708  0.9243858413571799  0.08956991003420148 
C  3.983285628894043  1.765678489192847  0.13573837784348558 
H  3.8720144385692796  2.7062175884335695  0.06735895840192652 
C  5.246026397359042  1.2218867621553495  0.28010887088829145 
H  6.003557206213416  1.7938538942255806  0.3293512659134441 
C  5.415649156524581  -0.15568961885051852  0.3544589486814227 
H  6.287018907702117  -0.5226302468436209  0.4573744589841722 
C  4.320731857006268  -0.9863313612625224  0.2799015011775919 
H  4.441810683809998  -1.9274907989361703  0.3107637852417708 
C  3.03239798238892  -0.4586189321801222  0.1588157476577799 
C  2.1209095536791547  -2.699049230944899  -1.2367754033262865 
C  2.2811713368273385  -4.058855925042989  -0.9761622310801675 
H  2.052875945922886  -4.406920659081399  -0.12144426273588772 
C  2.7745455717708714  -4.907803681841674  -1.9624371575182253 
H  2.8768042712340405  -5.833946741284802  -1.7771042898734817 
C  3.116398031406953  -4.421353986146533  -3.1989179901229767 
H  3.452581853304123  -5.0080344783997015  -3.866514558346785 
C  2.9699300259329364  -3.076473190675028  -3.4683430434205342 
H  3.208532710209143  -2.737001601331359  -4.323227046712601 
C  2.4739643553331527  -2.218429693759503  -2.497677472285192 
H  2.3748689921879667  -1.2948135843898079  -2.6931896241221436 
C  1.4696358741025022  -2.4633953062320484  1.5422349946496292 
C  2.321987654055781  -2.2286777601899534  2.6201286474700383 
H  2.99201935921611  -1.5573790791632842  2.5519300776036866 
C  2.201641011111744  -2.9634076792750053  3.7877283226727188 
H  2.7960262669298555  -2.8026784354504413  4.510248251702817 
C  1.230426168630182  -3.9200894794267898  3.904624357686243 
H  1.1411961050655446  -4.408541329230945  4.71461635528817 
C  0.3834908188378785  -4.174308704425771  2.8508791375722153 
H  -0.2837992809347496  -4.845076772796531  2.9311769656847404 
C  0.5028420769126583  -3.4488260922068386  1.6680641239018363 
H  -0.08122545393009051  -3.6316263084026006  0.9426344302897971 
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


