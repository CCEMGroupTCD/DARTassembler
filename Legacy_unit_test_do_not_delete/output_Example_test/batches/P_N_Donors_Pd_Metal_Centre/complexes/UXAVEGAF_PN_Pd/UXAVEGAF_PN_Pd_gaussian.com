%chk=UXAVEGAF_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5297833784559176  -1.4669229069722785  6.60006832410065e-17 
N  1.5297833784559174  1.4669229069722785  4.2398360663078176e-17 
C  3.00618642436987  -1.0417979075357635  -0.9562953993772892 
C  3.75128905929369  -2.0220192969255013  -1.6075960835566114 
H  3.459297634870894  -2.925797609813389  -1.5810280723952195 
C  4.921170228607021  -1.6959510572609904  -2.2979096085351003 
H  5.413854298453439  -2.374203684083459  -2.744278415160284 
C  5.362472132652235  -0.384813743324187  -2.331013286407079 
H  6.161528841398197  -0.16635782210510575  -2.7955418287943026 
C  4.642416038413719  0.6114981371798507  -1.688711318547336 
H  4.958405124225074  1.5072900590008886  -1.7016947868158732 
C  3.4500832354280053  0.3006987639932849  -1.0214335453890122 
C  2.770516738847412  1.4139022747152428  -0.34374896504098196 
H  3.3013812415712334  2.174226641250667  -0.13750712301977525 
C  1.0110860464940135  2.624729139956611  0.805909067251545 
H  0.28094781979108485  3.0445839680956057  0.2668614689863303 
C  0.36021645435627825  2.0661756494455488  2.0640449948318698 
H  0.01317753974712188  2.8042553261987737  2.6069602299383825 
H  -0.3764548706995954  1.4703176729659997  1.8132352448098181 
H  1.024952573528659  1.5646111831038254  2.58024769530041 
C  2.033490325180043  3.711470142106717  1.0733679709550847 
C  2.7587309872718517  3.7858612259744002  2.2571322374751 
H  2.682759106020831  3.0944889458228255  2.9039186732459927 
C  3.595238220454525  4.871755769190571  2.494649331273451 
H  4.094835607388948  4.91232701169904  3.3015641441562154 
C  3.7055343025759484  5.889463215099459  1.569048192115617 
H  4.245876268654529  6.647916852786793  1.756642566160268 
C  3.0257993711534508  5.8003328154182245  0.36189373436932637 
H  3.135463032411361  6.471967681379601  -0.3007529801208052 
C  2.1803219344312135  4.719033397799969  0.13246615657616448 
H  1.695798882460353  4.670110898424646  -0.683369470963929 
C  2.0387279461131387  -1.506375948249991  1.7353755197611345 
C  3.3427744520052043  -1.7838634958677162  2.1136642473064127 
H  4.009951528518519  -1.9243701741807662  1.4525517667672734 
C  3.673688116481137  -1.8560215073739057  3.4637151085848434 
H  4.565843412773813  -2.047855313814913  3.727670953383391 
C  2.683963732772849  -1.6443903415338692  4.43122281084797 
H  2.9141867990082986  -1.6721777015963906  5.3522574151891655 
C  1.3831059055837447  -1.400986084702839  4.064139957533595 
H  0.7107533167068414  -1.283404609098233  4.725081696891323 
C  1.0606871280351835  -1.3273058186941216  2.704581223868513 
H  0.16473082353999802  -1.1528624967010355  2.4425730410814666 
C  1.1769076925811128  -3.1868178741913558  -0.44117513761870314 
C  1.3041146177659826  -4.191728369482885  0.4743965660932254 
H  1.5748617739791062  -3.9959328551890807  1.363478448703002 
C  1.0340320807251226  -5.513733817610126  0.09439679200450733 
H  1.119232655273062  -6.211334908394194  0.7340186335902931 
C  0.6563606081429849  -5.808615990997696  -1.163899668538828 
H  0.48494505562772217  -6.709763613491658  -1.4106753578385873 
C  0.5206019844049374  -4.816466452611966  -2.075740898996067 
H  0.2533936182516512  -5.030213588545077  -2.9622486321479538 
C  0.7684494935145791  -3.482863769265927  -1.734579477964315 
H  0.6584162631529437  -2.791408471800069  -2.377157403061798 
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
