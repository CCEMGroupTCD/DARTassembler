%chk=TEWIJAKI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.6137414035402338  -1.3740228100362821  -2.1778355799279935e-15 
N  1.6137414035402322  1.3740228100362821  -1.6826926362663845e-16 
N  1.6447694109503859  2.681785422763614  0.3352685811535827 
C  2.906204502132665  3.1122010284436863  0.3891061806640305 
C  3.7238628401074534  2.0690428075389313  0.10040510948635614 
C  2.8881013063398457  0.9864159296338109  -0.1363909249226496 
C  3.1284971838379043  -0.4517702457090998  -0.4888379523387874 
C  1.5490817383450757  -2.893107045577381  -0.9663324236655733 
C  1.06185541043756  -4.055754749335527  -0.38954959422309315 
C  0.975949949559129  -5.213019038097963  -1.1270593960199606 
C  1.3441499347601005  -5.2176506460694725  -2.4387914605459344 
C  1.8077960572022589  -4.076431261975991  -3.0233312885419887 
C  1.9127975761598073  -2.9020300524666753  -2.2909902886521403 
C  1.8589149963252716  -1.8812728133907355  1.7104749262990786 
C  0.9406721243721178  -1.5063211122713078  2.686085611939031 
C  1.1519819371543616  -1.8604686755036473  4.0046301908132556 
C  2.253436413487286  -2.6033209466191805  4.356564389016076 
C  3.159658546910481  -2.9762900605645797  3.3883029341720134 
C  2.969658850627824  -2.6196060620019503  2.0836311859602583 
H  0.9520950013499276  3.1666082729624474  0.4909390864195121 
H  3.177893644192098  3.9783058604068295  0.5895320276426723 
H  4.652918670646986  2.0778575811449684  0.06841552513832069 
H  3.9005452716367044  -0.7934436091866708  -0.010702822050166053 
H  3.287883031762876  -0.5439770664944517  -1.4410894914196697 
H  0.7927589335176154  -4.054254439453521  0.5004030141665615 
H  0.6650769664843524  -5.993904735309693  -0.7293917756579965 
H  1.2799271006665232  -6.000860744285809  -2.9360657607394565 
H  2.0557242417843136  -4.084586082928129  -3.9196526395010785 
H  2.2268725737331603  -2.125675411766182  -2.694580461676834 
H  0.18486940215682957  -1.0171118028696318  2.450895419357156 
H  0.5465292531988073  -1.5935758973624183  4.657677227761676 
H  2.3854519700406005  -2.8530817529852257  5.242321992670327 
H  3.907287950073014  -3.474325784572633  3.625440242790048 
H  3.592022900789403  -2.874348026577311  1.4406092885486939 
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
