%chk=EYOVAHEZ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.6220537213668353  -1.3642000311537898  1.8995910089389847e-16 
N  1.6220537213668353  1.3642000311537898  -5.604401769240518e-17 
N  3.732900038818635  1.921828599272656  0.08819526131247855 
C  1.7031131688485206  2.702387205549332  0.33868560975813955 
C  2.982627151780262  3.0455929813210956  0.38109554914846666 
C  5.202727934956048  1.8372156009117226  -0.031182317171544297 
C  2.8720908458699284  0.9301487875841816  -0.14737739080256151 
C  3.149126027035093  -0.4781699238404889  -0.48629634195156507 
C  1.7663838922101371  -1.6884025746956919  1.7806168741605763 
C  0.6210840463188392  -1.6363955937692989  2.559936338139472 
C  0.7232982806666721  -1.8685725933051487  3.9352522659114384 
C  1.9261272923395767  -2.152980038063711  4.488431038139643 
C  3.058371233910555  -2.203224938690218  3.724597000709482 
C  2.991409956843583  -1.9889558884595284  2.350166710921445 
C  1.7001235567766064  -2.9622982500463158  -0.8417146449787132 
C  1.686117441883219  -2.995998146652473  -2.230362255322818 
C  1.656473319730525  -4.210333773517188  -2.885955041546748 
C  1.6719420133053358  -5.394769657625336  -2.1699005701857823 
C  1.7074813439712502  -5.368971864635817  -0.8122136991242962 
C  1.7092616404348233  -4.159451338086299  -0.12471884357998886 
H  0.9470666844082688  3.2793068155968097  0.5000625502057033 
H  3.338968207212428  3.9085397032943687  0.5849939239558141 
H  5.615532853641433  2.6709521417605346  0.1959054402320919 
H  5.550710154766033  1.1668839493088565  0.5694947604296343 
H  5.462163561825902  1.6010101892116007  -0.9174672594977978 
H  3.901919589294058  -0.8056245087976428  0.00686295017625581 
H  3.323235240642095  -0.5816985829628241  -1.4154577641369106 
H  -0.21725349878446854  -1.4348018181355275  2.1624361996541013 
H  -0.050756232582802596  -1.8354023954114218  4.493573972895638 
H  1.9849146562034456  -2.339715728712059  5.429242415456057 
H  3.9018390797846947  -2.386949817209698  4.139753371236802 
H  3.7769295796859113  -2.05638668152453  1.8059974906633653 
H  1.677750706654677  -2.1760784798467303  -2.7342566467648846 
H  1.6179970002531885  -4.235248466541535  -3.8393547693162806 
H  1.6785053690055167  -6.224073809629465  -2.6381104011704335 
H  1.7353821232689717  -6.188825763960289  -0.3254877903449877 
H  1.718008970858915  -4.1524900239218585  0.8429379314779732 
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


