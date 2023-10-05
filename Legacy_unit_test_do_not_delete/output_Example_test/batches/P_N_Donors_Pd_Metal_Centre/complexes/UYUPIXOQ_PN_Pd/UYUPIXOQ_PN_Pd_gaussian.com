%chk=UYUPIXOQ_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.5872363056583587  1.4045571935667125  2.287758239613427e-15 
N  1.5872363056583607  -1.404557193566712  -1.6653345369377348e-16 
N  2.9303481390273594  -1.2693458964619357  -0.27125222247835606 
C  3.664752564429291  -2.0510842585800098  0.5675262185554051 
C  5.157423182819648  -2.0803802674491894  0.5635476687500824 
H  5.466968655325294  -2.510557267551038  -0.23626610272090137 
H  5.470144694505219  -2.5661126271742725  1.3299152899639604 
H  5.495262638122873  -1.1820808847120419  0.5943005777909489 
C  1.4860920465874425  -2.273417540151122  1.0090653333606454 
C  2.7680307936836632  -2.7017545511006342  1.3814943114816465 
H  2.973471083723782  -3.3117407118096573  2.0527627086013798 
C  0.1746063154832176  -2.643940187704158  1.5958488664973272 
H  -0.4891165708741285  -2.0054921101838405  1.3253780132874342 
H  0.24231269556882173  -2.6485736052892928  2.5537192113442675 
H  -0.07862470656259957  -3.5180222013261164  1.2894779753601013 
C  3.4263626887336875  -0.2925812077160397  -1.1919274742901889 
C  3.008129239040321  1.0492321400296203  -1.0755251917989341 
C  3.5881811145914693  2.000703576095798  -1.9068466567087443 
H  3.3368536116926903  2.89266711146357  -1.831922576061098 
C  4.539849082127993  1.6283575994453563  -2.8479773560471267 
H  4.944457569831162  2.2744795901743133  -3.381386293706595 
C  4.887401381178161  0.3020300931446911  -2.993241961010639 
H  5.4982554292214925  0.05502154013661431  -3.650279255789197 
C  4.332791474220141  -0.6746941771687909  -2.1649652615408117 
H  4.5685725361313025  -1.5690003796882253  -2.264693577375941 
C  2.1214409047159757  1.1015564528711455  1.6951959731794102 
C  3.4301028833981637  1.410972453864599  2.1037647114363986 
H  4.04630760623974  1.7396386926723912  1.4886425016149025 
C  3.794436274923336  1.2278474724269635  3.4157125728216173 
H  4.663433358410526  1.4252695220736946  3.683296472195434 
C  2.8876811680237715  0.7531040778153666  4.337948245187557 
H  3.146735783523864  0.63757456378544  5.222980972026573 
C  1.5798819886102102  0.44312958051888485  3.9505748757478645 
H  0.9686593636374157  0.12412062022002535  4.573579703955843 
C  1.2061232487690718  0.6171664735123699  2.6244520555424415 
H  0.340054915126915  0.4097754528161641  2.3557672290320846 
C  1.2307790616678442  3.177687388198199  -0.10590862142807447 
C  1.3031871128231003  3.9753044344326276  1.0467687756280932 
H  1.5772881622186432  3.6081324514317763  1.856359387986903 
C  0.9571707516775361  5.326761739762956  0.963798602138632 
H  1.0201889347400694  5.865966254515975  1.7197072076416189 
C  0.5248201558292049  5.868376565762469  -0.22503507947000484 
H  0.27875192009387595  6.764741353018248  -0.26483509869122124 
C  0.4553448186377309  5.078510061869686  -1.364256315029783 
H  0.172487716815382  5.450882822762287  -2.168833599629526 
C  0.8058351985263656  3.7383486582460796  -1.3122920403487435 
H  0.7574811616231046  3.2144890437212044  -2.078542157076323 
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