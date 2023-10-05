%chk=QASOHICA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.1242715630576094  1.6781279607050235  -9.963509747049376e-16 
P  2.525344281640245  -1.427349751382299  0.6944514648068633 
C  1.1381912777858398  2.9517434030590577  1.3627076122702357 
H  0.2735432907630958  3.4463482600398647  1.2789194031443214 
C  1.1101567862691115  2.299229118904016  2.743580334956555 
H  0.38503790337226895  1.640886825664439  2.7785397485458745 
H  0.961133132922843  2.9858721115610845  3.4256775137942146 
H  1.9668218876797479  1.8524391406922214  2.910743965554315 
C  2.231629465326616  4.014146105205194  1.2962404888169572 
H  3.107268149988056  3.590972168754215  1.417067031567339 
H  2.0889940082874974  4.675268728734982  2.0047212175967566 
H  2.201094122525001  4.459170124521483  0.42388832078412275 
C  0.7369122632444897  2.746986827252023  -1.4742743117385662 
H  1.2541738048941622  3.596459053673363  -1.3643870328161258 
C  -0.742619845404598  3.1264416244490585  -1.5107437011030496 
H  -1.287007319010412  2.3157217736050124  -1.5880743032505533 
H  -0.9120297090338674  3.7082690135772376  -2.281610276695373 
H  -0.9780798722725914  3.6006612034289978  -0.6865246849068442 
C  1.1203485427843343  2.1651920071539044  -2.835549988021566 
H  2.0087123055874594  1.7553897744856157  -2.7770073109527336 
H  1.1349777173796722  2.880551576334584  -3.5052185433773246 
H  0.46400298696324305  1.48536467770478  -3.0976587698049705 
C  2.90054473465051  1.2421132483986992  -0.27082398441017225 
C  3.726738532790862  2.2336891585580854  -0.8234613268926185 
H  3.352036985373652  3.077410128087224  -1.0459155602530188 
C  5.07975671068811  2.0134888468698433  -1.0533086740864286 
H  5.617691105113773  2.7040648649901344  -1.4229388277414223 
C  5.639014327762917  0.7948098200685716  -0.7455671471867731 
H  6.565668489201729  0.6435956779161516  -0.8939353321230051 
C  4.842440063556473  -0.21564838973779055  -0.21679614019138355 
H  5.230840137980075  -1.0591954273577828  -0.01506315451411333 
C  3.4752914300748916  -0.013342598433670139  0.02383343150396558 
N  1.124271563057609  -1.678127960705024  -1.3322676295501878e-15 
C  0.7789696907227154  -2.7743572251451045  -0.84111953772766 
C  0.6851170435264866  -2.563851784992861  -2.2242698116605597 
C  0.2816804413829739  -3.602647637168297  -3.057767253500691 
H  0.22092604272250127  -3.4468795008158226  -3.993901651298714 
C  -0.035502112333362756  -4.8617479736600995  -2.5604005575265605 
C  0.08971581179473254  -5.060723944395798  -1.1947828458125807 
H  -0.10104482121257896  -5.921789809952872  -0.8402461113729752 
C  0.48626614364825715  -4.0455491842515965  -0.3164304901642079 
C  1.0083142523648652  -1.2084209696667223  -2.791989565644348 
H  0.4101805010660362  -0.5389133181000432  -2.40149631339679 
H  0.891266243995784  -1.2248304036579376  -3.7650689736131566 
H  1.9373224470976016  -0.9793870556864527  -2.580080504095916 
C  -0.519704369933695  -5.969929875485221  -3.4696713786805273 
H  -0.3446205455579845  -6.83694759050424  -3.0488020090595214 
H  -0.045421769837549864  -5.920116147356359  -4.325413902887407 
H  -1.4824108081126084  -5.868603600623222  -3.622606004947687 
C  0.5533856475910351  -4.344852522114451  1.1561941511258529 
H  -0.16834435558182848  -3.8710187033477714  1.6190889066730705 
H  1.416882520192951  -4.0495929640818  1.5108122899697278 
H  0.45432852867614393  -5.309636872838102  1.2978833313444775 
C  3.6194292173459024  -2.86122552626009  0.5472545026077698 
C  3.9100295895968227  -3.3485815647450243  -0.7380068931401763 
H  3.5578780215444255  -2.9096163086685323  -1.5032227692666325 
C  4.709089970713356  -4.469231136193626  -0.8876772487285054 
H  4.90658408956874  -4.799560317269371  -1.7562652092410729 
C  5.223502002550879  -5.110475499077554  0.23135306300144548 
H  5.772504623003917  -5.878655143410764  0.12164665130253863 
C  4.94574690938455  -4.643208619935294  1.5022884515546304 
H  5.300176655100901  -5.088290134867091  2.2629555603594413 
C  4.141575535604414  -3.517041420977458  1.662884627050973 
H  3.948188982019321  -3.1943799537479167  2.5362750944009287 
C  2.3557273571485826  -1.08327633533754  2.4618886329396896 
C  3.413390248685128  -0.5028154391429722  3.15935118810347 
H  4.190761090814353  -0.217108343959022  2.6918971480429827 
C  3.3319426351477066  -0.3418358468081255  4.535966997652661 
H  4.056760123564699  0.04661930396960812  5.013531709832222 
C  2.1965491648541215  -0.7477191292211315  5.210701262880367 
H  2.1482887560631996  -0.6500520181081366  6.153949527310349 
C  1.1277115140022804  -1.2972402304681285  4.519041624470033 
H  0.34391249312822936  -1.556943581586181  4.987298043854732 
C  1.2023537807139324  -1.4689455930788233  3.141060125874785 
H  0.4715337137003919  -1.8470507375359997  2.665683202386301 
Br  -1.5839191898578664  1.5839191898578664  0.0 
C  -1.2937473829822488  -1.2398794889499034  -0.4598005675376139 
C  -1.601562203225019  -1.4938654694550033  -1.7919791134336869 
H  -1.1309417375171225  -1.0428217585005106  -2.480879663633267 
C  -2.60218398542082  -2.411889969420291  -2.121542410022959 
H  -2.8231586547736818  -2.5626138747145775  -3.03270050271447 
C  -3.2683208079706407  -3.09759503004543  -1.143064365570095 
H  -3.9596020908062854  -3.707191723557817  -1.3683280760515548 
C  -2.9347449854640226  -2.896087320001541  0.14198826503765638 
H  -3.3668350933519986  -3.405606940434718  0.8171022443470677 
C  -1.9675509241539992  -1.9516175659140573  0.5112421903797603 
H  -1.7755298364628027  -1.8022912010877177  1.4301324877763764 

-Ni 0
lanl2dz
F 1 1.0
3.130 1.0
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

-Ni 0
lanl2dz

