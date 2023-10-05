%chk=ARAREVAD_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.3164691375038005  1.531995107694538  -9.529489683048541e-18 
N  1.3164691375038005  -1.531995107694538  1.3877787807814457e-17 
C  0.8941840556925083  3.203317311753773  -0.586796191107462 
C  2.1253693045689737  1.8273345168920607  1.5989073434995746 
C  2.7003987965099223  1.0653188197552739  -1.0961626390269372 
C  3.09942740277911  -0.28020589481255165  -1.1576607192087338 
C  1.1966847321447862  4.37923270531598  0.10177553204806172 
H  1.6334927909887138  4.3354758619270495  0.9221804128146596 
C  1.0602203930620704  -2.7948842972765746  0.6516524885311757 
C  4.1934397111842205  -0.6230782364811098  -1.9700246748203007 
H  4.492315344654406  -1.503593450416743  -1.9814769124043006 
C  3.3815511109145198  2.0037233047170155  -1.862017897665918 
H  3.1318938943588073  2.8984913644289385  -1.8163070287065979 
C  -0.07538833678694257  4.539728870587695  -2.3591845546319017 
H  -0.49221182786933326  4.591132853279591  -3.1884001284396604 
C  4.833682990755115  0.3144490002484086  -2.7484314793606215 
H  5.532345565132779  0.06196377489929228  -3.3068654308216914 
C  1.3005441237288404  2.11491033080115  2.688623844387131 
H  0.3816674727975955  2.1905592473139492  2.5615776645790076 
C  0.9794906630011329  -2.82842601013317  2.0474669870169113 
C  1.8301569800640216  2.288598274492061  3.9522372586167096 
H  1.2730868998437455  2.496962079569987  4.6667374968650055 
C  0.25245636684395256  3.305117310472239  -1.8132549566703506 
H  0.03635552431388911  2.5301558278689145  -2.28031133050444 
C  0.8880883980143721  -3.946707864556424  -0.12516027278240632 
C  4.429780091229764  1.6281227906488391  -2.695740776480603 
H  4.858261525620316  2.267191783155838  -3.21808153996622 
C  4.025291308846002  1.8720570291611152  3.0826418529588877 
H  4.941913977524868  1.7915324639190984  3.2156853332826865 
C  0.22563756457086326  5.689154838827212  -1.6592222483264962 
H  0.005067708794018744  6.52031163796099  -2.0155178755467533 
C  0.7497755272515133  -4.048394801556135  2.6562637054390223 
H  0.6990101101995768  -4.096973507601742  3.583579742047082 
C  2.5150832708749444  -1.4116619688521932  -0.4408938477449405 
H  3.081680173941777  -2.1328360473142576  -0.2858327694903695 
C  3.4915455808087876  1.710856220647019  1.8067158099859948 
H  4.054093708683409  1.523847496819473  1.089917420592459 
C  0.8534880617417555  5.608767460707364  -0.42797775081947853 
H  1.0443004262845461  6.385582539308311  0.04531938121894888 
C  3.196551207538964  2.1515042400887365  4.152893534008939 
H  3.5533131414642276  2.2477021974188878  5.005902734139529 
C  1.1836661469032568  -1.5770300165225604  2.870318066464853 
H  1.0460811292818817  -1.7777334120201571  3.7985317756729753 
H  0.5559407585215995  -0.9056652681404294  2.5926292632860153 
H  2.0781366951011444  -1.2515851064860652  2.7414643938826 
C  0.5932088654056404  -5.204132460969753  1.899776796754153 
H  0.44952263117471336  -6.019618387605689  2.3247147917031508 
C  0.6508992907252673  -5.151424590301996  0.5357027902588495 
H  0.5291060035380162  -5.928754311538578  0.04077068412873809 
C  0.9821530302909067  -3.9033516360238436  -1.631102858806904 
H  0.4955947213505082  -3.143528961556075  -1.9601868076189186 
H  0.6072262324264218  -4.70701702230097  -1.9989230909083042 
H  1.9032754347821332  -3.831247136835047  -1.8929090221053446 
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
-P 0
6-31+g(d)
****
-N 0
6-31+g(d)
****
-C 0
6-31g(d)
****
-H 0
6-31g(d,p)
****
-Br 0
6-31+g(d)
****

-Ni 0
lanl2dz

