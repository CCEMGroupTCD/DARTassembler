%chk=ONOQOKIP_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.3164691375038005  1.531995107694538  -8.154592623070045e-18 
N  1.3164691375038005  -1.531995107694538  1.3877787807814457e-17 
C  0.9361455540259651  3.203317311753773  -0.614823903300102 
C  2.0118647262910656  1.8273345168920607  1.6514385088371104 
C  2.7734920537652186  1.0653188197552739  -0.9969543890090269 
C  3.1758385356201573  -0.28020589481255165  -1.0304678344620157 
C  1.1898770187067926  4.37923270531598  0.09317187423983066 
H  1.5683924835485228  4.3354758619270495  0.9420484798382549 
C  1.015387622417994  -2.7948842972765746  0.6321900870074383 
C  4.323853529967151  -0.6230782364811098  -1.7645384714399726 
H  4.622799985125987  -1.503593450416743  -1.7551143016698831 
C  3.506408477607581  2.0037233047170155  -1.7134292791942063 
H  3.2541707844532812  2.8984913644289385  -1.6852449666859504 
C  0.0925705535086474  4.539728870587695  -2.4505287689899493 
H  -0.2653944220964022  4.591132853279591  -3.3068005522477657 
C  5.016836122910683  0.3144490002484086  -2.4963880026057748 
H  5.752751173644761  0.06196377489929228  -3.004725399502296 
C  1.1130339976531856  2.11491033080115  2.6809636190266977 
H  0.20525797756711106  2.1905592473139492  2.490129322352489 
C  0.8374874485067257  -2.82842601013317  2.0189730302328917 
C  1.553211527736148  2.288598274492061  3.9784428597957047 
H  0.9476574251235081  2.496962079569987  4.652343366949756 
C  0.38153452025256795  3.305117310472239  -1.8830597376297238 
H  0.19854029419209884  2.5301558278689145  -2.3640528183593683 
C  0.8978626512475827  -3.946707864556424  -0.15473771845512585 
C  4.610241592916491  1.6281227906488391  -2.472000493635666 
H  5.074115917821846  2.267191783155838  -2.963179507346274 
C  3.80365852870306  1.8720570291611152  3.264090575010857 
H  4.708767706588295  1.7915324639190984  3.460750333224288 
C  0.3440362686585804  5.689154838827212  -1.7312730302983532 
H  0.148857636529941  6.52031163796099  -2.1020869146313017 
C  0.5658643749782336  -4.048394801556135  2.6102626326113914 
H  0.4505363231706154  -4.096973507601742  3.531778557479322 
C  2.542918707202352  -1.4116619688521932  -0.3562087571672384 
H  3.097318894656378  -2.1328360473142576  -0.16200159791367097 
C  3.36021707994947  1.710856220647019  1.9540404038869603 
H  3.971396196887729  1.523847496819473  1.278229473055723 
C  0.8844700791110742  5.608767460707364  -0.45923114578399216 
H  1.0418020953598242  6.385582539308311  0.026223455925901487 
C  2.902280213341707  2.1515042400887365  4.273925189639 
H  3.198670179254307  2.2477021974188878  5.149742956706925 
C  0.9837663814730626  -1.5770300165225604  2.854062247589658 
H  0.7817675988042911  -1.7777334120201571  3.7704174291847257 
H  0.37694069229611404  -0.9056652681404294  2.5332619708056145 
H  1.8850464221997851  -1.2515851064860652  2.7879175673879235 
C  0.46244896083505693  -5.204132460969753  1.8446969497789238 
H  0.2894705629615255  -6.019618387605689  2.258576772104012 
C  0.6151518477972493  -5.151424590301996  0.4879700396444181 
H  0.5281799614026652  -5.928754311538578  -0.01425230699024975 
C  1.0967473912697443  -3.9033516360238436  -1.6504502869768256 
H  0.6343300497505043  -3.143528961556075  -2.012673195731831 
H  0.7483917386042378  -4.70701702230097  -2.0435280988106155 
H  2.0338886626352606  -3.831247136835047  -1.8473644527957886 
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

-Ni 0
lanl2dz
