%chk=EGURAGIS_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.2813461261501506  1.5614903473925146  -1.5245706772643693e-15 
N  1.2813461261501513  -1.5614903473925155  -1.1102230246251565e-16 
C  1.4026187490750983  2.7401838071845805  1.3835485314120177 
C  0.8098089045165268  2.443479587543613  2.6134840192493685 
H  0.349626172188763  1.620191398701005  2.7302886466389693 
C  0.8902649005366898  3.3496497904928813  3.6678959719468516 
H  0.4804162577680444  3.146765683313471  4.499724168796885 
C  1.5706621412276927  4.552412995616874  3.5044153086868155 
H  1.6275230684647857  5.1683019677171  4.224886278954954 
C  2.167065775243128  4.851745047704486  2.2843448689472736 
H  2.627380435467899  5.675176487512669  2.1716111228124944 
C  2.0929362568113863  3.951312792488909  1.2284524983599638 
H  2.5097857172978335  4.156413781768028  0.40060349552858393 
C  0.890739113147254  2.6027525349654788  -1.450069280456543 
C  0.026032436776950885  3.6901349692119236  -1.3321119342125898 
H  -0.3401619886007736  3.9107107163339903  -0.483208869897324 
C  -0.3057587500953063  4.45467668219371  -2.447093811834285 
H  -0.9049827006872659  5.185698386494429  -2.3606469791667277 
C  0.23885861744415293  4.148303651583009  -3.6820308381999456 
H  0.01198191176060348  4.67226707944936  -4.44183361963521 
C  1.107292917426763  3.0876033771523934  -3.8136913654529634 
H  1.4871643615000558  2.8894969780212305  -4.661841063679058 
C  1.4320709589683025  2.300950003309138  -2.703987768471545 
H  2.020408147342953  1.5619314260076658  -2.8023504909772123 
C  2.9705609842266276  1.047391770772747  -0.36504612964833494 
C  4.1604615667521  1.9539771828208627  -0.5972514576818463 
H  4.5268919914517145  2.2869282515662857  0.259769866251079 
H  3.9135299586277013  2.728061330710992  -1.1639717792641244 
C  5.170376439762223  1.0484685561674236  -1.3130205449783285 
H  6.096262851124919  1.2669589653907083  -1.0379521456094123 
H  5.09651602171413  1.1452584948643725  -2.295844840777346 
C  4.790213631857769  -0.3801659555509107  -0.8721461394912935 
H  4.913668456201499  -1.0265745111244367  -1.6126561904883085 
H  5.331816719108772  -0.6701735279679867  -0.09653980002457552 
C  3.324926384198481  -0.24849678427276212  -0.4983031227288366 
C  2.544256484146625  -1.4630109329940948  -0.2860585475125813 
H  3.018067546465743  -2.282782276529801  -0.3673937277342381 
C  0.8543751212959498  -2.9494818541086243  0.1685164629872712 
C  0.421568170559803  -3.654957130878354  -0.9589987818148478 
C  0.08346818543874268  -5.0038324880121134  -0.7827635952645623 
H  -0.21976468814735872  -5.512990470530303  -1.525059320171228 
C  0.18747010075186243  -5.607795014309391  0.47152439337275986 
H  -0.02078776786504255  -6.529352424128124  0.5714238957105113 
C  0.5907538493473007  -4.873440955115357  1.566999406261723 
H  0.6459711654044147  -5.2939589372447315  2.416435851626122 
C  0.9212058445578635  -3.5125936731911467  1.444775619129901 
C  0.3173054098122098  -2.991198592933668  -2.294070611353181 
H  0.5947751168948628  -2.054485709986639  -2.21569716832219 
H  -0.6101663156041481  -3.03128338763286  -2.6068968375349875 
H  0.898514895738304  -3.4530115437496414  -2.9353353475801103 
C  1.332792084040331  -2.6957889570237974  2.6426764969639125 
H  1.52128185549638  -1.775626185989503  2.362088291832362 
H  2.1375298532107374  -3.086470568356646  3.043992097753385 
H  0.6072551266707653  -2.6944459548057544  3.301235666542449 
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


