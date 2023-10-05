%chk=JOBINITA_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
P  1.2813461261501506  -1.5614903473925146  1.7157980928477425e-15 
N  1.2813461261501513  1.5614903473925155  -8.020511312085761e-17 
C  1.4026187490750983  -2.7401838071845805  -1.3835485314120173 
C  0.8098089045165268  -2.4434795875436133  -2.613484019249368 
H  0.349626172188763  -1.6201913987010055  -2.7302886466389693 
C  0.8902649005366898  -3.3496497904928817  -3.667895971946851 
H  0.4804162577680444  -3.1467656833134714  -4.499724168796885 
C  1.5706621412276927  -4.552412995616874  -3.504415308686815 
H  1.6275230684647857  -5.168301967717101  -4.224886278954953 
C  2.167065775243128  -4.851745047704486  -2.284344868947273 
H  2.627380435467899  -5.675176487512669  -2.1716111228124935 
C  2.0929362568113863  -3.951312792488909  -1.2284524983599634 
H  2.5097857172978335  -4.156413781768028  -0.40060349552858343 
C  0.890739113147254  -2.6027525349654788  1.4500692804565432 
C  0.026032436776950885  -3.6901349692119236  1.3321119342125902 
H  -0.3401619886007736  -3.9107107163339903  0.4832088698973245 
C  -0.3057587500953063  -4.45467668219371  2.4470938118342853 
H  -0.9049827006872659  -5.185698386494429  2.360646979166728 
C  0.23885861744415293  -4.148303651583008  3.682030838199946 
H  0.01198191176060348  -4.672267079449359  4.441833619635211 
C  1.107292917426763  -3.087603377152393  3.813691365452964 
H  1.4871643615000558  -2.88949697802123  4.661841063679058 
C  1.4320709589683025  -2.3009500033091377  2.7039877684715456 
H  2.020408147342953  -1.5619314260076653  2.8023504909772123 
C  2.9705609842266276  -1.047391770772747  0.36504612964833505 
C  4.1604615667521  -1.9539771828208627  0.5972514576818465 
H  4.5268919914517145  -2.2869282515662857  -0.25976986625107873 
H  3.9135299586277013  -2.728061330710992  1.1639717792641249 
C  5.170376439762223  -1.0484685561674234  1.3130205449783288 
H  6.096262851124919  -1.266958965390708  1.0379521456094125 
H  5.09651602171413  -1.1452584948643723  2.295844840777346 
C  4.790213631857769  0.3801659555509108  0.8721461394912935 
H  4.913668456201499  1.0265745111244369  1.6126561904883083 
H  5.331816719108772  0.6701735279679867  0.09653980002457543 
C  3.324926384198481  0.24849678427276217  0.49830312272883653 
C  2.544256484146625  1.4630109329940948  0.28605854751258114 
H  3.018067546465743  2.282782276529801  0.36739372773423784 
C  0.8543751212959498  2.9494818541086243  -0.16851646298727155 
C  0.421568170559803  3.654957130878354  0.9589987818148473 
C  0.08346818543874268  5.0038324880121134  0.7827635952645616 
H  -0.21976468814735872  5.512990470530303  1.5250593201712273 
C  0.18747010075186243  5.607795014309391  -0.47152439337276053 
H  -0.02078776786504255  6.529352424128124  -0.5714238957105121 
C  0.5907538493473007  4.873440955115357  -1.5669994062617236 
H  0.6459711654044147  5.2939589372447315  -2.4164358516261224 
C  0.9212058445578635  3.5125936731911467  -1.4447756191299015 
C  0.3173054098122098  2.9911985929336686  2.2940706113531806 
H  0.5947751168948628  2.0544857099866394  2.2156971683221895 
H  -0.6101663156041481  3.0312833876328606  2.606896837534987 
H  0.898514895738304  3.453011543749642  2.93533534758011 
C  1.332792084040331  2.695788957023797  -2.642676496963913 
H  1.52128185549638  1.7756261859895028  -2.362088291832362 
H  2.1375298532107374  3.0864705683566456  -3.0439920977533856 
H  0.6072551266707653  2.694445954805754  -3.3012356665424494 
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

