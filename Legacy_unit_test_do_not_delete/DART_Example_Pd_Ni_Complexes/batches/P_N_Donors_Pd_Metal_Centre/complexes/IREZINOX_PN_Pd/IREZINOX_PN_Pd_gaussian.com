%chk=IREZINOX_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  2.9153215570795754  0.9302449660989522  0.3529258943163101 
P  1.4713292748395919  -1.5255458580455712  -7.6251812536861e-17 
N  1.471329274839592  1.5255458580455716  1.4624142218674578e-16 
C  3.0627720779120136  -0.7051175448623357  -0.39766801163326754 
H  3.8246084364870034  -1.204181896802162  -0.009932863141261258 
H  3.1888154350723372  -0.6337913999115694  -1.3769759423487615 
C  4.301240007537612  1.9192383267658784  -0.24492254905107613 
C  5.282833503405342  1.3246080308518884  -1.0437736178063741 
H  5.232844276937905  0.39797197298140136  -1.2501098082612638 
C  6.331134867662716  2.090782439071536  -1.5360125088800436 
H  7.002208908899989  1.6881297906700279  -2.073721813163213 
C  6.398187891100333  3.441111997593704  -1.2428584957244944 
H  7.098248172492498  3.972641830488297  -1.6049318711060097 
C  5.44964032574412  4.015950885834979  -0.42426745080502587 
H  5.513609665821573  4.9388483081410754  -0.20685114177826797 
C  4.407646246062684  3.265111781943694  0.08064608434390734 
H  3.764448312900581  3.668442971075092  0.6511840678733848 
C  3.0322848600941557  0.7514652863131135  2.1388316622233754 
C  2.261839343211199  1.572462871141262  2.9557296553787813 
H  1.6740276052424496  2.2141362465022  2.572130671987189 
C  2.358144179595262  1.449004451974964  4.340498357085329 
H  1.8364771197968266  2.0087087935369423  4.90205780964393 
C  3.2080767828092256  0.5157539485623887  4.896884291373886 
H  3.264720412837002  0.4304547624948452  5.842029686555719 
C  3.975273338833712  -0.29589466787862023  4.081140006776318 
H  4.562158954323935  -0.9358290689661793  4.469959935066231 
C  3.8972900203348058  -0.18416435567524336  2.702717444798002 
H  4.429197235032283  -0.7398904228589116  2.1457233575714274 
C  1.9361934027254741  -2.704572667049086  1.2961156531173927 
C  1.535289316882158  -2.5191073278671623  2.616682024761151 
H  0.9431273841758464  -1.808907143266985  2.837855719506753 
C  1.9979470373888617  -3.365146972098849  3.6049339582497066 
H  1.7311689945642301  -3.231525977407053  4.5061668724534885 
C  2.834074254690168  -4.389891204961694  3.291255119782176 
H  3.143936246445657  -4.969391578770576  3.9776244393330296 
C  3.2424171988490462  -4.598421432566623  1.9885281514074726 
H  3.8276130901461567  -5.318376718871225  1.7835090924615924 
C  2.799841588724406  -3.7611508233196576  0.9804249492430949 
H  3.0790260706804427  -3.900718823768215  0.08304604471850732 
C  1.0820608885725798  -2.5185061823852593  -1.4748552022962602 
C  -0.1639341234197429  -3.1481945230757704  -1.5594629069616301 
H  -0.7834295596030263  -3.0677844492959148  -0.8433874887389174 
C  -0.49873922344217747  -3.8804014633712574  -2.6682950455570746 
H  -1.3498043323038142  -4.298897499811564  -2.716904385833327 
C  0.3902433563204186  -4.013664543551918  -3.711173736651796 
H  0.15253425972573065  -4.518697873711556  -4.479267989217228 
C  1.6306977073490816  -3.4069522083450443  -3.6330316524595103 
H  2.2492003629945803  -3.5040083926735117  -4.348605682202842 
C  1.979501482923967  -2.6588543527603536  -2.5212286732341846 
H  2.8321484933333347  -2.2434763740138965  -2.4761556168423153 
C  1.2815519755822278  2.269443492466011  -1.2874347463814548 
H  2.1313409935226453  2.769250054149339  -1.4560080713357881 
C  0.19239097906576985  3.303624645840995  -1.151679799247682 
H  -0.6623714699427976  2.856194706976705  -0.9807539510771839 
H  0.12969639731408988  3.8222718899255868  -1.9806230576330095 
H  0.4015426007810392  3.903352104163437  -0.40505069797450266 
C  1.0637775730761243  1.3595688153597554  -2.4924019001138604 
C  -0.10011093603582455  0.6255843670426658  -2.6533791282732184 
H  -0.7779558210841409  0.6720107721398797  -1.9885496082279652 
C  -0.2909712676031688  -0.1737189010755276  -3.765676296034427 
H  -1.0840557294006155  -0.690319851380397  -3.8474598799591178 
C  0.664341237319805  -0.22310645609568605  -4.753503217790968 
H  0.518933877945626  -0.751466417055916  -5.529420359948508 
C  1.8428434647950853  0.49835426840199826  -4.617044055976476 
H  2.5077545124731246  0.45922079385900183  -5.294262240792832 
C  2.0447664023061636  1.2771971453524111  -3.4855864826090452 
H  2.857903747760038  1.759412990336879  -3.385075208245686 
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


