%chk=OKOWOTIZ_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Ni  0.0  0.0  0.0 
C  2.96149194606992  -0.4779613527704972  -0.20550781716447686 
H  3.612138880228538  -0.7730962753971785  0.4634960725971087 
C  3.4679173440584323  -0.8643134374329423  -1.5932988291504944 
C  4.5856506283546805  -0.3812138788570034  -2.1560304807217214 
H  5.046320599311398  0.29721591911311007  -1.719081688430254 
C  5.1204519104481  -0.8920791178795546  -3.4623814313585433 
H  6.088973613388927  -0.9331317084675584  -3.4185640177122663 
H  4.878064547826265  -0.2785385416777739  -4.171616553017946 
C  4.572797641811704  -2.274427479099748  -3.782062929534744 
H  4.956907998541157  -2.9258148094431395  -3.1725526099163828 
H  4.825359008902218  -2.522757751303537  -4.685198761932404 
C  3.0498002931175816  -2.2939850976284757  -3.6536070404980974 
H  2.6604252698992994  -1.7051220854901055  -4.319353940024551 
H  2.721296802714522  -3.1914775614493878  -3.818688457313211 
C  2.6250972171301807  -1.846512039822763  -2.273677245260334 
C  1.5394627416456934  -2.284555714325719  -1.5877968877143234 
C  0.5505596754557117  -3.2860531293218793  -2.0136789763680687 
C  0.0834538317267024  -3.3599212693066733  -3.334430852629999 
H  0.34789494170252344  -2.7288537456662327  -3.965391358494734 
C  -0.7807779593132922  -4.3912934927666445  -3.6836875605556845 
H  -1.0886038327874819  -4.475697405044665  -4.557136878611401 
C  -1.1750882078267397  -5.288088526236737  -2.710600824513782 
H  -1.7463748013165954  -5.993479672058642  -2.91669490374663 
C  -0.7036777253386488  -5.1134629513322185  -1.4271195719862348 
H  -1.0050603293297946  -5.701291561519234  -0.7707609877831191 
C  2.7104894637664163  0.9926521644580752  -0.025780727325223962 
C  3.759663737084715  1.877110549988679  0.11496911626835471 
H  4.632608204543653  1.557173960028622  0.16139362636331245 
C  3.520870869514578  3.2427488304728453  0.18718478820922824 
H  4.224639512679803  3.843037976899829  0.2850289887493127 
C  2.206146669215019  3.694262076768027  0.11050385468134966 
H  2.019444215630503  4.6050250367700585  0.10407232967866181 
C  1.1814103182820261  2.7642542025568724  0.04434484182740605 
H  0.3019619246760912  3.0659296896181973  0.030923832540647844 
C  1.609132606796159  -2.579207416757804  1.4021453967025503 
H  0.7691909274891351  -3.0752645522257764  1.4969514287586698 
C  2.736552454156977  -3.6111383944479694  1.212549075376575 
H  2.5770766830232645  -4.124870903381496  0.4053062013379264 
H  3.585645865620873  -3.15329927582208  1.1179132995596137 
C  2.7861594047806513  -4.550386776999513  2.4222549430667852 
H  3.5147091743123235  -5.179437512417846  2.313102088843886 
H  1.9579626228588414  -5.053046981013204  2.4749929721849115 
C  2.986773378130017  -3.7489751077792257  3.721351656223516 
H  2.982562285449814  -4.356757358825837  4.478366413022327 
H  3.8475211167557104  -3.306356657287299  3.6971918189952797 
C  1.8925987077551825  -2.7176636860107046  3.89621046223548 
H  2.0647339884874905  -2.2000328719143454  4.698572759097149 
H  1.0402307367670325  -3.1680214990678346  4.008399509169495 
C  1.8135426151172505  -1.774577044127335  2.688039456078628 
H  1.0762837286772071  -1.154800761424351  2.8053700734129965 
H  2.633174925110361  -1.2595407556756049  2.6230208587910733 
N  1.4144073378981026  1.4420651450263953  -1.2645833150701269e-15 
N  0.15377407627182005  -4.156975091018008  -1.0623756467640615 
P  1.4144073378981037  -1.4420651450263953  2.220446049250313e-16 
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

