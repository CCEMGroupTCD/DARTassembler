%chk=AQAPAKUP_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.2724242000999517  1.5687691528711292  -7.095935914785308e-16 
N  1.272424200099951  -1.5687691528711292  1.6653345369377348e-16 
C  2.5471073234460135  -1.4779289424174835  0.21406165285145595 
H  3.0104173593315173  -2.299706797776616  0.325828003284515 
C  3.355189105632089  -0.27324353568240634  0.30431637072560436 
C  2.9730396655674247  1.0215183203131506  0.2655394224545703 
C  4.1928023447465685  1.921283029108361  0.3360787397782307 
H  4.394752647379547  2.3127984530800196  -0.5504680260828594 
H  4.053679592156989  2.6543211137082214  0.9867896404384822 
C  5.299375050012927  1.0107700818172956  0.7839437091326207 
H  5.445620006390735  1.0995843664078628  1.7590753646472688 
H  6.142243065095107  1.2364872654780887  0.3163751877887497 
C  4.863956996347764  -0.4148467407132026  0.43834746514191086 
H  5.27291588888463  -0.7196956976962867  -0.41001772605902215 
H  5.102551541578031  -1.0483817564808438  1.1606749440478847 
C  0.8035924284182944  2.5953181177412277  1.4632890051476553 
H  -0.1446004578897584  2.8774423481864697  1.316936147469711 
C  0.7992552748583582  1.695890918769583  2.705144372269411 
H  0.22251408332660838  0.9081883541580829  2.538823364290811 
H  1.7177376269509381  1.3714885080843664  2.8801379288254703 
C  0.28366977501100665  2.454699292724007  3.932134911745482 
H  -0.6721422124483927  2.677261082239689  3.801942652081176 
H  0.35589412744131854  1.8737499648181002  4.730414443815361 
C  1.0696559172522166  3.732744861972081  4.16788680915071 
H  2.0038855047064814  3.505022676495536  4.402890520114692 
H  0.6738466190162027  4.225949453982285  4.929441442680027 
C  1.0619189065796637  4.622001048866647  2.930316890828684 
H  1.615305647392558  5.425160451060414  3.1002842399290986 
H  0.1361732086431262  4.921927382284711  2.7477242916511213 
C  1.6045153436916118  3.8831604311258583  1.7091149591196293 
H  2.557681063855999  3.6584408225366647  1.8538124277189958 
H  1.5437205645679963  4.468597353645105  0.913086282134612 
C  1.4839979262098564  2.7217379256760243  -1.4227365804521914 
H  2.187767091203863  3.3838218614157536  -1.1651747151866947 
C  1.987783804651753  1.986038936511357  -2.6722789550936867 
H  2.8282441447887683  1.508326490018879  -2.4588504070431125 
H  1.3163999716922319  1.3157317003731017  -2.955096861212249 
C  2.236963801781707  2.975735608681671  -3.8150351084900738 
H  2.975854645980108  3.5855841407292623  -3.5654588710009993 
H  2.5120087227501173  2.4782009404537084  -4.625585819649326 
C  0.9846352195421955  3.7980588502108685  -4.126060024769988 
H  1.1966427295866584  4.4665914689716075  -4.82477171304898 
H  0.279242598354033  3.2000230504356355  -4.479366240623172 
C  0.4650780283392716  4.513098811920358  -2.884812846232255 
H  -0.37612321227746115  4.986487339888823  -3.1043047922639104 
H  1.1282287695714728  5.186420953082303  -2.589763564544799 
C  0.20911365900705237  3.5165341694178376  -1.749288101919648 
H  -0.5110364128972824  2.891167028430592  -2.0148395148506513 
H  -0.09109505119481676  4.005363324807915  -0.9423812565158072 
C  0.7714778980212841  -2.943643847501984  0.02127314028123059 
C  0.5253628985751174  -3.522790427696525  1.2655489359538017 
C  -0.012256259400294178  -4.816037787441533  1.288814521920662 
H  -0.18833035944615273  -5.2375264900909  2.1216901157712846 
C  -0.2881828536997648  -5.481703987015118  0.1156440546382676 
H  -0.6517827718725373  -6.358605739884432  0.14813459663870382 
C  -0.042962785410197935  -4.891941115717762  -1.1021764303870054 
H  -0.24550448326586194  -5.362774181658774  -1.9019644302009266 
C  0.5052085420367667  -3.597343394770074  -1.1755503813647623 
C  0.7959724072542591  -2.774959695014866  2.545741354721122 
H  0.47878296801966347  -3.3047042488767495  3.306918191791013 
H  0.3256476578087114  -1.9153542058758428  2.5276194590960612 
H  1.7591841127163796  -2.618484575937761  2.6355491359574215 
C  0.7984800541431365  -2.956420749676097  -2.500331028563637 
H  0.34529193881074827  -3.454413152554621  -3.211858483267088 
H  1.7656408652283608  -2.965210685387647  -2.6603012211970967 
H  0.4774314534925759  -2.0304261714799687  -2.4935246218979144 
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
