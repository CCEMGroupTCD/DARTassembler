%chk=OBEFARUX_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.4655122014162827  -1.3900985531608905  -9.948735041333616e-16 
N  1.465512201416283  1.3900985531608905  -3.9228257928781644e-16 
C  2.9694900156872013  -0.4190055706306354  0.4784644328071103 
C  2.6931820460954397  0.9682111670712782  0.10777511736930871 
C  1.7780241923041686  -2.0677100274395626  -1.6153310265940064 
C  0.8308660682469942  -1.8439538525318697  -2.6359150712582076 
C  1.0932547673077846  -2.407070196491837  -3.8895970818854697 
C  2.2670971216553455  -3.081893429612221  -4.153276226265063 
C  3.2048092166564306  -3.309061038902906  -3.117859964958606 
C  2.932881751420961  -2.758468134349745  -1.8962295596778773 
C  1.4064815754209379  -2.7815268663562733  1.1387366440195412 
C  1.0519162145055492  -4.084417449820091  0.6654310953934313 
C  1.0074398579671444  -5.125954862191188  1.5280640256368827 
C  1.3471156989000241  -4.9438474609394  2.8961849075387147 
C  1.6429295198115328  -3.6664105475909445  3.3769969350219067 
C  1.6497874635365593  -2.60781528709384  2.4955702602151386 
C  3.909444866419215  1.888230857108304  -0.05410218753237107 
C  4.7135452792651895  2.176293494170595  1.0550821673832702 
C  5.906410667632823  2.894367099621265  0.815304932455338 
C  6.208556301158502  3.312982440725536  -0.392882301516928 
C  5.420643134936076  3.059235072096216  -1.5181067748586232 
C  4.297067449044949  2.3198848924539197  -1.3387986603260191 
C  1.2486396543557423  2.7519276196651377  -0.46283195865848487 
C  1.2697155631062524  3.7502534878882234  0.5044190723086497 
C  1.1113844203019896  5.155530067421433  0.012930044726400365 
C  0.8629503426527227  5.256063795538926  -1.4082372566753008 
C  0.6982292666858025  4.2309865002926  -2.20965111723847 
C  0.9338648202553335  2.915724049409099  -1.7836744146009975 
C  1.4955445170333816  3.5570649082804318  1.9484970017422931 
C  0.7568275025892258  1.803870805992972  -2.718754744443621 
H  3.1276350489681857  -0.48927767096328917  1.4635136855672959 
H  3.771982265677659  -0.7577201127969877  -0.012116139199457218 
H  0.011084368200944716  -1.300840802537479  -2.4727902442225558 
H  0.41094490932396766  -2.323174180820996  -4.61623782960933 
H  2.4571281334288773  -3.4112059693304153  -5.080329814656114 
H  4.031132110826884  -3.8544576212359085  -3.2727757567595064 
H  3.615645600559003  -2.85909098564467  -1.16997203778397 
H  0.8370506587399702  -4.225977078489976  -0.3003765119882971 
H  0.7379805857145584  -6.031753779927538  1.2004438320822568 
H  1.3680045717479246  -5.726441578116422  3.515653908176425 
H  1.8498813449790121  -3.5232325588451627  4.34441506365019 
H  1.8344614653937694  -1.6871288853349062  2.842102903200214 
H  4.455933822300933  1.8852339293073752  1.9765549596796055 
H  6.52608281546805  3.0934430423813977  1.5764882337773338 
H  7.047485977577565  3.840665599111718  -0.5163938085880915 
H  5.68226398825468  3.3917123874268986  -2.42449100568976 
H  3.731413823956877  2.083421201262698  -2.129044058397373 
H  1.1680962378252833  5.948772030369058  0.6114511287094205 
H  0.8190714193958026  6.1654183172925325  -1.81036112299026 
H  0.3946504652346694  4.395402249435205  -3.147263543960395 
H  1.4577272328975477  4.44288760567121  2.410523886033561 
H  2.3966863798033997  3.144801393800195  2.0923664140292035 
H  0.7908900588807092  2.9552005526219256  2.321486615529777 
H  0.5315808433464261  2.166391336536732  -3.6240207481321893 
H  1.6047217061531593  1.2765458200535336  -2.7704349027895923 
H  0.015167911999441008  1.216676367049098  -2.3996792634683946 
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


