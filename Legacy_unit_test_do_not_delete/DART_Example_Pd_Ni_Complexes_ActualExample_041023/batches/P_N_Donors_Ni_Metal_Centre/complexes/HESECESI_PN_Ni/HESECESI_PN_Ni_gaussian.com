%chk=HESECESI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
P  1.280604340731357  -1.5620987556809587  3.8833603903735035e-15 
C  2.8468496920088207  -0.6824961601324941  0.38479701597881455 
H  3.522094363618753  -1.3403371969189528  0.6889345239947018 
H  3.188921870021417  -0.25612819678970117  -0.44015697962618733 
N  2.6767668576762116  0.34282456890284835  1.4203045109182986 
C  2.2278389393666638  1.5612003713791363  0.9737884103099091 
N  1.280604340731359  1.5620987556809591  3.638095882029033e-16 
C  0.994836346126271  2.7408879604628895  -0.6145926325209091 
H  0.4016304945516296  2.7325757843262792  -1.3561312595562751 
C  1.5254431437818026  3.9461649048000327  -0.21235128241991233 
H  1.3183235095883419  4.748284306946148  -0.6770363928251556 
C  2.3722880478131803  3.9589740949534917  0.8898367837302051 
H  2.6936455321254793  4.783795883789967  1.2357319758579275 
C  2.74292785073828  2.77346639048527  1.4779417226094291 
H  3.339673609866889  2.770603543316526  2.2165007676885167 
C  3.672545834566029  0.320030248204168  2.492663866791234 
H  4.531146901655754  0.6447423562929318  2.1488750443246527 
H  3.778837803432103  -0.5979133675788927  2.8204430845547894 
H  3.3754499205021764  0.8964585640999049  3.2273099002523664 
C  1.5993828546263968  -2.3031145498427037  -1.6930646094823947 
C  1.8943445000296055  -1.1290280951098883  -2.6433410506056965 
H  2.071097843652612  -1.473979138083546  -3.5439249575624543 
H  2.679179671878967  -0.6383303897232314  -2.321639659628803 
H  1.121300081942353  -0.5281145859055344  -2.669545479999481 
C  0.36480110528719367  -3.0448540722552333  -2.2138486771053065 
H  0.5589051533128724  -3.4238219370250493  -3.0961083135003094 
H  -0.3854522723923228  -2.4196296416946645  -2.2863370948745967 
H  0.13056461238611128  -3.7657425569970977  -1.5922379565092042 
C  2.8112318770792513  -3.2488249760417265  -1.6898294033780852 
H  2.9430968151272827  -3.6131130518951218  -2.5903122937789225 
H  2.6512891890923767  -3.9829678276568146  -1.059999053218616 
H  3.612059209260522  -2.753306551111817  -1.4170210739846474 
C  1.1061098568322651  -2.8371384692435946  1.3704285475452402 
C  0.14069112484651392  -3.978778995179862  1.0093196422673234 
H  0.0913946958811398  -4.613268848267248  1.7550464343436405 
H  0.4656842502602556  -4.4407072614791625  0.20758125274893585 
H  -0.7499007178131105  -3.6105091096122095  0.8354773540513868 
C  2.4556050639674023  -3.4502268425087754  1.764603017147714 
H  2.3191525129614945  -4.109783087754531  2.477077796126821 
H  3.0540333183874386  -2.7433632681128937  2.0848908275683944 
H  2.8552319925818606  -3.888908131068714  0.9844452542602429 
C  0.5480275792231377  -2.0805707390012063  2.587273180591934 
H  0.4390447191413469  -2.702676611927514  3.336669710527184 
H  -0.3202931713705248  -1.689759023173501  2.3585480685900584 
H  1.1700836364608254  -1.367338444341118  2.842487983273555 
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

-Ni 0
lanl2dz

