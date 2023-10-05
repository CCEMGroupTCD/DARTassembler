%chk=BAXIKINO_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
P  1.3793296750595954  1.4756522786551038  3.81417668715752e-15 
N  1.3793296750595914  -1.4756522786551038  0.0 
C  -0.1188271411809041  -3.3516814964075845  -0.6341193011501962 
H  -0.7741195385436865  -2.9265978848594116  -0.04188021332189251 
H  -0.10336941635077079  -4.316673188961208  -0.46487778128710644 
H  -0.3680902346775905  -3.187510926887785  -1.567692717826136 
C  1.2408867591813357  -2.775460555133428  -0.3708215115680071 
C  2.3478240715348777  -3.6104242303990803  -0.5295200544296843 
H  2.2204560123910664  -4.51624116595671  -0.7855111221779998 
C  3.621435932866966  -3.133200862681438  -0.3188907497770441 
H  4.37613478599528  -3.706441718427582  -0.3931802671413813 
C  3.772781472971573  -1.797040085837047  0.0031254546491275015 
H  4.640437929758327  -1.4302571432502706  0.12635205901537083 
C  2.648097740483731  -0.9931165381560205  0.14604473732705492 
C  2.830784222901423  0.4583919343901868  0.4868256541466838 
H  3.6370249576836153  0.8012749680958933  0.02599353419992534 
H  2.977268554133829  0.5471928510393502  1.4619496810815478 
C  1.6765751381049632  1.969317063738434  -1.8020946769433768 
C  2.8195651424844286  2.954819296232798  -2.0673039332623437 
H  2.6066228903520354  3.8196497501681286  -1.6580708661687937 
H  2.9338591143587847  3.0700936074377587  -3.0336127457904194 
H  3.6485793114797573  2.6053074820606543  -1.6792992047356403 
C  0.3620801360442263  2.548334910530037  -2.3632475550917027 
H  -0.3746077811645483  1.9306317610924555  -2.172760306664118 
H  0.445235153808097  2.669308772450065  -3.3317858649963146 
H  0.17965251322155829  3.413799924526396  -1.9402611345300633 
C  1.975152586332576  0.6599511054105254  -2.551055927330816 
H  2.839302865766318  0.3045204273391193  -2.2554337230898382 
H  2.006113848424439  0.8344753505360198  -3.514865777930546 
H  1.2704588340353749  0.006210678919799095  -2.3599457228805294 
C  1.5661696511865952  2.994059895686741  1.1225623331740573 
C  3.025113008159848  3.476787117610062  1.2229069411897895 
H  3.5779530202549736  2.766834813700897  1.6113896012928346 
H  3.067776207316325  4.272919680134017  1.7921386113617697 
H  3.3585738158400993  3.6968318765934853  0.32764931693564425 
C  0.6773876186302856  4.1392030465786975  0.6054791458914568 
H  1.0222562633595804  4.458183379417672  -0.254487920422439 
H  0.6829297228672186  4.875192341421129  1.2528685518582001 
H  -0.2400827493421429  3.8142013460174393  0.48914979222329796 
C  1.1105758088972622  2.6493852915424974  2.55251932487907 
H  0.16239178197766724  2.4025866247644956  2.5432156067435967 
H  1.240244065268571  3.4277412237725526  3.1334283849322926 
H  1.6403286749875465  1.8972454439277524  2.8913672708654845 
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


