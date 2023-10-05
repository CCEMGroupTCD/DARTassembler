%chk=WIFUKIDI_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Ni  0.0  0.0  0.0 
N  1.387300615403886  1.468161095554571  2.6559142876140813e-15 
C  1.4857275655250504  2.282179767360179  -1.1116520082479555 
C  0.22386930068386568  2.5753464071133036  -1.8873131912317782 
H  -0.49076513740100713  1.9366961201392965  -1.5999733213022154 
C  -0.2158450967186003  4.009045969357125  -1.5374643325243975 
H  -0.30745119309475655  4.094763814204706  -0.5650530027974083 
H  -1.07645300635182  4.199277164588806  -1.9641936357570784 
H  0.457122847573046  4.6444301774147805  -1.859833847377818 
C  0.44135675566536625  2.4348925763828904  -3.392203877746541 
H  1.1242850923234888  3.0749953844776767  -3.68484217335831 
H  -0.3983737133354375  2.6168188144563143  -3.8624678816126337 
H  0.7387728248839167  1.5237632772546423  -3.5949808291997476 
N  2.583240746607911  2.8567341316967223  -1.5266001493104648 
C  3.7259472865383687  2.6466930906054396  -0.8114591389793127 
C  4.936404211061494  3.2009522434081665  -1.2753944743532148 
H  4.9507976157993365  3.703924275725277  -2.081232190402086 
C  6.082959292503054  3.015227692590395  -0.5676390600721997 
H  6.88949845140587  3.416111248068127  -0.8705366832123078 
C  6.096117051472644  2.2379207073723135  0.6059385871715173 
H  6.91113951646952  2.1035075386829067  1.0744073818988584 
C  4.942005135004147  1.6764386687213577  1.0741106650850285 
H  4.954464404207023  1.151279037169646  1.865623076863523 
C  3.721538736319353  1.8824043498142624  0.3719475516325034 
C  2.480518663924903  1.290121341748628  0.735484727545491 
C  2.394843446138271  0.3436097230933932  1.8815079789709892 
C  2.6846699533532385  0.7892881108270215  3.200659941567698 
C  2.9431702262396744  2.1546189393508666  3.506982251561663 
H  2.8907686362468032  2.81130162759212  2.8224875325677896 
C  3.2659284776850726  2.5245317858664302  4.779403678254369 
H  3.448777622081047  3.4374619738807874  4.967324944229596 
C  3.3330600530835373  1.576973316831728  5.822658996642513 
H  3.5960334665462597  1.8461538558043398  6.694404487740793 
C  3.018024780398365  0.2683951261841291  5.572507870469476 
H  3.0200291684131524  -0.36085069894859445  6.2855585096121365 
C  2.6890990948138596  -0.15749203081332186  4.262583119627888 
C  2.394889041289847  -1.5007074490471697  3.978212807140309 
H  2.404286332069555  -2.1390039396962557  4.683237211647768 
C  2.09511073172994  -1.9088288853216187  2.7091804717088133 
H  1.898587624489518  -2.8245371301938715  2.5474524327249943 
C  2.0732206214259286  -0.9907275890256239  1.6330894287154853 
P  1.3873006154038856  -1.468161095554571  5.551115123125783e-17 
C  2.663136985516785  -1.1125987757378784  -1.2312063935115718 
C  4.016254190681402  -1.0665615661941472  -0.8902412128480173 
H  4.285528217152803  -1.2238692354269718  0.006994679138268872 
C  4.975261762670029  -0.7882919939809716  -1.8667744749793276 
C  6.440139296945035  -0.6996970057376221  -1.5197732103784118 
H  6.955440249094558  -1.2691073841567695  -2.129090967925616 
H  6.575738787951239  -1.0018975451357925  -0.5971999005399349 
H  6.7406332365227755  0.22948073172633365  -1.6079184568662557 
C  4.551883027102677  -0.5841154660048163  -3.169142491974927 
H  5.203129089526826  -0.42143379181852336  -3.841133637051257 
C  3.210074710008125  -0.6089945546215888  -3.5290318098568134 
C  2.7926551432140574  -0.33252727781592006  -4.941985160152989 
H  2.821749832735102  0.6323040777411018  -5.107766716351713 
H  1.8792417319717294  -0.6603228218179987  -5.082910358379022 
H  3.4033650159281263  -0.7880202833347707  -5.558210162317768 
C  2.259397642195564  -0.8772996174602179  -2.5432437515042383 
H  1.3360856846276872  -0.8994700659883501  -2.7661874220226927 
C  1.260284087874155  -3.2771266396230154  0.06754614798376471 
C  0.3107857580711346  -3.8692744615700443  0.9021720154625388 
H  -0.3436379269417029  -3.3306887693089933  1.3335437027415153 
C  0.3189219603334188  -5.249398414680666  1.1058728163030171 
C  -0.6589927957714317  -5.88382779904476  2.0598755081481745 
H  -0.48974836313895764  -5.555858547946931  2.9684569257202407 
H  -0.5514736270787055  -6.857549929626079  2.0397156276970327 
H  -1.5734795946867692  -5.6515412563536245  1.7939176981591194 
C  1.251359529029425  -6.025074206287058  0.42296136025610875 
H  1.2507887727170255  -6.965804170581813  0.5505532820819923 
C  2.181479468280141  -5.461020358580931  -0.4392432416199164 
C  3.1714871007600896  -6.332798063272482  -1.1677082931043326 
H  3.651340022597892  -5.796761876655829  -1.8327023382214367 
H  2.69599850910093  -7.0624487557227935  -1.6180012046727668 
H  3.811399151443595  -6.707515811662479  -0.527214611299768 
C  2.180255832698586  -4.080319016515807  -0.6145586052807468 
H  2.810183506145077  -3.6806810566375887  -1.2021610961286604 
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

