%chk=ICARIPOR_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Ni  0.0  0.0  0.0 
P  1.382708702149517  1.4724865517212702  0.0 
N  1.382708702149517  -1.4724865517212702  0.0 
H  1.1642805901279032  -1.7692414385412583  -0.8874016548480969 
H  1.4263315748552876  -2.124045032088028  0.5600179565327013 
C  2.7342451842745406  -0.9243858413571802  -0.0895699100342032 
C  3.8415654691414725  -1.7656784891928472  -0.13573837784348952 
H  3.7302942788167153  -2.706217588433568  -0.06735895840192807 
C  5.104306237606481  -1.2218867621553489  -0.2801088708882896 
H  5.861837046460866  -1.7938538942255808  -0.3293512659134502 
C  5.2739289967720255  0.15568961885051946  -0.3544589486814249 
H  6.145298747949579  0.5226302468436215  -0.45737445898417906 
C  4.179011697253692  0.9863313612625226  -0.2799015011775849 
H  4.300090524057432  1.9274907989361705  -0.31076378524175385 
C  2.8906778226363787  0.4586189321801226  -0.1588157476577818 
C  1.9791893939266163  2.6990492309448966  1.2367754033262959 
C  2.1394511770747964  4.058855925042987  0.9761622310801727 
H  1.9111557861703468  4.406920659081397  0.1214442627358932 
C  2.6328254120183225  4.907803681841671  1.962437157518212 
H  2.7350841114814903  5.833946741284798  1.7771042898734721 
C  2.9746778716544044  4.421353986146528  3.1989179901229665 
H  3.310861693551587  5.0080344783996935  3.8665145583467826 
C  2.8282098661803827  3.0764731906750225  3.4683430434205156 
H  3.0668125504565964  2.737001601331353  4.323227046712598 
C  2.3322441955806013  2.218429693759499  2.4976774722851625 
H  2.2331488324354245  1.2948135843898039  2.6931896241221134 
C  1.3279157143499645  2.4633953062320506  -1.5422349946496259 
C  2.1802674943032323  2.2286777601899574  -2.6201286474700085 
H  2.850299199463546  1.5573790791632878  -2.5519300776036484 
C  2.0599208513591973  2.9634076792750106  -3.7877283226726965 
H  2.6543061071773106  2.8026784354504475  -4.510248251702791 
C  1.0887060088776446  3.920089479426795  -3.904624357686207 
H  0.9994759453129981  4.408541329230951  -4.714616355288143 
C  0.24177065908535567  4.174308704425774  -2.8508791375721834 
H  -0.42551944068726844  4.845076772796535  -2.931176965684712 
C  0.36112191716012343  3.4488260922068408  -1.6680641239018306 
H  -0.22294561368262933  3.631626308402602  -0.942634430289791 
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
