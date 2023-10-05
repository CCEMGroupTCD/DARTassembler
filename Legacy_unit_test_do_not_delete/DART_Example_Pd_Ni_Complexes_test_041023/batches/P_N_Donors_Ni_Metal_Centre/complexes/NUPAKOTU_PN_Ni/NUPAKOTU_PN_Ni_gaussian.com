%chk=NUPAKOTU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
C  2.8421262813342887  -0.7294581554360186  -0.48178673482717227 
H  3.257519119156756  -0.32044875530931555  0.3183080601151555 
H  3.4726819656799517  -1.4009563643084868  -0.8443888506792825 
P  1.2864440631057361  -1.557293059285887  -1.9146603087536327e-15 
N  2.604154620544173  0.30955105696101676  -1.4928056244701486 
C  2.187166048120188  1.5289436489184876  -1.0185588773573626 
N  1.2864440631057386  1.5572930592858873  -1.3009364206640417e-15 
C  1.051763565343032  2.7437506396254943  0.6080261211695918 
H  0.4913468897341401  2.7532851809319387  1.3745164421923468 
C  1.5828793198617017  3.933951404310041  0.1731322892446107 
H  1.4068701258580563  4.743463732437981  0.638279668292452 
C  2.379778785951024  3.9235594039710575  -0.9571216013815369 
H  2.7045295458229726  4.7401391366150385  -1.318244481019037 
C  2.699615888232536  2.7301998471948945  -1.5550020127140858 
H  3.260020522265817  2.7122335563717477  -2.32212908729042 
C  3.492310799166054  0.2688295693631212  -2.6550942559311874 
H  3.0888391485949973  0.7738796884207779  -3.3923931060763293 
H  3.6267073546522655  -0.662334498538751  -2.9314274868714536 
H  4.356149776835953  0.6674864928613132  -2.4204466520541494 
C  1.7639865416636367  -2.513550782339553  1.495924482810893 
H  2.7126417006754147  -2.8014783835873316  1.3654189054585613 
C  1.740681246605381  -1.6028659619433654  2.740981277146888 
H  2.3071840336144476  -0.8070988758380256  2.5776365588901253 
H  0.8157785284236552  -1.291306472798472  2.9043972697052487 
C  2.249321507587923  -2.3453697110644405  3.9652395372574127 
H  2.176060061856969  -1.7588567297384754  4.759662487143376 
H  3.205426051875233  -2.5713321963537297  3.838936466688827 
C  1.455240570761879  -3.625252176420049  4.206055381098191 
H  1.838954709824593  -4.111402432780381  4.977714468859217 
H  0.5162678782488828  -3.3967119617971164  4.423082621358564 
C  1.489342927646774  -4.518542611836175  2.9733658349232335 
H  0.9467425993287368  -5.330064227799401  3.139223968680955 
H  2.4215299335999827  -4.803859230535899  2.7981697361488638 
C  0.9447134749793459  -3.7844144570980727  1.7523965717114234 
H  -0.003927571495137228  -3.543720386855517  1.9035555591844615 
H  0.9914877255002035  -4.375565438289851  0.9598195490037099 
C  0.8811041986117115  -2.771124867453587  -1.3108577137631459 
H  0.1482927685601816  -3.3504195060707276  -0.9550191131936032 
C  2.066798025688727  -3.692977367495103  -1.63148136635084 
H  2.80708203346335  -3.1611654197941537  -2.0171082976142762 
H  2.393211010676194  -4.117243337256104  -0.7984371597601302 
C  1.6392246313890912  -4.77126440215897  -2.6235059191532923 
H  2.413313406382124  -5.350205969787499  -2.83680816879361 
H  0.9372905859374941  -5.3369235404613935  -2.214904260944573 
C  1.1050351575592912  -4.145048572591739  -3.9050418144168293 
H  0.7874649300204108  -4.859593865377071  -4.5122139688666625 
H  1.8359425231292292  -3.6587726237685576  -4.362170433682516 
C  -0.040794556012786254  -3.1791086995341997  -3.6171682151880495 
H  -0.8171074304352832  -3.68882497950272  -3.273820264937943 
H  -0.31507013932080974  -2.7386181727771097  -4.460482041849449 
C  0.3490373196575449  -2.1122595184546618  -2.5979030075304244 
H  -0.44008634361527443  -1.554256061867963  -2.382722556512726 
H  1.0466064848963732  -1.5241599154037306  -2.9835755715006247 
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

