%chk=FEBEKALI_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.4667760352214645  1.5299242015537895  -7.166602905332491e-17 
N  1.4667760352214645  -1.5299242015537895  0.0 
C  1.425728308527328  2.8928886921358985  1.2138838253908544 
C  0.524840012291349  2.8470930932893386  2.2825331321474747 
H  -0.06586481088430096  2.133145837673226  2.3643711922128214 
C  0.5133250603116637  3.8697182294961237  3.219815025005223 
H  -0.07949448453060737  3.8321995439679744  3.9361472911284365 
C  1.3627759397067305  4.9303801226487085  3.1026196671238666 
H  1.3239577524062007  5.6257668084201455  3.718186499877633 
C  2.274906009202317  4.972603162152492  2.0767348240458903 
H  2.864713321511073  5.689440815865131  2.0122947544704606 
C  2.3253315152015808  3.954945124964262  1.1380398801525295 
H  2.958666378890431  3.980877622116543  0.4572928522282878 
C  1.5810171400201192  2.258164952872363  -1.6708040543682596 
C  1.3997439333228345  3.5907482259060512  -1.9546351572525444 
H  1.2001853822516249  4.185759144796984  -1.2680963520165063 
C  1.512395252429461  4.053415424896112  -3.260193786070996 
H  1.4161183417738918  4.961834798606594  -3.4356573376766875 
C  1.7608658250944973  3.1970324820879914  -4.278827365817326 
H  1.8062310730419966  3.5118476699978527  -5.152267055629945 
C  1.9442392347445905  1.8755179763715355  -4.017890182161381 
H  2.1364078801990276  1.2903177526671707  -4.714474362200354 
C  1.8486589176356372  1.3904951633155838  -2.7179438087302863 
H  1.9641359645824346  0.48252536089713294  -2.5512690370200755 
C  3.1099938543071968  0.7853266039453248  0.3204524873907201 
H  3.4354415385323405  0.34979603038927753  -0.4839938901450743 
H  3.740949248826194  1.478388024192896  0.5674560660931748 
C  2.9955428449291874  -0.23973432286874274  1.4534467339654618 
H  3.8487065701762386  -0.30527865173036695  1.9115038252507093 
H  2.3351121408871904  0.06607609213406196  2.095458949763911 
C  2.596338114207585  -1.615601018004484  0.9387000387980545 
H  3.3542433375748306  -2.025539153738949  0.4917514430253266 
H  2.3487851760717646  -2.1806857484257547  1.6875450496414641 
C  1.4252663193131738  -2.27222366243337  -1.0476344657671004 
C  0.31942900468299973  -2.131016375967416  -2.0401902605819453 
H  -0.25405328820934425  -2.899983837757475  -1.992316024848052 
H  0.6891222592965016  -2.062537241580498  -2.9244031962750148 
H  -0.18926330184479245  -1.3410415008902428  -1.8440232806229346 
C  2.465344262904153  -3.3098777294699606  -1.3928884424538488 
H  3.236485689561425  -2.8769984539297755  -1.7673021212916105 
H  2.1015083075303553  -3.926002834911571  -2.033277911466121 
H  2.7181296732086366  -3.7865597632983916  -0.5995175928744092 
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

