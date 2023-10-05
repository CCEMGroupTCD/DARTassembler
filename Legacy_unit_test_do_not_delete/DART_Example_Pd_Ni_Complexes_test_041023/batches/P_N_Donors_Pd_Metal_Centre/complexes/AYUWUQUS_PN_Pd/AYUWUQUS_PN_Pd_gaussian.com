%chk=AYUWUQUS_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.4667760352214645  -1.5299242015537895  2.5902770669041676e-16 
N  1.4667760352214645  1.5299242015537895  -1.8736167763709184e-16 
C  1.425728308527328  -2.8928886921358985  -1.213883825390854 
C  0.524840012291349  -2.847093093289339  -2.282533132147474 
H  -0.06586481088430096  -2.1331458376732266  -2.364371192212821 
C  0.5133250603116637  -3.869718229496124  -3.2198150250052224 
H  -0.07949448453060737  -3.832199543967975  -3.936147291128436 
C  1.3627759397067305  -4.9303801226487085  -3.102619667123866 
H  1.3239577524062007  -5.625766808420146  -3.718186499877632 
C  2.274906009202317  -4.972603162152492  -2.07673482404589 
H  2.864713321511073  -5.689440815865131  -2.0122947544704597 
C  2.3253315152015808  -3.954945124964262  -1.138039880152529 
H  2.958666378890431  -3.980877622116543  -0.4572928522282873 
C  1.5810171400201192  -2.258164952872363  1.6708040543682598 
C  1.3997439333228345  -3.590748225906051  1.9546351572525449 
H  1.2001853822516249  -4.185759144796984  1.2680963520165067 
C  1.512395252429461  -4.053415424896112  3.2601937860709964 
H  1.4161183417738918  -4.961834798606594  3.435657337676688 
C  1.7608658250944973  -3.197032482087991  4.278827365817326 
H  1.8062310730419966  -3.5118476699978523  5.152267055629945 
C  1.9442392347445905  -1.875517976371535  4.017890182161381 
H  2.1364078801990276  -1.29031775266717  4.714474362200354 
C  1.8486589176356372  -1.3904951633155835  2.7179438087302863 
H  1.9641359645824346  -0.4825253608971326  2.5512690370200755 
C  3.1099938543071968  -0.7853266039453248  -0.32045248739071996 
H  3.4354415385323405  -0.3497960303892775  0.4839938901450744 
H  3.740949248826194  -1.478388024192896  -0.5674560660931746 
C  2.9955428449291874  0.23973432286874258  -1.4534467339654618 
H  3.8487065701762386  0.3052786517303667  -1.9115038252507093 
H  2.3351121408871904  -0.06607609213406221  -2.095458949763911 
C  2.596338114207585  1.6156010180044837  -0.9387000387980547 
H  3.3542433375748306  2.025539153738949  -0.4917514430253268 
H  2.3487851760717646  2.1806857484257547  -1.6875450496414643 
C  1.4252663193131738  2.27222366243337  1.0476344657671002 
C  0.31942900468299973  2.1310163759674166  2.040190260581945 
H  -0.25405328820934425  2.8999838377574756  1.9923160248480516 
H  0.6891222592965016  2.0625372415804986  2.9244031962750143 
H  -0.18926330184479245  1.341041500890243  1.8440232806229344 
C  2.465344262904153  3.3098777294699606  1.3928884424538484 
H  3.236485689561425  2.8769984539297755  1.76730212129161 
H  2.1015083075303553  3.9260028349115714  2.0332779114661204 
H  2.7181296732086366  3.7865597632983916  0.5995175928744088 
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

-Pd 0
lanl2dz


