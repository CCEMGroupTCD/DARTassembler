%chk=ODIGOMOC_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-13 at 16:51:13

0 1
Pd  0.0  0.0  0.0 
P  1.5709416952898032  1.4227586548673679  2.055286847288396e-15 
N  1.570941695289803  -1.4227586548673679  -7.771561172376096e-16 
N  2.8407469968662604  -1.2706913177544759  -0.4896859616854047 
C  1.6555825478134083  -2.233176033848267  1.0733667117671506 
C  2.9746202110101563  -2.6122952387497365  1.2492492906184296 
H  3.305087198432122  -3.1998779796032317  1.9187819207981414 
C  3.7161592668552546  -1.9831301115950901  0.28567943771263826 
C  0.45184803220724357  -2.6028624934600932  1.8881956554836425 
H  -0.35218121003077907  -2.5449257522440045  1.3298619276438934 
H  0.3683145999895028  -1.9867955962380126  2.645626757201219 
H  0.5523887120076476  -3.5194736316447153  2.2234997608285747 
C  5.195527343928054  -2.027513162637493  0.07730298194536998 
H  5.559044752185508  -1.118811398232817  0.12574363905916783 
H  5.390175034127285  -2.4094239653875578  -0.8040362917982529 
H  5.606509317964599  -2.5817310177195325  0.7726197415055839 
C  2.9819670094412993  -0.3643870506261857  -1.6326070715759315 
H  2.1681540423625427  -0.4917940492621544  -2.200027428377732 
C  4.186375035108499  -0.6434620776110076  -2.528921596183177 
H  4.201613871084557  -1.6007483505431996  -2.7812709166373666 
H  5.022930876840689  -0.4375569305941305  -2.041329569899224 
C  4.086834483824532  0.23003027876889837  -3.791582671846794 
H  4.872273071836267  0.06749352019155719  -4.372271544284816 
H  3.272377709039105  -0.010630404494996304  -4.300966071200474 
C  4.032700406311221  1.6934135608717171  -3.4078548238688855 
H  3.9190853688801717  2.239298393371579  -4.226839306576677 
H  4.891379739363981  1.9514137626235875  -2.988159530966928 
C  2.889145020173686  2.003439121771306  -2.4375098952868797 
H  2.0220599917862847  1.852256080807798  -2.890924299195126 
H  2.9321514344623827  2.9564116694072182  -2.171250413235355 
C  2.9732173411617273  1.113247425753776  -1.1772143866123685 
H  3.8338514269315405  1.3086430320633475  -0.7066184692415403 
C  1.2869488833962281  3.204732794561742  -0.12918474466159913 
C  0.3823584966018563  3.7369166455681415  -1.0583854231198329 
H  -0.15988518868310542  3.1595200033550728  -1.5815864870805947 
C  0.2828311394113072  5.1072070341244125  -1.2072598876897258 
H  -0.33351126143643217  5.47176464073293  -1.833192317254702 
C  1.0727868518715193  5.953776876366758  -0.4522714211835999 
H  0.9901605975277123  6.894578006214591  -0.5545114630955476 
C  1.9733720476196366  5.440848476667069  0.4404347913974949 
H  2.525192135974805  6.024878562976728  0.9471137512654246 
C  2.081974387034234  4.075001596486848  0.6032942454541622 
H  2.7090262916159995  3.724822317453436  1.2260833923538539 
C  2.139380524535293  1.1990367685790182  1.7111859703146481 
C  3.464494555388221  1.0601568824338732  2.090398816169158 
H  4.146179839857496  1.0280622753522137  1.43099401702501 
C  3.7933269682432593  0.9711112037682876  3.4257524777011628 
H  4.703022717543131  0.8749084363264958  3.67905299337829 
C  2.817493492968814  1.019962103224934  4.399556717862485 
H  3.0559924224686212  0.971403185826226  5.318573383744306 
C  1.4810826584261743  1.1399789658290258  4.0290443644817 
H  0.8029819911899515  1.1601958591872896  4.693560957439897 
C  1.143259121162912  1.2310992409062398  2.7036639537539506 
H  0.23021142467639977  1.3153852085858238  2.454667428491022 
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


