%chk=FAFAGIWO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5987812373805217  1.391401651213624  -1.79385735910041e-15 
N  1.5987812373805241  -1.391401651213624  1.6653345369377348e-16 
N  1.9634508628760396  -2.6245125279352854  0.2047343155206251 
N  2.574194730511592  -3.0100628825234725  -0.9213033136119586 
C  1.1296066863601024  2.9162986054103857  -0.8725072376526448 
C  1.1126314045203518  2.6064283715897383  -2.189935902893555 
C  0.8084058251665549  3.3705772024270533  -3.450781484714306 
H  1.5920773693744488  3.9017067503583847  -3.7423326879584073 
H  0.03646156871643513  3.9773047560875354  -3.3231927682740885 
C  0.48104869869905875  2.2614616184014404  -4.4690008240623875 
H  -0.4862149040462913  2.0495502664616034  -4.453424183758289 
H  0.7305500909643297  2.5433673964044146  -5.384531362922917 
C  1.3108690320194427  1.0409297802232031  -4.031965417714917 
H  0.8453436902182794  0.19524757985089547  -4.253984903603806 
H  2.2009657072940323  1.0429696081891429  -4.4650556251394296 
C  1.429928419826489  1.2172638422002011  -2.5541389154617944 
C  1.7459684886545617  0.43719587688386596  -1.5131680684240525 
C  1.9764002759260442  -0.9650313292562323  -1.2551170583981188 
C  2.6026358887394325  -2.0191856805326247  -1.850825024543702 
H  2.977483944079766  -2.0537583037274025  -2.723463190757543 
C  1.0419302942279927  4.217669175454359  -0.21474319205108652 
C  0.874581646163428  4.309306076193112  1.1666897011485802 
H  0.7477713600359499  3.5189157679447614  1.678127138142931 
C  0.88795987370505  5.539735972913961  1.8018103214116692 
H  0.7723245621663843  5.585324890351764  2.74379850357335 
C  1.0692794650010384  6.697814546843011  1.0797353167073642 
H  1.089688820202567  7.538177061132194  1.5212902881410768 
C  1.2197254345797495  6.629081701142718  -0.29329437931508484 
H  1.3252966612309112  7.426886105091617  -0.7974226756277589 
C  1.2182390560643457  5.409989871317391  -0.934761822863126 
H  1.3375981586111743  5.375703538569121  -1.8771561575372602 
C  3.1685368957759836  -4.307699624168242  -0.9764135424482292 
C  3.1533479412856185  -5.018971376028519  -2.1495928043972103 
H  2.761233855980702  -4.648214345123721  -2.9314128208589083 
C  3.7132081194447526  -6.280429031673146  -2.1781281145794744 
H  3.701998415147413  -6.7825021996837425  -2.9837789176532397 
C  4.289263127355053  -6.817062630137492  -1.0423927673877487 
H  4.668438583631296  -7.6877527999696955  -1.065091699417053 
C  4.311247668250852  -6.081798571691344  0.12266250354104419 
H  4.709285232339515  -6.447452856028461  0.9041920435587558 
C  3.7581273207221084  -4.814482625600613  0.16619598414915066 
H  3.7832707278887003  -4.302702910966172  0.9663938176658162 
C  3.2522361683355587  1.6365131613965587  0.6720239273094248 
C  4.324399055727852  1.8612586256032775  -0.18431655358423688 
H  4.189906518359448  1.8830105510617536  -1.1243553565764801 
C  5.583515583482438  2.049844344019759  0.3411827388459429 
H  6.317323815154901  2.2159485417533267  -0.2381925839156141 
C  5.782994748259778  2.001040808433345  1.708138768765732 
H  6.654716996112287  2.1319841955143386  2.0631044177780433 
C  4.729265119959269  1.7635478496519794  2.5577378594967577 
H  4.873394452132078  1.7209503294116506  3.4952069613218715 
C  3.4574520351497284  1.5835024725450224  2.042674562911099 
H  2.7251099915942456  1.4248867846827133  2.6266523320127853 
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


