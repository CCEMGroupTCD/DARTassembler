%chk=UFUVUBUR_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Pd  0.0  0.0  0.0 
P  1.5906298838824824  -1.40071288010784  7.159128449528954e-17 
N  1.5906298838824826  1.40071288010784  5.050675041017767e-17 
C  2.1147902219249772  0.9173410844919722  -1.0992611096669944 
C  3.2352423698554795  1.416994750092674  -1.9477999824652716 
C  3.8180711656385293  2.683171220079483  -2.0239983529060566 
C  4.901469461176502  2.827666760796162  -2.849824049974448 
C  5.468828156562674  1.744115580140816  -3.5208711270637263 
C  4.920285650580356  0.4873316453315158  -3.377218767751722 
C  3.7705818107697935  0.3248288524804387  -2.596741539999088 
C  3.088901280024344  -0.9932132646577031  -2.2989653296453127 
C  3.8843510096185723  -1.9042909474028986  -1.328491376521308 
C  3.424437105468161  -1.6135105699363597  0.13612203078097537 
C  1.8274618490608685  -0.5139774380061197  -1.5726561912860877 
C  1.884341302252049  2.777521651832367  0.34167794605701646 
C  2.83860432897185  3.031356809664731  1.3171659308542751 
C  3.121079097733965  4.335247889654935  1.605133561908136 
C  2.4771857669464095  5.341386945712946  0.9886390619797034 
C  1.487124309228897  5.102904920064491  0.005870715224003752 
C  1.1757298321825895  3.756537147780704  -0.3256689376392017 
C  3.5654920351160695  1.875661322038909  1.9998620730271561 
C  0.1449398501581376  3.4668598372038484  -1.353808752399344 
C  1.0107361225468545  -3.0625535064471845  -0.27710737668291396 
C  0.8866596864081081  -3.947907601930711  0.8036181868186095 
C  0.47286641505882954  -5.267246969223258  0.6118431297904289 
C  0.1726891115785929  -5.697616132355098  -0.6532283345547941 
C  0.30216740881900694  -4.876784728340677  -1.7411408433198459 
C  0.731073983872378  -3.5551864898369536  -1.5538974032215587 
H  3.4763812554035765  3.4163860991494395  -1.5226266281806973 
H  5.278451848570873  3.6927405901356267  -2.9690729467675774 
H  6.229280782963879  1.8711274771016722  -4.075380704104511 
H  5.319370028042541  -0.260794012597269  -3.8047231532429837 
H  2.8580881143637207  -1.4794332028562907  -3.142478307764752 
H  4.853949380645406  -1.7259145124116213  -1.4178912442053266 
H  3.7217604998722615  -2.8539837780261057  -1.5499782239938298 
H  3.8513008889454596  -0.7912712185631284  0.484270346027198 
H  3.6488963228211913  -2.369471962908256  0.7326510259532608 
H  1.0162643340033584  -0.5594736855990552  -2.155910383743808 
H  3.7850934831988243  4.536019735187153  2.254588019077651 
H  2.6920773290033995  6.2370173763429655  1.2196813704130116 
H  1.0436557828461313  5.825104850804127  -0.4227771427521763 
H  4.514963396791323  2.0978345342792655  2.0957858298413528 
H  3.4729074870591736  1.0653597312427883  1.4567412914192555 
H  3.175222123425167  1.7224589362769693  2.8857595134939773 
H  -0.4239610719887148  4.256692347286018  -1.4758805214443076 
H  -0.4042810984835068  2.709992566219391  -1.0624599173480929 
H  0.5831368009302293  3.245509288206162  -2.2016206577135384 
H  1.0881954093699893  -3.644844668028569  1.6815213252484054 
H  0.3999803133254689  -5.861319392499565  1.3512673058062992 
H  -0.1325254166901102  -6.589381091997529  -0.7794684964988196 
H  0.10488463150293881  -5.1988206976003255  -2.6127220846722943 
H  0.831274953484218  -2.983831755854366  -2.306003242611554 
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


