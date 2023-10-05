%chk=XIFUREGO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-08-11 at 15:55:49

0 0
Pd  0.0  0.0  0.0 
P  1.563374083512965  1.4310700454554983  3.872929588258118e-16 
N  1.563374083512964  -1.4310700454554974  -8.604228440844963e-16 
C  2.4076035846554635  -0.9457075978200874  1.0759798435965635 
C  3.1214987104706724  -1.766337649248636  1.9344558642457406 
H  3.0262397655387785  -2.689799533030407  1.8836346815294716 
C  3.977009352756177  -1.2087400372141983  2.867399080089976 
H  4.460095427899295  -1.7592441075787273  3.4382035058917744 
C  4.114787461727306  0.1661110724488215  2.954871012426616 
H  4.708667865847474  0.5338574428457308  3.569409117560859 
C  3.380779060547802  0.9839917930444779  2.141741075776949 
H  3.4561074308350292  1.9082021202748451  2.2245594775401516 
C  2.515563718078451  0.43738010028766516  1.184939143035816 
C  1.769732694527314  -2.605494791705464  -0.4783238746455567 
H  2.372163517337314  -3.1500315860256634  -0.025489287450678882 
C  1.136317039256769  -3.1581704853318056  -1.6768527668016375 
C  0.8033588291254103  -2.3723958451798746  -2.7716012288696192 
H  0.9422752699131731  -1.4529514167408764  -2.742348613749846 
C  0.2640392146184176  -2.9512269880931785  -3.9037420378689904 
H  0.040543763799975485  -2.4228953020230195  -4.637088867791732 
C  0.0583831060311677  -4.309599313174768  -3.9486374881194664 
H  -0.33963499680515374  -4.690730416248341  -4.696939220099608 
C  0.432556043413993  -5.103958218673467  -2.9053185569266655 
H  0.32119423497463795  -6.024456033705794  -2.962298642411204 
C  0.9804636079293785  -4.5424431359199655  -1.7588385995845623 
H  1.240276080790748  -5.0856815137119895  -1.0520082636201447 
C  2.7490480538176794  2.0085842192896273  -1.2303107677178944 
C  2.4919391541304323  1.8318098847255369  -2.5791730691931294 
H  1.7139664857227104  1.4017522666834312  -2.854157545799576 
C  3.389763613807303  2.291798857376053  -3.5138404938131824 
H  3.2231450893136135  2.162087690833596  -4.419119837756133 
C  4.541738367588565  2.950847331624428  -3.1084090628146335 
H  5.130039788401674  3.2870225220197193  -3.745753286630942 
C  4.814239729689532  3.106917197463112  -1.7907117589431336 
H  5.599127419532993  3.530819383052437  -1.5249637930340616 
C  3.9295327499474233  2.6400504482626257  -0.843096442220888 
H  4.121941274098519  2.7445365183015005  0.06086744588628726 
C  0.9389098350305852  2.859130468835769  0.9307992101105045 
C  0.028367317428392935  2.6438297496442553  1.9373934004147595 
H  -0.23994075245578372  1.7762655699993466  2.135868421090116 
C  -0.4963220107901256  3.705302068783255  2.6618272376664973 
H  -1.1196740426278728  3.5504366932526947  3.334966792191101 
C  -0.09209329534121102  4.970581634063842  2.3804793617372644 
H  -0.4248052323953806  5.681914236708033  2.8797783056146744 
C  0.7965050876409221  5.209417459106361  1.3752576685581888 
H  1.0501014831616664  6.081717780089369  1.1778897172798142 
C  1.322666461978493  4.159497289390444  0.6508542275262457 
H  1.9379929512766316  4.32571107169646  -0.025723979438456474 
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

-Pd 0
lanl2dz