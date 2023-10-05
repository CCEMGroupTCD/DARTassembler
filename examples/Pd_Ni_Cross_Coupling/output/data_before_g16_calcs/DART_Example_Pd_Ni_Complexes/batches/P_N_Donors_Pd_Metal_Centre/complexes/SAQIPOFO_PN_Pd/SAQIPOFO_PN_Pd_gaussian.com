%chk=SAQIPOFO_PN_Pd_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 16:41:16

0 1
Pd  0.0  0.0  0.0 
P  1.477483463528447  -1.519586330222801  2.7986545896824825e-15 
N  2.9603557155655698  -1.0219936926994464  0.7079610308632506 
N  1.4774834635284422  1.519586330222801  -1.860956535335426e-16 
N  1.3318249649532221  -3.0099822754587247  0.7118442162560158 
N  1.7847438634749553  -1.6920390134221674  -1.5920577780602405 
C  3.999369121671484  -1.8190164290269177  1.1042959898427256 
H  3.9720349905379377  -2.7675899824548797  1.1613864357353874 
C  5.0766603323868305  -1.0461695057278597  1.4009039642174437 
H  5.92340696253892  -1.352446898952213  1.7057475843340948 
C  4.710813840396343  0.27424898915673895  1.1791994025482135 
H  5.2673276228134345  1.032544156125667  1.3031909334340122 
C  3.4213257425608408  0.29631770241969757  0.7527376956675782 
C  2.714984413950017  1.4306720893450287  0.3003915203111687 
H  3.223593841292877  2.227768862269017  0.20642650015929237 
C  1.0785880302445015  2.7665760486168645  -0.6091098304285706 
C  1.5630258984203749  3.9590406779441873  -0.10468026182609663 
H  2.0942049665860814  3.964626425085488  0.6835458522880624 
C  1.268529682378472  5.141704021419056  -0.7544379818555424 
C  0.4926481858619085  5.101456801643767  -1.886629434156627 
H  0.30036132198390186  5.912597863842766  -2.3438237027062705 
C  -0.014876771426268753  3.914889843163932  -2.3818506481343387 
C  0.28287495734685253  2.745103023786999  -1.7248474630396062 
H  -0.06411774865014674  1.921130930134104  -2.0442699276004865 
C  1.8168810276530247  6.430112271783441  -0.21968896359690543 
H  1.5493340108585423  7.16833032342617  -0.8069197036243606 
H  2.7935445642928007  6.37942136764198  -0.18318446475584 
H  1.4627966306402045  6.58696145189225  0.6817540709267905 
C  -0.8357491472117742  3.8954950615619692  -3.630961930254143 
H  -0.3097177616930322  3.508242723016833  -4.3614268149273245 
H  -1.0959555409150963  4.8109867071047585  -3.8664064451982534 
H  -1.6407094656242416  3.3558810437919018  -3.486209114712025 
C  3.0341781751552315  -2.378924012477772  -2.0125343158937716 
H  3.3308339958887214  -2.9097704599328034  -1.2177281550914474 
C  4.16974399808979  -1.4264840515393276  -2.3164754822245888 
H  4.008036314332671  -0.9914465470211917  -3.1786422901432063 
H  5.012169496043123  -1.9250655042225155  -2.3544538370923216 
H  4.22312749346308  -0.7480240525941946  -1.6126086111096667 
C  2.866288599032149  -3.375808105122837  -3.130073525971179 
H  2.222113377354694  -4.063464810909858  -2.8592822433143246 
H  3.729414512392299  -3.7952952953679215  -3.3254371738346653 
H  2.536761184579269  -2.9161466162638368  -3.9300265300807733 
C  0.8359611518567862  -1.076767434963762  -2.5619672382037173 
H  0.24134003845138285  -0.4810778337896837  -2.0220218978723015 
C  1.5271183228430594  -0.17645825029281914  -3.55950484500929 
H  2.0994405201185504  0.4609644750471636  -3.0819000693569727 
H  0.8554422465353724  0.3148450209668913  -4.077262968234287 
H  2.0752919529266576  -0.7188879652792309  -4.164484221170788 
C  -0.07722145153703508  -2.074866091706498  -3.2304917510432016 
H  0.43597646723857886  -2.6071424875341744  -3.8715468770730124 
H  -0.7938619803115183  -1.5985565589049773  -3.698832814244258 
H  -0.46709003788389203  -2.6659209775900212  -2.553397239728028 
C  0.9196472941170515  -4.1696501205502425  -0.146978417775534 
H  0.8217638887944905  -3.805486279382452  -1.073595790669223 
C  -0.41749565682868206  -4.7619038392764725  0.20875929654659278 
H  -0.34921808674624755  -5.225865083845333  1.0698578721822576 
H  -0.6871827668481818  -5.398129103857688  -0.4836448740054791 
H  -1.0841909462899384  -4.047634670481791  0.27472895259790925 
C  1.9891095627738284  -5.253739660075423  -0.22868096493046394 
H  2.86889598745614  -4.8382589022908835  -0.3421704219540653 
H  1.803433056227965  -5.8386296214570095  -0.9936060645264889 
H  1.9822040188026069  -5.783886430270364  0.5952438893580512 
C  0.9681428897464018  -3.123376106117056  2.1609216632345922 
H  -0.0059238642724606105  -3.345317461663832  2.204949193408729 
C  1.1796029977087668  -1.874585775973872  2.981249978046071 
H  2.1341606836726514  -1.6538401241287302  3.00141136090436 
H  0.860995103973841  -2.026270846716912  3.895790102569577 
H  0.6808939208167645  -1.132441099774628  2.580546913153059 
C  1.7169044919582226  -4.26206388948656  2.822433065296777 
H  1.4301749550777605  -5.113119758625883  2.4293462297453505 
H  1.5234130727643478  -4.266616763424946  3.7830213652657414 
H  2.67965479606925  -4.146516423371716  2.6830858191503864 
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
-H 0
6-31g(d,p)
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
-Br 0
6-31+g(d)
****

-Pd 0
lanl2dz

