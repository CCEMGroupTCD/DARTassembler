%chk=GIWIWIBU_PN_Ni_gaussain.chk
%nprocshared=32
%mem=64GB
#p opt rwb97xd/gen pseudo=read

This Gaussian input file was generated using the DART python program. This file was generated on the 2023-10-04 at 09:55:46

0 1
Ni  0.0  0.0  0.0 
N  1.394792696962528  1.4610452876280045  2.1835932222132284e-16 
P  1.3947926969625275  -1.4610452876280045  -4.440892098500626e-16 
C  2.634640420548883  -1.0513278715636152  -1.2700790689131634 
C  3.0357397154322023  -2.033585698992024  -2.1567110817007342 
H  2.806171861603793  -2.9187462590042306  -1.9903210775313256 
C  3.778887350928243  -1.7223680975594182  -3.2973528337049296 
C  4.028400967864501  -0.38348267999934427  -3.50997334919766 
H  4.472898609164898  -0.15600282512120955  -4.294564078497341 
C  3.6719566750257613  0.6661698482188405  -2.650635018789935 
C  2.9996586002222188  0.3029091508987367  -1.4547039737885963 
C  3.940992182976625  2.1067248452627636  -3.204258587610323 
C  3.5643849505925673  3.275725176370135  -2.295331369189123 
H  2.640051513631554  3.2026632761649028  -2.0466117759635867 
H  3.706268333164456  4.101893331608953  -2.7619225553673843 
H  4.110377317360722  3.2562100181190923  -1.5056104475786052 
C  3.104706456743208  2.2468892650340835  -4.487979530622654 
H  2.173870045637421  2.1507980136616305  -4.275706340677672 
H  3.3625130724456467  1.5663195295029495  -5.113561982405898 
H  3.256296523491809  3.1123131370434054  -4.875565280385896 
C  5.431982086471491  2.261673303682924  -3.5838936135361235 
H  5.586417514839797  3.1477324326930374  -3.9203247927727434 
H  5.658415633554462  1.618105586814587  -4.260696873398455 
H  5.977341478048193  2.1148684531993225  -2.807665786224061 
C  4.167554888304525  -2.794581808705418  -4.303894647550937 
C  5.6312766786365005  -2.642088931252104  -4.702948482056355 
H  6.1881841228228485  -2.781637305779912  -3.9345101372449234 
H  5.778154630284235  -1.7592370163096724  -5.049156589416972 
H  5.848386792894391  -3.290373267829876  -5.376988211771772 
C  3.331346185204901  -2.5774514139898446  -5.569773427905183 
H  2.399481128738245  -2.6750031505492164  -5.358458815938453 
H  3.578508002055907  -3.225970791194286  -6.232086694933191 
H  3.4907580422636544  -1.6937729504043095  -5.911506266280748 
C  3.9046735039404163  -4.207334809915287  -3.8040432443045913 
H  4.420418917053696  -4.366536329518849  -3.010020898895773 
H  4.157197688984188  -4.83916639844784  -4.481867240413031 
H  2.97167969216132  -4.307702800724164  -3.6036049655180658 
C  2.7082543110963915  1.2537969206101738  -0.3334867202233921 
C  3.7304336134574845  1.8200243893993786  0.42691377781246254 
C  3.3673075205267926  2.744431162712126  1.4180856943266544 
H  4.028132374143087  3.1669103460933057  1.9182755607300117 
C  2.046937892368753  3.031023579487788  1.6565306810311913 
H  1.8124150616174255  3.677494341207587  2.283699116220503 
C  1.086345931634741  2.363452640202616  0.9693940235931566 
H  0.19290959463636548  2.5277770331786726  1.1674747283677338 
C  5.185268510165702  1.4727149361627738  0.23435891027457445 
H  5.269713939450633  0.830196652883513  -0.4745724228774939 
H  5.53602349002618  1.1003373277429298  1.0466449380078 
H  5.6759885658123315  2.2666879463240797  0.006894130316080105 
C  0.9589266224328836  -3.193424558437053  -0.31540899978478154 
C  1.3908157706942141  -4.2268565867154475  0.5144258371210473 
H  1.8647968921190528  -4.031448478960275  1.289762730320652 
C  1.1172305116970267  -5.545996627288039  0.18578976268979205 
H  1.4031075501733805  -6.233349845223122  0.7425239206739058 
C  0.4243732605755728  -5.834124023861378  -0.9544928570638983 
H  0.23756102626287134  -6.720183474558984  -1.1677877209039567 
C  -0.001870979234008141  -4.827394696266091  -1.796061439526862 
H  -0.4598517486049214  -5.03342683472821  -2.5788494202771846 
C  0.25515426282229003  -3.514638810982173  -1.4676670263613159 
H  -0.04636516857176587  -2.834377125450454  -2.024476761988623 
C  2.19580199321239  -1.4642325365376951  1.6277344865689096 
C  1.3963085176678198  -1.4841307493078126  2.754389098102668 
H  0.47339118715458917  -1.4098365806607713  2.669522864868279 
C  1.969605062068998  -1.6139630218638281  4.009653987749717 
H  1.4292401101385792  -1.6384670261566203  4.766912815693539 
C  3.314298564326217  -1.706311343399466  4.135445134086922 
H  3.690504094640341  -1.823206939618304  4.977945412155554 
C  4.134754721000025  -1.6288685830975629  3.0302847322207933 
H  5.059886978486042  -1.6451340850647083  3.1294550887019232 
C  3.5642887733972466  -1.5258367563125153  1.7593925002455462 
H  4.106796773499741  -1.4994410892332204  1.0052140872793287 
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


