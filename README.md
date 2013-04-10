C++ version of low rank representation
====
**It's the C++ version of low rank representation.**
          
Just need install [Eigen library](http://eigen.tuxfamily.org/index.php?title=Main_Page')    

Install eigen3
#yaourt eigen3
run makefile.sh
$sh makefile.sh
run lrr
$./lrr       


**Result:**      
         
         initial, rank=0
         iter 1, mu=1.0e-06, stopALM=2.232e+00
         iter 50, mu=1.1e-04, stopALM=2.232e+00
         iter 100, mu=1.3e-02, stopALM=1.203e-01
         iter 150, mu=1.5e+00, stopALM=8.300e-06
         iter 186, mu=4.5e+01, stopALM=8.415e-09
         LRR done.
         Z = 
         0.0658979 0.0778617   0.08988  0.101959
         0.0678977 0.0802245 0.0926075  0.105053
         0.0698975 0.0825873  0.095335  0.108147
         0.0718972 0.0849501 0.0980626  0.111241
         E = 
           -3.00606   -2.73336   -2.46397   -2.19828
           -0.10842 -0.0358526  0.0324873  0.0961225
            2.78922    2.66165    2.52895    2.39052

