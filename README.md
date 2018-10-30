### Requirements

I assume you have GNU Make and a C++11 compiler.

#### qmake

`sudo apt install qt4-qmake`

Later versions also work.

#### NTL

Needed for calculating in Finite Fields.  
Tested with [ntl-11.1.0](http://www.shoup.net/ntl/ntl-11.1.0.tar.gz) https-rehosted [here](https://tasty.limo//random/ntl-11.1.0.tar.gz), [other versions](http://www.shoup.net/ntl/download.html).  
Build it with the following commands, extracted from [here](http://www.shoup.net/ntl/doc/tour-unix.html):  

`gunzip ntl-11.1.0.tar.gz`  
`tar xf ntl-11.1.0.tar.gz`  
`cd ntl-11.1.0/src`  
`./configure`  
`make`  
`make check`  
`sudo make install`  

If you dont have [GNU GMP](https://gmplib.org/) see the linked build instructions above.

#### Build the Project

`git clone https://gitlab.com/ggwpez/Niederreiter/`  
`mkdir Niederreiter/bin`  
`cd Niederreiter/bin`  
`qmake ../`  
`make` or `make -f Makefile.debug`  

Usage will follow...
