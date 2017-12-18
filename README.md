# Idhim

Download from:
https://github.com/piotrmp/idhim/raw/master/Idhim/target/Idhim-1.0-full.jar

Run as:
```
$ java -jar Idhim-1.0-full.jar <options>
```

Options:
```
 -b,--bins <file>    file with newline-separated bin names (possibly as
                     multiple tab-separated columns, defaults to
                     './bins.tsv')
 -h,--help           print this message
 -o,--out <file>     path to TSV file to which the computed moment values
                     will be written (none if undefined)
 -r,--rhos <dir>     directory containing files with mass densities named
                     'rho_<type>_<bin>.tsv' (defaults to './rho')
 -t,--types <file>   file with newline-separated <k> particle types
                     (defaults to './types.tsv')
 -v,--verbose        be extra verbose
 -W,--meanW <file>   file with mean sums of identity variables with moment
                     provided as tab-separated <k> columns, followed by a
                     single moment value (defaults to './meanW.tsv')
```

To see the computations on a simple example, download and unzip SimpleExample.zip and having bin.tsv, types.tsv, meanW.tsv and rho in current directory, run the following:
```
user@localhost:~/SimpleExample$ java -jar /path/to/Idhim-1.0-full.jar -v
```
