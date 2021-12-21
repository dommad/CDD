# CDD - Common Decoy Distribution Project

The source code of PeptideProphet with CDD negative models (PP-CDD) is located in "PeptideProphet + CDD" directory.

To install PP-CDD, follow the instructions:


# Compile TPP 5.1 with CDD

## Compile TPP 5.1

First, following the [official wiki page](http://tools.proteomecenter.org/wiki/index.php?title=TPP_5.1.0:_Installing_on_Ubuntu_16.04_LTS) to compile TPP 5.1. The following steps are included. Please refer to the official wiki page for detailed explanation for each step. 

http://tools.proteomecenter.org/wiki/index.php?title=TPP_5.1.0:_Installing_on_Ubuntu_16.04_LTS

- Installing prerequisite packages

    ```bash
    sudo apt update
    sudo apt --yes upgrade
    sudo apt --yes install subversion
    sudo apt --yes install make
    sudo apt --yes install g++
    sudo apt --yes install build-essential
    sudo apt --yes install zlib1g-dev
    sudo apt --yes install libghc-bzlib-dev
    sudo apt --yes install gnuplot
    sudo apt --yes install unzip
    sudo apt --yes install expat
    sudo apt --yes install libexpat1-dev
    ```

- Creating a suitable place to compile and install

  ```bash
  sudo mkdir /local
  cd /local
  sudo mkdir tpp data svn
  sudo chown tpp.tpp tpp data svn
  ```

  

- Pulling the TPP 5.1.0 source code from SourceForge

  ```bash
  cd /local/svn
  svn checkout svn://svn.code.sf.net/p/sashimi/code/tags/release_5-1-0
  ```

  

- Compiling the source code.

    ```bash
    cd /local/svn/release_5-1-0
    cat > site.mk
    INSTALL_DIR = /local/tpp
    TPP_BASEURL = /tpp
    TPP_DATADIR = /local/data
    ```

    hit CTRL-D here. Then we can compile the source code with following commands.

    ```bash
    make libgd
    make all
    make install
    ```

    Then make sure Perl modules are installed.

    ```bash
    sudo cpan
          (answer yes)
      make install
      install Bundle::CPAN
          (this takes a long time. At one point you need to hit [ENTER] to accept the 'exit')
      install CGI
      install XML::Parser
      install FindBin::libs
      install JSON
           (you may need to answer 'y' to petulant question)
      quit
    ```

    Now test and make sure all the modules are installed.

    ```bash
    cd /local/tpp/cgi-bin
    ./tpp_gui.pl
    ```

    

    

## Compile with CDD 

Now, get the source code file **TPP_5.1.0-src.zip** , unzip it .

```bash
unzip TPP_5.1.0-src.zip
```

Then copy the subfolder *src* to TPP5.1 source code folder */local/svn/release_5-1-0* as follows.

```bash
cp TPP_5.1.0-src/src /local/svn/release_5-1-0/
```

Then run the following command to compile TPP 5.1 with CDD support. 

```bash
make all
```

To install the TPP to */local/tpp*, run the following command:

```bash
make install
```



## Run with  command line Option -OC

In xinteract, use *-OC</path/to/param.txt*> to specify a parameter file for fixed Gumbel distribution. For example.

```bash
/local/tpp/bin/xinteract  -p0 -OMNEA -OC</path/to/param.txt> -N<outputName>  <input.pep.xml>
```

(To be added) We also provide a docker image named as *tpp5.1cdd_option:1.1*. Here is the command line to use this image. 

```bash
docker run -v `pwd`:`pwd` -v </path/to/fasta>:</path/to/fasta> -v --workdir `pwd`  -it tpp5.1cdd_option:1.2 /local/tpp/bin/xinteract  -p0 -OMNEA -OC/data/wulong/param.txt -N<outputName>  <input.pep.xml> 
```

 

