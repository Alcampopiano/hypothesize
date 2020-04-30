# Hypothesize: robust statistics in Python

![Screenshot](img/dist_overlay.png)

Hypothesize is a robust statistics library for 
Python based on Rand R. Wilcox's R package [WRS](https://dornsife.usc.edu/labs/rwilcox/software/). 
With Hypothesize you can compare groups and 
measure associations using methods that outperform 
traditional statistcal approaches in terms of power 
and accuracy. 

For more information on robust methods please see Wilcox's book 
[Introduction to Robust Estimation and Hypothesis Testing](https://play.google.com/store/books/details?id=8f8nBb4__EYC&gl=ca&hl=en-CA&source=productsearch&utm_source=HA_Desktop_US&utm_medium=SEM&utm_campaign=PLA&pcampaignid=MKTAD0930BO1&gclid=CjwKCAiA44LzBRB-EiwA-jJipJzyqx9kwNMq5MMU7fG2RrwBK9F7sirX4pfhS8wO7k9Uz_Sqf2P28BoCYzcQAvD_BwE&gclsrc=aw.ds).

## Getting Started

- [Overview](overview.md#overview)
- [Installation](install_dep.md#installation)
- [Dependencies](install_dep.md#dependencies)
- [Basic Tutorial](basic_tutorial.md#basic-tutorial)

## Function Guide

### [Comparing groups with a single factor](function_guide.md#Comparing-groups-with-a-single-factor)

#### independent groups
- [l2drmci](function_guide.md#l2drmci)
- [linconb](function_guide.md#linconb)
- [pb2gen](function_guide.md#pb2gen)
- [tmcppb](function_guide.md#tmcppb)
- [yuenbt](function_guide.md#yuenbt)

#### dependent groups
- [bootdpci](function_guide.md#bootdpci)
- [rmmcppb](function_guide.md#rmmcppb)
- [l2drmci](function_guide.md#l2drmci)
- [lindepbt](function_guide.md#lindepbt)
- [ydbt](function_guide.md#ydbt)

### [Comparing groups with two factors](function_guide.md#comparing-groups-with-two-factors)

#### dependent groups
- [wwmcppb](function_guide.md#wwmcppb)
- [wwmcpbt](function_guide.md#wwmcpbt)

#### mixed designs
- [bwamcp](function_guide.md#bwamcp)  
- [bwbmcp](function_guide.md#bwbmcp) 
- [bwcmp](function_guide.md#bwcmp)
- [bwimcp](function_guide.md#bwimcp)
- [bwmcppb](function_guide.md#bwmcppb)  
- [spmcpa](function_guide.md#spmcpa)
- [spmcpb](function_guide.md#spmcpb)
- [spmcpi](function_guide.md#spmcpi)

### [Measuring associations](function_guide.md#measuring-associations)
- [corb](function_guide.md#corb)
- [pball](function_guide.md#pball)
- [pbcor](function_guide.md#pbcor)
- [winall](function_guide.md#winall)
- [wincor](function_guide.md#wincor)

## Bug reports and Questions
Hypothesize is BSD-licenced and the source code is available
on [GitHub](https://github.com/Alcampopiano/hypothesize).
For issues and questions, 
please use [GitHub Issues](https://github.com/Alcampopiano/hypothesize/issues)

