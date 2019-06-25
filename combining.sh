#!/bin/bash

echo $1;

matlab2016b2 -nodisplay -nosplash -r "runOSort_channelmerge_2(pwd,'$1')";








exit




