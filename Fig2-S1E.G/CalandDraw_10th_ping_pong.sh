#!/usr/bin/env bash
#CalandDraw_10th_ping_pong.sh
#This is a ping-pong analysis pipeline to calculate 10th seld ping-pong from BED2, and draw the figure.
#Version: Yu Sun, 2016-12-02
#Version: Yu Sun, 2017-02-15, add all-zero tolerant

if [ ! -n "$1" ]
then
  echo "    This is a ping-pong analysis pipeline to calculate 10th seld ping-pong from BED2, and draw the figure."
  echo "    Usage: `basename $0` [BED2]"
  echo "    [Input]: any bed2 files that can be calculated for self ping-pong score"
  echo "    [Output]: The 10th ping-pong score will be outputted, and another [BED2].pp.pdf file will be generated"
#  echo "    [Warning]: Please make sure the ping-pong score can be calculated, or (e.g. all 0 values) the program won't stop"
else
  module load r
  Data=$1
  echo "1. Calculating self ping-pong score"
  ./piPipes_local_ping_pong -a $Data -b $Data > $Data.pp
  echo "2. Drawing the figure using R"
  echo "   The 10th Z-score is:"
  Rscript --vanilla Draw_10th_ping-pong.R $Data.pp
  echo "3. Done"
fi
