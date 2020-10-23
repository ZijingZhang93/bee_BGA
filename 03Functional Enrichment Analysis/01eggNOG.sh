#!/bin/bash
nohup /home/share/biosoft/eggnog-mapper/emapper.py -i ./genelist.fa --output EGGNOG -m diamond --cpu 8 --translate &