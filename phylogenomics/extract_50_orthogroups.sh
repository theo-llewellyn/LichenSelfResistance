#!/bin/bash

#CHANGE TAXON $NUMBER IF 79 taxon tree or 117 taxon tree
NUMBER=79
mkdir Single_Copy_Orthologues_50percent

cat Orthologues_Leca${NUMBER}T_50p.txt | while read line; do cp Orthogroup_Sequences/${line}.fa Single_Copy_Orthologues_50percent; done

