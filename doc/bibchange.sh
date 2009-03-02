#!/bin/sh
#
# Process the .bbl file that BibTeX gives to make the 
# following changes:
#  - {\it -> \textit{
#  - {\tt -> \texttt{
#  - {\bf -> \textbf{
#  - {\em -> \emph{
#
# so that the JHEP bib style can be used with memoir.
#
# Alberto Ramos <alberto@martin.ft.uam.es> 2006
#

cat $1 | sed "s/\{\\\it/\\\textit\{/g" | sed "s/\{\\\tt/\\\texttt\{/g" | sed "s/\{\\\bf/\\\textbf\{/g" | sed "s/\{\\\em/\\\emph\{/g" 

