# -*- coding: utf-8 -*-
"""
@author: Feng Ju
@email: richieju520@gmail.com
The script was written and tested in python 2.7 and biopython 1.58

Copyright (C) 2018 Feng Ju

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

"""

import sys
import os
try:  
    from Bio import SeqIO
    from Bio.Seq import reverse_complement
except:
    print "Please install biopython before continue!"

### N50 calculation
def main(fn1, ft1):

    f=open(ft1,'w')
    f.write(';'.join(["","no.of.seq","N50","N80","Max.len","Avg.len","Min.len","no.of.bases","A%","T%","C%","G%","N%"]) + '\n')
     
    for root,dirs,files in os.walk(fn1):
        for file in files:
            if file.endswith((".fa",".fasta","fna")):
                print '------','Processing',file,'in prograss','------'
                BaseSum,Length= 0,[]
                ValueSum,N = 0,[]
                no_c,no_g,no_a,no_t,no_n = 0,0,0,0,0
                
                for record in SeqIO.parse(os.path.join(root, file), 'fasta'):    
                    BaseSum += len(record.seq)
                    Length.append(len(record.seq))
                    seq =record.seq.lower()
                    no_c+=seq.count('c')
                    no_g+=seq.count('g')
                    no_a+=seq.count('a')
                    no_t+=seq.count('t')
                    no_n+=seq.count('n')
        
                #N50 calcuation
                N50_pos = float(BaseSum / 2.0)
                N80_pos = float(BaseSum / 1.25)
                
                Length.sort()   
                Length.reverse()
                average=sum(Length)/len(Length)
        
                for value in Length:
                    ValueSum += value
                    if N50_pos <= float(ValueSum):
                        N.append(value)              
                        break
                    else:
                        continue
        
                per_A = float(no_a*100)/BaseSum
                per_T = float(no_t*100)/BaseSum
                per_C = float(no_c*100)/BaseSum
                per_G = float(no_g*100)/BaseSum
                per_N = float(no_n*100)/BaseSum
        
                f.write(';'.join([file, str(len(Length)),str(N[0]),str(N[-1]),str(max(Length)),str(average),str(min(Length)),str(BaseSum),str(per_A),str(per_T),str(per_C),str(per_G),str(per_N)]) + '\n')
            else:
                print 'File without .fa, .fna, .fasta ingored: ',file
    f.close()
            
if __name__=="__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Calculate N50|N80 length, GC content of fasta files, e.g., from DNA assembly")
    parser.add_argument("-i", dest="input", required=True, help="The folder contains fasta files")
    parser.add_argument("-o", dest="output", required=True, help="The statistics for N50, length, GC content")

    args = parser.parse_args()
    main(args.input,args.output)
